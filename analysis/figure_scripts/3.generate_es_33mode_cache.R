#!/usr/bin/env Rscript
# ==============================================================================
# Generate .es_33mode_cache.rds (global effect-size forest data, hop0)
#
# This is the generator for the cache consumed by plot_global_es_forest.R
# (Fig 5b / Fig 6b). It computes the rank-biserial effect size of each
# looplook assignment mode's target genes versus the non-assigned background,
# plus a single ChIPseeker nearest-TSS comparison, with bootstrap CIs.
#
# NOTE on naming: "33Mode" is a historical name. The hop0 cache actually
# contains 16 looplook modes + 1 ChIPseeker = 17 hop0 rows (the ChIPseeker row
# is duplicated for hop1 so the older 33-mode renderer could reuse the file).
# plot_global_es_forest.R filters Hop == "hop0".
#
# Inputs (environment variables):
#   LOOPLOOK_OUT_BASE      : res_v13 output root (size200 etc. live under it)
#   LOOPLOOK_RDATA_FILE    : tss5000 looplook annotation RData (tss_5000_res_chromatin.RData)
#   LOOPLOOK_DIFF_FILE     : DESeq2 differential table (log2FoldChange column)
#   LOOPLOOK_PEAK_FILE     : factor ChIP-seq narrowPeak (BRD4/FOSL2)
#   LOOPLOOK_R_LIB         : optional non-default R library path
#
# Output: $LOOPLOOK_OUT_BASE/size200/.es_33mode_cache.rds
#
# Usage:
#   Rscript 3.generate_es_33mode_cache.R
# ==============================================================================

.libPaths(c("~/R/library", .libPaths()))
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(boot)
  library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db); library(parallel)
})

out_base  <- looplook_env_path("LOOPLOOK_OUT_BASE")
rdata_file <- looplook_env_path("LOOPLOOK_RDATA_FILE")
diff_file  <- looplook_env_path("LOOPLOOK_DIFF_FILE")
peak_file  <- looplook_env_path("LOOPLOOK_PEAK_FILE")
n_cores    <- max(1, detectCores() - 1)

# One cache per size directory; the manuscript primary analysis uses size200.
SIZE <- "size200"

rank_biserial <- function(x, y, R = 2000L, seed = 42L) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(est = NA_real_, lo = NA_real_, hi = NA_real_))
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(est = NA_real_, lo = NA_real_, hi = NA_real_))
  wt <- wilcox.test(x, y, exact = FALSE)
  U <- as.numeric(wt$statistic)
  r <- 2 * U / (n1 * n2) - 1
  boot_fn <- function(d, i) {
    xb <- d$v[i[1:n1]]; yb <- d$v[i[(n1 + 1):(n1 + n2)]]
    if (length(xb) < 2 || length(yb) < 2) return(NA_real_)
    wtb <- wilcox.test(xb, yb, exact = FALSE)
    2 * as.numeric(wtb$statistic) / (length(xb) * length(yb)) - 1
  }
  set.seed(seed)
  bd <- data.frame(v = c(x, y), stringsAsFactors = FALSE)
  br <- tryCatch(boot::boot(bd, boot_fn, R = R, strata = rep(c(1, 2), c(n1, n2))), error = function(e) NULL)
  ci_obj <- if (!is.null(br)) {
    tryCatch(boot::boot.ci(br, type = "bca", conf = 0.95),
             error = function(e) boot::boot.ci(br, type = "perc", conf = 0.95))
  } else NULL
  if (!is.null(ci_obj)) {
    if (!is.null(ci_obj$bca)) list(est = r, lo = ci_obj$bca[4], hi = ci_obj$bca[5])
    else if (!is.null(ci_obj$percent)) list(est = r, lo = ci_obj$percent[4], hi = ci_obj$percent[5])
    else list(est = r, lo = NA_real_, hi = NA_real_)
  } else list(est = r, lo = NA_real_, hi = NA_real_)
}

get_target_col <- function(md) {
  bc <- if (identical(md$map, "promoter")) "Regulated_promoter_genes" else "Assigned_Target_Genes"
  if (md$fill) paste0(bc, "_Filled") else bc
}
get_obj <- function(nm) {
  if (grepl("^anno_", nm)) "res"
  else if (grepl("^refined_", nm)) "refined_res"
  else if (grepl("^chrom_only_", nm)) "cr_only"
  else if (grepl("^chrom_", nm)) "cr"
  else NA_character_
}

# ---- load inputs ----
message("Loading...")
tmp_env <- new.env(parent = emptyenv())
load(rdata_file, envir = tmp_env)
res <- tmp_env$res; refined_res <- tmp_env$refined_res
cr <- tmp_env$cr; cr_only <- tmp_env$cr_only
rm(tmp_env)

diff_df <- read.csv(diff_file, stringsAsFactors = FALSE, row.names = 1)
gs <- diff_df[!is.na(diff_df$log2FoldChange) & is.finite(diff_df$log2FoldChange), ]
gs$gene <- toupper(rownames(gs))
gs$signal_lfc <- -gs$log2FoldChange

peak_anno <- ChIPseeker::annotatePeak(peak_file, tssRegion = c(-5000, 5000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db")
padf <- as.data.frame(peak_anno)
padf$SYMBOL <- toupper(padf$SYMBOL)
chip_genes <- unique(na.omit(padf$SYMBOL))

# ---- build 16 hop0 modes ----
mode_names <- c(
  "anno_all_F", "anno_all_T", "anno_promoter_F", "anno_promoter_T",
  "refined_all_F", "refined_all_T", "refined_promoter_F", "refined_promoter_T",
  "chrom_all_F", "chrom_all_T", "chrom_promoter_F", "chrom_promoter_T",
  "chrom_only_all_F", "chrom_only_all_T", "chrom_only_promoter_F", "chrom_only_promoter_T"
)
modes <- list(); mode_genes <- list()
for (nm in mode_names) {
  ann_obj <- get(get_obj(nm))
  md <- list(name = nm, ann = ann_obj,
             map = if (grepl("_promoter_", nm)) "promoter" else "all",
             fill = grepl("_T(_hop1)?$", nm), near = FALSE)
  modes[[length(modes) + 1]] <- md
  bi <- ann_obj$target_annotation
  dc <- get_target_col(md)
  mode_genes[[nm]] <- if (!is.null(bi) && dc %in% colnames(bi))
    unique(toupper(looplook:::clean_gene_names(bi[[dc]], "[;,]"))) else character(0)
}

# ---- compute ES ----
out_dir <- file.path(out_base, SIZE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_fn <- file.path(out_dir, ".es_33mode_cache.rds")
if (file.exists(cache_fn)) {
  message("Cache exists, loading: ", cache_fn)
  df <- readRDS(cache_fn)
} else {
  lo_es <- parallel::mclapply(mode_names, function(nm) {
    loop_g <- intersect(unique(toupper(mode_genes[[nm]])), gs$gene)
    bg_g <- setdiff(gs$gene, union(loop_g, intersect(chip_genes, gs$gene)))
    if (length(loop_g) >= 5 && length(bg_g) >= 5) {
      x <- gs$signal_lfc[gs$gene %in% loop_g]
      y <- gs$signal_lfc[gs$gene %in% bg_g]
      rb <- rank_biserial(x, y)
      data.frame(Mode = nm, Type = "looplook", ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi,
                 stringsAsFactors = FALSE)
    } else NULL
  }, mc.cores = min(8L, n_cores))
  lo_es <- lo_es[!vapply(lo_es, is.null, logical(1))]

  # ChIPseeker global ES (background excludes all looplook targets)
  chip_g <- intersect(chip_genes, gs$gene)
  all_loop <- unique(unlist(lapply(mode_genes, function(g) intersect(g, gs$gene))))
  chip_bg <- setdiff(gs$gene, union(all_loop, chip_g))
  ch_row <- NULL
  if (length(chip_g) >= 5 && length(chip_bg) >= 5) {
    x <- gs$signal_lfc[gs$gene %in% chip_g]
    y <- gs$signal_lfc[gs$gene %in% chip_bg]
    ch_rb <- rank_biserial(x, y)
    ch_row <- data.frame(Mode = "ChIPseeker", Type = "chipseeker",
                         ES = ch_rb$est, ES_lo = ch_rb$lo, ES_hi = ch_rb$hi,
                         stringsAsFactors = FALSE)
  }

  if (length(lo_es) > 0) {
    df <- bind_rows(c(lo_es, if (!is.null(ch_row)) list(ch_row) else list())) %>%
      mutate(Pipeline = ifelse(Type == "chipseeker", "ChIPseeker",
        as.character(factor(dplyr::case_when(
          grepl("^anno_", Mode) ~ "Basic",
          grepl("^refined_", Mode) ~ "E-refined",
          grepl("^chrom_only_", Mode) ~ "C-refined",
          grepl("^chrom_", Mode) ~ "I-refined",
          TRUE ~ "Other"),
          levels = c("Basic", "E-refined", "C-refined", "I-refined")))),
        Hop = "hop0")
    # duplicate ChIPseeker row for hop1 so the historical 33-mode renderer works
    ch_rows <- df %>% filter(Type == "chipseeker")
    if (nrow(ch_rows) > 0) {
      ch_rows$Hop <- "hop1"
      df <- bind_rows(df, ch_rows)
    }
    saveRDS(df, cache_fn)
    message("ES cached: ", cache_fn)
  } else {
    stop("No looplook modes had enough genes to compute ES.")
  }
}

message(sprintf("size200: %d looplook + ChIPseeker (hop0 rows = %d)",
                sum(df$Type == "looplook" & df$Hop == "hop0"),
                sum(df$Hop == "hop0")))
message("Next: run plot_global_es_forest.R to render Global_ES_Forest_33Mode.pdf")