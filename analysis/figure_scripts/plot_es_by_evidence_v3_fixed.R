#!/usr/bin/env Rscript
# ES by Evidence Type: looplook target genes stratified by assignment evidence
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(openxlsx); library(boot)
  library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
})

# ═══ CONFIG ──────────────────────────────────────────────────

out_base   = looplook_env_path("LOOPLOOK_OUT_BASE")
rdata = looplook_env_path("LOOPLOOK_RDATA_FILE")
diff_csv  = looplook_env_path("LOOPLOOK_DIFF_FILE")
peak_file  = looplook_env_path("LOOPLOOK_PEAK_FILE")
expr_csv = looplook_env_path("LOOPLOOK_EXPR_FILE")

# ═══════════════════════════════════════════════════════════════

rank_biserial_point <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  wt <- wilcox.test(x, y, exact = FALSE)
  2 * as.numeric(wt$statistic) / (length(x) * length(y)) - 1
}

rank_biserial <- function(x, y, R = 500L, seed = 42L) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(est = NA_real_, lo = NA_real_, hi = NA_real_))
  wt <- wilcox.test(x, y, exact = FALSE)
  U <- as.numeric(wt$statistic); r <- 2 * U / (n1 * n2) - 1
  boot_fn <- function(d, i) {
    xb <- d$v[i[1:n1]]; yb <- d$v[i[(n1+1):(n1+n2)]]
    if (length(xb) < 2 || length(yb) < 2) return(NA_real_)
    2 * as.numeric(wilcox.test(xb, yb, exact = FALSE)$statistic) / (length(xb) * length(yb)) - 1
  }
  set.seed(seed)
  bd <- data.frame(v = c(x, y), stringsAsFactors = FALSE)
  br <- tryCatch(boot::boot(bd, boot_fn, R = R, strata = rep(c(1, 2), c(n1, n2))),
                 error = function(e) NULL)
  ci_obj <- if (!is.null(br))
    tryCatch(boot::boot.ci(br, type = "perc", conf = 0.95), error = function(e) NULL)
  else NULL
  if (!is.null(ci_obj) && !is.null(ci_obj$percent))
    list(est = r, lo = ci_obj$percent[4], hi = ci_obj$percent[5])
  else list(est = r, lo = NA_real_, hi = NA_real_)
}

message("Loading...")
load(rdata, verbose = FALSE)
diff_df <- read.csv(diff_csv, row.names = 1)
gs <- diff_df[!is.na(diff_df$log2FoldChange) & is.finite(diff_df$log2FoldChange), ]
gs$gene <- toupper(rownames(gs))
gs$signal_lfc <- -gs$log2FoldChange

padf <- as.data.frame(annotatePeak(peak_file, tssRegion = c(-5000, 5000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db"))
padf$SYMBOL <- toupper(padf$SYMBOL)
chip_genes <- unique(na.omit(padf$SYMBOL))

# ── Expression adjustment model ──
library(splines)
expr_mat <- as.matrix(read.table(expr_csv, header = TRUE, row.names = 1, check.names = FALSE))
dmso_cols <- grep("dmso|DMSO|Dmso", colnames(expr_mat), ignore.case = TRUE, value = TRUE)
expr_lookup <- setNames(rowMeans(expr_mat[, dmso_cols, drop = FALSE], na.rm = TRUE), toupper(rownames(expr_mat)))
expr_lookup <- expr_lookup[is.finite(expr_lookup) & expr_lookup >= 0]

analysis_universe <- intersect(gs$gene, names(expr_lookup))
model_df <- data.frame(
  gene = analysis_universe,
  signal_lfc = gs$signal_lfc[match(analysis_universe, gs$gene)],
  baseline_tpm = expr_lookup[analysis_universe]) %>%
  mutate(log_expr = log2(baseline_tpm + 0.1)) %>%
  filter(is.finite(signal_lfc), is.finite(log_expr))
fit <- lm(signal_lfc ~ ns(log_expr, df = 4), data = model_df)
model_df$residual_signal <- residuals(fit)
resid_lookup <- setNames(model_df$residual_signal, model_df$gene)
message("Spline model R2 = ", round(summary(fit)$r.squared, 3))

# Mode helpers
get_obj <- function(nm) {
  if (grepl("^anno_", nm))      return(if (grepl("_hop1$", nm)) "res2" else "res")
  if (grepl("^refined_", nm))   return(if (grepl("_hop1$", nm)) "refined_res2" else "refined_res")
  if (grepl("^chrom_only_", nm)) return(if (grepl("_hop1$", nm)) "cr2_only" else "cr_only")
  if (grepl("^chrom_", nm))     return(if (grepl("_hop1$", nm)) "cr2" else "cr")
  NA_character_
}
get_target_col <- function(md) {
  bc <- if (identical(md$map, "promoter")) "Regulated_promoter_genes" else "Assigned_Target_Genes"
  if (md$fill) paste0(bc, "_Filled") else bc
}

mode_names <- c(
  "anno_all_F","anno_all_T","refined_all_F","refined_all_T",
  "chrom_all_F","chrom_all_T",
  "chrom_only_all_F","chrom_only_all_T",
  "anno_promoter_F","anno_promoter_T","refined_promoter_F","refined_promoter_T",
  "chrom_promoter_F","chrom_promoter_T",
  "chrom_only_promoter_F","chrom_only_promoter_T")

# Extract evidence per mode
results <- list()
message("Computing ES by evidence...")
pb <- txtProgressBar(max = length(mode_names), style = 3)
for (ki in seq_along(mode_names)) {
  nm <- mode_names[[ki]]
  ao <- get_obj(nm); ann_obj <- if (!is.na(ao)) get(ao) else NULL
  if (is.null(ann_obj)) next

  md <- list(name = nm, ann = ann_obj,
    map = if (grepl("_promoter_", nm)) "promoter" else "all",
    fill = grepl("_T($|_hop1)", nm))

  bi <- ann_obj$target_annotation; tgl <- ann_obj$target_gene_links
  if (is.null(bi) || is.null(tgl) || nrow(tgl) == 0) next

  dc <- get_target_col(md)
  if (!dc %in% colnames(bi) || !all(c("gene","evidence","source") %in% colnames(tgl))) next

  # Build peak-gene pairs from target_annotation
  genes_assigned <- sort(unique(toupper(looplook:::clean_gene_names(bi[[dc]], "[;,]"))))
  peak_gene_pairs <- lapply(seq_len(nrow(bi)), function(i) {
    gs_i <- unique(toupper(looplook:::clean_gene_names(bi[[dc]][i], "[;,]")))
    gs_i <- gs_i[!is.na(gs_i) & gs_i != "" & gs_i != "NA"]
    if (length(gs_i) == 0) return(NULL)
    data.frame(Input_ID = as.character(bi$input_id[i]), Gene = gs_i, stringsAsFactors = FALSE)
  }) %>% bind_rows() %>% distinct(Input_ID, Gene)

  if (nrow(peak_gene_pairs) == 0) next

  # Normalize target_gene_links (keep all evidence, filter after join)
  tgl_norm <- tgl %>%
    transmute(
      Input_ID = as.character(input_id),
      Gene = toupper(trimws(as.character(gene))),
      Evidence = as.character(evidence),
      Source = as.character(source),
      GeneRole = if ("gene_role" %in% colnames(tgl)) as.character(gene_role) else NA_character_) %>%
    filter(!is.na(Gene), Gene != "", Gene != "NA")

  # Join: only keep pairs in the actual target column
  out <- inner_join(tgl_norm, peak_gene_pairs, by = c("Input_ID", "Gene"))

  # F (strict) modes: restrict to loop_anchor source
  if (!md$fill) out <- out %>% filter(Source == "loop_anchor")
  # Promoter strict modes only: restrict to promoter gene_role
  if (md$map == "promoter" && !md$fill)
    out <- out %>% filter(!is.na(GeneRole), GeneRole == "promoter")

  # Evidence priority: one primary evidence per peak-gene pair
  evidence_priority <- c(
    "local_promoter_overlap"=1, "direct_opposite_promoter"=2,
    "local_promoter"=3, "local_enhancer_candidate"=4,
    "gene_body_context"=5, "distal_promoter"=6,
    "distal_gene_body_context"=7, "distal_enhancer_candidate"=8,
    "linear_annotation"=9, "linear_fallback"=10, "positional_candidate"=11)
  out$EvPriority <- evidence_priority[out$Evidence]
  out$EvPriority[is.na(out$EvPriority)] <- 99L
  out <- out[order(out$Input_ID, out$Gene, out$EvPriority), ]
  out <- out[!duplicated(paste(out$Input_ID, out$Gene, sep = "\r")), ]

  # Drop none/other (after priority assignment)
  out <- out %>% filter(!Evidence %in% c("none", "other"))

  # Gene-level multi-label allowed
  ev_groups <- out %>%
    distinct(Gene, Evidence) %>%
    filter(Gene %in% gs$gene)

  if (ki == 1) {
    cat(sprintf("\n  [%s] pairs=%d tgl=%d ev_rows=%d ev_genes=%d\n",
      nm, nrow(peak_gene_pairs), nrow(tgl_norm), nrow(out), nrow(ev_groups)))
    cat("  evidence counts:\n"); print(sort(table(out$Evidence), decreasing=TRUE))
    flush.console()
  }

  # ── Multi-evidence QC ──
  gene_n_ev <- ev_groups %>% count(Gene, name = "n_evidence")
  message(sprintf("  %s: %d genes, multi-ev: %d (%.1f%%)",
    nm, nrow(gene_n_ev), sum(gene_n_ev$n_evidence > 1),
    100 * mean(gene_n_ev$n_evidence > 1)))

  if (nrow(ev_groups) == 0) next

  pipeline <- case_when(
    grepl("^anno_", nm) ~ "Basic", grepl("^refined_", nm) ~ "E-refined",
    grepl("^chrom_only_", nm) ~ "C-refined", grepl("^chrom_", nm) ~ "I-refined")

  mode_targets <- intersect(genes_assigned, analysis_universe)
  bg_full <- setdiff(analysis_universe, mode_targets)

  ev_size_qc <- ev_groups %>%
    filter(Gene %in% analysis_universe) %>%
    distinct(Gene, Evidence) %>%
    count(Evidence, name = "N_genes") %>%
    arrange(desc(N_genes))
  message("  ", nm, " evidence sizes: ",
          paste0(ev_size_qc$Evidence, "=", ev_size_qc$N_genes, collapse = "; "))

  for (ev in unique(ev_groups$Evidence)) {
    ev_genes <- ev_groups$Gene[ev_groups$Evidence == ev]
    ev_genes <- intersect(ev_genes, analysis_universe)
    if (length(ev_genes) < 10) next

    # Full ES
    x_full <- gs$signal_lfc[gs$gene %in% ev_genes]
    y_full <- gs$signal_lfc[gs$gene %in% bg_full]
    rb_full <- rank_biserial(x_full, y_full)

    # Expression-adjusted ES
    ra <- resid_lookup[ev_genes]; ra <- ra[is.finite(ra)]
    rb <- resid_lookup[intersect(bg_full, names(resid_lookup))]; rb <- rb[is.finite(rb)]
    rb_adj <- if (length(ra) >= 10 && length(rb) >= 10)
      rank_biserial(ra, rb, seed = 4242) else list(est = NA_real_, lo = NA_real_, hi = NA_real_)

    results[[length(results)+1]] <- data.frame(
      Mode = nm, Pipeline = pipeline, Evidence = ev,
      N_genes = length(ev_genes),
      ES = rb_full$est, ES_lo = rb_full$lo, ES_hi = rb_full$hi,
      ES_adj = rb_adj$est, ES_adj_lo = rb_adj$lo, ES_adj_hi = rb_adj$hi,
      stringsAsFactors = FALSE)
  }
  setTxtProgressBar(pb, ki)
}
close(pb)

res <- bind_rows(results)
if (nrow(res) == 0) stop("No results")

# Color palette
ev_colors <- c(
  "local_promoter_overlap"    = "#2166AC",
  "direct_opposite_promoter"  = "#053061",
  "local_promoter"            = "#74ADD1",
  "distal_promoter"           = "#4393C3",
  "gene_body_context"         = "#92C5DE",
  "distal_gene_body_context"  = "#7B8FD4",
  "local_enhancer_candidate"  = "#F4A582",
  "distal_enhancer_candidate" = "#E64B35",
  "linear_annotation"         = "#D6604D",
  "linear_fallback"           = "#CA0020",
  "positional_candidate"      = "#FFA500")

unexpected_ev <- setdiff(unique(as.character(res$Evidence)), names(ev_colors))
if (length(unexpected_ev) > 0L) {
  stop("Unexpected evidence labels not present in ev_colors: ",
       paste(unexpected_ev, collapse = ", "))
}
res$Evidence <- factor(res$Evidence, levels = names(ev_colors))
res$Pipeline <- factor(res$Pipeline, levels = c("Basic","E-refined","C-refined","I-refined"))

# Plot: Full ES (diamond) + Expression-adjusted ES (circle)
p <- ggplot(res, aes(x = ES, y = reorder(Mode, ES), color = Evidence)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  # Full ES
  geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.25, linewidth = 0.5, alpha = 0.5) +
  geom_point(size = 2, shape = 18, alpha = 0.8) +
  # Expression-adjusted ES
  geom_errorbarh(aes(x = ES_adj, xmin = ES_adj_lo, xmax = ES_adj_hi),
    height = 0.25, linewidth = 0.35, alpha = 0.4, linetype = "dotted") +
  geom_point(aes(x = ES_adj), size = 1.6, shape = 1, stroke = 0.7, alpha = 0.7) +
  scale_color_manual(values = ev_colors) +
  facet_wrap(~ Pipeline, ncol = 4, scales = "free_y") +
  labs(title = "Looplook Target ES by Assignment Evidence (hop0, vs BG)",
       subtitle = "filled diamond = Full ES, open circle = Expression-adjusted ES | rank-biserial, perc bootstrap CI",
       x = "Rank-Biserial Effect Size (vs BG)", y = NULL, color = "Evidence") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5))

out_dir <- file.path(out_base, "expression_sensitivity_v5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "ES_by_Evidence.pdf"), p, width = 12, height = 8)

# Excel
wb <- createWorkbook()
addWorksheet(wb, "ES_by_Evidence"); writeData(wb, "ES_by_Evidence", as.data.frame(res))
saveWorkbook(wb, file.path(out_dir, "ES_by_Evidence.xlsx"), overwrite = TRUE)
message("Done: ", out_dir)
