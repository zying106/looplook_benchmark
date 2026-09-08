#!/usr/bin/env Rscript
# =============================================================================
# looplook BRD4 expression-sensitivity integrated analysis v4
#
# Four-layer analysis for looplook target validation:
#   1) Full ES                : original target sets, no expression adjustment
#   2) Common-support ES      : remove non-overlapping baseline-expression tails only
#   3) Expression-adjusted ES : compare residual BRD4 response after spline adjustment
#   4) Strict-matched ES      : exact equal-count matching within expression bins
#
# Interpretation:
#   - Full ES is the descriptive/full-set effect.
#   - Expression-adjusted ES is the main expression sensitivity analysis.
#   - Common-support ES diagnoses whether extreme expression tails drive the result.
#   - Strict-matched ES is a stress test, NOT the primary benchmark.
#
# Compatible with the user's current looplook RData / DESeq2 CSV / TPM / narrowPeak
# structure used by plot_expression_matched_es.R v1-v3.
# =============================================================================

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(openxlsx)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(splines)
})

if (!requireNamespace("looplook", quietly = TRUE)) {
  stop("Package 'looplook' is required.")
}

# =============================================================================
# CONFIG — paths are intentionally kept compatible with the current scripts
# =============================================================================
cfg <- list(

  out_base   = looplook_env_path("LOOPLOOK_OUT_BASE"),
  rdata = looplook_env_path("LOOPLOOK_RDATA_FILE"),
  diff_csv  = looplook_env_path("LOOPLOOK_DIFF_FILE"),
  peak_file  = looplook_env_path("LOOPLOOK_PEAK_FILE"),
  expr_file = looplook_env_path("LOOPLOOK_EXPR_FILE"),
  # Optional metadata. If NULL, DMSO/control columns are identified by regex.
  meta_file = NULL,
  dmso_regex = "dmso|vehicle|control",

  # Current looplook benchmark comparator definition:
  # all looplook modes, including promoter modes, are compared with ChIPseeker-all.
  chip_reference_policy = "all",

  # TRUE = include hop1 correctly as primary + eligible expanded targets.
  # FALSE = legacy 16 hop0 modes, matching the old v1/v2 expression scripts.
  include_hop1 = FALSE,

  # Expression universe.
  # 0 keeps exact compatibility with prior scripts.
  # For a formal expressed-gene sensitivity analysis, 0.5 is also reasonable.
  min_baseline_tpm = 0,

  # Common-support trimming.
  support_quantiles = c(0.005, 0.995),

  # Strict expression matching.
  strict_n_iter = 500L,
  strict_n_bins = 20L,
  strict_max_per_bin = 50L,
  min_valid_matched_iter = 50L,
  smd_target = 0.10,

  # Effect-size bootstrap for Full / Common-support / Expression-adjusted ES.
  # Increase to 1000-2000 for final figures if runtime allows.
  bootstrap_R = 300L,
  bootstrap_seed = 42L,

  # Nonlinear expression adjustment.
  residual_spline_df = 4L
)

# =============================================================================
# Utility functions
# =============================================================================

check_file <- function(path, label) {
  path2 <- path.expand(path)
  if (!file.exists(path2)) stop(label, " not found: ", path2)
  path2
}

cfg$rdata     <- check_file(cfg$rdata, "RData")
cfg$diff_csv  <- check_file(cfg$diff_csv, "DESeq2 CSV")
cfg$expr_file <- check_file(cfg$expr_file, "TPM file")
cfg$peak_file <- check_file(cfg$peak_file, "Peak file")
cfg$out_base  <- path.expand(cfg$out_base)

out_dir <- file.path(cfg$out_base, "expression_sensitivity_v5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_ids <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | !nzchar(x) | x == "NA"] <- NA_character_
  x
}

clean_gene_vector <- function(x) {
  if (length(x) == 0L) return(character(0))
  g <- looplook:::clean_gene_names(x, "[;,]")
  g <- normalize_ids(g)
  sort(unique(g[!is.na(g)]))
}

rank_biserial_point <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  wt <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
  2 * as.numeric(wt$statistic) / (length(x) * length(y)) - 1
}

rank_biserial_boot <- function(x, y, R = cfg$bootstrap_R, seed = cfg$bootstrap_seed) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x)
  n2 <- length(y)

  if (n1 < 2L || n2 < 2L) {
    return(list(est = NA_real_, lo = NA_real_, hi = NA_real_, n1 = n1, n2 = n2))
  }

  est <- rank_biserial_point(x, y)

  if (R <= 0L) {
    return(list(est = est, lo = NA_real_, hi = NA_real_, n1 = n1, n2 = n2))
  }

  set.seed(seed)
  boots <- vapply(seq_len(R), function(i) {
    xb <- sample(x, n1, replace = TRUE)
    yb <- sample(y, n2, replace = TRUE)
    rank_biserial_point(xb, yb)
  }, numeric(1))

  boots <- boots[is.finite(boots)]
  if (length(boots) < max(20L, floor(R * 0.5))) {
    return(list(est = est, lo = NA_real_, hi = NA_real_, n1 = n1, n2 = n2))
  }

  list(
    est = est,
    lo = unname(stats::quantile(boots, 0.025, na.rm = TRUE)),
    hi = unname(stats::quantile(boots, 0.975, na.rm = TRUE)),
    n1 = n1,
    n2 = n2
  )
}

calc_smd <- function(genes_a, genes_b, expr_lookup) {
  ea <- log2(expr_lookup[genes_a] + 0.1)
  eb <- log2(expr_lookup[genes_b] + 0.1)
  ea <- ea[is.finite(ea)]
  eb <- eb[is.finite(eb)]

  if (length(ea) < 2L || length(eb) < 2L) return(NA_real_)

  psd <- sqrt((stats::var(ea) + stats::var(eb)) / 2)
  if (!is.finite(psd) || psd == 0) return(NA_real_)

  (mean(ea) - mean(eb)) / psd
}

common_support_sets <- function(genes_a, genes_b, expr_lookup,
                                q = cfg$support_quantiles) {
  ea <- log2(expr_lookup[genes_a] + 0.1)
  eb <- log2(expr_lookup[genes_b] + 0.1)

  ea <- ea[is.finite(ea)]
  eb <- eb[is.finite(eb)]

  if (length(ea) < 2L || length(eb) < 2L) return(NULL)

  qa <- stats::quantile(ea, q, na.rm = TRUE, names = FALSE)
  qb <- stats::quantile(eb, q, na.rm = TRUE, names = FALSE)

  overlap_min <- max(qa[1], qb[1])
  overlap_max <- min(qa[2], qb[2])

  if (!is.finite(overlap_min) || !is.finite(overlap_max) ||
      overlap_min >= overlap_max) {
    return(NULL)
  }

  ea_all <- log2(expr_lookup[genes_a] + 0.1)
  eb_all <- log2(expr_lookup[genes_b] + 0.1)

  keep_a <- is.finite(ea_all) & ea_all >= overlap_min & ea_all <= overlap_max
  keep_b <- is.finite(eb_all) & eb_all >= overlap_min & eb_all <= overlap_max

  list(
    a = genes_a[keep_a],
    b = genes_b[keep_b],
    overlap_min = overlap_min,
    overlap_max = overlap_max
  )
}

strict_match_once <- function(genes_a, genes_b, expr_lookup,
                              seed,
                              n_bins = cfg$strict_n_bins,
                              max_per_bin = cfg$strict_max_per_bin) {
  match_df <- bind_rows(
    data.frame(
      gene = genes_a,
      group = "A",
      expr = log2(expr_lookup[genes_a] + 0.1),
      stringsAsFactors = FALSE
    ),
    data.frame(
      gene = genes_b,
      group = "B",
      expr = log2(expr_lookup[genes_b] + 0.1),
      stringsAsFactors = FALSE
    )
  ) %>%
    filter(is.finite(expr))

  if (nrow(match_df) < 20L) return(NULL)

  probs <- seq(0, 1, length.out = n_bins + 1L)
  breaks <- unique(as.numeric(stats::quantile(
    match_df$expr,
    probs = probs,
    type = 8,
    na.rm = TRUE,
    names = FALSE
  )))

  if (length(breaks) < 3L) return(NULL)

  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf

  match_df$expr_bin <- cut(
    match_df$expr,
    breaks = breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = FALSE
  )

  bins <- sort(unique(match_df$expr_bin))
  bins <- bins[is.finite(bins)]

  bin_info <- lapply(bins, function(b) {
    ga <- match_df$gene[match_df$group == "A" & match_df$expr_bin == b]
    gb <- match_df$gene[match_df$group == "B" & match_df$expr_bin == b]
    n_match <- min(length(ga), length(gb), max_per_bin)

    if (n_match < 1L) return(NULL)

    list(a = ga, b = gb, n = n_match, bin = b)
  })

  bin_info <- bin_info[!vapply(bin_info, is.null, logical(1))]
  if (length(bin_info) == 0L) return(NULL)

  total_n <- sum(vapply(bin_info, function(z) z$n, integer(1)))
  if (total_n < 10L) return(NULL)

  set.seed(seed)

  sa <- unlist(lapply(bin_info, function(z) {
    sample(z$a, size = z$n, replace = FALSE)
  }), use.names = FALSE)

  sb <- unlist(lapply(bin_info, function(z) {
    sample(z$b, size = z$n, replace = FALSE)
  }), use.names = FALSE)

  if (length(sa) != length(sb)) {
    stop("Internal strict-matching error: N_A != N_B.")
  }

  list(
    a = sa,
    b = sb,
    n = length(sa),
    eligible_bins = length(bin_info)
  )
}

safe_ratio <- function(num, den) {
  ifelse(is.finite(num) & is.finite(den) & den != 0, num / den, NA_real_)
}

# =============================================================================
# Load and validate data
# =============================================================================
message(">>> Loading looplook RData...")
loaded_objects <- load(cfg$rdata, verbose = FALSE)

required_hop0 <- c("res", "refined_res", "cr", "cr_only")
required_hop1 <- c("res2", "refined_res2", "cr2", "cr2_only")
required_objects <- if (isTRUE(cfg$include_hop1)) {
  c(required_hop0, required_hop1)
} else {
  required_hop0
}

missing_objects <- setdiff(required_objects, loaded_objects)
if (length(missing_objects) > 0L) {
  stop("RData missing object(s): ", paste(missing_objects, collapse = ", "))
}

message(">>> Loading differential expression...")
diff_df <- read.csv(cfg$diff_csv, row.names = 1, check.names = FALSE)
if (!"log2FoldChange" %in% colnames(diff_df)) {
  stop("DESeq2 CSV requires column: log2FoldChange")
}

diff_gene <- normalize_ids(rownames(diff_df))
if (anyNA(diff_gene)) stop("DESeq2 CSV contains blank/NA gene row names after normalization.")
if (anyDuplicated(diff_gene)) {
  stop("DESeq2 CSV has duplicated/case-colliding gene IDs after toupper().")
}
rownames(diff_df) <- diff_gene

gs <- diff_df %>%
  mutate(
    gene = rownames(diff_df),
    signal_lfc = -log2FoldChange
  ) %>%
  filter(is.finite(signal_lfc), !is.na(gene), nzchar(gene))

message("    DE genes with finite LFC: ", nrow(gs))

message(">>> Loading baseline TPM...")
expr_mat <- as.matrix(read.table(
  cfg$expr_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "",
  quote = "",
  comment.char = ""
))

expr_gene <- normalize_ids(rownames(expr_mat))
if (anyNA(expr_gene)) stop("TPM file contains blank/NA gene IDs after normalization.")
if (anyDuplicated(expr_gene)) {
  dup <- unique(expr_gene[duplicated(expr_gene)])
  stop(
    "TPM file has duplicated/case-colliding gene IDs after toupper(): ",
    paste(head(dup, 20), collapse = ", ")
  )
}
rownames(expr_mat) <- expr_gene

# DMSO/control columns: metadata first, regex fallback.
dmso_cols <- character(0)

if (!is.null(cfg$meta_file)) {
  meta_path <- check_file(cfg$meta_file, "Metadata file")
  meta_df <- read.table(
    meta_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("sample", "condition") %in% colnames(meta_df))) {
    stop("Metadata requires columns: sample, condition")
  }

  dmso_cols <- as.character(meta_df$sample[
    grepl(cfg$dmso_regex, as.character(meta_df$condition), ignore.case = TRUE)
  ])
  dmso_cols <- intersect(dmso_cols, colnames(expr_mat))
}

if (length(dmso_cols) == 0L) {
  dmso_cols <- grep(cfg$dmso_regex, colnames(expr_mat),
                    ignore.case = TRUE, value = TRUE)
}

if (length(dmso_cols) == 0L) {
  stop("No DMSO/vehicle/control columns found in TPM matrix.")
}

message("    DMSO/control samples: ", paste(dmso_cols, collapse = ", "))

expr_lookup <- setNames(
  rowMeans(expr_mat[, dmso_cols, drop = FALSE], na.rm = TRUE),
  rownames(expr_mat)
)

expr_lookup <- expr_lookup[
  is.finite(expr_lookup) &
  expr_lookup >= cfg$min_baseline_tpm
]

analysis_universe <- intersect(gs$gene, names(expr_lookup))
if (length(analysis_universe) < 100L) {
  stop("Expression + DE analysis universe is unexpectedly small: ",
       length(analysis_universe))
}

message(
  "    Analysis universe: ", length(analysis_universe),
  " genes; min baseline TPM=", cfg$min_baseline_tpm
)

# =============================================================================
# ChIPseeker comparator
# =============================================================================
message(">>> ChIPseeker annotation...")
peak_anno <- as.data.frame(annotatePeak(
  cfg$peak_file,
  tssRegion = c(-5000, 5000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db"
))

peak_anno$SYMBOL <- normalize_ids(peak_anno$SYMBOL)

chip_genes <- sort(unique(na.omit(peak_anno$SYMBOL)))
chip_prom  <- sort(unique(na.omit(
  peak_anno$SYMBOL[grepl("Promoter", peak_anno$annotation)]
)))

message("    ChIPseeker-all genes: ", length(chip_genes))
message("    ChIPseeker-promoter genes: ", length(chip_prom))

# =============================================================================
# looplook mode definitions and correct hop1 expansion
# =============================================================================

get_obj_name <- function(nm) {
  is_hop1 <- grepl("_hop1$", nm)

  if (grepl("^anno_", nm)) {
    return(if (is_hop1) "res2" else "res")
  }
  if (grepl("^refined_", nm)) {
    return(if (is_hop1) "refined_res2" else "refined_res")
  }
  if (grepl("^chrom_only_", nm)) {
    return(if (is_hop1) "cr2_only" else "cr_only")
  }
  if (grepl("^chrom_", nm)) {
    return(if (is_hop1) "cr2" else "cr")
  }

  NA_character_
}

get_target_col <- function(md) {
  base_col <- if (identical(md$map, "promoter")) {
    "Regulated_promoter_genes"
  } else {
    "Assigned_Target_Genes"
  }

  if (isTRUE(md$fill)) paste0(base_col, "_Filled") else base_col
}

get_primary_mode_genes <- function(md) {
  bed_info <- md$ann$target_annotation
  if (is.null(bed_info) || nrow(bed_info) == 0L) return(character(0))

  col <- get_target_col(md)
  if (!col %in% colnames(bed_info)) {
    stop("[", md$name, "] missing target column: ", col)
  }

  clean_gene_vector(bed_info[[col]])
}

get_expanded_mode_genes <- function(md) {
  if (!grepl("_hop1$", md$name)) return(character(0))

  tgl <- md$ann$target_gene_links
  required <- c(
    "gene", "gene_role", "source", "path_length", "in_expanded_target"
  )

  if (is.null(tgl) || nrow(tgl) == 0L) return(character(0))

  missing <- setdiff(required, colnames(tgl))
  if (length(missing) > 0L) {
    stop(
      "[", md$name, "] target_gene_links missing expanded-target columns: ",
      paste(missing, collapse = ", ")
    )
  }

  keep <- as.character(tgl$source) == "loop_anchor" &
    is.finite(tgl$path_length) &
    tgl$path_length > 1L &
    tgl$in_expanded_target %in% TRUE

  if ("anchor_role" %in% colnames(tgl)) {
    keep <- keep & as.character(tgl$anchor_role) == "expanded_anchor"
  }

  if (identical(md$map, "promoter")) {
    keep <- keep & as.character(tgl$gene_role) == "promoter"
  }

  clean_gene_vector(tgl$gene[keep])
}

get_mode_genes <- function(md) {
  sort(unique(c(
    get_primary_mode_genes(md),
    get_expanded_mode_genes(md)
  )))
}

base_modes <- c(
  "anno_all_F", "anno_all_T",
  "refined_all_F", "refined_all_T",
  "chrom_all_F", "chrom_all_T",
  "chrom_only_all_F", "chrom_only_all_T",
  "anno_promoter_F", "anno_promoter_T",
  "refined_promoter_F", "refined_promoter_T",
  "chrom_promoter_F", "chrom_promoter_T",
  "chrom_only_promoter_F", "chrom_only_promoter_T"
)

mode_names <- if (isTRUE(cfg$include_hop1)) {
  c(base_modes, paste0(base_modes, "_hop1"))
} else {
  base_modes
}

modes <- vector("list", length(mode_names))
names(modes) <- mode_names
mode_genes <- vector("list", length(mode_names))
names(mode_genes) <- mode_names

for (nm in mode_names) {
  obj_name <- get_obj_name(nm)
  ann_obj <- get(obj_name)

  md <- list(
    name = nm,
    ann = ann_obj,
    map = if (grepl("_promoter_", nm)) "promoter" else "all",
    fill = grepl("_T($|_hop1$)", nm)
  )

  modes[[nm]] <- md
  mode_genes[[nm]] <- get_mode_genes(md)

  message(sprintf(
    "    %-36s -> %d genes",
    nm, length(mode_genes[[nm]])
  ))
}

# =============================================================================
# Expression-response residual model
# =============================================================================
message(">>> Fitting expression-adjusted BRD4 response model...")

model_df <- data.frame(
  gene = analysis_universe,
  signal_lfc = gs$signal_lfc[match(analysis_universe, gs$gene)],
  baseline_tpm = unname(expr_lookup[analysis_universe]),
  stringsAsFactors = FALSE
) %>%
  mutate(log_expr = log2(baseline_tpm + 0.1)) %>%
  filter(is.finite(signal_lfc), is.finite(log_expr))

if (nrow(model_df) < 100L) {
  stop("Too few genes for expression-adjustment model.")
}

spline_formula <- stats::as.formula(
  paste0(
    "signal_lfc ~ splines::ns(log_expr, df = ",
    as.integer(cfg$residual_spline_df),
    ")"
  )
)

residual_model <- tryCatch(
  stats::lm(spline_formula, data = model_df),
  error = function(e) {
    warning(
      "Spline model failed: ", conditionMessage(e),
      ". Falling back to linear signal_lfc ~ log_expr."
    )
    stats::lm(signal_lfc ~ log_expr, data = model_df)
  }
)

model_df$expected_signal <- stats::fitted(residual_model)
model_df$residual_signal <- stats::residuals(residual_model)

resid_lookup <- setNames(model_df$residual_signal, model_df$gene)

model_qc <- data.frame(
  N_Genes = nrow(model_df),
  Model = paste(deparse(stats::formula(residual_model)), collapse = ""),
  R2 = summary(residual_model)$r.squared,
  Adjusted_R2 = summary(residual_model)$adj.r.squared,
  Residual_SD = stats::sd(model_df$residual_signal, na.rm = TRUE),
  BaselineTPM_Median = stats::median(model_df$baseline_tpm, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(
  model_qc,
  file.path(out_dir, "Residual_Model_QC.csv"),
  row.names = FALSE
)

# Diagnostic plot: baseline expression vs response + fitted mean.
p_model <- ggplot(model_df, aes(x = log_expr, y = signal_lfc)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_line(
    data = model_df[order(model_df$log_expr), ],
    aes(y = expected_signal),
    linewidth = 0.8
  ) +
  labs(
    title = "Baseline expression vs BRD4 perturbation response",
    subtitle = "Spline fit defines the expression-adjusted residual response",
    x = "log2(DMSO TPM + 0.1)",
    y = "BRD4 response (-log2FoldChange)"
  ) +
  theme_classic(base_size = 10)

ggsave(
  file.path(out_dir, "Expression_Response_Residual_Model.pdf"),
  p_model, width = 6.5, height = 5
)

# =============================================================================
# Four-layer comparison
# =============================================================================

get_chip_reference <- function(md) {
  if (identical(cfg$chip_reference_policy, "promoter_specific") &&
      identical(md$map, "promoter")) {
    chip_prom
  } else {
    chip_genes
  }
}

make_comparison_sets <- function(md, nm) {
  loop_raw <- intersect(mode_genes[[nm]], analysis_universe)
  chip_ref <- intersect(get_chip_reference(md), analysis_universe)

  list(
    all_l_vs_bg = list(
      a = loop_raw,
      b = setdiff(analysis_universe, loop_raw),
      label = "All_looplook vs BG",
      scope = "Global", priority = "Primary"),
    l_vs_c = list(
      a = setdiff(loop_raw, chip_ref),
      b = setdiff(chip_ref, loop_raw),
      label = "Only_looplook vs Only_ChIPseeker",
      scope = "Unique", priority = "Supplementary"),
    l_vs_bg = list(
      a = setdiff(loop_raw, chip_ref),
      b = setdiff(analysis_universe, union(loop_raw, chip_ref)),
      label = "Only_looplook vs BG",
      scope = "Unique", priority = "Supplementary"),
    c_vs_bg = list(
      a = setdiff(chip_ref, loop_raw),
      b = setdiff(analysis_universe, union(loop_raw, chip_ref)),
      label = "Only_ChIPseeker vs BG",
      scope = "Unique", priority = "Supplementary")
  )
}

analyse_pair <- function(genes_a, genes_b, mode_name, comparison_key, seed_offset) {
  genes_a <- unique(intersect(genes_a, analysis_universe))
  genes_b <- unique(intersect(genes_b, analysis_universe))

  if (length(genes_a) < 10L || length(genes_b) < 10L) {
    return(NULL)
  }

  # ----- 1) Full ES -----
  xa <- gs$signal_lfc[match(genes_a, gs$gene)]
  xb <- gs$signal_lfc[match(genes_b, gs$gene)]

  full <- rank_biserial_boot(
    xa, xb,
    R = cfg$bootstrap_R,
    seed = cfg$bootstrap_seed + seed_offset
  )

  smd_before <- calc_smd(genes_a, genes_b, expr_lookup)

  # ----- 2) Common-support ES -----
  cs <- common_support_sets(genes_a, genes_b, expr_lookup)

  if (is.null(cs) || length(cs$a) < 10L || length(cs$b) < 10L) {
    common <- list(est = NA_real_, lo = NA_real_, hi = NA_real_)
    n_a_common <- 0L
    n_b_common <- 0L
    smd_common <- NA_real_
    overlap_min <- NA_real_
    overlap_max <- NA_real_
  } else {
    xcs <- gs$signal_lfc[match(cs$a, gs$gene)]
    ycs <- gs$signal_lfc[match(cs$b, gs$gene)]

    common <- rank_biserial_boot(
      xcs, ycs,
      R = cfg$bootstrap_R,
      seed = cfg$bootstrap_seed + 10000L + seed_offset
    )

    n_a_common <- length(cs$a)
    n_b_common <- length(cs$b)
    smd_common <- calc_smd(cs$a, cs$b, expr_lookup)
    overlap_min <- cs$overlap_min
    overlap_max <- cs$overlap_max
  }

  # ----- 3) Expression-adjusted residual ES -----
  ra <- resid_lookup[genes_a]
  rb <- resid_lookup[genes_b]
  ra <- ra[is.finite(ra)]
  rb <- rb[is.finite(rb)]

  adjusted <- if (length(ra) >= 10L && length(rb) >= 10L) {
    rank_biserial_boot(
      ra, rb,
      R = cfg$bootstrap_R,
      seed = cfg$bootstrap_seed + 20000L + seed_offset
    )
  } else {
    list(est = NA_real_, lo = NA_real_, hi = NA_real_)
  }

  # ----- 4) Strict expression-matched stress test -----
  # Strict matching is intentionally performed only inside the same
  # common-support population used above. If common support is unavailable,
  # the strict stress test is reported as unavailable rather than silently
  # falling back to the full sets.
  if (is.null(cs) || length(cs$a) < 10L || length(cs$b) < 10L) {
    iter_df <- data.frame(
      Iteration = seq_len(cfg$strict_n_iter),
      ES = NA_real_,
      SMD = NA_real_,
      Matched_N = NA_integer_,
      Eligible_Bins = NA_integer_
    )
  } else {
    strict_source_a <- cs$a
    strict_source_b <- cs$b

    iter_df <- lapply(seq_len(cfg$strict_n_iter), function(i) {
      mt <- strict_match_once(
        strict_source_a,
        strict_source_b,
        expr_lookup,
        seed = cfg$bootstrap_seed + 30000L + seed_offset * 1000L + i
      )

      if (is.null(mt)) {
        return(data.frame(
          Iteration = i,
          ES = NA_real_,
          SMD = NA_real_,
          Matched_N = NA_integer_,
          Eligible_Bins = NA_integer_
        ))
      }

      xm <- gs$signal_lfc[match(mt$a, gs$gene)]
      ym <- gs$signal_lfc[match(mt$b, gs$gene)]

      data.frame(
        Iteration = i,
        ES = rank_biserial_point(xm, ym),
        SMD = calc_smd(mt$a, mt$b, expr_lookup),
        Matched_N = mt$n,
        Eligible_Bins = mt$eligible_bins
      )
    })

    iter_df <- bind_rows(iter_df)
  }
  valid_iter <- iter_df %>%
    filter(is.finite(ES), is.finite(SMD), is.finite(Matched_N))

  if (nrow(valid_iter) >= cfg$min_valid_matched_iter) {
    strict_es <- stats::median(valid_iter$ES, na.rm = TRUE)
    strict_lo <- unname(stats::quantile(valid_iter$ES, 0.025, na.rm = TRUE))
    strict_hi <- unname(stats::quantile(valid_iter$ES, 0.975, na.rm = TRUE))

    strict_smd_med <- stats::median(valid_iter$SMD, na.rm = TRUE)
    strict_abs_smd_med <- stats::median(abs(valid_iter$SMD), na.rm = TRUE)
    strict_abs_smd_p95 <- unname(stats::quantile(abs(valid_iter$SMD), 0.95, na.rm = TRUE))
    strict_pass_rate <- mean(abs(valid_iter$SMD) < cfg$smd_target, na.rm = TRUE)
    matched_n_median <- stats::median(valid_iter$Matched_N, na.rm = TRUE)
    eligible_bins_median <- stats::median(valid_iter$Eligible_Bins, na.rm = TRUE)
  } else {
    strict_es <- strict_lo <- strict_hi <- NA_real_
    strict_smd_med <- strict_abs_smd_med <- strict_abs_smd_p95 <- NA_real_
    strict_pass_rate <- NA_real_
    matched_n_median <- eligible_bins_median <- NA_real_
  }

  # A positive signal_lfc means stronger downregulation after BRD4 perturbation.
  # Positive rank-biserial ES therefore means group A is more downregulated than group B.
  summary_row <- data.frame(
    Mode = mode_name,
    Comparison_Key = comparison_key,

    Full_ES = full$est,
    Full_ES_Lo = full$lo,
    Full_ES_Hi = full$hi,

    CommonSupport_ES = common$est,
    CommonSupport_ES_Lo = common$lo,
    CommonSupport_ES_Hi = common$hi,

    ExpressionAdjusted_ES = adjusted$est,
    ExpressionAdjusted_ES_Lo = adjusted$lo,
    ExpressionAdjusted_ES_Hi = adjusted$hi,

    StrictMatched_ES = strict_es,
    StrictMatched_ES_Lo = strict_lo,
    StrictMatched_ES_Hi = strict_hi,

    N_A_Original = length(genes_a),
    N_B_Original = length(genes_b),

    N_A_Common = n_a_common,
    N_B_Common = n_b_common,

    Matched_N_Median = matched_n_median,

    CommonRetention_A = safe_ratio(n_a_common, length(genes_a)),
    CommonRetention_B = safe_ratio(n_b_common, length(genes_b)),
    StrictRetention_A = safe_ratio(matched_n_median, length(genes_a)),
    StrictRetention_B = safe_ratio(matched_n_median, length(genes_b)),
    StrictRetention_FromCommon_A = safe_ratio(matched_n_median, n_a_common),
    StrictRetention_FromCommon_B = safe_ratio(matched_n_median, n_b_common),

    SMD_Before = smd_before,
    SMD_CommonSupport = smd_common,
    Strict_SMD_Median = strict_smd_med,
    Strict_AbsSMD_Median = strict_abs_smd_med,
    Strict_AbsSMD_P95 = strict_abs_smd_p95,
    Strict_P_absSMD_lt_0.10 = strict_pass_rate,
    Strict_Balance_Pass = ifelse(
      is.finite(strict_abs_smd_med),
      strict_abs_smd_med < cfg$smd_target,
      NA
    ),

    Eligible_Bins_Median = eligible_bins_median,
    CommonSupport_Min_LogExpr = overlap_min,
    CommonSupport_Max_LogExpr = overlap_max,

    Full_to_Common_AbsRetention = safe_ratio(abs(common$est), abs(full$est)),
    Full_to_Adjusted_AbsRetention = safe_ratio(abs(adjusted$est), abs(full$est)),
    Full_to_Strict_AbsRetention = safe_ratio(abs(strict_es), abs(full$est)),

    stringsAsFactors = FALSE
  )

  iter_df$Mode <- mode_name
  iter_df$Comparison_Key <- comparison_key

  list(summary = summary_row, strict_iterations = iter_df)
}

# =============================================================================
# Run all modes/comparisons
# =============================================================================
message(">>> Running integrated expression-sensitivity analysis...")

summary_list <- list()
strict_iter_list <- list()
mode_count_list <- list()

counter <- 0L

for (nm in mode_names) {
  md <- modes[[nm]]
  sets <- make_comparison_sets(md, nm)

  pipeline <- case_when(
    grepl("^anno_", nm) ~ "Basic",
    grepl("^refined_", nm) ~ "E-refined",
    grepl("^chrom_only_", nm) ~ "C-refined",
    grepl("^chrom_", nm) ~ "I-refined",
    TRUE ~ "Other"
  )

  for (cmp_key in names(sets)) {
    counter <- counter + 1L
    s <- sets[[cmp_key]]

    message(sprintf(
      "    %-36s | %-8s | A=%d B=%d",
      nm, cmp_key, length(s$a), length(s$b)
    ))

    ans <- analyse_pair(
      s$a, s$b,
      mode_name = nm,
      comparison_key = cmp_key,
      seed_offset = counter
    )

    if (is.null(ans)) next

    ans$summary$Comparison <- s$label
    ans$summary$Pipeline <- pipeline
    ans$summary$Hop <- ifelse(grepl("_hop1$", nm), "hop1", "hop0")
    ans$summary$Map <- ifelse(grepl("_promoter_", nm), "promoter", "all")
    ans$summary$Filled <- grepl("_T($|_hop1$)", nm)
    ans$summary$Scope <- s$scope
    ans$summary$Priority <- s$priority

    summary_list[[length(summary_list) + 1L]] <- ans$summary
    strict_iter_list[[length(strict_iter_list) + 1L]] <- ans$strict_iterations
  }

  mode_count_list[[nm]] <- data.frame(
    Mode = nm,
    Pipeline = pipeline,
    Hop = ifelse(grepl("_hop1$", nm), "hop1", "hop0"),
    Map = ifelse(grepl("_promoter_", nm), "promoter", "all"),
    Filled = grepl("_T($|_hop1$)", nm),
    Looplook_N = length(intersect(mode_genes[[nm]], analysis_universe)),
    ChIPseeker_N = length(intersect(get_chip_reference(md), analysis_universe)),
    stringsAsFactors = FALSE
  )
}

# ── ChIPseeker global (computed once) ──
message("    Computing ChIPseeker global...")
chip_global <- intersect(chip_genes, analysis_universe)
chip_bg_global <- setdiff(analysis_universe, chip_global)
ans_cg <- analyse_pair(chip_global, chip_bg_global, "ChIPseeker", "c_vs_bg", 99999)
if (!is.null(ans_cg)) {
  ans_cg$summary$Comparison <- "All_ChIPseeker vs BG"
  ans_cg$summary$Pipeline <- "ChIPseeker"
  ans_cg$summary$Hop <- "global"
  ans_cg$summary$Map <- "all"
  ans_cg$summary$Filled <- FALSE
  ans_cg$summary$Scope <- "Global"
  ans_cg$summary$Priority <- "Primary"
  summary_list[[length(summary_list) + 1L]] <- ans_cg$summary
  strict_iter_list[[length(strict_iter_list) + 1L]] <- ans_cg$strict_iterations
}

res_df <- bind_rows(summary_list)
strict_iter_df <- bind_rows(strict_iter_list)
mode_counts <- bind_rows(mode_count_list)

if (nrow(res_df) == 0L) {
  stop("No valid comparison results were produced.")
}

# =============================================================================
# Long-format result table
# =============================================================================
method_rows <- bind_rows(
  res_df %>% transmute(
    Mode, Pipeline, Hop, Map, Filled, Comparison,
    Method = "Full",
    ES = Full_ES, ES_Lo = Full_ES_Lo, ES_Hi = Full_ES_Hi
  ),
  res_df %>% transmute(
    Mode, Pipeline, Hop, Map, Filled, Comparison,
    Method = "Common support",
    ES = CommonSupport_ES,
    ES_Lo = CommonSupport_ES_Lo,
    ES_Hi = CommonSupport_ES_Hi
  ),
  res_df %>% transmute(
    Mode, Pipeline, Hop, Map, Filled, Comparison,
    Method = "Expression adjusted",
    ES = ExpressionAdjusted_ES,
    ES_Lo = ExpressionAdjusted_ES_Lo,
    ES_Hi = ExpressionAdjusted_ES_Hi
  ),
  res_df %>% transmute(
    Mode, Pipeline, Hop, Map, Filled, Comparison,
    Method = "Strict matched",
    ES = StrictMatched_ES,
    ES_Lo = StrictMatched_ES_Lo,
    ES_Hi = StrictMatched_ES_Hi
  )
)

method_rows$Method <- factor(
  method_rows$Method,
  levels = c(
    "Full",
    "Common support",
    "Expression adjusted",
    "Strict matched"
  )
)

# =============================================================================
# Output tables (Global + Unique split)
# =============================================================================
# Global
res_global <- res_df %>% filter(Scope == "Global")
write.csv(res_global, file.path(out_dir, "ExpressionSensitivity_v5_Global_ES.csv"), row.names = FALSE)

# Unique
res_unique <- res_df %>% filter(Scope == "Unique")
write.csv(res_unique, file.path(out_dir, "ExpressionSensitivity_v5_Unique_ES.csv"), row.names = FALSE)

# Full
write.csv(res_df, file.path(out_dir, "ExpressionSensitivity_v5_All_ES.csv"), row.names = FALSE)

write.csv(method_rows, file.path(out_dir, "ExpressionSensitivity_v5_ES_Long.csv"), row.names = FALSE)
write.csv(mode_counts, file.path(out_dir, "ExpressionSensitivity_v5_Mode_Gene_Counts.csv"), row.names = FALSE)
write.csv(strict_iter_df, file.path(out_dir, "ExpressionSensitivity_v5_StrictMatching_Iterations.csv"), row.names = FALSE)

# =============================================================================
# Figures
# =============================================================================
pipeline_colors <- c(
  "Basic" = "#1F77B4",
  "E-refined" = "#9467BD",
  "C-refined" = "#FFA500",
  "I-refined" = "#E64B35",
  "ChIPseeker" = "#94A3B8",
  "Other" = "grey50"
)

# 1) Global ES (Primary: All_looplook vs BG + All_ChIPseeker vs BG)
method_global <- method_rows %>% filter(Comparison %in% c("All_looplook vs BG", "All_ChIPseeker vs BG"))
p_es_global <- ggplot(method_global,
  aes(x = ES, y = reorder(Mode, ES, FUN = median, na.rm = TRUE), color = Pipeline, shape = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ES_Lo, xmax = ES_Hi), height = 0.20, linewidth = 0.35,
    alpha = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.55)) +
  facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
  scale_color_manual(values = pipeline_colors) +
  labs(title = "Global Target Validation: all assigned genes vs background",
    subtitle = "Full = full-set ES; Expression adjusted = spline-residual ES",
    x = "Rank-biserial effect size", y = NULL, color = "Pipeline", shape = "Analysis") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey35", hjust = 0.5))
ggsave(file.path(out_dir, "ExpressionSensitivity_v5_Global_ES.pdf"), p_es_global, width = 8, height = 6)

# 2) Unique ES (Supplementary)
method_unique <- method_rows %>% filter(!Comparison %in% c("All_looplook vs BG", "All_ChIPseeker vs BG"))
p_es_unique <- ggplot(method_unique,
  aes(x = ES, y = reorder(Mode, ES, FUN = median, na.rm = TRUE), color = Pipeline, shape = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ES_Lo, xmax = ES_Hi), height = 0.20, linewidth = 0.35,
    alpha = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.55)) +
  facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
  scale_color_manual(values = pipeline_colors) +
  labs(title = "Method-unique decomposition (Supplementary)",
    subtitle = "Full = full-set ES; Expression adjusted = spline-residual ES",
    x = "Rank-biserial effect size", y = NULL, color = "Pipeline", shape = "Analysis") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey35", hjust = 0.5))
ggsave(file.path(out_dir, "ExpressionSensitivity_v5_Unique_ES.pdf"), p_es_unique, width = 11, height = 12)

# 2) Retention after common-support and strict matching
retention_df <- bind_rows(
  res_df %>% transmute(
    Mode, Pipeline, Comparison,
    Stage = "Common support",
    Retention_A = CommonRetention_A,
    Retention_B = CommonRetention_B
  ),
  res_df %>% transmute(
    Mode, Pipeline, Comparison,
    Stage = "Strict matched",
    Retention_A = StrictRetention_A,
    Retention_B = StrictRetention_B
  )
) %>%
  tidyr::pivot_longer(
    c(Retention_A, Retention_B),
    names_to = "Group",
    values_to = "Retention"
  ) %>%
  mutate(
    Group = recode(
      Group,
      Retention_A = "Group A",
      Retention_B = "Group B"
    )
  )

p_ret <- ggplot(
  retention_df,
  aes(x = Retention, y = reorder(Mode, Retention, median, na.rm = TRUE),
      color = Stage, shape = Group)
) +
  geom_vline(
    xintercept = c(0.25, 0.50, 0.75),
    linetype = "dotted",
    color = "grey80",
    linewidth = 0.3
  ) +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 1.05),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Gene retention through expression conditioning",
    subtitle = "Low strict-match retention means the stress test represents only a subset of the original target set",
    x = "Fraction of original genes retained",
    y = NULL,
    color = "Stage",
    shape = "Group"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey35", hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "ExpressionSensitivity_Gene_Retention.pdf"),
  p_ret, width = 10, height = 11
)

# 3) SMD QC
smd_plot_df <- bind_rows(
  res_df %>% transmute(
    Mode, Pipeline, Comparison,
    Stage = "Before",
    SMD = SMD_Before
  ),
  res_df %>% transmute(
    Mode, Pipeline, Comparison,
    Stage = "Common support",
    SMD = SMD_CommonSupport
  ),
  res_df %>% transmute(
    Mode, Pipeline, Comparison,
    Stage = "Strict matched",
    SMD = Strict_SMD_Median
  )
)

p_smd <- ggplot(
  smd_plot_df,
  aes(
    x = abs(SMD),
    y = reorder(Mode, abs(SMD), median, na.rm = TRUE),
    color = Stage
  )
) +
  geom_vline(
    xintercept = cfg$smd_target,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
  labs(
    title = "Baseline-expression balance QC",
    subtitle = paste0(
      "Dashed line = |SMD| ",
      cfg$smd_target
    ),
    x = "|SMD| of log2(DMSO TPM + 0.1)",
    y = NULL,
    color = "Stage"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "ExpressionSensitivity_SMD_QC.pdf"),
  p_smd, width = 9, height = 11
)

# =============================================================================
# Excel workbook
# =============================================================================
wb <- createWorkbook()

addWorksheet(wb, "ES_4Layer")
writeData(wb, "ES_4Layer", as.data.frame(res_df))

addWorksheet(wb, "ES_Long")
writeData(wb, "ES_Long", as.data.frame(method_rows))

balance_cols <- c(
  "Mode", "Pipeline", "Hop", "Map", "Filled", "Comparison",
  "N_A_Original", "N_B_Original",
  "N_A_Common", "N_B_Common", "Matched_N_Median",
  "CommonRetention_A", "CommonRetention_B",
  "StrictRetention_A", "StrictRetention_B",
  "StrictRetention_FromCommon_A", "StrictRetention_FromCommon_B",
  "SMD_Before", "SMD_CommonSupport",
  "Strict_SMD_Median", "Strict_AbsSMD_Median",
  "Strict_AbsSMD_P95", "Strict_P_absSMD_lt_0.10",
  "Strict_Balance_Pass", "Eligible_Bins_Median"
)

addWorksheet(wb, "Balance_Retention_QC")
writeData(
  wb, "Balance_Retention_QC",
  as.data.frame(res_df[, intersect(balance_cols, colnames(res_df)), drop = FALSE])
)

addWorksheet(wb, "Residual_Model_QC")
writeData(wb, "Residual_Model_QC", model_qc)

addWorksheet(wb, "Mode_Gene_Counts")
writeData(wb, "Mode_Gene_Counts", mode_counts)

addWorksheet(wb, "Strict_Iterations")
writeData(wb, "Strict_Iterations", strict_iter_df)

saveWorkbook(
  wb,
  file.path(out_dir, "ExpressionSensitivity_Integrated_v5.xlsx"),
  overwrite = TRUE
)

# =============================================================================
# Provenance / session information
# =============================================================================
config_df <- data.frame(
  Parameter = names(cfg),
  Value = vapply(cfg, function(x) paste(x, collapse = ","), character(1)),
  stringsAsFactors = FALSE
)

write.csv(
  config_df,
  file.path(out_dir, "Analysis_Config.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("\n>>> DONE")
message("Output directory: ", out_dir)
message("Primary output:")
message("  - ExpressionSensitivity_Integrated_v5.xlsx")
message("  - ExpressionSensitivity_4Layer_ES.pdf")
message("  - ExpressionSensitivity_Gene_Retention.pdf")
message("  - ExpressionSensitivity_SMD_QC.pdf")
message("\nInterpretation priority:")
message("  Primary benchmark: matched-N raw-LFC benchmark from the main pipeline")
message("  Main expression sensitivity: ExpressionAdjusted_ES")
message("  Diagnostic: CommonSupport_ES")
message("  Stress test: StrictMatched_ES")
