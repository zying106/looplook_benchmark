#!/usr/bin/env Rscript
# Expression-Matched ES Forest v3 (balance-corrected)
# Common-support trimming + equal-count bin matching + balance QC + empirical CI
# NOTE: expression matching is a sensitivity analysis; primary benchmark remains matched-N.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
library(dplyr); library(ggplot2); library(stringr); library(openxlsx); library(boot)
library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)

# ═══ CONFIG ══════════════════════════════════════════════════════════

out_base   = looplook_env_path("LOOPLOOK_OUT_BASE")
rdata = looplook_env_path("LOOPLOOK_RDATA_FILE")
diff_csv  = looplook_env_path("LOOPLOOK_DIFF_FILE")
peak_file  = looplook_env_path("LOOPLOOK_PEAK_FILE")
expr_csv = looplook_env_path("LOOPLOOK_EXPR_FILE")

n_iter     <- 500L
n_bins     <- 20L
max_per_bin <- 50L
smd_target <- 0.10          # balance target for |SMD| after matching
min_success_iter <- 50L     # minimum valid matched iterations required
# Keep comparator consistent with the primary looplook benchmark.
# 'all' = ChIPseeker-all for every looplook mode (recommended/current benchmark definition).
# 'promoter_specific' = legacy behavior; promoter modes use ChIPseeker promoter genes.
chip_reference_policy <- "all"
# ══════════════════════════════════════════════════════════════════════

rank_biserial <- function(x, y, R = 2000L, seed = 42L) {
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

expr_mat <- as.matrix(read.table(expr_csv, header = TRUE, row.names = 1, check.names = FALSE))
expr_mat <- expr_mat[, grep("dmso|DMSO|Dmso", colnames(expr_mat)), drop = FALSE]
expr_lookup <- setNames(rowMeans(expr_mat, na.rm = TRUE), toupper(rownames(expr_mat)))
expr_lookup <- expr_lookup[is.finite(expr_lookup) & expr_lookup >= 0]
message("TPM: ", length(expr_lookup), " genes")

padf <- as.data.frame(annotatePeak(peak_file, tssRegion = c(-5000, 5000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db"))
padf$SYMBOL <- toupper(padf$SYMBOL)
chip_genes <- unique(na.omit(padf$SYMBOL))
chip_prom <- unique(na.omit(padf$SYMBOL[grepl("Promoter", padf$annotation)]))

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

mode_genes <- list()
for (nm in mode_names) {
  ao <- get_obj(nm); ann_obj <- if (!is.na(ao)) get(ao) else NULL
  md <- list(name = nm, ann = ann_obj,
    map = if (grepl("_promoter_", nm)) "promoter" else "all",
    fill = grepl("_T($|_hop1)", nm))
  if (!is.null(ann_obj)) {
    bi <- ann_obj$target_annotation; dc <- get_target_col(md)
    g <- if (!is.null(bi) && dc %in% colnames(bi))
      unique(toupper(looplook:::clean_gene_names(bi[[dc]], "[;,]"))) else character(0)
  } else g <- character(0)
  mode_genes[[nm]] <- g
}

# ═══ Core: expression-matched ES computation ═══
# Balance-corrected implementation:
#   1) trim A and B independently to their common expression support;
#   2) construct common pooled-expression bins;
#   3) sample EXACTLY the same number from A and B within every eligible bin;
#   4) report SMD distribution and balance-pass rate.
#
# The previous v2 implementation calculated bin_avail correctly but did not
# use it in matched_iter(), so A/B could contribute different counts per bin.
# It also used one recycled logical vector to trim A and B simultaneously.
do_expr_match <- function(genes_a, genes_b, label_a, label_b, ki) {
  valid_genes <- intersect(gs$gene[!is.na(gs$signal_lfc)], names(expr_lookup))
  genes_a <- unique(intersect(genes_a, valid_genes))
  genes_b <- unique(intersect(genes_b, valid_genes))

  if (length(genes_a) < 10L || length(genes_b) < 10L) return(NULL)

  calc_smd <- function(ga, gb) {
    ea <- log2(expr_lookup[ga] + 0.1)
    eb <- log2(expr_lookup[gb] + 0.1)
    ea <- ea[is.finite(ea)]
    eb <- eb[is.finite(eb)]
    if (length(ea) < 2L || length(eb) < 2L) return(NA_real_)
    psd <- sqrt((stats::var(ea) + stats::var(eb)) / 2)
    if (!is.finite(psd) || psd == 0) return(NA_real_)
    (mean(ea) - mean(eb)) / psd
  }

  # 1. Full ES before expression matching.
  xa <- gs$signal_lfc[gs$gene %in% genes_a]
  xb <- gs$signal_lfc[gs$gene %in% genes_b]
  xa <- xa[is.finite(xa)]
  xb <- xb[is.finite(xb)]
  if (length(xa) < 10L || length(xb) < 10L) return(NULL)
  rb_full <- rank_biserial(xa, xb, R = 500L, seed = 42 + ki)

  # 2. Baseline expression imbalance before trimming/matching.
  smd_before <- calc_smd(genes_a, genes_b)

  # 3. Common-support trimming.
  ea_all <- log2(expr_lookup[genes_a] + 0.1)
  eb_all <- log2(expr_lookup[genes_b] + 0.1)

  q_a <- stats::quantile(ea_all, c(0.005, 0.995), na.rm = TRUE, names = FALSE)
  q_b <- stats::quantile(eb_all, c(0.005, 0.995), na.rm = TRUE, names = FALSE)
  overlap_min <- max(q_a[1], q_b[1])
  overlap_max <- min(q_a[2], q_b[2])

  # No genuine common support = matching is not interpretable.
  if (!is.finite(overlap_min) || !is.finite(overlap_max) || overlap_min >= overlap_max) {
    warning("[", label_b, " | ", label_a, "] no usable common expression support; skipped.")
    return(NULL)
  }

  keep_a <- is.finite(ea_all) & ea_all >= overlap_min & ea_all <= overlap_max
  keep_b <- is.finite(eb_all) & eb_all >= overlap_min & eb_all <= overlap_max
  genes_a_trim <- genes_a[keep_a]
  genes_b_trim <- genes_b[keep_b]

  if (length(genes_a_trim) < 10L || length(genes_b_trim) < 10L) return(NULL)

  # 4. Build COMMON pooled-expression bins using shared numeric cut points.
  #    Unlike separate ranks, these bins represent the same expression ranges
  #    in both groups.
  match_df <- dplyr::bind_rows(
    data.frame(
      gene = genes_a_trim, group = "A",
      expr = log2(expr_lookup[genes_a_trim] + 0.1),
      stringsAsFactors = FALSE
    ),
    data.frame(
      gene = genes_b_trim, group = "B",
      expr = log2(expr_lookup[genes_b_trim] + 0.1),
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::filter(is.finite(expr))

  probs <- seq(0, 1, length.out = n_bins + 1L)
  bin_breaks <- unique(as.numeric(stats::quantile(
    match_df$expr, probs = probs, na.rm = TRUE, type = 8, names = FALSE
  )))

  if (length(bin_breaks) < 3L) {
    warning("[", label_b, " | ", label_a, "] too few distinct expression values for bin matching; skipped.")
    return(NULL)
  }

  # Ensure boundary values are included.
  bin_breaks[1] <- -Inf
  bin_breaks[length(bin_breaks)] <- Inf
  match_df$expr_bin <- cut(
    match_df$expr,
    breaks = bin_breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = FALSE
  )

  bin_ids <- sort(unique(match_df$expr_bin))
  bin_ids <- bin_ids[is.finite(bin_ids)]

  # Number sampled from EACH group in each bin.
  bin_avail <- vapply(bin_ids, function(b) {
    na <- sum(match_df$group == "A" & match_df$expr_bin == b)
    nb <- sum(match_df$group == "B" & match_df$expr_bin == b)
    min(na, nb, max_per_bin)
  }, integer(1))

  eligible_bins <- bin_ids[bin_avail > 0L]
  bin_avail <- bin_avail[bin_avail > 0L]

  if (length(eligible_bins) == 0L || sum(bin_avail) < 10L) {
    warning("[", label_b, " | ", label_a, "] insufficient overlap after exact bin-count balancing; skipped.")
    return(NULL)
  }

  a_bins <- lapply(eligible_bins, function(b) {
    match_df$gene[match_df$group == "A" & match_df$expr_bin == b]
  })
  b_bins <- lapply(eligible_bins, function(b) {
    match_df$gene[match_df$group == "B" & match_df$expr_bin == b]
  })

  # 5. Repeated equal-count matched sampling.
  matched_iter <- function(ri) {
    set.seed(42 + ki * 100L + ri)

    sa <- unlist(lapply(seq_along(eligible_bins), function(j) {
      sample(a_bins[[j]], size = bin_avail[[j]], replace = FALSE)
    }), use.names = FALSE)

    sb <- unlist(lapply(seq_along(eligible_bins), function(j) {
      sample(b_bins[[j]], size = bin_avail[[j]], replace = FALSE)
    }), use.names = FALSE)

    # This should always hold; fail loudly if a future edit breaks matching.
    if (length(sa) != length(sb)) {
      stop("Internal matching error: unequal matched N (A=", length(sa), ", B=", length(sb), ").")
    }
    if (length(sa) < 10L) return(c(es = NA_real_, smd = NA_real_, n = length(sa)))

    xm <- gs$signal_lfc[gs$gene %in% sa]
    ym <- gs$signal_lfc[gs$gene %in% sb]
    xm <- xm[is.finite(xm)]
    ym <- ym[is.finite(ym)]
    if (length(xm) < 10L || length(ym) < 10L) {
      return(c(es = NA_real_, smd = NA_real_, n = min(length(xm), length(ym))))
    }

    wt <- suppressWarnings(stats::wilcox.test(xm, ym, exact = FALSE))
    es_val <- 2 * as.numeric(wt$statistic) / (length(xm) * length(ym)) - 1
    smd_val <- calc_smd(sa, sb)

    c(es = es_val, smd = smd_val, n = length(sa))
  }

  mes <- t(vapply(seq_len(n_iter), matched_iter, numeric(3)))
  keep_iter <- is.finite(mes[, "es"]) & is.finite(mes[, "smd"])
  mes_es  <- mes[keep_iter, "es"]
  mes_smd <- mes[keep_iter, "smd"]
  mes_n   <- mes[keep_iter, "n"]

  if (length(mes_es) < min_success_iter) return(NULL)

  abs_smd <- abs(mes_smd)

  data.frame(
    Full_ES = rb_full$est,
    Full_ES_lo = rb_full$lo,
    Full_ES_hi = rb_full$hi,
    Matched_ES_Median = median(mes_es),
    Matched_ES_Lo = stats::quantile(mes_es, 0.025, names = FALSE),
    Matched_ES_Hi = stats::quantile(mes_es, 0.975, names = FALSE),
    SMD_Before = smd_before,
    SMD_After_Median = median(mes_smd),
    Abs_SMD_After_Median = median(abs_smd),
    Abs_SMD_After_P95 = stats::quantile(abs_smd, 0.95, names = FALSE),
    P_absSMD_lt_0.10 = mean(abs_smd < smd_target),
    Balance_Pass = median(abs_smd) < smd_target,
    N_A = length(genes_a),
    N_B = length(genes_b),
    N_A_Trim = length(genes_a_trim),
    N_B_Trim = length(genes_b_trim),
    Matched_N_Median = median(mes_n),
    Eligible_Bins = length(eligible_bins),
    Bin_N = sum(bin_avail),
    CommonSupport_Min = overlap_min,
    CommonSupport_Max = overlap_max,
    stringsAsFactors = FALSE
  )
}

# ═══ Compute ═══
results <- list()
message("Computing...")
pb <- txtProgressBar(max = length(mode_names), style = 3)
for (ki in seq_along(mode_names)) {
  nm <- mode_names[[ki]]
  md_map <- if (grepl("_promoter_", nm)) "promoter" else "all"
  loop_raw <- intersect(toupper(mode_genes[[nm]]), gs$gene)
  chip_ref_raw <- if (identical(chip_reference_policy, "promoter_specific") && md_map == "promoter") {
    chip_prom
  } else {
    chip_genes
  }
  chip_ref <- intersect(chip_ref_raw, gs$gene)
  only_loop_full <- intersect(setdiff(loop_raw, chip_ref), gs$gene)
  only_chip_full <- intersect(setdiff(chip_ref, loop_raw), gs$gene)
  bg_full  <- setdiff(setdiff(gs$gene, loop_raw), chip_ref)

  pipeline <- case_when(
    grepl("^anno_", nm) ~ "Basic", grepl("^refined_", nm) ~ "E-refined",
    grepl("^chrom_only_", nm) ~ "C-refined", grepl("^chrom_", nm) ~ "I-refined")

  message(sprintf("[%s] loop=%d chip=%d bg=%d", nm,
    length(only_loop_full), length(only_chip_full), length(bg_full)))
  comparisons <- list(
    list(label = "l_vs_c",  a = only_loop_full, b = only_chip_full, offset = 0),
    list(label = "l_vs_bg", a = only_loop_full, b = bg_full,        offset = 1000),
    list(label = "c_vs_bg", a = only_chip_full,  b = bg_full,        offset = 2000))

  for (cmp in comparisons) {
    cm <- do_expr_match(cmp$a, cmp$b, cmp$label, nm, ki + cmp$offset)
    if (!is.null(cm)) {
      cm$Comparison <- cmp$label
      cm$Mode <- nm
      cm$Pipeline <- pipeline
      results[[length(results) + 1]] <- cm
    }
  }
  setTxtProgressBar(pb, ki)
}
close(pb)

res_df <- bind_rows(results)
if (nrow(res_df) == 0) stop("No results.")

res_df$Pipeline <- factor(res_df$Pipeline, levels = c("Basic","E-refined","C-refined","I-refined"))
res_df <- res_df %>% mutate(
  Comparison = factor(Comparison,
    levels = c("l_vs_bg", "c_vs_bg", "l_vs_c"),
    labels = c("Only_looplook vs BG", "Only_ChIPseeker vs BG", "Only_looplook vs Only_ChIPseeker")))

pipeline_colors <- c("Basic"="#1F77B4","E-refined"="#9467BD","C-refined"="#FFA500","I-refined"="#E64B35")

p <- ggplot(res_df, aes(x = Full_ES, y = reorder(Mode, Full_ES), color = Pipeline)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = Full_ES_lo, xmax = Full_ES_hi), height = 0.3, linewidth = 0.5, alpha = 0.4) +
  geom_point(size = 2.5, shape = 18) +
  geom_point(aes(x = Matched_ES_Median), shape = 1, size = 2.2, stroke = 0.7) +
  geom_errorbarh(aes(xmin = Matched_ES_Lo, xmax = Matched_ES_Hi),
    height = 0.3, linewidth = 0.5, linetype = "dotted", alpha = 0.4) +
  facet_wrap(~ Comparison, ncol = 3, scales = "free_y") +
  scale_color_manual(values = pipeline_colors) +
  labs(title = "Expression-Matched Effect Size (Trimmed to Common Support)",
       subtitle = "filled diamond = full-set ES, open circle = expression-matched ES (median, 95% CI)",
       x = "Rank-Biserial ES", y = NULL, color = "Pipeline") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7, color = "grey40", hjust = 0.5),
        legend.position = "bottom", axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold", size = 8))

out_dir <- file.path(out_base, "expression_matched_es")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "ExpressionMatched_ES.pdf"), p, width = 14, height = 7.5)

wb <- createWorkbook()
addWorksheet(wb, "ES_Results")
writeData(wb, "ES_Results", as.data.frame(res_df))

balance_cols <- intersect(c(
  "Mode", "Pipeline", "Comparison",
  "SMD_Before", "SMD_After_Median", "Abs_SMD_After_Median",
  "Abs_SMD_After_P95", "P_absSMD_lt_0.10", "Balance_Pass",
  "N_A", "N_B", "N_A_Trim", "N_B_Trim",
  "Matched_N_Median", "Eligible_Bins", "Bin_N",
  "CommonSupport_Min", "CommonSupport_Max"
), colnames(res_df))
addWorksheet(wb, "Balance_QC")
writeData(wb, "Balance_QC", as.data.frame(res_df[, balance_cols, drop = FALSE]))

saveWorkbook(wb, file.path(out_dir, "ExpressionMatched_ES.xlsx"), overwrite = TRUE)
message("Done. Output: ", out_dir)
