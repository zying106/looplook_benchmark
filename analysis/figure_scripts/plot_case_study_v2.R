#!/usr/bin/env Rscript
# Case Study v2: Loop-supported distal-intergenic targets
#   - peak 位置和 peak name 记录
#   - 三层比较：looplook vs ChIPseeker（basic）、pipeline 间、优化层级（综合 vs 单一）
#   - 基因名 + 表达变化记录
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(openxlsx)
  library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
})

# ═══ CONFIG ──────────────────────────────────────────────────

out_base   = looplook_env_path("LOOPLOOK_OUT_BASE")
rdata = looplook_env_path("LOOPLOOK_RDATA_FILE")
diff_csv  = looplook_env_path("LOOPLOOK_DIFF_FILE")
peak_file  = looplook_env_path("LOOPLOOK_PEAK_FILE")
expr_csv = looplook_env_path("LOOPLOOK_EXPR_FILE")

dd_csv   = looplook_env_path("LOOPLOOK_DD_FILE")
# ═══════════════════════════════════════════════════════════════

# ── 工具函数 ──
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

message("Loading...")
load(rdata, verbose = FALSE)
diff_df <- read.csv(diff_csv, row.names = 1)
gs <- diff_df[!is.na(diff_df$log2FoldChange) & is.finite(diff_df$log2FoldChange), ]
gs$gene <- toupper(rownames(gs))
gs$signal_lfc <- -gs$log2FoldChange
# pvalue 列（若存在）；不存在则置 1
gs$pvalue <- if ("pvalue" %in% colnames(diff_df)) diff_df$pvalue[match(rownames(gs), rownames(diff_df))] else 1
gs$significant_down <- gs$log2FoldChange <= -1 & gs$pvalue < 0.05

# ChIPseeker 注释
padf <- as.data.frame(annotatePeak(peak_file, tssRegion = c(-5000, 5000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db"))
padf$SYMBOL <- toupper(padf$SYMBOL)
padf$peak_id <- paste0(padf$seqnames, ":", padf$start, "-", padf$end)
chip_genes <- unique(na.omit(padf$SYMBOL))

# ChIPseeker distal intergenic peaks
is_distal_anno <- grepl("Distal Intergenic", padf$annotation)
distal_coords <- padf$peak_id[is_distal_anno]
message("ChIPseeker distal intergenic peaks: ", length(distal_coords))

# TPM 基础表达（DMSO 平均）
expr_mat <- as.matrix(read.table(expr_csv, header = TRUE, row.names = 1, check.names = FALSE))
dmso_cols <- grep("dmso|DMSO|Dmso", colnames(expr_mat), ignore.case = TRUE, value = TRUE)
expr_lookup <- setNames(rowMeans(expr_mat[, dmso_cols, drop = FALSE], na.rm = TRUE), toupper(rownames(expr_mat)))
message("TPM lookup: ", sum(is.finite(expr_lookup) & expr_lookup >= 0), " genes")

# DD_vs_N 差异文件：标记在 DD 中显著上调的基因
dd_df <- read.csv(dd_csv, stringsAsFactors = FALSE)
dd_df$gene <- toupper(trimws(dd_df$gene))
dd_up_genes <- unique(dd_df$gene[dd_df$sig %in% TRUE & dd_df$dir == "up"])
message("DD_vs_N significant up-regulated genes: ", length(dd_up_genes))

# 8 个 mode（4 pipeline × all/promoter，F only，hop0）
cs_modes <- c(
  "anno_all_F", "anno_promoter_F",
  "refined_all_F", "refined_promoter_F",
  "chrom_all_F", "chrom_promoter_F",
  "chrom_only_all_F", "chrom_only_promoter_F")

# ── 1. 逐 mode 提取 distal peak → gene 明细 ──
peak_rows <- list()
gene_rows <- list()
summary_rows <- list()

for (nm in cs_modes) {
  ao <- get_obj(nm); ann_obj <- if (!is.na(ao)) get(ao) else NULL
  if (is.null(ann_obj)) next
  md <- list(name = nm, ann = ann_obj,
    map = if (grepl("_promoter_", nm)) "promoter" else "all",
    fill = grepl("_T($|_hop1)", nm))
  bi <- ann_obj$target_annotation
  if (is.null(bi) || nrow(bi) == 0) next
  dc <- get_target_col(md)
  if (!dc %in% colnames(bi)) next

  # 标记 distal peak
  loop_peak_ids <- paste0(bi$seqnames, ":", bi$start, "-", bi$end)
  is_distal <- loop_peak_ids %in% distal_coords
  if (sum(is_distal) == 0) next

  # distal peak → gene 明细
  for (i in which(is_distal)) {
    genes_i <- unique(toupper(looplook:::clean_gene_names(bi[[dc]][i], "[;,]")))
    genes_i <- genes_i[!is.na(genes_i) & genes_i != "" & genes_i != "NA"]
    if (length(genes_i) == 0) next
    for (g in genes_i) {
      if (!g %in% gs$gene) next
      # 只保留显著下调基因：log2FC <= -1 & p < 0.05
      if (!gs$significant_down[gs$gene == g]) next
      peak_rows[[length(peak_rows)+1]] <- data.frame(
        Mode = nm,
        Peak_ID = loop_peak_ids[i],
        Peak_Chr = as.character(bi$seqnames[i]),
        Peak_Start = bi$start[i], Peak_End = bi$end[i],
        Gene = g,
        Gene_TPM = expr_lookup[g],
        Gene_signal_lfc = gs$signal_lfc[gs$gene == g],
        Gene_log2FC = gs$log2FoldChange[gs$gene == g],
        Gene_pvalue = gs$pvalue[gs$gene == g],
        significant_down = TRUE,
        DD_up_regulated = g %in% dd_up_genes,
        stringsAsFactors = FALSE)
    }
  }

  # 汇总（只统计显著下调）
  pd <- bind_rows(peak_rows[grep(paste0("^", nm, "$"), sapply(peak_rows, function(x) x$Mode[1]))])
  distal_genes <- unique(pd$Gene)
  if (length(distal_genes) == 0) next

  gene_rows[[length(gene_rows)+1]] <- data.frame(
    Mode = nm, Gene = distal_genes,
    TPM = expr_lookup[distal_genes],
    signal_lfc = gs$signal_lfc[match(distal_genes, gs$gene)],
    log2FC = gs$log2FoldChange[match(distal_genes, gs$gene)],
    pvalue = gs$pvalue[match(distal_genes, gs$gene)],
    significant_down = gs$significant_down[match(distal_genes, gs$gene)],
    DD_up_regulated = distal_genes %in% dd_up_genes,
    stringsAsFactors = FALSE)

  summary_rows[[length(summary_rows)+1]] <- data.frame(
    Mode = nm,
    Pipeline = case_when(grepl("^anno_",nm)~"Basic", grepl("^refined_",nm)~"E-refined",
      grepl("^chrom_only_",nm)~"C-refined", grepl("^chrom_",nm)~"I-refined", TRUE~"Other"),
    Map = ifelse(grepl("_promoter_", nm), "promoter", "all"),
    Distal_Peaks = sum(is_distal),
    Assigned_Genes = length(distal_genes),
    N_Downregulated = sum(pd$Gene_signal_lfc > 0, na.rm = TRUE),
    Median_signal_lfc = median(pd$Gene_signal_lfc, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

peak_df <- bind_rows(peak_rows)
gene_df <- bind_rows(gene_rows)
summary_df <- bind_rows(summary_rows) %>% arrange(desc(Assigned_Genes))

out_dir <- file.path(out_base, "case_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 输出：peak 明细、基因明细、汇总
write.csv(peak_df, file.path(out_dir, "CaseStudy2_Peak_Gene_Detail.csv"), row.names = FALSE)
write.csv(gene_df, file.path(out_dir, "CaseStudy2_Genes.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "CaseStudy2_Summary.csv"), row.names = FALSE)

# ── 1.5 逐级对照表：每个 distal peak 在 4 个 pipeline 下的基因保留情况 ──
peak_df <- peak_df %>%
  mutate(Pipeline = case_when(
    grepl("^anno_", Mode) ~ "Basic",
    grepl("^refined_", Mode) ~ "E-refined",
    grepl("^chrom_only_", Mode) ~ "C-refined",
    grepl("^chrom_", Mode) ~ "I-refined",
    TRUE ~ "Other"))

# 只保留 all map 的 4 个 mode（递进比较基线一致）
hier_peak <- peak_df %>% filter(grepl("_all_F$", Mode))

if (nrow(hier_peak) > 0) {
  hier_levels <- c("Basic", "E-refined", "C-refined", "I-refined")

  # 每 peak × 每 pipeline 的显著下调基因集合（- 换 _ 避免列名解析错误）
  peak_hier <- hier_peak %>%
    mutate(Pipeline = str_replace_all(Pipeline, "-", "_")) %>%
    distinct(Peak_ID, Pipeline, Gene) %>%
    arrange(Peak_ID, match(Pipeline, str_replace_all(hier_levels, "-", "_")))

  # 转宽：每个 peak 一行，4 列各 pipeline 的基因列表
  hier_wide <- peak_hier %>%
    group_by(Peak_ID, Pipeline) %>%
    summarise(Genes = paste(sort(Gene), collapse = "; "), N = n(), .groups = "drop") %>%
    tidyr::pivot_wider(id_cols = Peak_ID, names_from = Pipeline,
      values_from = c(Genes, N), names_sep = "_")

  # 递进指标：每级新增/保留/滤掉的基因
  hier_compare <- hier_wide %>%
    mutate(
      Basic_genes   = ifelse(!is.na(N_Basic), N_Basic, 0),
      Eref_genes    = ifelse(!is.na(N_E_refined), N_E_refined, 0),
      Cref_genes    = ifelse(!is.na(N_C_refined), N_C_refined, 0),
      Iref_genes    = ifelse(!is.na(N_I_refined), N_I_refined, 0),
      Retained_B_to_E = pmin(Basic_genes, Eref_genes),
      Retained_E_to_C = pmin(Eref_genes, Cref_genes),
      Retained_C_to_I = pmin(Cref_genes, Iref_genes),
      Iref_only = Iref_genes)

  # 对照表：peak 坐标 + 各级基因数 + 各级基因列表
  hier_final <- hier_compare %>%
    dplyr::select(Peak_ID, Basic_genes, Eref_genes, Cref_genes, Iref_genes,
      Retained_B_to_E, Retained_E_to_C, Retained_C_to_I,
      Genes_Basic, Genes_E_refined, Genes_C_refined, Genes_I_refined)

  write.csv(hier_final, file.path(out_dir, "CaseStudy2_Pipeline_Progression.csv"), row.names = FALSE)
  message("Pipeline progression table saved: ", nrow(hier_final), " distal peaks")
}

# ── 2. 三层比较 ──
pipeline_colors <- c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")

# 图 A：pipeline 间 Assigned_Genes + 中位 signal_lfc 对比
p_bar <- ggplot(summary_df, aes(x = reorder(Mode, Assigned_Genes), y = Assigned_Genes, fill = Pipeline)) +
  geom_col(alpha = 0.85) + coord_flip() +
  scale_fill_manual(values = pipeline_colors) +
  labs(title = "Loop-supported distal-intergenic candidate genes",
       subtitle = "ChIPseeker-distal peaks assigned by looplook | 8 modes",
       x = NULL, y = "Assigned genes") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "CaseStudy2_Gene_Counts.pdf"), p_bar, width = 8, height = 5)

# 图 B：信号分布（含 ChIPseeker 最近基因对照）
chip_distal_genes <- unique(na.omit(padf$SYMBOL[is_distal_anno]))
chip_ref_df <- data.frame(Mode = "ChIPseeker", Gene = chip_distal_genes,
  TPM = expr_lookup[chip_distal_genes],
  signal_lfc = gs$signal_lfc[match(chip_distal_genes, gs$gene)],
  log2FC = gs$log2FoldChange[match(chip_distal_genes, gs$gene)],
  pvalue = gs$pvalue[match(chip_distal_genes, gs$gene)],
  significant_down = gs$significant_down[match(chip_distal_genes, gs$gene)],
  DD_up_regulated = chip_distal_genes %in% dd_up_genes,
  stringsAsFactors = FALSE) %>%
  filter(significant_down %in% TRUE)

gene_df$Pipeline <- case_when(grepl("^anno_",gene_df$Mode)~"Basic", grepl("^refined_",gene_df$Mode)~"E-refined",
  grepl("^chrom_only_",gene_df$Mode)~"C-refined", grepl("^chrom_",gene_df$Mode)~"I-refined", TRUE~"Other")
gene_df$ModeLabel <- paste0(gene_df$Mode, " (n=", ave(seq_along(gene_df$Mode), gene_df$Mode, FUN=length), ")")

p_dist <- ggplot(gene_df, aes(x = signal_lfc, y = reorder(ModeLabel, signal_lfc, median), fill = Pipeline)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, alpha = 0.7, width = 0.6) +
  scale_fill_manual(values = pipeline_colors) +
  labs(title = "Signal distribution of distal-intergenic loop targets",
       subtitle = "signal_lfc = -log2FoldChange (>0 = downregulated)",
       x = "signal_lfc", y = NULL) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "CaseStudy2_SignalDistribution.pdf"), p_dist, width = 10, height = 6)

# 图 C：优化层级比较（Basic vs E-refined vs I-refined vs C-refined，all map）
hier_df <- summary_df %>% filter(Map == "all") %>%
  mutate(Hierarchy = case_when(
    Pipeline == "Basic" ~ "1_Basic",
    Pipeline == "E-refined" ~ "2_Expr-refined",
    Pipeline == "C-refined" ~ "3_Chrom-only",
    Pipeline == "I-refined" ~ "4_Chrom-refined"))
p_hier <- ggplot(hier_df, aes(x = Hierarchy, y = Assigned_Genes, fill = Pipeline)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = pipeline_colors) +
  labs(title = "Optimization hierarchy (all map)",
       subtitle = "More filtering = fewer, higher-confidence distal targets",
       x = NULL, y = "Assigned genes") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "CaseStudy2_Optimization_Hierarchy.pdf"), p_hier, width = 7, height = 5)

# Excel 汇总
wb <- createWorkbook()
addWorksheet(wb, "Peak_Gene_Detail"); writeData(wb, "Peak_Gene_Detail", peak_df)
addWorksheet(wb, "Genes"); writeData(wb, "Genes", gene_df)
addWorksheet(wb, "Summary"); writeData(wb, "Summary", summary_df)
addWorksheet(wb, "ChIPseeker_Distal_Genes"); writeData(wb, "ChIPseeker_Distal_Genes", chip_ref_df)
saveWorkbook(wb, file.path(out_dir, "CaseStudy2_Detailed.xlsx"), overwrite = TRUE)

message("Done: ", out_dir)
