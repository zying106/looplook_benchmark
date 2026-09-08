#!/usr/bin/env Rscript
# Pipeline UpSet: anno / refined / chrom / chrom_only (all_T, hop0)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr)
})

# ═══ CONFIG ──────────────────────────────────────────────────
out_base   = looplook_env_path("LOOPLOOK_OUT_BASE")
rdata = looplook_env_path("LOOPLOOK_RDATA_FILE")
# ═══════════════════════════════════════════════════════════════

message("Loading...")
load(rdata, verbose = FALSE)

get_target_col <- function(md) {
  bc <- if (identical(md$map, "promoter")) "Regulated_promoter_genes" else "Assigned_Target_Genes"
  if (md$fill) paste0(bc, "_Filled") else bc
}

# all_T modes
mode_map <- list(
  "anno_all_T"      = "res",
  "refined_all_T"   = "refined_res",
  "chrom_all_T"     = "cr",
  "chrom_only_all_T" = "cr_only")

pipeline_names <- c("anno_all_T", "refined_all_T", "chrom_all_T", "chrom_only_all_T")
gene_sets <- lapply(pipeline_names, function(nm) {
  ann_obj <- get(mode_map[[nm]])
  bi <- ann_obj$target_annotation
  md <- list(map = "all", fill = TRUE)
  dc <- get_target_col(md)
  genes <- unique(toupper(looplook:::clean_gene_names(bi[[dc]], "[;,]")))
  genes[!is.na(genes) & genes != "" & genes != "NA"]
})
# 按位置映射：anno→Basic, refined→E-refined, chrom→I-refined, chrom_only→C-refined
names(gene_sets) <- c("Basic", "E-refined", "I-refined", "C-refined")

message("Gene set sizes:")
for (nm in names(gene_sets)) message("  ", nm, ": ", length(gene_sets[[nm]]))

# ----------------------------------------------------------------------
# UpSet 图（单页 PDF + 统一 Nature/Cell 规范配色）
# ----------------------------------------------------------------------
out_pdf <- file.path(out_base, "Pipeline_UpSet_AllT.pdf")

# 1. 明确集合顺序与 Pipeline 专属配色映射（与 master 脚本命名一致）
pipeline_cols <- c(
  "Basic"     = "#1F77B4",  # anno
  "E-refined" = "#9467BD",  # refined
  "C-refined" = "#FFA500",  # chrom_only
  "I-refined" = "#E64B35"   # chrom
)

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  message("UpSetR not installed. Falling back to ComplexHeatmap...")
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    m <- ComplexHeatmap::make_comb_mat(gene_sets)
    ht <- ComplexHeatmap::UpSet(
      m, 
      comb_order = order(ComplexHeatmap::comb_size(m), decreasing = TRUE)
    )
    
    cairo_pdf(out_pdf, width = 10, height = 6)
    ComplexHeatmap::draw(ht)
    dev.off()
    
    message("Saved: ", out_pdf)
  } else {
    stop("Neither UpSetR nor ComplexHeatmap available.")
  }
} else {
  set.seed(42)
  
  # 🌟 关键 1：必须在设备打开【之前】构造对象！
  # 若设备已打开再构造，UpSetR 1.4.x 会在 PDF 里残留一页空白页（实测 2 页：1 白 + 1 图）
  # 在 null device 上构造 → 无任何页面副作用
  p_upset <- UpSetR::upset(
    UpSetR::fromList(gene_sets),
    sets = c("Basic", "C-refined", "E-refined", "I-refined"), # UpSetR sets 顺序=从下到上；倒序后视觉上→下: I-refined, E-refined, C-refined, Basic（Basic 最底）
    keep.order = TRUE,
    order.by = "freq", 
    decreasing = TRUE,
    nsets = 4, 
    nintersects = 20,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size (Target Genes)",
    text.scale = c(1.3, 1.2, 1.0, 1.0, 1.2, 1.1),
    # 🌟 统一视觉语言配色
    main.bar.color = "#334155",   # 柱状图用高级深板岩灰
    matrix.color   = "#E64B35",   # 交集圆点连线用 I-refined 主色调点睛
    sets.bar.color = unname(pipeline_cols[c("Basic", "C-refined", "E-refined", "I-refined")]) # 与 sets 同序（从下到上）
  )
  # 🌟 关键 2：构造完成后再打开设备（优先 cairo_pdf，无 cairo 退回 base pdf()）
  if (capabilities("cairo")) {
    cairo_pdf(out_pdf, width = 6, height = 5)
  } else {
    pdf(out_pdf, width = 6, height = 5, useDingbats = FALSE)
  }
  print(p_upset)
  
  # 🌟 直接 close 设备，产出完美的单页 PDF
  dev.off()
  
  message("Saved: ", out_pdf)
}

# 逐级保留率表
sizes <- sapply(gene_sets, length)
retention <- data.frame(
  Pipeline = names(gene_sets),
  N_genes = sizes,
  stringsAsFactors = FALSE)
retention$Retained_from_Basic <- sapply(gene_sets, function(g) length(intersect(g, gene_sets[["Basic"]])))
retention$Retained_from_Refined <- sapply(gene_sets, function(g) length(intersect(g, gene_sets[["E-refined"]])))
write.csv(retention, file.path(out_base, "Pipeline_Retention_AllT.csv"), row.names = FALSE)
message("Retention table saved.")

# ── 交集组合 gene list 明细 ──
library(openxlsx)
comb_names <- names(gene_sets)
comb_list <- list()
# 全 15 种非空组合
for (k in 1:15) {
  # 用位运算枚举组合
  sets_in <- comb_names[as.logical(intToBits(k)[1:4])]
  if (length(sets_in) == 0) next
  sets_out <- setdiff(comb_names, sets_in)
  genes_combo <- Reduce(intersect, gene_sets[sets_in])
  if (length(sets_out) > 0) genes_combo <- setdiff(genes_combo, Reduce(union, gene_sets[sets_out]))
  comb_list[[length(comb_list)+1]] <- data.frame(
    Combination = paste(sets_in, collapse = "+"),
    Excluded = if (length(sets_out) > 0) paste(sets_out, collapse = "+") else "",
    N = length(genes_combo),
    Genes = paste(genes_combo, collapse = "; "),
    stringsAsFactors = FALSE)
}
comb_df <- bind_rows(comb_list) %>% arrange(desc(N))

wb <- createWorkbook()
addWorksheet(wb, "Combinations"); writeData(wb, "Combinations", comb_df)
addWorksheet(wb, "Gene_Sets"); writeData(wb, "Gene_Sets", data.frame(
  Pipeline = names(gene_sets),
  N = sizes,
  Genes = sapply(gene_sets, paste, collapse = "; "),
  stringsAsFactors = FALSE))
saveWorkbook(wb, file.path(out_base, "Pipeline_UpSet_AllT_Genes.xlsx"), overwrite = TRUE)
message("UpSet gene tables saved.")
message("Done.")
