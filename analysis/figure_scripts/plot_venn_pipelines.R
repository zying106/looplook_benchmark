#!/usr/bin/env Rscript
# Custom Venn: 4 pipelines (all_T) vs ChIPseeker, 2x2 panel + full mode gene tables
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(ggvenn); library(patchwork)
  library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
  library(openxlsx)
})

# ═══ CONFIG ──────────────────────────────────────────────────
out_base   = looplook_env_path("LOOPLOOK_OUT_BASE")
rdata = looplook_env_path("LOOPLOOK_RDATA_FILE")
peak_file  = looplook_env_path("LOOPLOOK_PEAK_FILE")


# ═══════════════════════════════════════════════════════════════

message("Loading...")
load(rdata, verbose = FALSE)

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

# 全部 16 个 hop0 mode
all_modes <- c(
  "anno_all_F","anno_all_T","anno_promoter_F","anno_promoter_T",
  "refined_all_F","refined_all_T","refined_promoter_F","refined_promoter_T",
  "chrom_all_F","chrom_all_T","chrom_promoter_F","chrom_promoter_T",
  "chrom_only_all_F","chrom_only_all_T","chrom_only_promoter_F","chrom_only_promoter_T")

loop_gene_sets <- list()
for (nm in all_modes) {
  ao <- get_obj(nm); ann_obj <- if (!is.na(ao)) get(ao) else NULL
  if (is.null(ann_obj)) next
  md <- list(map = if (grepl("_promoter_", nm)) "promoter" else "all",
             fill = grepl("_T($|_hop1)", nm))
  bi <- ann_obj$target_annotation; dc <- get_target_col(md)
  genes <- if (!is.null(bi) && dc %in% colnames(bi))
    unique(toupper(looplook:::clean_gene_names(bi[[dc]], "[;,]"))) else character(0)
  loop_gene_sets[[nm]] <- genes[!is.na(genes) & genes != "" & genes != "NA"]
}

# ChIPseeker genes
padf <- as.data.frame(annotatePeak(peak_file, tssRegion = c(-5000, 5000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db"))
chip_genes <- unique(na.omit(toupper(padf$SYMBOL)))
message("ChIPseeker genes: ", length(chip_genes))

# ── 全 mode 三集合表格 ──
venn_rows <- list()
for (nm in all_modes) {
  g1 <- loop_gene_sets[[nm]]
  shared <- intersect(g1, chip_genes)
  only_loop <- setdiff(g1, chip_genes)
  only_chip <- setdiff(chip_genes, g1)
  max_n <- max(length(only_loop), length(shared), length(only_chip))
  venn_rows[[length(venn_rows)+1]] <- data.frame(
    Mode = nm,
    Pipeline = case_when(grepl("^anno_",nm)~"Basic", grepl("^refined_",nm)~"E-refined",
      grepl("^chrom_only_",nm)~"C-refined", grepl("^chrom_",nm)~"I-refined", TRUE~"Other"),
    Map = ifelse(grepl("_promoter_", nm), "promoter", "all"),
    Filled = grepl("_T($|_hop1)", nm),
    Only_Looplook_N = length(only_loop),
    Shared_N = length(shared),
    Only_ChIPseeker_N = length(only_chip),
    Only_Looplook_Genes = paste(only_loop, collapse = "; "),
    Shared_Genes = paste(shared, collapse = "; "),
    Only_ChIPseeker_Genes = paste(only_chip, collapse = "; "),
    stringsAsFactors = FALSE)
}
venn_df <- bind_rows(venn_rows)

out_dir <- file.path(out_base, "venn_custom")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 汇总表 + 每 mode 基因明细 sheet
write.csv(venn_df, file.path(out_dir, "Venn_AllModes_vs_ChIPseeker.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", venn_df %>% dplyr::select(Mode, Pipeline, Map, Filled,
  Only_Looplook_N, Shared_N, Only_ChIPseeker_N))
for (nm in all_modes) {
  row <- venn_df %>% filter(Mode == nm)
  gl <- loop_gene_sets[[nm]]
  shared <- intersect(gl, chip_genes)
  only_loop <- setdiff(gl, chip_genes)
  only_chip <- setdiff(chip_genes, gl)
  max_n <- max(length(only_loop), length(shared), length(only_chip))
  detail <- data.frame(
    Only_Looplook = c(only_loop, rep(NA, max_n - length(only_loop))),
    Shared = c(shared, rep(NA, max_n - length(shared))),
    Only_ChIPseeker = c(only_chip, rep(NA, max_n - length(only_chip))),
    stringsAsFactors = FALSE)
  addWorksheet(wb, substr(nm, 1, 31))
  writeData(wb, substr(nm, 1, 31), detail)
}
saveWorkbook(wb, file.path(out_dir, "Venn_AllModes_vs_ChIPseeker_Genes.xlsx"), overwrite = TRUE)
message("Venn tables saved: ", nrow(venn_df), " modes")

# ── 4 pipelines all_T 可视化 2x2 ──
pipeline_colors <- c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")
title_map <- c("anno_all_T"="Basic", "refined_all_T"="E-refined",
               "chrom_all_T"="I-refined", "chrom_only_all_T"="C-refined")

plots <- lapply(c("anno_all_T","refined_all_T","chrom_all_T","chrom_only_all_T"), function(nm) {
  g1 <- loop_gene_sets[[nm]]
  inter_n <- length(intersect(g1, chip_genes))
  total_union <- length(union(g1, chip_genes))
  pct <- round(100 * inter_n / total_union, 1)

  # unname() 去掉名字属性，保证 fill_color 是纯颜色值
  current_color <- unname(pipeline_colors[title_map[[nm]]])

  ggvenn(list("looplook" = g1, "ChIPseeker" = chip_genes),
    fill_color = c(current_color, "#94A3B8"),
    stroke_size = 0.5, stroke_color = "white",
    set_name_size = 4.5, text_size = 4, fill_alpha = 0.6,
    show_percentage = TRUE) +
    labs(title = sprintf("%s (%s)", title_map[[nm]], nm)) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
})

p_final <- wrap_plots(plots, ncol = 2, nrow = 2)
ggsave(file.path(out_dir, "Venn_4Pipelines_allT_vs_ChIPseeker.pdf"), p_final,
  width = 11, height = 9, device = cairo_pdf)
message("Done: ", out_dir)
