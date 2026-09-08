#!/usr/bin/env Rscript
# I-refined GO top terms: Nature-style Nested Sub-Tile Heatmap (No Clustering)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(openxlsx)
  library(tidyr)
  library(scales)
  library(ggtext)
})

out_base <- looplook_env_path("LOOPLOOK_OUT_BASE")
go_file  <- file.path(out_base, "go_enrichment", "GO_Enrichment_300iter.xlsx")
stopifnot(file.exists(go_file))

go <- read.xlsx(go_file, sheet = 1)

# ----------------------------------------------------------------------
# 1. 数据筛选与 Panel 提取
# ----------------------------------------------------------------------
chrom_modes <- c("chrom_all_F","chrom_all_T","chrom_promoter_F","chrom_promoter_T")
gene_sets   <- c("only_looplook","intersection")

go_chrom <- go %>%
  filter(Mode %in% chrom_modes, GeneSet %in% gene_sets) %>%
  mutate(
    ModeShort = case_when(
      Mode == "chrom_all_F" ~ "all (F)",
      Mode == "chrom_all_T" ~ "all (T)",
      Mode == "chrom_promoter_F" ~ "promoter (F)",
      Mode == "chrom_promoter_T" ~ "promoter (T)"),
    Panel = paste0(ModeShort, "\n", GeneSet))

# 面板逻辑顺序（only_looplook 左半边，intersection 右半边；每半边先 F 再 T，先 promoter 再 all）
panel_order <- c(
  "promoter (F)\nonly_looplook", "promoter (T)\nonly_looplook",
  "all (F)\nonly_looplook", "all (T)\nonly_looplook",
  "promoter (F)\nintersection", "promoter (T)\nintersection",
  "all (F)\nintersection", "all (T)\nintersection")

# ----------------------------------------------------------------------
# 2. 提取 Top 10 通路 (按 pvalue)
# ----------------------------------------------------------------------
go_top_only <- go_chrom %>%
  filter(GeneSet == "only_looplook", !is.na(pvalue)) %>%
  group_by(Description) %>%
  summarise(min_p = min(pvalue), .groups = "drop") %>%
  slice_min(min_p, n = 10, with_ties = FALSE)

go_top_inter <- go_chrom %>%
  filter(GeneSet == "intersection", !is.na(pvalue)) %>%
  group_by(Description) %>%
  summarise(min_p = min(pvalue), .groups = "drop") %>%
  slice_min(min_p, n = 10, with_ties = FALSE)

go_top_terms <- c(go_top_only$Description, go_top_inter$Description)

# ----------------------------------------------------------------------
# 3. 补全网格数据 & 解析 GeneRatio 为纯数值 (支持 "15/300" 或 0.05 格式)
# ----------------------------------------------------------------------
grid_df <- expand.grid(
  Panel = factor(panel_order, levels = panel_order),
  Description = go_top_terms,
  stringsAsFactors = FALSE
)

parse_ratio <- function(x) {
  sapply(x, function(v) {
    if (is.na(v) || v == "") return(NA_real_)
    if (grepl("/", v)) {
      pts <- as.numeric(unlist(strsplit(v, "/")))
      return(pts[1] / pts[2])
    } else {
      return(as.numeric(v))
    }
  })
}

plot_data <- grid_df %>%
  left_join(go_chrom %>% distinct(Panel, Description, .keep_all = TRUE), by = c("Panel", "Description")) %>%
  mutate(
    p_safe = ifelse(!is.na(pvalue) & pvalue == 0, 1e-300, pvalue),
    neg_log10_p = ifelse(is.na(p_safe), 0, -log10(p_safe)),
    # 统一填色：p > 0.01 → 0（灰色段），p <= 0.01 → -log10(p)（红色渐变）
    fill_val = ifelse(!is.na(p_safe) & p_safe <= 0.01, neg_log10_p,
      ifelse(!is.na(p_safe), 0, NA)),
    gene_ratio = parse_ratio(GeneRatio),
    # 有 p 值的格子都显示方块（灰或红）
    gene_ratio = ifelse(is.na(p_safe), NA_real_, gene_ratio)
  )

# 计算 GeneRatio 的实际极值范围（用于 Size Scale 的自适应映射）
min_ratio <- min(plot_data$gene_ratio, na.rm = TRUE)
max_ratio <- max(plot_data$gene_ratio, na.rm = TRUE)

# ----------------------------------------------------------------------
# 4. 设置 GO 通路顺序（only_looplook 组与 intersection 组各自按 p 升序，互不影响）
# ----------------------------------------------------------------------
plot_data <- plot_data %>%
  mutate(Description_wrapped = str_wrap(Description, width = 50))
plot_data$Description_wrapped <- factor(plot_data$Description_wrapped,
  levels = str_wrap(rev(go_top_terms), width = 50))

# ----------------------------------------------------------------------
# 5. X 轴 Panel 文本着色与格式化
# ----------------------------------------------------------------------
panel_color_map <- c(
  "all (F)\nonly_looplook"        = "<span style='color:#E64B35; font-weight:bold;'>all (F)<br>only_looplook</span>",
  "all (F)\nintersection"         = "<span style='color:#1F77B4; font-weight:bold;'>all (F)<br>intersection</span>",
  "all (T)\nonly_looplook"        = "<span style='color:#E64B35; font-weight:bold;'>all (T)<br>only_looplook</span>",
  "all (T)\nintersection"         = "<span style='color:#1F77B4; font-weight:bold;'>all (T)<br>intersection</span>",
  "promoter (F)\nonly_looplook"   = "<span style='color:#E64B35; font-weight:bold;'>promoter (F)<br>only_looplook</span>",
  "promoter (F)\nintersection"    = "<span style='color:#1F77B4; font-weight:bold;'>promoter (F)<br>intersection</span>",
  "promoter (T)\nonly_looplook"   = "<span style='color:#E64B35; font-weight:bold;'>promoter (T)<br>only_looplook</span>",
  "promoter (T)\nintersection"    = "<span style='color:#1F77B4; font-weight:bold;'>promoter (T)<br>intersection</span>"
)

plot_data <- plot_data %>%
  mutate(Panel_colored = factor(panel_color_map[as.character(Panel)],
    levels = unname(panel_color_map[panel_order])))

# ----------------------------------------------------------------------
# 6. 核心绘图逻辑：外层网格 + 内层嵌套缩放正方形
# ----------------------------------------------------------------------
# 动态色带上限与刻度（避免超出 8 的尾部无标注；步长 4 防止文字重叠）
fill_max <- max(8, ceiling(max(plot_data$neg_log10_p, na.rm = TRUE)))
fill_breaks <- sort(unique(c(0, seq(4, fill_max, by = 4))))
fill_labels <- c(">0.01", paste0("10^-", seq(4, fill_max, by = 4)))

p_main <- ggplot(plot_data, aes(x = Panel_colored, y = Description_wrapped)) +
  # 图层 1：固定尺寸的背景底框（浅灰细线，保证网格可见）
  geom_tile(fill = "white", color = "#CBD5E1", linewidth = 0.5, width = 0.9, height = 0.9) +

  # 图层 2：内嵌动态正方形 (shape = 22)，fill 统一灰→红渐变
  geom_point(
    aes(fill = fill_val, size = gene_ratio),
    shape = 22,             # 22 代表带填充色和边框的正方形
    color = "#64748B",      # 内嵌小方块描深灰细边，防白底色块融合
    stroke = 0.35,          # 细边框厚度
    na.rm = TRUE
  ) +

# 颜色渐变：灰色段（p>0.01）→ 玫红渐变（p<=0.01）
scale_fill_gradientn(
    colors = c("#D1D5DB", "#FDA4AF", "#BE123C", "#881337"),
    values = c(0, 0.2, 0.55, 1),
    limits = c(0, fill_max),
    oob = scales::squish,
    breaks = fill_breaks,
    labels = fill_labels,
    na.value = "transparent"
  ) +

  # 动态尺寸映射 (GeneRatio)
  scale_size_continuous(
    limits = c(min_ratio, max_ratio),
    range = c(2.2, 7.5),     # 内嵌方块大小范围（7.5pt 刚好最大化嵌套在 0.9 网格内）
    name = "Gene Ratio",
    labels = scales::percent_format(accuracy = 0.1), # 图例显示精准百分比 (如 2.5%, 5.0%)
    breaks = scales::breaks_pretty(n = 4)
  ) +

  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +

  labs(fill = "p-value") +
  theme_classic(base_size = 9.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),

    axis.text.x.top = element_markdown(angle = 90, hjust = 0, vjust = 0.5, size = 8.5),
    axis.text.y.right = element_text(hjust = 0, size = 8.5, color = "#0F172A"),

    legend.position = "top",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8.5, face = "bold", vjust = 0.85),
    legend.text = element_text(size = 7.5),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.background = element_blank()
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      ticks = TRUE,
      frame.colour = "black",
      frame.linewidth = 0.4
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(fill = "white", color = "#9CA3AF")
    )
  )

# ----------------------------------------------------------------------
# 7. 导出 PDF
# ----------------------------------------------------------------------
out_dir <- file.path(out_base, "go_enrichment")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, "GO_Nested_Square_Heatmap.pdf")
ggsave(out_file, p_main, width = 11, height = 9, device = cairo_pdf)

message("Done! Rendered nested-tile GO heatmap saved to: ", out_file)
