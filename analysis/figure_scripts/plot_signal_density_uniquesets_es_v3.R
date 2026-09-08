#!/usr/bin/env Rscript
# Signal Density UniqueSets ES Forest — Native Legend Integration
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(purrr)
library(cowplot)
library(tidytext)

# 1. 数据加载与预处理 ----------------------------------------------------
out_base <- looplook_env_path("LOOPLOOK_OUT_BASE")
cache_file <- file.path(out_base, ".signal_density_cache.rds")
stopifnot(file.exists(cache_file))

sd <- readRDS(cache_file)

master_comp <- map_dfr(names(sd), ~ sd[[.x]]$comparisons %>% mutate(Mode = .x)) %>%
  mutate(
    Hop = if_else(str_detect(Mode, "_hop1$"), "hop1", "hop0"),
    Pipeline = factor(case_when(
      str_detect(Mode, "^anno_")       ~ "Basic",
      str_detect(Mode, "^refined_")    ~ "E-refined",
      str_detect(Mode, "^chrom_only_") ~ "C-refined",
      str_detect(Mode, "^chrom_")      ~ "I-refined",
      TRUE ~ "Other"
    ), levels = c("Basic", "E-refined", "C-refined", "I-refined"))
  )

# 2. 数据筛选与分面内部排序 ----------------------------------------------
df <- master_comp %>%
  filter(Hop == "hop0") %>%
  filter((Set_A == "Only_looplook" & Set_B == "Background") |
           (Set_A == "Only_ChIPseeker" & Set_B == "Background") |
           (Set_A == "Only_looplook" & Set_B == "Only_ChIPseeker")) %>%
  mutate(
    Comparison = factor(case_when(
      Set_A == "Only_looplook" & Set_B == "Background"        ~ "looplook-unique vs. Background",
      Set_A == "Only_ChIPseeker" & Set_B == "Background"      ~ "ChIPseeker-unique vs. Background",
      Set_A == "Only_looplook" & Set_B == "Only_ChIPseeker"   ~ "looplook-unique vs. ChIPseeker-unique"
    ), levels = c("looplook-unique vs. Background",
                  "ChIPseeker-unique vs. Background",
                  "looplook-unique vs. ChIPseeker-unique")),
    is_promoter = str_detect(Mode, "_promoter_"),
    is_all      = !is_promoter,
    is_F        = !str_detect(Mode, "_T($|_hop1)"),
    is_T        = !is_F,
    ModeRow     = reorder_within(Mode, ES, Comparison)
  )

pipeline_colors <- c(
  "Basic"     = "#1F77B4",
  "E-refined" = "#9467BD",
  "C-refined" = "#FFA500",
  "I-refined" = "#E64B35"
)

# 3. Header 设置 ---------------------------------------------------------
first_comp <- levels(df$Comparison)[1]
header_df <- data.frame(
  x = c(1.5, 3.5),
  label = c("Mode", "Fallback"),
  Comparison = factor(first_comp, levels = levels(df$Comparison))
)

ind <- df %>%
  distinct(ModeRow, Comparison, is_all, is_promoter, is_F, is_T) %>%
  pivot_longer(
    cols = c(is_all, is_promoter, is_F, is_T),
    names_to = "Prop",
    values_to = "match"
  ) %>%
  mutate(
    Prop = factor(Prop,
                  levels = c("is_all", "is_promoter", "is_F", "is_T"),
                  labels = c("All", "Promoter", "Strict (F)", "Filled (T)"))
  )

y_scale_aligned <- scale_y_reordered(expand = expansion(add = c(0.6, 0.6)))

# 4. 构建面板 ------------------------------------------------------------

# --- 左图：特征矩阵面板 ---
p_left <- ggplot(ind, aes(y = ModeRow)) +
  geom_point(data = filter(ind, !match), aes(x = Prop), color = "#CBD5E1", fill = "white", shape = 21, size = 1.8, stroke = 0.35) +
  geom_point(data = filter(ind, match),  aes(x = Prop), color = "#334155", fill = "#334155", shape = 21, size = 1.8) +
  geom_vline(xintercept = 2.5, color = "#CBD5E1", linewidth = 0.4) +
  geom_text(data = header_df, aes(x = x, y = Inf, label = label),
            vjust = -0.8, fontface = "bold", size = 3.0, color = "#0F172A") +
  facet_wrap(~Comparison, ncol = 1, scales = "free_y") +
  scale_x_discrete(position = "bottom", expand = expansion(add = c(1.0, 1.0))) +
  y_scale_aligned +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8, color = "#1E293B"),
    axis.text.y      = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(face = "bold", size = 8.5, color = "transparent"),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    panel.spacing    = unit(0.4, "lines"),
    plot.margin      = margin(t = 16, r = -5, b = 35, l = 5) # b 留出与右图图例同等高度的下边距
  )

# --- 右图：森林图面板（自带原生图例） ---
p_right <- ggplot(df, aes(y = ModeRow)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi, color = Pipeline), height = 0.25, linewidth = 0.6) +
  geom_point(aes(x = ES, color = Pipeline), size = 2.5, shape = 18) +
  scale_color_manual(values = pipeline_colors, name = "Pipeline") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.05))) +
  facet_wrap(~Comparison, ncol = 1, scales = "free_y") +
  y_scale_aligned +
  labs(x = "Effect Size (Rank-biserial Correlation)", y = NULL) +
  theme_bw(base_size = 10, base_family = "sans") +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.x     = element_text(size = 8.5, face = "bold", color = "#1E293B", margin = margin(t = 6)),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "#CBD5E1", linewidth = 0.5),
    strip.text       = element_text(face = "bold", size = 8.5, color = "#0F172A"),
    strip.background = element_rect(fill = "#F1F5F9", color = "#CBD5E1", linewidth = 0.5),
    panel.spacing    = unit(0.4, "lines"),
    # 图例直接内置在右图底部
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 8.5, color = "#0F172A"),
    legend.text      = element_text(size = 8, color = "#1E293B"),
    legend.key       = element_blank(),
    legend.background= element_blank(),
    legend.box.margin= margin(t = 4, b = 0),
    plot.margin      = margin(t = 16, r = 5, b = 2, l = -5)
  )

# 5. 精确对齐 Panel Box 顶级区 -------------------------------------------
# align = "v", axis = "tb" 将两幅图的顶边和底边Panel无缝对齐
aligned_plots <- align_plots(p_left, p_right, align = "v", axis = "tb")

p_main <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], rel_widths = c(0.85, 2.2), nrow = 1)

title_gg <- ggdraw() +
  draw_label("Signal Density in Unique Gene Sets (hop0)", fontface = "bold", size = 11, color = "#0F172A")

p_final <- plot_grid(
  title_gg,
  p_main,
  ncol = 1,
  rel_heights = c(0.04, 1)
)

# 6. 保存 PDF ------------------------------------------------------------
output_pdf <- file.path(out_base, "SignalDensity_UniqueSets_ES_hop0_NativeLegend.pdf")
ggsave(output_pdf, p_final, width = 6.0, height = 8.5, device = cairo_pdf)

message("Saved PDF with native legend to: ", output_pdf)
