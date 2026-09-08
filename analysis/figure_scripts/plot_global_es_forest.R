#!/usr/bin/env Rscript
# Global ES Forest (33-mode & 32-mode) — Option 1 with Explicit Axis Break (//)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(cowplot)
library(scales) # 提供 trans_new 支持

# 1. Configuration & Data Preparation ----------------------------------
out_base  <- looplook_env_path("LOOPLOOK_OUT_BASE")
size_dir  <- file.path(out_base, "size200")
cache_rds <- file.path(size_dir, ".es_33mode_cache.rds")

stopifnot(file.exists(cache_rds))

pipeline_colors <- c(
  "Basic"      = "#1F77B4",
  "E-refined"  = "#9467BD",
  "C-refined"  = "#FFA500",
  "I-refined"  = "#E64B35",
  "ChIPseeker" = "#94A3B8"
)

# Load global dataset
raw_es_df <- readRDS(cache_rds) %>%
  filter(Hop == "hop0")

sorted_looplook_modes <- raw_es_df %>%
  filter(Mode != "ChIPseeker") %>%
  arrange(desc(ES)) %>%
  pull(Mode)

mode_levels_33 <- c(sorted_looplook_modes, "ChIPseeker")

es_df <- raw_es_df %>%
  mutate(
    is_chip     = Mode == "ChIPseeker",
    is_promoter = grepl("_promoter_", Mode),
    is_all      = !is_promoter & !is_chip,
    is_F        = !grepl("_T($|_hop1)", Mode) & !is_chip,
    is_T        = !is_F & !is_chip,
    Pipeline = case_when(
      is_chip ~ "ChIPseeker",
      grepl("^anno_", Mode)       ~ "Basic",
      grepl("^refined_", Mode)    ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined",
      grepl("^chrom_", Mode)      ~ "I-refined",
      TRUE ~ "Other"
    ),
    Pipeline = factor(Pipeline, levels = c("Basic", "E-refined", "C-refined", "I-refined", "ChIPseeker")),
    ModeRow  = factor(Mode, levels = mode_levels_33)
  )

# 定义分段物理空间压缩变换：
# 0 ~ 0.01 占用最左侧 8% 物理宽度（保留 0 虚线）；
# 0.11 ~ 0.22 占用右侧 92% 物理宽度（让数据密集区彻底铺满）
piecewise_compress <- trans_new(
  name = "piecewise_compress",
  transform = function(x) {
    ifelse(x <= 0.01, x * 8, 0.08 + (x - 0.11) * (0.92 / 0.11))
  },
  inverse = function(x) {
    ifelse(x <= 0.08, x / 8, 0.11 + (x - 0.08) * (0.11 / 0.92))
  }
)

# 2. Plot Generator Function ---------------------------------------------
plot_global_forest <- function(df_subset, title_text, output_filename, is_33mode = TRUE) {
  
  # Prepare indicator matrix dataset
  ind_base <- df_subset %>%
    filter(!is_chip) %>%
    dplyr::select(ModeRow, is_all, is_promoter, is_F, is_T) %>%
    pivot_longer(-ModeRow, names_to = "Prop", values_to = "match") %>%
    mutate(Prop = factor(Prop,
                         levels = c("is_all", "is_promoter", "is_F", "is_T"),
                         labels = c("All", "Promoter", "Strict(F)", "Filled(T)")
    ))
  
  if (is_33mode) {
    ind_chip <- data.frame(
      ModeRow = factor("ChIPseeker", levels = levels(df_subset$ModeRow)),
      Prop    = factor(c("All", "Promoter", "Strict(F)", "Filled(T)"), levels = levels(ind_base$Prop)),
      match   = FALSE
    )
    ind <- bind_rows(ind_base, ind_chip)
  } else {
    ind <- ind_base
  }
  
  y_scale_aligned <- scale_y_discrete(limits = rev(levels(df_subset$ModeRow)), expand = expansion(add = c(0.6, 0.6)))
  
  # --- Left Panel: Indicator Matrix ---
  p_left <- ggplot(ind, aes(y = ModeRow)) +
    geom_point(data = filter(ind, !match), aes(x = Prop), color = "#CBD5E1", fill = "white", shape = 21, size = 1.8, stroke = 0.35) +
    geom_point(data = filter(ind, match),  aes(x = Prop), color = "#334155", fill = "#334155", shape = 21, size = 1.8) +
    geom_vline(xintercept = 2.5, color = "#CBD5E1", linewidth = 0.4) +
    
    annotate("text", x = 1.5, y = Inf, label = "Target", vjust = -0.8, fontface = "bold", size = 3.0, color = "#0F172A") +
    annotate("text", x = 3.5, y = Inf, label = "Fallback", vjust = -0.8, fontface = "bold", size = 3.0, color = "#0F172A") +
    
    scale_x_discrete(position = "bottom", expand = expansion(add = c(0.8, 0.8))) +
    y_scale_aligned +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10, base_family = "sans") +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8, color = "#1E293B"),
      axis.text.y      = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank(),
      strip.text       = element_blank(),
      strip.background = element_blank(),
      plot.margin      = margin(t = 16, r = -5, b = 35, l = 5)
    )
  
  # --- Right Panel: Forest Plot ---
  active_colors <- if (is_33mode) pipeline_colors else pipeline_colors[names(pipeline_colors) != "ChIPseeker"]
  
  # 计算 Y 轴底层边框在离散 scale 下的真实纵坐标位置 (y_min_border = 0.4)
  y_border <- 0.4
  y_slash_bottom <- y_border - 0.25
  y_slash_top    <- y_border + 0.25
  
  p_right <- ggplot(df_subset, aes(y = ModeRow)) +
    # 0 虚线保留在最左侧边缘
    geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi, color = Pipeline), height = 0.25, linewidth = 0.6) +
    geom_point(aes(x = ES, color = Pipeline), size = 2.5, shape = 18) +
    
    # ------------------ X 轴 // 截断符 ------------------
  # 1. 白色方块擦除 x = 0.03~0.07 处的底层边框
  annotate("rect", xmin = 0.03, xmax = 0.07, ymin = y_slash_bottom, ymax = y_slash_top, fill = "white", color = "white", clip = "off") +
    # 2. 第一条倾斜线段 /
    annotate("segment", x = 0.032, xend = 0.048, y = y_slash_bottom, yend = y_slash_top, color = "#0F172A", linewidth = 0.5, clip = "off") +
    # 3. 第二条倾斜线段 /
    annotate("segment", x = 0.052, xend = 0.068, y = y_slash_bottom, yend = y_slash_top, color = "#0F172A", linewidth = 0.5, clip = "off") +
    # ----------------------------------------------------
  
  # 强制显示完整图例
  scale_color_manual(values = active_colors, name = "Pipeline", drop = FALSE) +
    
    # 分段坐标映射：彻底挤掉中间空白，只在 0 和数据区打刻度
    scale_x_continuous(
      trans  = piecewise_compress,
      limits = c(0, 0.22),
      breaks = c(0, 0.12, 0.15, 0.18, 0.21),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    
    y_scale_aligned +
    coord_cartesian(clip = "off") + # 允许截断符微幅超出底边界而不被裁剪
    labs(x = "Effect Size (Rank-biserial Correlation)", y = NULL) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title.x     = element_text(size = 8.5, face = "bold", color = "#1E293B", margin = margin(t = 6)),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "#CBD5E1", linewidth = 0.5),
      strip.text       = element_blank(),
      strip.background = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 8.5, color = "#0F172A"),
      legend.text      = element_text(size = 8, color = "#1E293B"),
      legend.key       = element_blank(),
      legend.background= element_blank(),
      legend.box.margin= margin(t = 4, b = 0),
      plot.margin      = margin(t = 16, r = 5, b = 2, l = -5)
    )
  
  # --- Standard Alignment & Layout Composition ---
  aligned_plots <- align_plots(p_left, p_right, align = "v", axis = "tb")
  p_main <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], rel_widths = c(0.85, 2.0), nrow = 1)
  
  title_gg <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 10, color = "#0F172A")
  
  p_final <- plot_grid(
    title_gg,
    p_main,
    ncol = 1,
    rel_heights = c(0.04, 1)
  ) 
  
  out_path <- file.path(out_base, output_filename)
  ggsave(out_path, p_final, width = 3.5, height = 5, device = cairo_pdf)
  message("Saved: ", out_path)
}

# 3. Render Outputs ------------------------------------------------------

# Output 1: 33-Mode Global Forest Plot
plot_global_forest(
  df_subset = es_df,
  title_text = "Global Effect Size: 32 Looplook Modes + ChIPseeker",
  output_filename = "Global_ES_Forest_33Mode.pdf",
  is_33mode = TRUE
)

# Output 2: 32-Mode Global Forest Plot
es_32 <- es_df %>%
  filter(!is_chip) %>%
  mutate(
    Pipeline = factor(Pipeline, levels = c("Basic", "E-refined", "C-refined", "I-refined")),
    ModeRow  = factor(Mode, levels = sorted_looplook_modes)
  )

plot_global_forest(
  df_subset = es_32,
  title_text = "Global Effect Size: 32 Looplook Modes",
  output_filename = "Global_ES_Forest_32Mode.pdf",
  is_33mode = FALSE
)

message("Done.")
