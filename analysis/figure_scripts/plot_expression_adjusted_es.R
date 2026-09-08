#!/usr/bin/env Rscript
# Expression-adjusted ES figures only (no Common-support / Strict-matched)
# Reads v5 Global + Unique + Evidence data from CSV/xlsx
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(tidyr); library(openxlsx)
  library(cowplot)
})

out_base <- looplook_env_path("LOOPLOOK_OUT_BASE")
v5_dir   <- file.path(out_base, "expression_sensitivity_v5")
out_dir  <- v5_dir

pipeline_colors <- c("Basic"="#1F77B4","E-refined"="#9467BD","C-refined"="#FFA500",
                     "I-refined"="#E64B35","ChIPseeker"="#94A3B8")
ev_colors <- c(
  "local_promoter_overlap"="#2166AC","distal_promoter"="#4393C3",
  "gene_body_context"="#92C5DE","distal_gene_body_context"="#7B8FD4",
  "local_enhancer_candidate"="#F4A582","distal_enhancer_candidate"="#E64B35",
  "linear_fallback"="#CA0020","positional_candidate"="#FFA500",
  "direct_opposite_promoter"="#053061","local_promoter"="#74ADD1","linear_annotation"="#D6604D")

# ═══ 1. Global (v3-style: left indicators + right forest) ─────
res <- read.csv(file.path(v5_dir, "ExpressionSensitivity_v5_All_ES.csv"),
                stringsAsFactors = FALSE) %>%
  filter(Hop == "hop0" | Hop == "global") %>%
  mutate(Pipeline = factor(Pipeline, levels = names(pipeline_colors)))

make_global_figure <- function(df, chip_row = NULL, title_text, out_file) {
  # Sort by ES descending, ChIPseeker at bottom
  if (!is.null(chip_row)) {
    df <- bind_rows(df, chip_row)
    loop_order <- df %>% filter(Mode != "ChIPseeker") %>%
      arrange(desc(ExpressionAdjusted_ES)) %>% pull(Mode)
    df <- df %>% mutate(ModeRow = factor(Mode, levels = c(loop_order, "ChIPseeker")))
  } else {
    loop_order <- df %>% arrange(desc(ExpressionAdjusted_ES)) %>% pull(Mode)
    df <- df %>% mutate(ModeRow = factor(Mode, levels = loop_order))
  }

  df <- df %>% arrange(desc(ExpressionAdjusted_ES))

  # Indicator data
  ind <- df %>%
    mutate(is_chip = Mode == "ChIPseeker",
           is_promoter = grepl("_promoter_", Mode),
           is_all = !is_promoter & !is_chip,
           is_F = !grepl("_T($|_hop1)", Mode) & !is_chip,
           is_T = !is_F & !is_chip) %>%
    filter(!is_chip) %>%
    dplyr::select(ModeRow, is_all, is_promoter, is_F, is_T) %>%
    pivot_longer(-ModeRow, names_to = "Prop", values_to = "match") %>%
    mutate(Prop = factor(Prop, levels = c("is_all","is_promoter","is_F","is_T"),
      labels = c("All","Promoter","Strict(F)","Filled(T)")))

  if (!is.null(chip_row)) {
    ind_chip <- data.frame(ModeRow = factor("ChIPseeker", levels=levels(df$ModeRow)),
      Prop = factor(c("All","Promoter","Strict(F)","Filled(T)"), levels=levels(ind$Prop)),
      match = FALSE)
    ind <- bind_rows(ind, ind_chip)
  }

  # Left panel
  p_left <- ggplot(ind, aes(y = ModeRow)) +
    geom_point(data = filter(ind, !match), aes(x = Prop),
      color = "#CBD5E1", fill = "white", shape = 21, size = 1.8, stroke = 0.35) +
    geom_point(data = filter(ind, match), aes(x = Prop),
      color = "#334155", fill = "#334155", shape = 21, size = 1.8) +
    geom_vline(xintercept = 2.5, color = "#CBD5E1", linewidth = 0.4) +
    annotate("text", x = 1.5, y = Inf, label = "Target", vjust = -0.8,
      fontface = "bold", size = 3.0, color = "#0F172A") +
    annotate("text", x = 3.5, y = Inf, label = "Fallback", vjust = -0.8,
      fontface = "bold", size = 3.0, color = "#0F172A") +
    scale_x_discrete(position = "bottom", expand = expansion(add = c(0.8, 0.8))) +
    scale_y_discrete(limits = rev(levels(df$ModeRow)), expand = expansion(add = c(0.6, 0.6))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
      face = "bold", size = 8, color = "#1E293B"),
      axis.text.y = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), strip.text = element_blank(),
      strip.background = element_blank(),
      plot.margin = margin(t = 16, r = -5, b = 35, l = 5))

  # Right panel
  active_colors <- pipeline_colors[names(pipeline_colors) %in% unique(df$Pipeline)]
  p_right <- ggplot(df, aes(y = ModeRow)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ExpressionAdjusted_ES_Lo, xmax = ExpressionAdjusted_ES_Hi,
      color = Pipeline), height = 0.25, linewidth = 0.6) +
    geom_point(aes(x = ExpressionAdjusted_ES, color = Pipeline), size = 2.5, shape = 18) +
    scale_color_manual(values = active_colors, name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_discrete(limits = rev(levels(df$ModeRow)), expand = expansion(add = c(0.6, 0.6))) +
    labs(x = "Expression-adjusted ES (rank-biserial)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 8.5, face = "bold", color = "#1E293B", margin = margin(t = 6)),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#CBD5E1", linewidth = 0.5),
      strip.text = element_blank(), strip.background = element_blank(),
      legend.position = "bottom", legend.text = element_text(size = 8, color = "#1E293B"),
      legend.key = element_blank(), legend.background = element_blank(),
      legend.box.margin = margin(t = 4, b = 0),
      plot.margin = margin(t = 16, r = 5, b = 2, l = -5))

  aligned <- align_plots(p_left, p_right, align = "v", axis = "tb")
  p_main <- plot_grid(aligned[[1]], aligned[[2]], rel_widths = c(0.85, 2.0), nrow = 1)
  title_gg <- ggdraw() + draw_label(title_text, fontface = "bold", size = 10, color = "#0F172A")
  p_final <- plot_grid(title_gg, p_main, ncol = 1, rel_heights = c(0.04, 1))
  ggsave(file.path(out_dir, out_file), p_final, width = 5.2, height = 7, device = cairo_pdf)
  message("Saved: ", out_file)
}

# 33-mode: looplook global + ChIPseeker
df_loop <- res %>%
  filter(Comparison == "All_looplook vs BG", is.finite(ExpressionAdjusted_ES))
df_chip <- res %>%
  filter(Comparison == "All_ChIPseeker vs BG", is.finite(ExpressionAdjusted_ES))

make_global_figure(df_loop, df_chip,
  "Expression-adjusted ES: 32 Looplook Modes + ChIPseeker",
  "ExpressionAdjusted_Global_ES_33Mode.pdf")

# 32-mode: looplook only
make_global_figure(df_loop, NULL,
  "Expression-adjusted ES: 32 Looplook Modes",
  "ExpressionAdjusted_Global_ES_32Mode.pdf")

# ═══ 2. V5 Unique (Fig B: only_looplook/chipseeker vs BG) ─────
unique_df <- res %>%
  filter(Comparison %in% c("Only_looplook vs BG","Only_ChIPseeker vs BG","Only_looplook vs Only_ChIPseeker"),
         is.finite(ExpressionAdjusted_ES)) %>%
  mutate(Mode = factor(Mode, levels = rev(sort(unique(Mode)))))

p_unique <- ggplot(unique_df, aes(x = ExpressionAdjusted_ES, y = Mode, color = Pipeline)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ExpressionAdjusted_ES_Lo, xmax = ExpressionAdjusted_ES_Hi),
    height = 0.4, linewidth = 0.6) +
  geom_point(size = 2, shape = 18) +
  scale_color_manual(values = pipeline_colors) +
  facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
  labs(title = "Expression-adjusted ES: Method-unique Decomposition (Supplementary)",
       subtitle = "Only_looplook / Only_ChIPseeker vs background | expression-adjusted (spline residual)",
       x = "Rank-biserial ES (expression-adjusted)", y = NULL, color = "Pipeline") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5))

ggsave(file.path(out_dir, "ExpressionAdjusted_Unique_ES.pdf"), p_unique, width = 11, height = 12)

# ═══ 3. Evidence ES (1x4 Compact Forest) ─────
ev_file <- file.path(v5_dir, "ES_by_Evidence.xlsx")
if (file.exists(ev_file)) {
  ev_clean <- read.xlsx(ev_file) %>%
    filter(is.finite(ES_adj), grepl("_all_T($|_hop1)", Mode)) %>%
    mutate(
      Pipeline = factor(Pipeline, levels = names(pipeline_colors)[1:4]),
      Evidence = factor(Evidence, levels = rev(names(ev_colors))) %>% droplevels(),
      Facet_Label = paste0(Pipeline, "\n(", Mode, ")"))

  facet_levels <- ev_clean %>% arrange(Pipeline) %>% pull(Facet_Label) %>% unique()
  ev_clean$Facet_Label <- factor(ev_clean$Facet_Label, levels = facet_levels)

  plot_compact_forest <- function(df, es_var, lo_var, hi_var, title, xlab) {
    ggplot(df, aes(x = .data[[es_var]], y = Evidence, color = Evidence)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      geom_errorbarh(aes(xmin = .data[[lo_var]], xmax = .data[[hi_var]]),
        height = 0.35, linewidth = 0.5) +
      geom_point(size = 2.0, shape = 18) +
      scale_color_manual(values = ev_colors, drop = TRUE) +
      facet_wrap(~ Facet_Label, ncol = 4, scales = "free_x") +
      labs(title = title,
        subtitle = "Rows = Evidence | Points = rank-biserial ES (95% CI)",
        x = xlab, y = NULL) +
      theme_bw(base_size = 9.5) +
      theme(legend.position = "none",
        axis.text.y = element_text(size = 8, color = "black", face = "bold"),
        axis.text.x = element_text(size = 7.5, color = "black"),
        axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 6)),
        strip.text = element_text(face = "bold", size = 8.5, color = "white"),
        strip.background = element_rect(fill = "#2C3E50", color = NA),
        panel.grid = element_blank(),
        panel.spacing = unit(0.8, "lines"),
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle = element_text(size = 8.5, color = "grey35", hjust = 0, margin = margin(b = 8)))
  }

  p_ev <- plot_compact_forest(ev_clean, "ES_adj", "ES_adj_lo", "ES_adj_hi",
    "Expression-adjusted ES by Assignment Evidence (hop0)",
    "Rank-biserial ES (Adjusted)")
  ggsave(file.path(out_dir, "ExpressionAdjusted_Evidence_ES.pdf"), p_ev, width = 8, height = 3, device = cairo_pdf)

  p_ev_full <- plot_compact_forest(ev_clean, "ES", "ES_lo", "ES_hi",
    "Full ES by Assignment Evidence (hop0)",
    "Rank-biserial ES (Full)")
  ggsave(file.path(out_dir, "Full_Evidence_ES.pdf"), p_ev_full, width = 8, height = 3, device = cairo_pdf)
}

message("Done: ", out_dir)
