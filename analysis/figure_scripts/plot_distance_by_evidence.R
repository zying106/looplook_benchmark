#!/usr/bin/env Rscript
# Distance by Evidence: Auto-sorted per panel + Exact Color Mapping
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(scales); library(openxlsx)
})

out_base <- looplook_env_path("LOOPLOOK_OUT_BASE")
xlsx_path <- file.path(out_base, "distance_analysis", "Peak_Gene_Distance_Analysis.xlsx")
stopifnot(file.exists(xlsx_path))
out_dir <- file.path(out_base, "distance_analysis")

# ── 1. Load & Clean ──
ev_raw <- read.xlsx(xlsx_path, sheet = "Distance_by_Evidence")

ev <- ev_raw %>%
  filter(!Evidence %in% c("none", "other"),
         grepl("_all_T$", Mode),   # only all_T hop0, consistent with ES analysis
         is.finite(Distance), Distance >= 0) %>%
  mutate(
    # 避免 log10(0) 产生的绘图异常/剔除
    Dist_kb_plot = pmax(Distance, 1) / 1000,

    # 统一命名规范
    Evidence_Clean = case_when(
      Evidence == "local_promoter_overlap"     ~ "local promoter overlap",
      Evidence == "distal_promoter"            ~ "distal promoter",
      Evidence == "gene_body_context"          ~ "gene body context",
      Evidence == "distal_gene_body_context"   ~ "distal gene body",
      Evidence == "local_enhancer_candidate"   ~ "local enhancer",
      Evidence == "distal_enhancer_candidate"  ~ "distal enhancer",
      Evidence == "linear_fallback"            ~ "linear fallback",
      Evidence == "positional_candidate"       ~ "positional",
      TRUE ~ Evidence),

    Pipeline = case_when(
      grepl("^anno_", Mode) ~ "Basic",
      grepl("^refined_", Mode) ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined",
      grepl("^chrom_", Mode) ~ "I-refined",
      TRUE ~ "Other"),
    Pipeline = factor(Pipeline, levels = c("Basic","E-refined","C-refined","I-refined"))
  )

# ── 2. 按每个 Panel 内部的中位数距离排序 ──
ev_order <- ev %>%
  group_by(Pipeline, Evidence_Clean) %>%
  summarise(Med = median(Dist_kb_plot, na.rm = TRUE), .groups = "drop") %>%
  group_by(Pipeline) %>%
  arrange(Med) %>%
  mutate(Evidence_Sorted = factor(Evidence_Clean, levels = Evidence_Clean)) %>%
  ungroup() %>%
  dplyr::select(Pipeline, Evidence_Clean, Evidence_Sorted)

ev <- ev %>% left_join(ev_order, by = c("Pipeline", "Evidence_Clean"))

# ── 3. Color Palette（覆盖所有可能的 clean 名称） ──
ev_colors <- c(
  "local promoter overlap" = "#2166AC",
  "distal promoter"        = "#4393C3",
  "gene body context"      = "#92C5DE",
  "distal gene body"       = "#7B8FD4",
  "local enhancer"         = "#F4A582",
  "distal enhancer"        = "#E64B35",
  "linear fallback"        = "#CA0020",
  "positional"             = "#FFA500"
)

# ── 4. Plot ──
p <- ggplot(ev, aes(x = Evidence_Sorted, y = Dist_kb_plot, fill = Evidence_Clean)) +
  geom_boxplot(
    outlier.shape = 16, outlier.size = 0.2, outlier.alpha = 0.1,
    coef = 1.5, alpha = 0.75, width = 0.55, color = "grey30"
  ) +

  scale_y_log10(
    breaks = c(0.1, 1, 10, 100, 1000, 10000),
    labels = c("0.1", "1", "10", "100", "1,000", "10,000"),
    expand = expansion(mult = c(0.05, 0.08))
  ) +

  scale_x_discrete(expand = expansion(add = c(1.5, 1.0))) +

  annotation_logticks(
    sides = "l", scaled = TRUE,
    short = unit(0.07, "cm"), mid = unit(0.12, "cm"), long = unit(0.18, "cm"),
    color = "black", linewidth = 0.35
  ) +
  coord_cartesian(clip = "on") +

  # 填充颜色，若未定义则显式设为 grey70 防错
  scale_fill_manual(values = ev_colors, na.value = "grey70") +

  facet_wrap(~ Pipeline, ncol = 4, scales = "free_x") +

  labs(
    title = "Peak\u2013Gene Distance by Evidence Type (hop0)",
    subtitle = "Sorted by median distance per pipeline | Unmapped evidence omitted",
    x = NULL, y = "Linear distance to TSS (kb, log scale)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 7.5, color = "black"),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(10, 12, 25, 25),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5)
  )

ggsave(file.path(out_dir, "Distance_by_Evidence_custom.pdf"), p, width = 7, height = 3.5)
message("Done: ", out_dir)
