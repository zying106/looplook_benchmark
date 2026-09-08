#!/usr/bin/env Rscript
# Distance analysis v2 final: Density + Boxplot + Mode-level stats + Descriptive summary
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(scales)
  library(openxlsx)
})

out_base <- looplook_env_path("LOOPLOOK_OUT_BASE")
xlsx_path <- file.path(out_base, "distance_analysis", "Peak_Gene_Distance_Analysis.xlsx")
stopifnot(file.exists(xlsx_path))
out_dir <- file.path(out_base, "distance_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ═══ 1. LOAD & QC ────────────────────────────────────────────
message(">>> Loading...")
dist_raw <- read.xlsx(xlsx_path, sheet = "Peak_Gene_Pairs_Raw")

# QC: Distance must be absolute
if (any(dist_raw$Distance < 0, na.rm = TRUE))
  stop("Negative Distance detected. This analysis expects absolute distance.")

dist_df <- dist_raw %>%
  filter(is.finite(Distance)) %>%
  mutate(
    Dist_kb = Distance / 1000,
    Dist_kb_plot = pmax(Distance, 1) / 1000,
    GeneSet_Clean = case_when(
      GeneSet == "Only_ChIPseeker" ~ "Only ChIPseeker",
      GeneSet == "Intersection (ChIPseeker)" ~ "Overlap (ChIPseeker)",
      GeneSet == "Intersection (looplook)" ~ "Overlap (looplook)",
      GeneSet == "Only_looplook" ~ "Only looplook",
      TRUE ~ "OTHER_CHK"))

if (any(dist_df$GeneSet_Clean == "OTHER_CHK"))
  stop("Unexpected GeneSet: ", paste(unique(dist_df$GeneSet[dist_df$GeneSet_Clean == "OTHER_CHK"]), collapse=", "))

dist_df$GeneSet_Clean <- factor(dist_df$GeneSet_Clean,
  levels = c("Only ChIPseeker","Overlap (ChIPseeker)","Overlap (looplook)","Only looplook"))

hop0_modes <- unique(dist_df$Mode)[!grepl("_hop1$", unique(dist_df$Mode))]
dist_hop0 <- dist_df %>% filter(Mode %in% hop0_modes) %>%
  mutate(Pipeline = case_when(
    grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
    grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
    TRUE ~ "OTHER_CHK"))

if (any(dist_hop0$Pipeline == "OTHER_CHK")) stop("Unknown Pipeline assignment.")
dist_hop0$Pipeline <- factor(dist_hop0$Pipeline, levels = c("Basic","E-refined","C-refined","I-refined"))

# duplicate pair QC: check conflicts, then deduplicate
if ("Peak_ID" %in% colnames(dist_hop0)) {
  nr <- nrow(dist_hop0)

  pair_conflict <- dist_hop0 %>%
    group_by(Mode, GeneSet_Clean, Peak_ID, Gene) %>%
    summarise(N_rows = n(), N_distance = n_distinct(Distance), .groups = "drop") %>%
    filter(N_rows > 1, N_distance > 1)

  if (nrow(pair_conflict) > 0)
    stop("Conflicting Distance values for duplicate peak-gene pairs.")

  dist_hop0 <- dist_hop0 %>%
    distinct(Mode, GeneSet_Clean, Peak_ID, Gene, Distance, .keep_all = TRUE)

  message(sprintf("  duplicate QC: %d -> %d rows", nr, nrow(dist_hop0)))
}

dist_colors <- c("Only ChIPseeker"="#4DBBD5","Overlap (ChIPseeker)"="#8EA600",
                 "Overlap (looplook)"="#00A087","Only looplook"="#E64B35")

# ═══ 2. DESCRIPTIVE SUMMARY ──────────────────────────────────
message(">>> Descriptive...")
desc <- dist_hop0 %>% group_by(Pipeline, Mode, GeneSet_Clean) %>%
  summarise(N=n(), Median_kb=median(Dist_kb), Q1_kb=quantile(Dist_kb,.25),
    Q3_kb=quantile(Dist_kb,.75), P_gt10kb=mean(Dist_kb>10),
    P_gt100kb=mean(Dist_kb>100), P_gt1Mb=mean(Dist_kb>1000), .groups="drop")
write.csv(desc, file.path(out_dir,"Distance_Descriptive_Summary.csv"), row.names=FALSE)

# ═══ 3. MODE-LEVEL STATISTICS (per Pipeline) ──────────────
# Statistical-unit correction. Raw peak-gene pairs are NOT independent:
# one gene can link to many peaks, and looplook pairs are many-to-one while
# ChIPseeker pairs are one-to-one. Inference on raw pairs therefore inflates
# the effective sample size (the previous LMM on log-distance pairs). Each
# benchmark mode instead contributes ONE median distance per gene-set
# category, and the two categories are compared with a paired Wilcoxon test
# across modes within each pipeline.
message(">>> Mode-level statistics (per Pipeline)...")
pipelines_lv <- c("Basic","E-refined","C-refined","I-refined")

# Per-mode median distance per gene-set category
mode_med <- dist_hop0 %>%
  group_by(Pipeline, Mode, GeneSet_Clean) %>%
  summarise(Med_kb = median(Dist_kb, na.rm = TRUE),
            N_pairs = n(), .groups = "drop")

# Wide table: one row per mode, keyed on the two compared categories
med_loop <- mode_med %>%
  filter(GeneSet_Clean == "Only looplook") %>%
  dplyr::select(Pipeline, Mode, Med_loop_kb = Med_kb)
med_chip <- mode_med %>%
  filter(GeneSet_Clean == "Only ChIPseeker") %>%
  dplyr::select(Pipeline, Mode, Med_chip_kb = Med_kb)
mode_wide <- inner_join(med_loop, med_chip, by = c("Pipeline", "Mode"))

mode_paired_stats <- function(df, pipe) {
  n_pair <- nrow(df)
  if (n_pair < 3L) {
    return(data.frame(Pipeline = pipe, N_PairedModes = n_pair,
      Delta_Median_kb = NA_real_, Ratio = NA_real_,
      Pct_looplook_further = NA_real_, P = NA_real_,
      stringsAsFactors = FALSE))
  }
  delta <- df$Med_loop_kb - df$Med_chip_kb
  wt <- suppressWarnings(wilcox.test(df$Med_loop_kb, df$Med_chip_kb,
    paired = TRUE, exact = FALSE))
  data.frame(Pipeline = pipe, N_PairedModes = n_pair,
    Delta_Median_kb = median(delta),
    Ratio = median(df$Med_loop_kb / pmax(df$Med_chip_kb, 0.1)),
    Pct_looplook_further = 100 * mean(delta > 10),
    P = wt$p.value, stringsAsFactors = FALSE)
}

mode_level_stats <- bind_rows(lapply(pipelines_lv, function(p)
  mode_paired_stats(mode_wide %>% filter(Pipeline == p), p)))
mode_level_stats <- bind_rows(mode_level_stats,
  mode_paired_stats(mode_wide %>% mutate(Pipeline = "Overall"), "Overall"))
mode_level_stats$P_BH <- p.adjust(mode_level_stats$P, method = "BH")
write.csv(mode_level_stats,
  file.path(out_dir, "Distance_ModeLevel_Stats.csv"), row.names = FALSE)

# Mode-level QC: how many modes per pipeline could be paired
mode_qc <- data.frame(
  Pipeline = pipelines_lv,
  N_modes = vapply(pipelines_lv, function(p)
    length(unique(dist_hop0$Mode[dist_hop0$Pipeline == p])), integer(1)),
  N_paired = vapply(pipelines_lv, function(p)
    sum(mode_wide$Pipeline == p), integer(1)),
  stringsAsFactors = FALSE)
write.csv(mode_qc, file.path(out_dir, "Distance_ModeLevel_QC.csv"), row.names = FALSE)

# ═══ 4. DENSITY ──────────────────────────────────────────────
message(">>> Density...")
p_d <- ggplot(dist_hop0, aes(x=Dist_kb_plot, color=GeneSet_Clean, fill=GeneSet_Clean)) +
  geom_density(alpha=.15, linewidth=.85) +
  geom_vline(xintercept=c(10,100), linetype="dashed", color="grey50", linewidth=.35) +
  annotate("text", x=10, y=Inf, label="10 kb", vjust=1.5, hjust=1.1, size=3, color="grey40") +
  annotate("text", x=100, y=Inf, label="100 kb", vjust=1.5, hjust=-.1, size=3, color="grey40") +
  scale_x_log10(breaks=c(.1,1,10,100,1000,10000), labels=c("0.1","1","10","100","1,000","10,000")) +
  scale_color_manual(values=dist_colors) + scale_fill_manual(values=dist_colors) +
  guides(color=guide_legend(override.aes=list(alpha=.6,linewidth=1),nrow=1), fill="none") +
  annotation_logticks(sides="b",scaled=TRUE,short=unit(.08,"cm"),mid=unit(.14,"cm"),long=unit(.22,"cm"),color="black",linewidth=.4) +
  coord_cartesian(clip="on") +
  labs(title="Linear TSS-distance of peak\u2013gene assignment categories (hop0)",
    subtitle="Pooled across hop0 modes; row-weighted | pair-level, descriptive (stats at mode level)",
    x="Distance to TSS (kb, log)", y="Density", color=NULL) +
  theme_classic(10) + theme(legend.position="top",legend.justification="center",
    legend.key.size=unit(.4,"cm"),legend.key.spacing.x=unit(.3,"cm"),legend.text=element_text(size=8.5,face="bold"),
    plot.title=element_text(face="bold",hjust=.5,size=11),plot.subtitle=element_text(size=8.5,color="grey35",hjust=.5),
    axis.ticks.x=element_line(color="black",linewidth=.4),axis.ticks.length.x=unit(.22,"cm"),
    plot.margin=margin(10,10,15,10))

ggsave(file.path(out_dir,"Distance_Density_Overlay_LogTicks.pdf"), p_d, width=4.5, height=4.5)

# ═══ 5. BOXPLOT ──────────────────────────────────────────────
message(">>> Boxplot...")
p_b <- ggplot(dist_hop0, aes(x=GeneSet_Clean, y=Dist_kb_plot, fill=GeneSet_Clean)) +
  geom_hline(yintercept=c(10,100), linetype="dashed", color="grey55", linewidth=.35) +
  geom_boxplot(outlier.shape=16,outlier.size=.15,outlier.alpha=.08,coef=1.5,alpha=.8,width=.5,color="grey20") +
  scale_y_log10(breaks=c(.1,1,10,100,1000,10000), labels=c("0.1","1","10","100","1,000","10,000"),
    expand=expansion(mult=c(.05,.08))) +
  scale_fill_manual(values=dist_colors) + facet_wrap(~Pipeline, ncol=4) +
  annotation_logticks(sides="l",scaled=TRUE,short=unit(.07,"cm"),mid=unit(.12,"cm"),long=unit(.18,"cm"),color="black",linewidth=.35) +
  coord_cartesian(clip="on") +
  labs(title="Peak\u2013Gene Linear Distance across Pipelines (hop0)",
    subtitle=sprintf("%d modes | pair-level, descriptive | mode-level stats in Distance_ModeLevel_Stats.csv",
      length(hop0_modes)), x=NULL, y="Distance to TSS (kb, log)") +
  theme_bw(10) + theme(legend.position="none",
    axis.text.x=element_text(angle=40,hjust=1,vjust=1,size=8.5,color="black"),
    strip.background=element_rect(fill="grey95",color="grey80"),strip.text=element_text(face="bold",size=9),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    axis.ticks.y=element_line(color="black",linewidth=.35),axis.ticks.length.y=unit(.18,"cm"),
    plot.title=element_text(face="bold",hjust=.5,size=11),plot.subtitle=element_text(size=8.5,color="grey35",hjust=.5),
    plot.margin=margin(10,12,25,12))

ggsave(file.path(out_dir,"Distance_Boxplot_Pipeline_LogTicks.pdf"), p_b, width=7.5, height=4.5)

message("Done: ", out_dir)
