#!/usr/bin/env Rscript
# Unique GSEA: density + boxplot + ES forest — v13 data, v8 styling
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = FALSE)), "figure_path_config.R"))
configure_looplook_library()
library(dplyr); library(ggplot2); library(stringr); library(rstatix); library(looplook); library(boot); library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(parallel)

# ═══ CONFIG ═══
cfg <- list(
  out_base   = looplook_env_path("LOOPLOOK_OUT_BASE"),
  rdata_path = looplook_env_path("LOOPLOOK_RDATA_FILE"),
  diff_file  = looplook_env_path("LOOPLOOK_DIFF_FILE"),
  peak_file  = looplook_env_path("LOOPLOOK_PEAK_FILE"),
  n_cores    = max(12, detectCores() - 1)
)
out_base <- cfg$out_base
n_iterations <- 300
recompute_cache <- TRUE  # TRUE=recompute ES; FALSE=load .es_cache.rds and re-plot only
uniq_colors <- c("only_looplook"="#E64B35","intersection"="#00A087","only_ChIPseeker"="#4DBBD5")

# ── Dependencies for ES forest ──
message("Loading RData + diff + peak...")
tmp_env <- new.env(parent=emptyenv()); load(cfg$rdata_path, envir=tmp_env)
res<-tmp_env$res; res2<-tmp_env$res2; refined_res<-tmp_env$refined_res; refined_res2<-tmp_env$refined_res2
cr<-tmp_env$cr; cr2<-tmp_env$cr2; cr_only<-tmp_env$cr_only; cr2_only<-tmp_env$cr2_only; rm(tmp_env)

diff_df <- read.csv(cfg$diff_file, stringsAsFactors=FALSE, row.names=1)
diff_df_for_stat <- diff_df[!is.na(diff_df$log2FoldChange) & is.finite(diff_df$log2FoldChange),]
diff_df_for_stat$gene <- toupper(rownames(diff_df_for_stat))
diff_df_for_stat$signal_lfc <- -diff_df_for_stat$log2FoldChange
if (!"pvalue" %in% colnames(diff_df_for_stat)) diff_df_for_stat$pvalue <- 1
diff_df_for_stat$pvalue[is.na(diff_df_for_stat$pvalue)] <- 1
gene_stat <- diff_df_for_stat

peak_anno <- ChIPseeker::annotatePeak(cfg$peak_file, tssRegion=c(-5000,5000),
  TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
peak_anno_df <- as.data.frame(peak_anno); peak_anno_df$SYMBOL <- toupper(peak_anno_df$SYMBOL)

# ── Helpers ──
rank_biserial_point <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  2 * unname(wt$statistic) / (length(x) * length(y)) - 1
}

rank_biserial <- function(x, y, R = 2000L, seed = 42L) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(est = NA_real_, lo = NA_real_, hi = NA_real_))
  wt <- wilcox.test(x, y, exact = FALSE)
  U <- as.numeric(wt$statistic); r <- 2 * U / (n1 * n2) - 1
  boot_fn <- function(d, i) { xb <- d$v[i[1:n1]]; yb <- d$v[i[(n1+1):(n1+n2)]]
    if (length(xb) < 2 || length(yb) < 2) return(NA_real_)
    wtb <- wilcox.test(xb, yb, exact = FALSE)
    2 * as.numeric(wtb$statistic) / (length(xb) * length(yb)) - 1 }
  set.seed(seed)
  bd <- data.frame(v = c(x, y), stringsAsFactors = FALSE)
  br <- tryCatch(boot::boot(bd, boot_fn, R = R, strata = rep(c(1, 2), c(n1, n2))), error = function(e) NULL)
  ci_obj <- if (!is.null(br)) tryCatch(boot::boot.ci(br, type = "bca", conf = 0.95), error = function(e) boot::boot.ci(br, type = "perc", conf = 0.95)) else NULL
  if (!is.null(ci_obj)) {
    if (!is.null(ci_obj$bca)) { ci_lo <- ci_obj$bca[4]; ci_hi <- ci_obj$bca[5] }
    else if (!is.null(ci_obj$percent)) { ci_lo <- ci_obj$percent[4]; ci_hi <- ci_obj$percent[5] }
    else { ci_lo <- NA_real_; ci_hi <- NA_real_ }
  } else { ci_lo <- NA_real_; ci_hi <- NA_real_ }
  list(est = r, lo = ci_lo, hi = ci_hi)
}
get_target_col <- function(md){bc<-if(identical(md$map,"promoter"))"Regulated_promoter_genes" else "Assigned_Target_Genes";if(md$fill)paste0(bc,"_Filled")else bc}
get_obj <- function(nm){if(grepl("^anno_",nm)){if(grepl("_hop1$",nm))"res2" else "res"}else if(grepl("^refined_",nm)){if(grepl("_hop1$",nm))"refined_res2" else "refined_res"}else if(grepl("^chrom_only_",nm)){if(grepl("_hop1$",nm))"cr2_only" else "cr_only"}else if(grepl("^chrom_",nm)){if(grepl("_hop1$",nm))"cr2" else "cr"}else NA_character_}
is_distal <- grepl("Distal Intergenic",peak_anno_df$annotation)
is_genic <- grepl("Exon|Intron|UTR|Downstream",peak_anno_df$annotation)&!is_distal
is_prom <- grepl("Promoter",peak_anno_df$annotation)
chipseeker_genes <- unique(peak_anno_df$SYMBOL[!is.na(peak_anno_df$SYMBOL)])
chipseeker_promoter_genes <- unique(peak_anno_df$SYMBOL[is_prom])
get_chip_reference <- function(md){if(identical(md$map,"promoter"))chipseeker_promoter_genes else chipseeker_genes}

# ══════════════════════════════════════════════════════════════════════════
for (sz in 200L) {
  out_dir <- file.path(out_base, paste0("size", sz))
  uniq_dir <- file.path(out_base, paste0("size", sz), "tmp_unique")
  uniq_files <- list.files(uniq_dir, pattern="uniq_.*\\.csv$", full.names=TRUE)
  if (length(uniq_files) == 0) next

  unique_df <- bind_rows(lapply(uniq_files, function(f) {
    df <- read.csv(f, stringsAsFactors=FALSE); df$ID <- as.character(df$ID); df$Mode <- as.character(df$Mode); df
  }))
  unique_df$ID <- factor(unique_df$ID, levels=c("only_looplook","intersection","only_ChIPseeker"))
  message(sprintf("size=%d: %d files, %d rows", sz, length(uniq_files), nrow(unique_df)))

  # Build modes + mode_genes from unique_df mode names
  all_modes <- unique(as.character(unique_df$Mode))
  modes <- list(); mode_genes <- list()
  for(nm in all_modes){
    ao <- get_obj(nm); ann_obj <- if(!is.na(ao)) get(ao) else NULL
    md <- list(name=nm, ann=ann_obj, map=if(grepl("_promoter_",nm))"promoter" else "all",
               fill=grepl("_T(_hop1)?$",nm), near=FALSE)
    modes[[length(modes)+1]] <- md
    if(!is.null(ann_obj)){bi<-ann_obj$target_annotation;dc<-get_target_col(md)
      mode_genes[[nm]] <- if(!is.null(bi)&&dc%in%colnames(bi)) unique(toupper(looplook:::clean_gene_names(bi[[dc]],"[;,]"))) else character(0)
    } else {mode_genes[[nm]] <- character(0)}
  }

  methods <- c("all"=NULL, "anno"="anno", "refined"="refined", "chromatin"="chromatin", "chromatin_only"="chromatin_only")
  for (mn in names(methods)) {
    m <- methods[[mn]]
    for (hop in 0L) {
      if (hop == 0) keep <- grep("_(promoter|all)_(F|T)$", unique(unique_df$Mode), value=TRUE)
      else          keep <- grep("_hop1$", unique(unique_df$Mode), value=TRUE)
      if (!is.null(m)) {
        if (m == "chromatin_only") keep <- keep[grepl("^chrom_only_", keep)]
        else if (m == "chromatin") keep <- keep[grepl("^chrom_", keep) & !grepl("^chrom_only_", keep)]
        else keep <- keep[grepl(paste0("^", m, "_"), keep)]
      }
      if (length(keep) == 0) next

      sub <- unique_df %>% filter(Mode %in% keep) %>%
        mutate(BaseMode = Mode %>% str_replace("^anno_","") %>% str_replace("^refined_","") %>%
                 str_replace("^chrom_only_","") %>% str_replace("^chrom_",""))
      sub_med <- sub %>% group_by(BaseMode, ID) %>% summarise(Median_NES = median(NES, na.rm=TRUE), .groups="drop")

      hlab <- paste0("hop = ", hop); if (!is.null(m)) hlab <- paste(hlab, "|", m)
      nm <- if (is.null(m)) paste0("hop", hop) else paste0("hop", hop, "_", m)

      # Density (v8 style)
      p <- ggplot(sub, aes(x=NES, fill=ID)) + geom_density(alpha=0.5, color=NA) +
        geom_vline(data=sub_med, aes(xintercept=Median_NES, color=ID), linetype="dashed", linewidth=0.4) +
        scale_fill_manual(values=uniq_colors) + scale_color_manual(values=uniq_colors, guide="none") +
        facet_wrap(~BaseMode, ncol=4, scales="free_y") +
        labs(title=paste("NES by Gene Set Origin (ChIPseeker) \u2014", hlab),
             subtitle=sprintf("%d subsamples per mode", n_iterations), x="NES", y="Density") +
        theme_minimal() + theme(panel.border=element_rect(color="grey10", fill=NA, linewidth=0.4),
          panel.grid=element_blank(), legend.position="bottom",
          strip.text=element_text(size=9, face="bold", color="black"),
          plot.title=element_text(face="bold", size=13.5, hjust=0.5),
          plot.subtitle=element_text(size=9.5, color="grey40", hjust=0.5))
      ggsave(file.path(out_dir, paste0("GSEA_Unique_Density_", nm, "_size", sz, "_newplot.pdf")), p, width=12, height=4)

      # Boxplot (v8 style)
      p2 <- ggplot(sub, aes(x=ID, y=NES)) +
        geom_hline(yintercept=0, linetype="dashed", color="grey80", linewidth=0.3) +
        geom_boxplot(aes(fill=ID), outlier.shape=16, outlier.size=0.8, outlier.alpha=0.3,
                     linewidth=0.3, alpha=0.6, width=0.65, fatten=1, notch=FALSE) +
        scale_fill_manual(values=uniq_colors) +
        scale_y_continuous(expand=expansion(mult=c(0.05,0.3)),
                           breaks=function(x)seq(floor(min(x)),ceiling(max(x)),by=1)) +
        facet_wrap(~BaseMode, ncol=4, scales="free_y") +
        labs(title=paste("NES by Gene Set Origin \u2014", hlab),
             subtitle=sprintf("%d subsamples per mode | Boxplot = descriptive only", n_iterations),
             x=NULL, y="Normalized Enrichment Score (NES)") +
        theme_minimal() + theme(panel.border=element_rect(color="grey10", fill=NA, linewidth=0.4),
          panel.grid=element_blank(), axis.ticks.x=element_blank(),
          legend.position="bottom", strip.text=element_text(size=9, face="bold"),
          plot.title=element_text(face="bold", size=13.5, hjust=0.5),
          plot.subtitle=element_text(size=9.5, color="grey40", hjust=0.5))
      ggsave(file.path(out_dir, paste0("GSEA_Unique_Boxplot_", nm, "_size", sz, "_newplot.pdf")), p2, width=8, height=4)
    }
    message(sprintf("  %s done", mn))
  }

  # ── Effect Size forest (cache-aware) ──
  cache_fn <- file.path(out_dir, ".es_cache.rds")
  if (file.exists(cache_fn) && !recompute_cache) {
    message("  Loading ES from cache for re-plotting...")
    es_cached <- readRDS(cache_fn)
    uniq_es_df <- es_cached$pairwise
    es3_df     <- es_cached$three_set
    # re-plot pairwise
    if (nrow(uniq_es_df) > 1) {
      for (hv in "hop0") {
        hv_label <- if (hv=="hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
        es_sub <- uniq_es_df %>% filter(Hop==hv_label) %>% arrange(desc(ES))
        if (nrow(es_sub) < 2) next
        p_es <- ggplot(es_sub, aes(x=ES, y=reorder(Mode,ES), color=Pipeline)) +
          geom_vline(xintercept=0, linetype="dashed", color="grey50") +
          geom_errorbarh(aes(xmin=ES_lo, xmax=ES_hi), height=0.5, linewidth=1.0) +
          geom_point(size=3, shape=18) +
          scale_color_manual(values=c("Basic"="#1F77B4","E-refined"="#9467BD","C-refined"="#FFA500","I-refined"="#E64B35")) +
          labs(title=sprintf("Unique GSEA: Effect Size %s size=%d", hv_label, sz),
               subtitle=sprintf("Only_looplook vs Only_ChIPseeker | %d modes", nrow(es_sub)),
               x="Rank-Biserial Effect Size", y=NULL) +
          theme_classic() + theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
            plot.subtitle=element_text(size=9,color="grey40",hjust=0.5),
            legend.position="bottom", axis.text.y=element_text(size=8))
        ggsave(file.path(out_dir, sprintf("GSEA_Unique_EffectSize_Forest_%s_size%d_newplot.pdf", hv, sz)),
          p_es, width=10, height=8)
      }
    }
    # re-plot three-set
    if (nrow(es3_df) > 0) {
      for (hv in "hop0") {
        es3_sub <- es3_df %>% filter(Hop==hv)
        if (nrow(es3_sub) < 3) next
        p <- ggplot(es3_sub, aes(x=ES, y=reorder(Mode,ES), color=Pipeline)) +
          geom_vline(xintercept=0, linetype="dashed", color="grey50") +
          geom_errorbarh(aes(xmin=ES_lo, xmax=ES_hi), height=0.5, linewidth=1.0) +
          geom_point(size=3, shape=18) +
          scale_color_manual(values=c("Basic"="#1F77B4","E-refined"="#9467BD","C-refined"="#FFA500","I-refined"="#E64B35")) +
          facet_wrap(~GeneSet, ncol=3, scales="free_x") +
          labs(title=sprintf("Three-Set ES vs BG %s size=%d", hv, sz),
               subtitle="rank-biserial, 2000 bootstrap, bca CI", x="Effect Size", y=NULL) +
          theme_classic()+theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
            plot.subtitle=element_text(size=9,color="grey40",hjust=0.5),
            legend.position="bottom", axis.text.y=element_text(size=7), strip.text=element_text(face="bold"))
        ggsave(file.path(out_dir, sprintf("GSEA_Unique_ThreeSet_ES_Forest_%s_size%d_newplot.pdf", hv, sz)),
          p, width=16, height=10)
      }
    }
    message(sprintf("  Re-plotted from cache: pairwise=%d three_set=%d", nrow(uniq_es_df), nrow(es3_df)))
    next  # skip to next sample_size
  }
  # ── Original ES computation (if cache not used) ──
  es_raw <- parallel::mclapply(all_modes, function(md_name) {
    md_genes <- mode_genes[[md_name]]
    if (length(md_genes) == 0) return(NULL)
    md <- modes[[which(sapply(modes, function(m) m$name == md_name))[1]]]
    loop_genes <- intersect(unique(toupper(md_genes)), gene_stat$gene)
    chip_g     <- intersect(get_chip_reference(md), gene_stat$gene)
    only_loop  <- intersect(setdiff(loop_genes, chip_g), gene_stat$gene)
    only_chip  <- intersect(setdiff(chip_g, loop_genes), gene_stat$gene)
    if (length(only_loop) >= 5 && length(only_chip) >= 5) {
      x <- gene_stat$signal_lfc[gene_stat$gene %in% only_loop]
      y <- gene_stat$signal_lfc[gene_stat$gene %in% only_chip]
      rb <- rank_biserial(x, y)
      data.frame(Mode=md_name, N_loop=length(only_loop), N_chip=length(only_chip),
        ES=rb$est, ES_lo=rb$lo, ES_hi=rb$hi, stringsAsFactors=FALSE)
    } else NULL
  }, mc.cores = min(8L, cfg$n_cores))
  es_raw <- es_raw[!vapply(es_raw, is.null, logical(1))]
  uniq_es_list <- if (length(es_raw) > 0) setNames(es_raw, sapply(es_raw, function(x) x$Mode[1])) else list()
  if (length(uniq_es_list) > 1) {
    uniq_es_df <- bind_rows(uniq_es_list) %>%
      mutate(Pipeline=factor(case_when(
        grepl("^anno_",Mode)~"Basic",grepl("^refined_",Mode)~"E-refined",
        grepl("^chrom_only_",Mode)~"C-refined",grepl("^chrom_",Mode)~"I-refined",TRUE~"Other"),
        levels=c("Basic","E-refined","C-refined","I-refined")),
        Hop=factor(ifelse(grepl("_hop1$",Mode),"hop1 (primary+expanded)","hop0 (primary)"),
                   levels=c("hop0 (primary)","hop1 (primary+expanded)")))
    for (hv in "hop0") {
      hv_label <- if (hv=="hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
      es_sub <- uniq_es_df %>% filter(Hop==hv_label) %>% arrange(desc(ES))
      if (nrow(es_sub) < 2) next
      p_es <- ggplot(es_sub, aes(x=ES, y=reorder(Mode,ES), color=Pipeline)) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        geom_errorbarh(aes(xmin=ES_lo, xmax=ES_hi), height=0.5, linewidth=1.0) +
        geom_point(size=3, shape=18) +
        scale_color_manual(values=c("Basic"="#1F77B4","E-refined"="#9467BD","C-refined"="#FFA500","I-refined"="#E64B35")) +
        labs(title=sprintf("Unique GSEA: Effect Size %s size=%d", hv_label, sz),
             subtitle=sprintf("Only_looplook vs Only_ChIPseeker (rank-biserial, 2000 bootstrap, bca CI) | %d modes", nrow(es_sub)),
             x="Rank-Biserial Effect Size", y=NULL) +
        theme_classic() + theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(size=9,color="grey40",hjust=0.5),
          legend.position="bottom", axis.text.y=element_text(size=8))
      ggsave(file.path(out_dir, sprintf("GSEA_Unique_EffectSize_Forest_%s_size%d_newplot.pdf", hv, sz)),
        p_es, width=10, height=8)
    }
    message(sprintf("  ES forest: %d modes", nrow(uniq_es_df)))
  }
  # Default empty DFs in case ES returns no data
  if (!exists("uniq_es_df")) uniq_es_df <- data.frame()
  es3_df <- data.frame()

  # ── Three-set ES forest (parallelized) ──
  es3_raw <- parallel::mclapply(all_modes, function(md_name) {
    md_genes <- mode_genes[[md_name]]
    if (length(md_genes) == 0) return(NULL)
    md <- modes[[which(sapply(modes, function(m) m$name == md_name))[1]]]
    loop_genes <- intersect(unique(toupper(md_genes)), gene_stat$gene)
    chip_g     <- intersect(get_chip_reference(md), gene_stat$gene)
    only_loop  <- intersect(setdiff(loop_genes, chip_g), gene_stat$gene)
    inter_set  <- intersect(loop_genes, chip_g)
    only_chip  <- intersect(setdiff(chip_g, loop_genes), gene_stat$gene)
    bg_genes   <- setdiff(setdiff(gene_stat$gene, loop_genes), chip_g)
    rows <- list()
    for (setname in c("only_loop","intersection","only_chip")) {
      g <- switch(setname, only_loop=only_loop, intersection=inter_set, only_chip=only_chip)
      if (length(g) >= 5 && length(bg_genes) >= 5) {
        x <- gene_stat$signal_lfc[gene_stat$gene %in% g]
        y <- gene_stat$signal_lfc[gene_stat$gene %in% bg_genes]
        rb <- rank_biserial(x, y)
        rows[[length(rows)+1]] <- data.frame(Mode=md_name, GeneSet=setname, N=length(g),
          ES=rb$est, ES_lo=rb$lo, ES_hi=rb$hi, stringsAsFactors=FALSE)
      }
    }
    if (length(rows) > 0) bind_rows(rows) else NULL
  }, mc.cores = min(8L, cfg$n_cores))
  es3_raw <- es3_raw[!vapply(es3_raw, is.data.frame, logical(1))]
  es3_list <- if (length(es3_raw) > 0) unname(do.call(rbind, es3_raw)) else list()
  if (length(es3_list) > 0) {
    es3_df <- bind_rows(es3_list) %>%
      mutate(Pipeline=factor(case_when(
        grepl("^anno_",Mode)~"Basic",grepl("^refined_",Mode)~"E-refined",
        grepl("^chrom_only_",Mode)~"C-refined",grepl("^chrom_",Mode)~"I-refined",TRUE~"Other"),
        levels=c("Basic","E-refined","C-refined","I-refined")),
        Hop=ifelse(grepl("_hop1$",Mode),"hop1","hop0"),
        GeneSet=factor(GeneSet, levels=c("only_loop","intersection","only_chip"),
                       labels=c("Only looplook","Intersection","Only ChIPseeker")))
    for (hv in "hop0") {
      es3_sub <- es3_df %>% filter(Hop==hv)
      if (nrow(es3_sub) < 3) next
      p <- ggplot(es3_sub, aes(x=ES, y=reorder(Mode,ES), color=Pipeline)) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        geom_errorbarh(aes(xmin=ES_lo, xmax=ES_hi), height=0.5, linewidth=1.0) +
        geom_point(size=3, shape=18) +
        scale_color_manual(values=c("Basic"="#1F77B4","E-refined"="#9467BD","C-refined"="#FFA500","I-refined"="#E64B35")) +
        facet_wrap(~GeneSet, ncol=3, scales="free_x") +
        labs(title=sprintf("Three-Set Effect Size vs Background %s size=%d", hv, sz),
             subtitle="rank-biserial, 2000 bootstrap, bca CI", x="Effect Size", y=NULL) +
        theme_classic() + theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(size=9,color="grey40",hjust=0.5),
          legend.position="bottom", axis.text.y=element_text(size=7),
          strip.text=element_text(face="bold"))
      ggsave(file.path(out_dir, sprintf("GSEA_Unique_ThreeSet_ES_Forest_%s_size%d_newplot.pdf", hv, sz)),
        p, width=16, height=10)
    }
    message(sprintf("  Three-set ES forest: %d rows", nrow(es3_df)))
    }
    saveRDS(list(pairwise=uniq_es_df, three_set=es3_df),
      file.path(out_dir, ".es_cache.rds"))
    message("    ES cache saved")

  message(sprintf("size=%d complete", sz))
}
message(sprintf("\nDone. Output: %s", out_base))
