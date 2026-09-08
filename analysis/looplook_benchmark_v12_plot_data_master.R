# ── brd4_r21_resampling_v4_density_lb_main_optimized.R ──
# BRD4 R21 多模式分析：looplook vs ChIPseeker (OPTIMIZED VERSION)
# 16 primary looplook assignment modes: annotation/refined/chromatin/chromatin_only × all/promoter × strict/filled, hop0
#   hop0 = primary targets (path_length <= 1)
#   hop1 = primary targets + eligible expanded targets (path_length > 1)
#   Display labels: chromatin_only = "Chromatin, no expression reclassification"
# 功能: GSEA重抽样、密度图、箱线图、ridge图、Venn图、GO富集、距离分析、
#       Signal Density、Validation、Convergence QC、hop可解释性、跨TSS汇总，
#       以及 expanded-only GSEA/effect-size/表达匹配/拓扑敏感性分析。

# ── 加载所有依赖包 ──────────────────────────────────────────────────────────────

# unlink("resampling_v4_lb/tmp_sim", recursive = TRUE)
# unlink("resampling_v4_lb/tmp_unique", recursive = TRUE)
.required_env_path <- function(name) {
  value <- trimws(Sys.getenv(name, ""))
  if (!nzchar(value)) stop("Set environment variable ", name, " to an absolute path.", call. = FALSE)
  path.expand(value)
}
.r_lib <- trimws(Sys.getenv("LOOPLOOK_R_LIB", ""))
if (nzchar(.r_lib)) .libPaths(c(path.expand(.r_lib), .libPaths()))

library(looplook)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(rstatix)
library(openxlsx)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)
library(boot)           # bootstrap CI for rank_biserial effect size
library(GenomicRanges)  # peak-gene distance computation
library(parallel)        # mclapply for parallel GSEA loops
library(ggvenn)          # Venn diagrams
library(ChIPseeker)      # peak annotation
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # hg38 gene models
if (isTRUE(capabilities("cairo"))) {
  options(bitmapType = "cairo", device = "cairo_pdf")
}

# Close graphics devices only in non-interactive batch runs. This avoids
# destroying user-owned RStudio devices while retaining the headless-server
# protection used after forked computations.
close_batch_devices <- function() {
  if (!interactive()) {
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
  }
  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 路径配置 — 修改此处即可适配不同环境
# ══════════════════════════════════════════════════════════════════════════════

# 基础目录（本地 Mac / 服务器分别设置）
cfg <- list(
  # ── 计算资源 ──
  n_cores     = min(10, max(1, detectCores() - 1)),   # 并行核心数，上限 8 核（降低内存压力）
  # ── 输入数据 ──
  data_base   = .required_env_path("LOOPLOOK_DATA_BASE"),
  rdata_base  = .required_env_path("LOOPLOOK_RDATA_BASE"),
  # ── 输出目录 ──
  out_base    = .required_env_path("LOOPLOOK_OUT_BASE"),
  # ── 单文件路径（相对于 data_base） ──
  expr_file   = "TPM_clean.txt",
  diff_file   = "clean_filter_count4_DESeq2_gene_lps141_ARV825_vs_lps141_dmso.csv",
  meta_file   = "ARV825_coldata.txt",
  target_file = "lps141_folsl2_peaks.narrowPeak",
  # ── 跨 TSS 窗口 RData（服务器脚本生成，见 服务器脚本/ 目录） ──
  #     主分析用 "5kb"；跨 TSS 汇总读取全部 4 个
  tss_subdirs = list(
    "1kb"  = "tss1000/tss_1000_res_chromatin.RData",
    "2kb"  = "tss2000/tss_2000_res_chromatin.RData",
    "5kb"  = "tss5000/tss_5000_res_chromatin.RData",
    "10kb" = "tss10000/tss_10000_res_chromatin.RData"),
  # ── Super-enhancer annotation (optional) ──
  # Provide BED files for SE and TE regions. Both relative to data_base.
  # If provided and valid, BRD4 peaks are classified by overlap with the
  # ROSE-called region sets. ROSE SE/TE labels are treated as authoritative;
  # no secondary ChIPseeker distal filter is applied.
  #   Super_Enhancer = peaks overlapping the ROSE SE BED
  #   Typical_Enhancer = peaks overlapping the ROSE TE BED (and not SE)
  # Set to NULL to skip SE/TE stratification (default: 4 peak categories).
  se_bed_file = "141_H3K27ac_peaks_Stitched_only_Super.bed",
  te_bed_file = "141_H3K27ac_peaks_Stitched_only_TE.bed"
 )

# Optional path overrides used by the modular entry-point scripts. Empty
# environment variables leave the cfg values above unchanged.
.env_path_override <- function(name, current) {
  value <- trimws(Sys.getenv(name, ""))
  if (nzchar(value)) value else current
}
cfg$data_base <- .env_path_override("LOOPLOOK_DATA_BASE", cfg$data_base)
cfg$rdata_base <- .env_path_override("LOOPLOOK_RDATA_BASE", cfg$rdata_base)
cfg$out_base <- .env_path_override("LOOPLOOK_OUT_BASE", cfg$out_base)
rm(.env_path_override)

# ── 组装完整路径 ──
setwd(cfg$rdata_base)
out_dir      <- cfg$out_base
base_out_dir <- out_dir  # 提前定义，供后续 gene_panel_dir 等使用
expr_path   <- file.path(cfg$data_base, cfg$expr_file)
diff_path   <- file.path(cfg$data_base, cfg$diff_file)
target_file <- file.path(cfg$data_base, cfg$target_file)

dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
unlink(file.path(base_out_dir, c(
  "RUN_COMPLETED.ok", "RUN_COMPUTE_COMPLETED.ok", "RUN_PLOTS_COMPLETED.ok",
  "RUN_FINISHED_WITH_ISSUES.txt"
)))
n_cores <- cfg$n_cores
options(mc.cores = n_cores)  # 减少 mclapply 端口协商信息
message("    Using ", n_cores, " parallel cores (detected ", detectCores(), " total)")

# ── 可调参数 ──────────────────────────────────────────────────────────────────
code_version <- "looplook_benchmark_v12_plot_data"
# Numerical cache compatibility is intentionally separated from the script
# file MD5. Plotting, logging, checkpointing, or cache-policy edits must not
# invalidate already-computed NES values.
gsea_cache_logic_version <- "looplook_benchmark_v9_bioinformatics_hardened"
mode_cache_logic_version <- "looplook_mode_gene_extraction_v1"
analysis_run_start <- Sys.time()

# Cache policy:
#   "resume"  = reuse validated completed outputs and resume missing iteration chunks
#   "rebuild" = delete managed caches/outputs and recompute from scratch
# Do not switch to rebuild merely because a previous run was interrupted.
cache_policy <- tolower(trimws(Sys.getenv("LOOPLOOK_CACHE_POLICY", "resume")))
if (!cache_policy %in% c("resume", "rebuild")) {
  stop("LOOPLOOK_CACHE_POLICY/cache_policy must be either 'resume' or 'rebuild'.")
}
force_recompute <- identical(cache_policy, "rebuild")

# Plot execution policy:
#   "queue"     = build plot objects and save them as independent render tasks;
#                 no graphics device is opened in the compute process.
#   "immediate" = render plots in the current process (legacy behaviour).
#   "off"       = neither render nor queue plots.
# Queue mode is the safe default for long server runs.
plot_policy <- tolower(trimws(Sys.getenv("LOOPLOOK_PLOT_POLICY", "queue")))
if (!plot_policy %in% c("queue", "immediate", "off")) {
  stop("LOOPLOOK_PLOT_POLICY/plot_policy must be queue, immediate, or off.")
}
plot_queue_dir <- file.path(base_out_dir, "plot_queue")
dir.create(plot_queue_dir, recursive = TRUE, showWarnings = FALSE)
message("    Plot policy: ", plot_policy,
  if (identical(plot_policy, "queue")) paste0("; queue: ", plot_queue_dir) else "")

# GSEA reliability controls. Four workers are usually faster in practice than
# eight when each child holds a ranked list and clusterProfiler/fgsea objects.
gsea_n_cores <- max(1L, n_cores)  # use all cores (was min(4L, n_cores))
gsea_chunk_size <- 300L  # one chunk = all iterations (no checkpoint I/O overhead)
gsea_iteration_timeout_sec <- 900L
message("    GSEA workers: ", gsea_n_cores, "; checkpoint chunk: ", gsea_chunk_size,
  " iterations; cache policy: ", cache_policy)

n_iterations <- 300L
sample_sizes <- 200L
primary_sample_size <- 200L
primary_hop <- 0L
bootstrap_R <- 2000L
bootstrap_seed <- 42L
# Expanded-only sensitivity settings. Lower these for smoke tests.
expanded_interval_R <- min(1000L, bootstrap_R)
expanded_interval_max_per_group <- 500L
expanded_expression_match_iterations <- 500L
expanded_expression_match_bins <- 10L
expanded_expression_match_max_per_bin <- 50L
gsea_min_size <- 10L  # Minimum gene set size for GSEA (must match minGSSize in clusterProfiler::GSEA)
gsea_max_size <- 50000L

# ── Script path and early run signature (for cache validation) ────────────────
script_path <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile) && length(ofile) == 1L && !is.na(ofile) && nzchar(ofile)) {
    normalizePath(ofile, mustWork = FALSE)
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0L) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    else NA_character_
  }
}, error = function(e) { warning("Unable to resolve script path: ", conditionMessage(e)); NA_character_ })
script_md5 <- if (length(script_path) == 1L && !is.na(script_path) && file.exists(script_path)) {
  unname(as.character(tools::md5sum(script_path)))
} else NA_character_
message("    Script: ", ifelse(is.na(script_path), "unknown", basename(script_path)),
  " MD5=", ifelse(is.na(script_md5), "NA", substr(script_md5, 1, 8)))

# Preliminary run signature (inputs only, before manifest is built)
main_rdata_path <- file.path(cfg$rdata_base, cfg$tss_subdirs[["5kb"]])
meta_path <- file.path(cfg$data_base, cfg$meta_file)
preliminary_signature <- list(
  ScriptMD5 = script_md5,
  DiffMD5 = unname(as.character(tryCatch(tools::md5sum(diff_path), error = function(e) NA))),
  ExprMD5 = unname(as.character(tryCatch(tools::md5sum(expr_path), error = function(e) NA))),
  MetadataMD5 = unname(as.character(tryCatch(tools::md5sum(meta_path), error = function(e) NA))),
  TargetMD5 = unname(as.character(tryCatch(tools::md5sum(target_file), error = function(e) NA))),
  MainRDataMD5 = unname(as.character(tryCatch(tools::md5sum(main_rdata_path), error = function(e) NA))),
  LooplookVersion = as.character(tryCatch(packageVersion("looplook"), error = function(e) "unknown")),
  Iterations = n_iterations,
  SampleSizes = sample_sizes,
  PrimarySampleSize = primary_sample_size,
  GSEAMinSize = gsea_min_size,
  GSEAMaxSize = gsea_max_size)

# ── Early input validation (before expensive computation) ────────────────────
required_input_files <- c(diff_path, expr_path, meta_path, target_file, main_rdata_path)
missing_inputs <- required_input_files[!file.exists(required_input_files)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input files: ", paste(missing_inputs, collapse = "; "))
}
# Validate metadata schema early
meta_df <- read.table(meta_path, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE)
required_meta_cols <- c("sample", "condition")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_df))
if (length(missing_meta_cols) > 0L) {
  stop("Metadata missing columns: ", paste(missing_meta_cols, collapse = ", "))
}
meta_df$sample <- trimws(as.character(meta_df$sample))
meta_df$condition <- trimws(as.character(meta_df$condition))
if (anyNA(meta_df$sample) || any(!nzchar(meta_df$sample))) {
  stop("Metadata contains NA or blank sample IDs after trim.")
}
if (anyNA(meta_df$condition) || any(!nzchar(meta_df$condition))) {
  stop("Metadata contains NA or blank condition labels after trim.")
}
if (anyDuplicated(meta_df$sample)) {
  dup_samples <- unique(meta_df$sample[duplicated(meta_df$sample)])
  stop("Metadata contains duplicated sample IDs: ", paste(dup_samples, collapse = ", "))
}
dmso_cols_early <- meta_df$sample[
  tolower(meta_df$condition) %in% c("dmso", "vehicle", "control")]
if (length(dmso_cols_early) == 0L) {
  stop("Metadata contains no DMSO/vehicle/control samples.")
}
message(sprintf("    Metadata validated: %d samples, %d DMSO/control",
  nrow(meta_df), length(dmso_cols_early)))

# Robust sampling: if a pool is smaller than the requested size, draw 80%
# while never dropping a statistically eligible pool below GSEA minGSSize.
sample_n <- function(pool, size, min_n = gsea_min_size) {
  n_pool <- length(pool)
  if (n_pool < min_n) return(n_pool)
  if (n_pool < size) {
    min(n_pool, max(min_n, round(0.8 * n_pool)))
  } else {
    size
  }
}


# ── 0. Literature-informed BRD4 gene panel ────────────────────────────────────
gene_panel_dir <- file.path(dirname(dirname(dirname(dirname(base_out_dir)))),
                            "looplook_ChIPseeker_analysis_package",
                            "brd4_literature_informed_gene_panel_package")
gene_panel_file <- file.path(gene_panel_dir, "brd4_literature_informed_gene_sets.R")
if (file.exists(gene_panel_file)) {
  source(gene_panel_file)
  gene_panel_md5 <- unname(as.character(tools::md5sum(gene_panel_file)))
  message("    BRD4 gene panel loaded: core=", length(brd4_core_panel),
          " context=", length(brd4_context_panel),
          " associated=", length(brd4_associated_panel),
          " MD5=", substr(gene_panel_md5, 1, 8))
} else {
  warning("Gene panel file not found: ", gene_panel_file)
  gene_panel_md5 <- NA_character_
  brd4_core_panel <- brd4_context_panel <- brd4_associated_panel <- character(0)
}

# ── Safe GSEA / GO wrappers ─────────────────────────────────────────────────
# Unified error tracking for all GSEA and GO enrichment calls.
# Returns list(result = data.frame_or_NULL, status = status_df)
safe_gsea_once <- function(ranked_list, t2g, expected_terms, sample_size_val,
    module, mode, iteration, gsea_args = list()) {
  # Prevent one pathological fgsea/GSEA call from blocking an entire mode.
  # setTimeLimit is best-effort for compiled code, but it reliably catches most
  # R-level stalls and returns a normal status row through tryCatch().
  setTimeLimit(elapsed = gsea_iteration_timeout_sec, transient = TRUE)
  on.exit(
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE),
    add = TRUE
  )
  args <- modifyList(list(
    pvalueCutoff = 1.1, pAdjustMethod = "none",
    minGSSize = gsea_min_size, maxGSSize = gsea_max_size,
    verbose = FALSE, seed = TRUE), gsea_args)
  tryCatch({
    gsea_obj <- do.call(clusterProfiler::GSEA,
      c(list(gene = ranked_list, TERM2GENE = t2g), args))
    res_df <- as.data.frame(gsea_obj@result)
    returned_terms <- if (nrow(res_df) > 0L) unique(as.character(res_df$ID)) else character(0)
    missing_terms <- setdiff(expected_terms, returned_terms)
    status <- if (nrow(res_df) == 0L) "no_result"
              else if (length(missing_terms) > 0L) "partial_terms"
              else "success"
    list(result = res_df, status = data.frame(
      SampleSize = sample_size_val, Module = module, Mode = mode, Iteration = iteration,
      Status = status, ExpectedTerms = paste(expected_terms, collapse = ";"),
      ReturnedTerms = paste(returned_terms, collapse = ";"),
      MissingTerms = paste(missing_terms, collapse = ";"),
      Error = NA_character_, stringsAsFactors = FALSE))
  }, error = function(e) {
    list(result = NULL, status = data.frame(
      SampleSize = sample_size_val, Module = module, Mode = mode, Iteration = iteration,
      Status = "error", ExpectedTerms = paste(expected_terms, collapse = ";"),
      ReturnedTerms = "", MissingTerms = paste(expected_terms, collapse = ";"),
      Error = conditionMessage(e), stringsAsFactors = FALSE))
  })
}

safe_go_once <- function(entrez, ontology, sample_size_val, mode, gene_set_label, universe) {
  tryCatch({
    ego <- clusterProfiler::enrichGO(gene = entrez, universe = universe,
      OrgDb = org.Hs.eg.db, ont = ontology, pAdjustMethod = "BH",
      pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
    res_df <- as.data.frame(ego@result)
    list(result = res_df, status = data.frame(
      SampleSize = sample_size_val, Module = "go", Mode = mode,
      Iteration = NA_integer_, GeneSet = gene_set_label, Ontology = ontology,
      N_Entrez = length(entrez),
      Status = if (nrow(res_df) > 0L) "success" else "no_result",
      ExpectedTerms = NA_character_, ReturnedTerms = NA_character_,
      MissingTerms = NA_character_, Error = NA_character_, stringsAsFactors = FALSE))
  }, error = function(e) {
    list(result = NULL, status = data.frame(
      SampleSize = sample_size_val, Module = "go", Mode = mode,
      Iteration = NA_integer_, GeneSet = gene_set_label, Ontology = ontology,
      N_Entrez = length(entrez), Status = "error",
      ExpectedTerms = NA_character_, ReturnedTerms = NA_character_,
      MissingTerms = NA_character_, Error = conditionMessage(e), stringsAsFactors = FALSE))
  })
}

# ── Gene ID normalisation (centralised trim + NA/blank/duplicate checks) ────
normalise_gene_ids <- function(ids, label) {
  ids <- trimws(as.character(ids))
  if (anyNA(ids) || any(!nzchar(ids))) {
    n_na <- sum(is.na(ids))
    n_blank <- sum(!is.na(ids) & !nzchar(ids))
    stop(label, " contains ", n_na, " NA and ", n_blank, " blank gene IDs after trim.")
  }
  ids_upper <- toupper(ids)
  dup <- unique(ids_upper[duplicated(ids_upper)])
  if (length(dup) > 0L) {
    stop(label, " contains ", length(dup), " duplicate/case-colliding gene IDs after trim+toupper: ",
      paste(head(dup, 20), collapse = ", "))
  }
  ids_upper
}

# ── Normalise rownames of a data.frame/matrix (trim → NA/blank check → toupper → duplicate check)
normalise_rownames <- function(x, label) {
  ids <- tryCatch(rownames(x), error = function(e) stop(label, ": cannot extract rownames: ", conditionMessage(e)))
  if (is.null(ids) || length(ids) == 0L) stop(label, ": rownames are NULL or empty.")
  ids_norm <- normalise_gene_ids(ids, label)
  rownames(x) <- ids_norm
  x
}

# Mode-level status helper (for insufficient_pool/terms/empty_mode_genes)
make_mode_status <- function(sample_size_val, module, mode, status, expected_terms = NA_character_, reason = NA_character_) {
  data.frame(SampleSize = sample_size_val, Module = module, Mode = mode, Iteration = NA_integer_,
    Status = status, ExpectedTerms = expected_terms, ReturnedTerms = "", MissingTerms = expected_terms,
    Error = reason, stringsAsFactors = FALSE)
}

# ── Object MD5 helper (for cache signatures) ────────────────────────────────
md5_object <- function(x) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf, version = 3)
  unname(as.character(tools::md5sum(tf)))
}

# Cache directories are reusable only when their biological inputs, ranked list,
# mode gene sets and analysis parameters match exactly. Structural CSV checks
# alone cannot protect against stale-but-complete results.
ensure_cache_signature <- function(cache_dir, expected_signature,
    force = force_recompute, signature_file = ".cache_signature.rds") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  sig_path <- file.path(cache_dir, signature_file)
  cache_files <- setdiff(list.files(cache_dir, full.names = TRUE), sig_path)

  if (file.exists(sig_path)) {
    old_signature <- tryCatch(readRDS(sig_path), error = function(e) NULL)
    if (is.null(old_signature) || !identical(old_signature, expected_signature)) {
      if (!force && length(cache_files) > 0L) {
        stop(
          "Cache signature mismatch in ", cache_dir,
          ". Set cache_policy='rebuild' or remove this cache directory.",
          call. = FALSE
        )
      }
      unlink(cache_files, recursive = TRUE, force = TRUE)
    }
  } else if (length(cache_files) > 0L) {
    if (!force) {
      stop(
        "Unsigned cache detected in ", cache_dir,
        ". Set cache_policy='rebuild' or remove this cache directory.",
        call. = FALSE
      )
    }
    unlink(cache_files, recursive = TRUE, force = TRUE)
  }

  saveRDS(expected_signature, sig_path, version = 3)
  invisible(sig_path)
}


# ── Atomic RDS writer and chunk-level resumable parallel execution ───────────
atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) {
    unlink(path)
    if (!file.rename(tmp, path)) {
      stop("Failed to atomically write checkpoint: ", path, call. = FALSE)
    }
  }
  invisible(path)
}


# ── Plot queue: decouple analysis from graphics rendering ───────────────────
# In queue mode, every plot is serialized to an independent task. Rendering is
# performed later by small stage scripts, one plot per fresh R process. This
# prevents a slow or broken graphics device from blocking GSEA computation.
infer_plot_group <- function(filename) {
  path <- tolower(gsub("\\\\", "/", as.character(filename)))
  if (grepl("/(go_enrichment|distance_analysis|hop_interpretability|background_sensitivity|signal_density)/", path) ||
      grepl("(go_|distance_|hop_|background|signaldensity|signal_density)", basename(path))) {
    return("functional_distance")
  }
  if (grepl("/(expanded_only|convergence_plots|cross_size_summary|cross_tss_summary|evidence_composition|peak_category_es)/", path) ||
      grepl("(expanded|convergence|cross.?size|cross.?tss|evidence|peakcategory|se_distance|se_vs_te|distancebinned)", basename(path))) {
    return("expanded_evidence")
  }
  if (grepl("/(case_study|expression_matched_es|mode_ranking|summary|lfc_pv)/", path) ||
      grepl("(casestudy|expressionmatched|ranking|summary|lfc_pv)", basename(path))) {
    return("summary_lfcpv")
  }
  "core_gsea"
}

plot_queue_records <- list()

queueable_ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  filename <- normalizePath(as.character(filename), mustWork = FALSE)
  args <- list(...)

  if (identical(plot_policy, "off")) {
    return(invisible(NULL))
  }
  if (identical(plot_policy, "immediate")) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    return(do.call(
      ggplot2::ggsave,
      c(list(filename = filename, plot = plot), args)
    ))
  }

  group <- infer_plot_group(filename)
  task_id <- paste0(group, "__", substr(md5_object(filename), 1L, 16L))
  task_path <- file.path(plot_queue_dir, paste0(task_id, ".plot.rds"))
  meta_path <- file.path(plot_queue_dir, paste0(task_id, ".meta.rds"))

  style_family <- if (exists("infer_plot_family", mode = "function")) infer_plot_family(filename) else group
  plot_ready_data <- if (exists("extract_plot_ready_data", mode = "function")) {
    extract_plot_ready_data(plot)
  } else {
    list(plot_data = plot$data, layer_data = lapply(plot$layers, function(z) z$data))
  }
  plot_data_path <- file.path(plot_queue_dir, paste0(task_id, ".data.rds"))
  atomic_save_rds(plot_ready_data, plot_data_path)

  task <- list(
    TaskID = task_id,
    Group = group,
    StyleFamily = style_family,
    Target = filename,
    Plot = plot,
    PlotDataRDS = plot_data_path,
    GgsaveArgs = args,
    CodeVersion = code_version,
    CreatedAt = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    RequiredPackages = c(
      "ggplot2", "ggpubr", "ggridges", "ggvenn", "scales"
    )
  )
  meta <- data.frame(
    TaskID = task_id,
    Group = group,
    StyleFamily = style_family,
    Target = filename,
    TaskRDS = task_path,
    PlotDataRDS = plot_data_path,
    MetaRDS = meta_path,
    CodeVersion = code_version,
    CreatedAt = task$CreatedAt,
    stringsAsFactors = FALSE
  )

  atomic_save_rds(task, task_path)
  atomic_save_rds(meta, meta_path)
  plot_queue_records[[task_id]] <<- meta
  message("    Plot queued [", group, "]: ", basename(filename))
  invisible(task_path)
}

queueable_print_plot <- function(x, ...) {
  if (identical(plot_policy, "immediate")) {
    return(base::print(x, ...))
  }
  invisible(x)
}

write_plot_queue_manifest <- function() {
  meta_files <- list.files(
    plot_queue_dir,
    pattern = "\\.meta\\.rds$",
    full.names = TRUE
  )
  rows <- lapply(meta_files, function(f) {
    tryCatch(readRDS(f), error = function(e) NULL)
  })
  rows <- rows[vapply(rows, is.data.frame, logical(1))]
  manifest <- if (length(rows) > 0L) {
    all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
    rows <- lapply(rows, function(x) {
      miss <- setdiff(all_names, names(x))
      for (nm in miss) x[[nm]] <- NA
      x[, all_names, drop = FALSE]
    })
    do.call(rbind, rows)
  } else {
    data.frame(
      TaskID = character(), Group = character(), StyleFamily = character(), Target = character(),
      TaskRDS = character(), PlotDataRDS = character(), MetaRDS = character(), CodeVersion = character(),
      CreatedAt = character(), stringsAsFactors = FALSE
    )
  }
  manifest <- manifest[!duplicated(manifest$TaskID, fromLast = TRUE), , drop = FALSE]
  manifest <- manifest[order(manifest$Group, manifest$Target), , drop = FALSE]
  utils::write.csv(
    manifest,
    file.path(plot_queue_dir, "plot_queue_manifest.csv"),
    row.names = FALSE
  )
  message("    Plot queue manifest: ", nrow(manifest), " task(s)")
  invisible(manifest)
}

checkpointed_lapply <- function(
    X,
    FUN,
    checkpoint_dir,
    label,
    workers = gsea_n_cores,
    chunk_size = gsea_chunk_size,
    ...) {
  X <- as.integer(X)
  if (length(X) == 0L) return(list())
  if (anyNA(X) || anyDuplicated(X)) {
    stop("[", label, "] checkpointed_lapply requires unique integer iteration IDs.")
  }

  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  chunk_size <- max(1L, as.integer(chunk_size))
  workers <- max(1L, min(as.integer(workers), length(X)))
  split_index <- split(seq_along(X), ceiling(seq_along(X) / chunk_size))
  collected <- vector("list", length(X))
  names(collected) <- as.character(X)

  for (chunk_no in seq_along(split_index)) {
    pos <- split_index[[chunk_no]]
    ids <- X[pos]
    checkpoint_file <- file.path(
      checkpoint_dir,
      sprintf("chunk_%05d_%05d.rds", min(ids), max(ids))
    )

    cached <- NULL
    if (file.exists(checkpoint_file) && !force_recompute) {
      cached <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
      valid_cached <- is.list(cached) &&
        identical(as.integer(cached$iterations), as.integer(ids)) &&
        is.list(cached$results) &&
        length(cached$results) == length(ids) &&
        identical(as.character(cached$CodeVersion), code_version)
      if (!valid_cached) {
        warning("[", label, "] damaged/stale checkpoint will be rebuilt: ",
          basename(checkpoint_file), call. = FALSE)
        unlink(checkpoint_file)
        cached <- NULL
      }
    }

    if (is.null(cached)) {
      worker_fun <- function(i) {
        tryCatch(
          FUN(i, ...),
          error = function(e) {
            structure(
              conditionMessage(e),
              class = c("checkpoint_worker_error", "character")
            )
          }
        )
      }

      chunk_results <- if (workers <= 1L || .Platform$OS.type == "windows") {
        lapply(ids, worker_fun)
      } else {
        parallel::mclapply(
          ids,
          worker_fun,
          mc.cores = min(workers, length(ids)),
          mc.set.seed = FALSE,
          mc.preschedule = TRUE,
          mc.cleanup = TRUE
        )
      }

      failed <- vapply(
        chunk_results,
        inherits,
        logical(1),
        what = "checkpoint_worker_error"
      )
      if (any(failed)) {
        stop(
          "[", label, "] iteration chunk ",
          min(ids), "-", max(ids),
          " failed and was not checkpointed. First error: ",
          as.character(chunk_results[[which(failed)[1L]]]),
          call. = FALSE
        )
      }

      cached <- list(
        CodeVersion = code_version,
        iterations = ids,
        results = chunk_results,
        completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
      atomic_save_rds(cached, checkpoint_file)
      message(sprintf(
        "    [%s] checkpoint %d/%d saved (iterations %d-%d)",
        label, chunk_no, length(split_index), min(ids), max(ids)
      ))
    } else {
      message(sprintf(
        "    [%s] checkpoint %d/%d reused (iterations %d-%d)",
        label, chunk_no, length(split_index), min(ids), max(ids)
      ))
    }

    collected[as.character(ids)] <- cached$results
    rm(cached)
    invisible(gc())
  }

  unname(collected[as.character(X)])
}


# ── CSV cache validator: skip only if file is readable, has expected columns,
#   contains the expected number of iterations, and every iteration has
#   the expected ID terms (looplook + Background for main GSEA,
#   only_looplook + intersection + only_ChIPseeker for Unique GSEA) ─────────
validate_csv_cache <- function(csv_path, expected_cols = c("ID", "NES", "pvalue", "Iteration", "Mode"), expected_iterations = NULL, expected_terms = NULL) {
  if (!file.exists(csv_path)) return(FALSE)
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(FALSE)
  if (!all(expected_cols %in% colnames(df))) return(FALSE)
  if (!is.null(expected_iterations)) {
    n_iter <- dplyr::n_distinct(df$Iteration)
    if (n_iter < expected_iterations) return(FALSE)
    # Check that iterations are strictly 1:expected_iterations
    iter_nums <- sort(unique(df$Iteration))
    if (!identical(iter_nums, seq_len(expected_iterations))) return(FALSE)
    # Per-iteration term completeness
    if (is.null(expected_terms)) {
      if (grepl("^(uniq_|expanded_)", basename(csv_path))) {
        # Unique analyses may legitimately contain only two eligible categories.
        # Require every iteration to contain the complete term set represented
        # by this signed cache, rather than hard-coding three categories.
        expected_terms <- sort(unique(as.character(df$ID)))
      } else {
        expected_terms <- if (any(grepl("^ctrl_", unique(df$Mode))))
          c("ChIPseeker", "Background") else c("looplook", "Background")
      }
    }
    iter_qc <- df %>%
      dplyr::filter(is.finite(NES)) %>%
      dplyr::distinct(Iteration, ID) %>%
      dplyr::group_by(Iteration) %>%
      dplyr::summarise(Complete = all(expected_terms %in% as.character(ID)), .groups = "drop")
    n_complete <- sum(iter_qc$Complete, na.rm = TRUE)
    if (n_complete < expected_iterations) return(FALSE)
  }
  return(TRUE)
}

# ── 1. 全局数据 ──────────────────────────────────────────────────────────────
message(">>> Loading global data...")
message("    Data base:  ", cfg$data_base)
message("    RData base: ", cfg$rdata_base)
message("    Output:     ", out_dir)

stopifnot(file.exists(diff_path), file.exists(expr_path), file.exists(target_file))
main_rdata_path <- file.path(cfg$rdata_base, cfg$tss_subdirs[["5kb"]])
stopifnot(file.exists(main_rdata_path))

# Load RData into isolated environment and verify required objects
main_env <- new.env(parent = emptyenv())
loaded_objects <- load(main_rdata_path, envir = main_env, verbose = FALSE)
required_objects <- c("res", "res2", "refined_res", "refined_res2", "cr", "cr2", "cr_only", "cr2_only")
missing_objects <- setdiff(required_objects, loaded_objects)
if (length(missing_objects) > 0L) {
  stop("Main RData is missing objects: ", paste(missing_objects, collapse = ", "))
}
# Copy to current environment
for (obj in required_objects) assign(obj, main_env[[obj]])
rm(main_env, loaded_objects, missing_objects)

# RData contains: res, res2, refined_res, refined_res2, cr, cr2, cr_only, cr2_only
#   res/res2           = annotation (hop0/hop1)
#   refined_res/2      = expression-refined (hop0/hop1)
#   cr/cr2             = refined → chromatin-refined (hop0/hop1)
#   cr_only/cr2_only   = annotation → chromatin-refined (hop0/hop1, no expression reclassification)

diff_df <- looplook:::read_robust_general(diff_path, header = TRUE, row_name = 1, desc = "Diff", min_cols = 1)
diff_df <- normalise_rownames(diff_df, "Differential expression matrix")
glist <- diff_df[["log2FoldChange"]]
names(glist) <- rownames(diff_df)
glist <- glist[!is.na(glist) & is.finite(glist)]
if (anyDuplicated(names(glist))) {
  dup <- unique(names(glist)[duplicated(names(glist))])
  stop("Duplicated gene names in ranked list: ", paste(head(dup, 20), collapse = ", "))
}
if (any(duplicated(glist))) { set.seed(42); glist <- glist + runif(length(glist), 0, 1e-6) }
glist <- sort(glist, decreasing = TRUE)
all_genes <- unique(names(glist))
# GO universe: all expressed/ranked genes converted to Entrez
universe_entrez <- tryCatch({
  e <- AnnotationDbi::select(org.Hs.eg.db, keys = all_genes,
    columns = "ENTREZID", keytype = "SYMBOL")
  unique(na.omit(e$ENTREZID))
}, error = function(cond) { warning("universe_entrez failed: ", cond$message); character(0) })
message("    Universe (all_genes -> entrez): ", length(universe_entrez))
if (length(universe_entrez) < 10) {
  stop("GO universe mapping from ranked genes to Entrez IDs yielded <10 IDs (", length(universe_entrez), "). ",
    "Cannot proceed with GO enrichment. Check gene symbol format in differential file.")
}
gouniv <- universe_entrez

stopifnot("log2FoldChange" %in% colnames(diff_df))
stopifnot(anyDuplicated(names(glist)) == 0)

# ── Early gene_stat and effect-size helpers (used by Unique/Hop/Signal Density) ──

# Point estimate only (no bootstrap CI) — used by expression-matched sampling
rank_biserial_point <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  2 * unname(wt$statistic) / (length(x) * length(y)) - 1
}

rank_biserial <- function(x, y, R = bootstrap_R, seed = bootstrap_seed) {
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

# Build gene_stat early so Unique/Hop effect-size code can access it
diff_df_for_stat <- diff_df
diff_df_for_stat$gene <- rownames(diff_df_for_stat)
diff_df_for_stat <- diff_df_for_stat[!is.na(diff_df_for_stat$log2FoldChange) & is.finite(diff_df_for_stat$log2FoldChange), ]
if (!"pvalue" %in% colnames(diff_df_for_stat)) diff_df_for_stat$pvalue <- 1
diff_df_for_stat$pvalue[is.na(diff_df_for_stat$pvalue)] <- 1
diff_df_for_stat$signal_lfc <- -diff_df_for_stat$log2FoldChange
neglogp <- pmin(-log10(pmax(diff_df_for_stat$pvalue, 1e-300)), 50)
diff_df_for_stat$signal_pv  <- -diff_df_for_stat$log2FoldChange * neglogp
gene_stat <- diff_df_for_stat
dup_genes <- unique(gene_stat$gene[duplicated(gene_stat$gene)])
if (length(dup_genes) > 0L) stop("Duplicated gene identifiers in DE input: ", paste(head(dup_genes, 20), collapse = ", "))

evidence_palette <- c(
  "local_promoter_overlap"                    = "#2166AC",
  "distal_promoter"                           = "#4393C3",
  "gene_body_context"                         = "#92C5DE",
  "distal_gene_body_context"                  = "#67A9CF",
  "local_enhancer_candidate"                  = "#1B9E77",
  "distal_enhancer_candidate"                 = "#66C2A5",
  "expanded_promoter_loop"                    = "#D1E5F0",
  "expanded_gene_body_context"                = "#A6D96A",
  "expanded_enhancer_candidate"               = "#B8E186",
  "positional_candidate"                      = "#FFD92F",
  "expanded_positional_candidate"             = "#FEE08B",
  "expanded_anchor"                           = "#F46D43",
  "linear_annotation"                         = "#F4A582",
  "linear_fallback"                           = "#CA0020",
  # Legacy labels retained for older benchmark objects
  "basic_gene_body_retained_after_chromatin"  = "#80B1D3",
  "positional_candidate_after_chromatin"      = "#E6AB02",
  "local_gene_body"                           = "#B2182B",
  "linear_nearest"                            = "#67001F",
  "local_promoter"                            = "#FDDBC7",
  "other"                                     = "grey60",
  "none"                                      = "grey80")

normalize_evidence <- function(x) {
  x <- as.character(x); x[is.na(x) | trimws(x) == ""] <- "none"
  x[!x %in% names(evidence_palette)] <- "other"
  factor(x, levels = names(evidence_palette))
}

get_target_gene_column <- function(md) {
  base_col <- if (identical(md$map, "promoter")) {
    "Regulated_promoter_genes"
  } else {
    "Assigned_Target_Genes"
  }
  if (isTRUE(md$fill)) paste0(base_col, "_Filled") else base_col
}

# In looplook, neighbor_hop=1 adds path_length > 1 genes to
# Expanded_Target_Genes / target_gene_links. It does not alter the primary
# Assigned_Target_Genes columns. For this benchmark:
#   hop0 modes = primary targets only
#   hop1 modes = primary targets + eligible expanded targets
is_hop1_mode <- function(md) {
  !is.null(md$name) &&
    length(md$name) == 1L &&
    !is.na(md$name) &&
    grepl("_hop1$", md$name)
}

clean_mode_genes <- function(x) {
  genes <- looplook:::clean_gene_names(x, "[;,]")
  genes <- toupper(trimws(as.character(genes)))
  unique(
    genes[
      !is.na(genes) &
        nzchar(genes) &
        genes != "NA"
    ]
  )
}

# Formal benchmark objects must expose the target membership columns expected
# by each mode. Interface/schema failures are structural errors, not biological
# empty sets, and therefore stop before any expensive computation.
validate_mode_schema <- function(md) {
  ta <- md$ann$target_annotation
  if (is.null(ta) || nrow(ta) == 0L) {
    stop("[", md$name, "] target_annotation is missing or empty.", call. = FALSE)
  }

  required_ta <- c("input_id", get_target_gene_column(md))
  missing_ta <- setdiff(required_ta, colnames(ta))
  if (length(missing_ta) > 0L) {
    stop(
      "[", md$name, "] target_annotation missing required column(s): ",
      paste(missing_ta, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(ta$input_id) || any(!nzchar(trimws(as.character(ta$input_id))))) {
    stop("[", md$name, "] target_annotation contains blank input_id values.", call. = FALSE)
  }

  tgl <- md$ann$target_gene_links
  required_tgl <- c("input_id", "gene", "gene_role", "source", "evidence")
  if (is_hop1_mode(md)) {
    required_tgl <- c(required_tgl, "path_length", "in_expanded_target")
  }
  if (is.null(tgl)) {
    stop("[", md$name, "] target_gene_links is missing.", call. = FALSE)
  }
  missing_tgl <- setdiff(required_tgl, colnames(tgl))
  if (length(missing_tgl) > 0L) {
    stop(
      "[", md$name, "] target_gene_links missing required column(s): ",
      paste(missing_tgl, collapse = ", "),
      call. = FALSE
    )
  }

  data.frame(
    Mode = md$name,
    TargetRows = nrow(ta),
    LinkRows = nrow(tgl),
    TargetColumn = get_target_gene_column(md),
    Hop1 = is_hop1_mode(md),
    Status = "ok",
    stringsAsFactors = FALSE
  )
}

# Primary target pairs are kept separate from expanded links so that the
# biological contribution of neighbor_hop can be evaluated explicitly.
get_primary_peak_gene_pairs <- function(md, row_index = NULL) {
  empty <- data.frame(
    Input_ID = character(),
    Gene = character(),
    stringsAsFactors = FALSE
  )

  bed_info <- md$ann$target_annotation
  if (is.null(bed_info) ||
      nrow(bed_info) == 0L ||
      !"input_id" %in% colnames(bed_info)) {
    return(empty)
  }

  desired_col <- get_target_gene_column(md)
  if (!desired_col %in% colnames(bed_info)) {
    warning(
      "[", md$name,
      "] missing primary target column: ",
      desired_col,
      call. = FALSE
    )
    return(empty)
  }

  if (is.null(row_index)) {
    row_index <- seq_len(nrow(bed_info))
  }

  row_index <- unique(as.integer(row_index))
  row_index <- row_index[
    is.finite(row_index) &
      row_index >= 1L &
      row_index <= nrow(bed_info)
  ]

  if (length(row_index) == 0L) {
    return(empty)
  }

  rows <- lapply(row_index, function(i) {
    genes <- clean_mode_genes(bed_info[[desired_col]][i])
    if (length(genes) == 0L) {
      return(NULL)
    }

    data.frame(
      Input_ID = as.character(bed_info$input_id[i]),
      Gene = genes,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows) %>%
    dplyr::filter(
      !is.na(Input_ID),
      nzchar(Input_ID),
      !is.na(Gene),
      nzchar(Gene),
      Gene != "NA"
    ) %>%
    dplyr::distinct(Input_ID, Gene)
}

get_primary_mode_genes <- function(md, row_index = NULL) {
  pairs <- get_primary_peak_gene_pairs(md, row_index = row_index)
  if (nrow(pairs) == 0L) {
    return(character())
  }
  sort(unique(pairs$Gene))
}

get_expanded_peak_gene_pairs <- function(md, input_ids = NULL) {
  empty <- data.frame(
    Input_ID = character(),
    Gene = character(),
    stringsAsFactors = FALSE
  )

  if (!is_hop1_mode(md)) {
    return(empty)
  }

  tgl <- md$ann$target_gene_links
  required_cols <- c(
    "input_id",
    "gene",
    "gene_role",
    "source",
    "path_length",
    "in_expanded_target"
  )

  if (is.null(tgl) || nrow(tgl) == 0L) {
    return(empty)
  }

  missing_cols <- setdiff(required_cols, colnames(tgl))
  if (length(missing_cols) > 0L) {
    stop(
      "[", md$name,
      "] target_gene_links lacks expanded-target provenance column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  keep <- as.character(tgl$source) == "loop_anchor" &
    is.finite(tgl$path_length) &
    tgl$path_length > 1L &
    tgl$in_expanded_target %in% TRUE

  if ("anchor_role" %in% colnames(tgl)) {
    keep <- keep &
      as.character(tgl$anchor_role) == "expanded_anchor"
  }

  if (identical(md$map, "promoter")) {
    keep <- keep &
      as.character(tgl$gene_role) == "promoter"
  }

  if (!is.null(input_ids)) {
    input_ids <- unique(as.character(input_ids))
    keep <- keep &
      as.character(tgl$input_id) %in% input_ids
  }

  out <- data.frame(
    Input_ID = as.character(tgl$input_id[keep]),
    Gene = toupper(trimws(as.character(tgl$gene[keep]))),
    stringsAsFactors = FALSE
  )

  out <- out[
    !is.na(out$Input_ID) &
      nzchar(out$Input_ID) &
      !is.na(out$Gene) &
      nzchar(out$Gene) &
      out$Gene != "NA",
    ,
    drop = FALSE
  ]

  unique(out)
}

# Expanded-only is defined at the mode level as eligible expanded genes that
# are absent from that mode's primary target set. A gene that is expanded for
# one peak but primary for another is therefore not counted as expanded-only.
get_expanded_only_peak_gene_pairs <- function(md, row_index = NULL) {
  empty <- data.frame(
    Input_ID = character(),
    Gene = character(),
    stringsAsFactors = FALSE
  )

  if (!is_hop1_mode(md)) {
    return(empty)
  }

  bed_info <- md$ann$target_annotation
  if (is.null(bed_info) ||
      nrow(bed_info) == 0L ||
      !"input_id" %in% colnames(bed_info)) {
    return(empty)
  }

  if (is.null(row_index)) {
    row_index <- seq_len(nrow(bed_info))
  }

  row_index <- unique(as.integer(row_index))
  row_index <- row_index[
    is.finite(row_index) &
      row_index >= 1L &
      row_index <= nrow(bed_info)
  ]

  if (length(row_index) == 0L) {
    return(empty)
  }

  input_ids <- as.character(bed_info$input_id[row_index])
  primary_pairs <- get_primary_peak_gene_pairs(md, row_index = row_index)
  expanded_pairs <- get_expanded_peak_gene_pairs(md, input_ids = input_ids)

  if (nrow(expanded_pairs) == 0L) {
    return(empty)
  }

  primary_genes <- unique(primary_pairs$Gene)
  expanded_only_genes <- setdiff(
    unique(expanded_pairs$Gene),
    primary_genes
  )

  expanded_pairs %>%
    dplyr::filter(Gene %in% expanded_only_genes) %>%
    dplyr::distinct(Input_ID, Gene)
}

get_expanded_only_genes <- function(md, row_index = NULL) {
  pairs <- get_expanded_only_peak_gene_pairs(md, row_index = row_index)
  if (nrow(pairs) == 0L) {
    return(character())
  }
  sort(unique(pairs$Gene))
}

get_mode_peak_gene_pairs <- function(md, row_index = NULL) {
  empty <- data.frame(
    Input_ID = character(),
    Gene = character(),
    stringsAsFactors = FALSE
  )

  bed_info <- md$ann$target_annotation
  if (is.null(bed_info) ||
      nrow(bed_info) == 0L ||
      !"input_id" %in% colnames(bed_info)) {
    return(empty)
  }

  desired_col <- get_target_gene_column(md)
  if (!desired_col %in% colnames(bed_info)) {
    warning(
      "[", md$name,
      "] missing target_annotation column: ",
      desired_col,
      call. = FALSE
    )
    return(empty)
  }

  if (is.null(row_index)) {
    row_index <- seq_len(nrow(bed_info))
  }

  row_index <- unique(as.integer(row_index))
  row_index <- row_index[
    is.finite(row_index) &
      row_index >= 1L &
      row_index <= nrow(bed_info)
  ]

  if (length(row_index) == 0L) {
    return(empty)
  }

  primary_pairs <- lapply(row_index, function(i) {
    genes <- clean_mode_genes(bed_info[[desired_col]][i])
    if (length(genes) == 0L) {
      return(NULL)
    }

    data.frame(
      Input_ID = as.character(bed_info$input_id[i]),
      Gene = genes,
      stringsAsFactors = FALSE
    )
  })

  primary_pairs <- dplyr::bind_rows(primary_pairs)

  expanded_pairs <- get_expanded_peak_gene_pairs(
    md,
    input_ids = bed_info$input_id[row_index]
  )

  dplyr::bind_rows(
    primary_pairs,
    expanded_pairs
  ) %>%
    dplyr::distinct(Input_ID, Gene)
}

get_mode_genes <- function(md, row_index = NULL) {
  pairs <- get_mode_peak_gene_pairs(
    md,
    row_index = row_index
  )

  if (nrow(pairs) == 0L) {
    return(character())
  }

  sort(unique(pairs$Gene))
}

build_mode_gene_evidence <- function(md, nm, valid_genes = NULL) {
  bed_info <- md$ann$target_annotation; tgl <- md$ann$target_gene_links
  if (is.null(bed_info) || is.null(tgl) || nrow(bed_info) == 0 || nrow(tgl) == 0)
    return(data.frame(Mode=character(),Input_ID=character(),Gene=character(),Evidence=character(),Source=character(),stringsAsFactors=FALSE))
  if (!"input_id" %in% colnames(bed_info))
    return(data.frame(Mode=character(),Input_ID=character(),Gene=character(),Evidence=character(),Source=character(),stringsAsFactors=FALSE))
  if (!all(c("input_id","gene","evidence","source","gene_role") %in% colnames(tgl)))
    return(data.frame(Mode=character(),Input_ID=character(),Gene=character(),Evidence=character(),Source=character(),stringsAsFactors=FALSE))

  peak_gene_pairs <- get_mode_peak_gene_pairs(md)
  if (nrow(peak_gene_pairs) == 0) return(data.frame(Mode=character(),Input_ID=character(),Gene=character(),Evidence=character(),Source=character(),stringsAsFactors=FALSE))

  tgl_clean <- tgl %>% transmute(Input_ID = as.character(input_id), Gene = as.character(gene),
    Evidence = as.character(normalize_evidence(evidence)), Source = as.character(source),
    GeneRole = as.character(gene_role)) %>%
    filter(!is.na(Input_ID), Input_ID != "", !is.na(Gene), Gene != "", Gene != "NA") %>%
    mutate(Gene = toupper(trimws(Gene)))

  evidence_priority <- c(
    "local_promoter_overlap"                    = 1,
    "distal_promoter"                           = 2,
    "gene_body_context"                         = 3,
    "distal_gene_body_context"                  = 4,
    "local_enhancer_candidate"                  = 5,
    "distal_enhancer_candidate"                 = 6,
    "expanded_promoter_loop"                    = 7,
    "expanded_gene_body_context"                = 8,
    "expanded_enhancer_candidate"               = 9,
    "positional_candidate"                      = 10,
    "expanded_positional_candidate"             = 11,
    "expanded_anchor"                           = 12,
    "linear_annotation"                         = 13,
    "linear_fallback"                           = 14,
    # Legacy labels retained for older benchmark objects
    "basic_gene_body_retained_after_chromatin"  = 15,
    "positional_candidate_after_chromatin"      = 16,
    "local_gene_body"                           = 17,
    "linear_nearest"                            = 18,
    "local_promoter"                            = 19,
    "other"                                     = 20,
    "none"                                      = 21)

  out <- inner_join(tgl_clean, peak_gene_pairs, by = c("Input_ID", "Gene"))
  if (!is.null(valid_genes)) out <- out %>% filter(Gene %in% unique(toupper(valid_genes)))
  if (md$map == "promoter" && !md$fill) out <- out %>% filter(GeneRole == "promoter")
  if (!md$fill) out <- out %>% filter(Source == "loop_anchor")
  # Select primary (highest-priority) evidence per Input_ID + Gene
  # Use base R order + !duplicated to avoid dplyr/S4Vectors dispatch conflict
  out$EvPriority <- evidence_priority[as.character(out$Evidence)]
  out <- out[order(out$Input_ID, out$Gene, out$EvPriority), ]
  out <- out[!duplicated(paste(out$Input_ID, out$Gene)), ]
  out$EvPriority <- NULL
  out %>% transmute(Mode = nm, Input_ID = Input_ID, Gene = Gene, Evidence = Evidence, Source = Source)
}


# Direction note: glist sorted decreasing (upregulated first), so negative NES = enrichment at bottom = downregulated.
# signal_lfc = -log2FoldChange gives positive values for downregulation. Both point to same biology.
# 第二套 ranked list: sign(LFC) × -log10(pvalue)
if ("pvalue" %in% colnames(diff_df)) {
  pvals <- diff_df[["pvalue"]]
  names(pvals) <- toupper(rownames(diff_df))
  pvals[is.na(pvals)] <- 1
  lfc_vals <- diff_df[["log2FoldChange"]]
  names(lfc_vals) <- toupper(rownames(diff_df))
  glist_pv <- sign(lfc_vals) * pmin(-log10(pmax(pvals, 1e-300)), 50)
  glist_pv <- glist_pv[!is.na(glist_pv) & is.finite(glist_pv)]
  if (any(duplicated(glist_pv))) { set.seed(42); glist_pv <- glist_pv + runif(length(glist_pv), 0, 1e-6) }
  glist_pv <- sort(glist_pv, decreasing = TRUE)
  use_glist_pv <- TRUE
  message("    glist_pv (sign(LFC)*-log10(pvalue)) prepared: ", length(glist_pv), " genes")
} else {
  use_glist_pv <- FALSE
  message("    No pvalue column found, skipping glist_pv")
}

# ── 2. ChIPseeker ──────────────────────────────────────────────────────────
message(">>> ChIPseeker annotation...")
peak_anno <- annotatePeak(target_file, tssRegion = c(-5000, 5000),
  # +/-5kb promoter definition: standard for ChIP-seq, balances sensitivity (captures
  # proximal regulatory elements) vs specificity (avoids misclassifying distal enhancers
  # as promoters). For BRD4, which binds both promoters and super-enhancers, this window
  # provides a reasonable promoter/distal boundary.
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db")
peak_anno_df <- as.data.frame(peak_anno)
chipseeker_genes <- unique(toupper(na.omit(peak_anno_df$SYMBOL)))
message("    ChIPseeker all: ", length(chipseeker_genes), " genes")
is_promoter_peak <- grepl("Promoter", peak_anno_df$annotation)
chipseeker_promoter_genes <- unique(toupper(na.omit(peak_anno_df$SYMBOL[is_promoter_peak])))
message("    ChIPseeker promoter: ", length(chipseeker_promoter_genes), " genes")

get_chip_reference <- function(md) {
  # Unified ChIPseeker-all reference for all modes.
  # looplook-promoter filters by target-side promoter role,
  # not peak-side annotation, so the external comparator
  # must always be the complete ChIPseeker nearest-gene set.
  chipseeker_genes
}

is_distal <- grepl("Distal Intergenic", peak_anno_df$annotation)
chipseeker_distal_genes <- unique(toupper(na.omit(peak_anno_df$SYMBOL[is_distal])))
message("    ChIPseeker Distal Intergenic: ", length(chipseeker_distal_genes), " genes")

chipseeker_promoter_only <- chipseeker_promoter_genes  # alias, same as above
message("    ChIPseeker Promoter: ", length(chipseeker_promoter_only), " genes")

is_genic <- grepl("Exon|Intron|UTR|Downstream", peak_anno_df$annotation) & !is_promoter_peak
chipseeker_genic_genes <- unique(toupper(na.omit(peak_anno_df$SYMBOL[is_genic])))
chipseeker_non_promoter_genes <- unique(c(chipseeker_distal_genes, chipseeker_genic_genes))
message("    ChIPseeker Non-Promoter (distal+genic): ", length(chipseeker_non_promoter_genes), " genes")
message("    ChIPseeker Genic (exon/intron/UTR): ", length(chipseeker_genic_genes), " genes")

# ── Super-enhancer / Typical enhancer classification (optional) ──────────────
# Each BED independently loaded; SE/TE only added if corresponding file exists
load_enhancer_bed <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  tryCatch({
    df <- read.table(path, header = FALSE, sep = "\t", comment.char = "#",
      stringsAsFactors = FALSE)
    if (ncol(df) < 3L || nrow(df) == 0L) {
      stop("BED must contain at least three columns and one region.")
    }
    chr <- trimws(as.character(df[[1L]]))
    start0 <- suppressWarnings(as.integer(df[[2L]]))
    end0 <- suppressWarnings(as.integer(df[[3L]]))
    invalid <- is.na(chr) | !nzchar(chr) | is.na(start0) | is.na(end0) |
      start0 < 0L | end0 <= start0
    if (any(invalid)) {
      stop("BED contains ", sum(invalid), " malformed interval(s).")
    }
    unique(GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = start0 + 1L, end = end0)
    ))
  }, error = function(e) {
    stop("Failed to read enhancer BED ", path, ": ", conditionMessage(e), call. = FALSE)
  })
}

use_se_stratification <- FALSE
se_gr <- NULL; te_gr <- NULL
if (!is.null(cfg$se_bed_file)) se_gr <- load_enhancer_bed(file.path(cfg$data_base, cfg$se_bed_file))
if (!is.null(cfg$te_bed_file)) te_gr <- load_enhancer_bed(file.path(cfg$data_base, cfg$te_bed_file))

# ROSE BED classifications are treated as the authoritative SE/TE calls.
# Do not re-cut them using ChIPseeker promoter/distal labels: ROSE stitching can
# legitimately span promoter-proximal or genic sequence. Genome build and
# chromosome naming concordance are checked through overlap/QC instead.
is_se_peak <- rep(FALSE, nrow(peak_anno_df))
is_te_peak <- rep(FALSE, nrow(peak_anno_df))
if (!is.null(se_gr) && length(se_gr) > 0) {
  peak_gr <- GenomicRanges::GRanges(seqnames = peak_anno_df$seqnames,
    ranges = IRanges::IRanges(start = peak_anno_df$start, end = peak_anno_df$end))
  if (length(intersect(GenomeInfoDb::seqlevels(peak_gr), GenomeInfoDb::seqlevels(se_gr))) == 0L) {
    stop("SE BED and BRD4 peaks share no chromosome names; check genome build/chr prefix.")
  }
  se_overlap <- GenomicRanges::findOverlaps(peak_gr, se_gr)
  if (length(se_overlap) == 0L) {
    warning("ROSE SE BED has zero overlap with the BRD4 peak set.", call. = FALSE)
  }
  is_se_peak[unique(queryHits(se_overlap))] <- TRUE
  use_se_stratification <- TRUE
  message(sprintf("    Super-Enhancer peaks: %d", sum(is_se_peak)))
}
if (!is.null(te_gr) && length(te_gr) > 0) {
  if (!exists("peak_gr")) peak_gr <- GenomicRanges::GRanges(seqnames = peak_anno_df$seqnames,
    ranges = IRanges::IRanges(start = peak_anno_df$start, end = peak_anno_df$end))
  if (length(intersect(GenomeInfoDb::seqlevels(peak_gr), GenomeInfoDb::seqlevels(te_gr))) == 0L) {
    stop("TE BED and BRD4 peaks share no chromosome names; check genome build/chr prefix.")
  }
  te_overlap <- GenomicRanges::findOverlaps(peak_gr, te_gr)
  if (length(te_overlap) == 0L) {
    warning("ROSE TE BED has zero overlap with the BRD4 peak set.", call. = FALSE)
  }
  is_te_peak[unique(queryHits(te_overlap))] <- TRUE
  # TE excludes SE regions (user-defined, not our heuristic)
  is_te_peak <- is_te_peak & !is_se_peak
  use_se_stratification <- TRUE
  message(sprintf("    Typical Enhancer peaks: %d", sum(is_te_peak)))
}

gc()

# ── Peak input concordance QC ────────────────────────────────────────────────
make_peak_id <- function(df) {
  paste0(as.character(df$seqnames), ":", as.integer(df$start), "-", as.integer(df$end))
}
target_peak_ids <- unique(make_peak_id(peak_anno_df))
if (!is.null(res$target_annotation)) {
  rdata_peak_ids <- unique(make_peak_id(res$target_annotation))
  n_shared <- length(intersect(target_peak_ids, rdata_peak_ids))
  n_union  <- length(union(target_peak_ids, rdata_peak_ids))
  peak_jaccard <- n_shared / max(1L, n_union)
  peak_recall_target <- n_shared / max(1L, length(target_peak_ids))
  peak_recall_rdata  <- n_shared / max(1L, length(rdata_peak_ids))
  message(sprintf("    Peak concordance: Target=%d RData=%d Shared=%d Jaccard=%.4f",
    length(target_peak_ids), length(rdata_peak_ids), n_shared, peak_jaccard))
  peak_qc <- data.frame(
    Target_N = length(target_peak_ids), RData_N = length(rdata_peak_ids),
    Shared_N = n_shared, Jaccard = peak_jaccard,
    TargetRecall = peak_recall_target, RDataRecall = peak_recall_rdata,
    stringsAsFactors = FALSE)
  write.csv(peak_qc, file.path(base_out_dir, "Peak_Input_Concordance.csv"), row.names = FALSE)
  if (peak_recall_target < 0.99 || peak_recall_rdata < 0.99) {
    stop("RData target peaks do not match current narrowPeak (target recall=",
      round(peak_recall_target, 4), ", RData recall=", round(peak_recall_rdata, 4),
      "). Check that both files use the same peak set.")
  }
} else {
  warning("RData target_annotation is NULL; cannot verify peak concordance.")
}

if (force_recompute) {
  message(">>> Force recompute: deleting old caches and all GSEA temp files...")
  for (sz in sample_sizes) {
    unlink(file.path(out_dir, paste0("size", sz), "tmp_sim"), recursive = TRUE)
    unlink(file.path(out_dir, paste0("size", sz), "tmp_unique"), recursive = TRUE)
  }
  unlink(file.path(base_out_dir, "go_enrichment"), recursive = TRUE)
  unlink(file.path(base_out_dir, "LFC_pv"), recursive = TRUE)
  unlink(file.path(out_dir, "mode_genes.rds"))
  unlink(file.path(out_dir, "modes.rds"))
  unlink(file.path(out_dir, "mode_cache.rds"))
  unlink(file.path(out_dir, "mode_ranking"), recursive = TRUE)
  unlink(file.path(base_out_dir, "evidence_composition"), recursive = TRUE)
  unlink(file.path(base_out_dir, "cross_tss_summary"), recursive = TRUE)
  unlink(file.path(base_out_dir, "cross_size_summary"), recursive = TRUE)
  unlink(file.path(base_out_dir, "hop_interpretability"), recursive = TRUE)
  unlink(file.path(base_out_dir, "background_sensitivity"), recursive = TRUE)
  unlink(file.path(base_out_dir, "convergence_plots"), recursive = TRUE)
  unlink(file.path(base_out_dir, "expanded_only"), recursive = TRUE)
  # Downstream modules are also regenerated. Removing them prevents stale PDFs
  # or workbooks from being mistaken for outputs of the current run.
  downstream_dirs <- c(
    "distance_analysis", "peak_category_es", "case_study",
    "expression_matched_es", "summary", "mode_ranking",
    "signal_density", "LFC_pv", "go_enrichment",
    "cross_tss_summary", "cross_size_summary",
    "hop_interpretability", "background_sensitivity",
    "convergence_plots", "evidence_composition"
  )
  unlink(file.path(base_out_dir, downstream_dirs), recursive = TRUE)
  # Delete all status/provenance files
  unlink(file.path(base_out_dir, c(
    "Analysis_Run_Status.csv", "GSEA_Run_Status.csv", "GSEA_Success_Summary.csv",
    "GSEA_Iteration_Summary.csv", "GSEA_Mode_Status.csv",
    "GO_Run_Status.csv", "GO_Call_Status.csv", "GO_Status_Summary.csv",
    "Output_Completeness_Report.csv", "analysis_manifest.csv", "run_signature.rds")))
  message("    All caches cleared. Starting fresh run.")
}

# ── 3. 预计算靶基因集 (single entry point for mode cache) ─────────────────────
mode_cache_file <- file.path(out_dir, "mode_cache.rds")

# Lightweight mode spec for cache validation (no looplook annotation)
mode_spec_defs <- list(
  list(name = "anno_all_F",              ann_source = "res", map = "all",      fill = FALSE, near = FALSE),
  list(name = "anno_all_T",              ann_source = "res", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "refined_all_F",           ann_source = "refined_res", map = "all",      fill = FALSE, near = FALSE),
  list(name = "refined_all_T",           ann_source = "refined_res", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "anno_all_F_hop1",         ann_source = "res2", map = "all",      fill = FALSE, near = FALSE),
  list(name = "anno_all_T_hop1",         ann_source = "res2", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "refined_all_F_hop1",      ann_source = "refined_res2", map = "all",      fill = FALSE, near = FALSE),
  list(name = "refined_all_T_hop1",      ann_source = "refined_res2", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "anno_promoter_F",         ann_source = "res", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "anno_promoter_T",         ann_source = "res", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "refined_promoter_F",      ann_source = "refined_res", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "refined_promoter_T",      ann_source = "refined_res", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "anno_promoter_F_hop1",    ann_source = "res2", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "anno_promoter_T_hop1",    ann_source = "res2", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "refined_promoter_F_hop1", ann_source = "refined_res2", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "refined_promoter_T_hop1", ann_source = "refined_res2", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "chrom_all_F",             ann_source = "cr", map = "all",      fill = FALSE, near = FALSE),
  list(name = "chrom_all_T",             ann_source = "cr", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "chrom_all_F_hop1",        ann_source = "cr2", map = "all",      fill = FALSE, near = FALSE),
  list(name = "chrom_all_T_hop1",        ann_source = "cr2", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "chrom_promoter_F",        ann_source = "cr", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "chrom_promoter_T",        ann_source = "cr", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "chrom_promoter_F_hop1",   ann_source = "cr2", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "chrom_promoter_T_hop1",   ann_source = "cr2", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "chrom_only_all_F",        ann_source = "cr_only", map = "all",      fill = FALSE, near = FALSE),
  list(name = "chrom_only_all_T",        ann_source = "cr_only", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "chrom_only_all_F_hop1",   ann_source = "cr2_only", map = "all",      fill = FALSE, near = FALSE),
  list(name = "chrom_only_all_T_hop1",   ann_source = "cr2_only", map = "all",      fill = TRUE,  near = FALSE),
  list(name = "chrom_only_promoter_F",   ann_source = "cr_only", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "chrom_only_promoter_T",   ann_source = "cr_only", map = "promoter", fill = TRUE,  near = FALSE),
  list(name = "chrom_only_promoter_F_hop1", ann_source = "cr2_only", map = "promoter", fill = FALSE, near = FALSE),
  list(name = "chrom_only_promoter_T_hop1", ann_source = "cr2_only", map = "promoter", fill = TRUE,  near = FALSE))

n_modes <- 16L  # primary manuscript modes (hop0)

if (file.exists(mode_cache_file) && !force_recompute) {
  mode_cache <- readRDS(mode_cache_file)
  if (is.null(mode_cache$signature)) {
    stop("Unsigned mode cache detected. Delete ", mode_cache_file, " and re-run with force_recompute=TRUE.")
  }
  # Validate input data keys (not DiffMD5, which does not affect mode gene assignment)
  input_keys <- c("TargetMD5", "MainRDataMD5")
  old_input_sig <- mode_cache$signature[input_keys]
  new_input_sig <- preliminary_signature[input_keys]
  if (!identical(old_input_sig, new_input_sig)) {
    stop("Mode cache signature mismatch (input files changed). ",
      "Set cache_policy='rebuild' or delete ", mode_cache_file)
  }
  # Validate mode definition
  new_mode_def_md5 <- md5_object(mode_spec_defs)
  if (!identical(mode_cache$signature$ModeDefinitionMD5, new_mode_def_md5)) {
    stop("Mode definition changed since cache was built. ",
      "Set cache_policy='rebuild' or delete ", mode_cache_file)
  }
  # Validate looplook package version
  if (!identical(mode_cache$signature$LooplookVersion,
      preliminary_signature$LooplookVersion)) {
    stop("looplook package version changed. ",
      "Set cache_policy='rebuild' or delete ", mode_cache_file)
  }
  # Do not invalidate mode genes solely because plotting/logging/runtime
  # configuration changed in the same script. Validate an explicit mode-logic
  # version when present; legacy caches are accepted after the input,
  # package-version and mode-definition checks above.
  if (!is.null(mode_cache$signature$ModeLogicVersion) &&
      !identical(mode_cache$signature$ModeLogicVersion, mode_cache_logic_version)) {
    stop(
      "Mode gene-extraction logic changed. Set cache_policy='rebuild' or delete ",
      mode_cache_file
    )
  }
  modes <- Filter(function(x) !grepl("_hop1$", x$name), mode_cache$modes)
  mode_genes <- mode_cache$mode_genes[vapply(modes, `[[`, character(1), "name")]
  message("    Mode cache loaded and validated: ", length(mode_genes), " modes")
} else {
  message(">>> Building mode cache from scratch...")
  modes <- list(
    list(name = "anno_all_F",             ann = res,          map = "all",      fill = FALSE, near = FALSE),
    list(name = "anno_all_T",             ann = res,          map = "all",      fill = TRUE,  near = FALSE),
    list(name = "refined_all_F",          ann = refined_res,  map = "all",      fill = FALSE, near = FALSE),
    list(name = "refined_all_T",          ann = refined_res,  map = "all",      fill = TRUE,  near = FALSE),
    list(name = "anno_all_F_hop1",        ann = res2,         map = "all",      fill = FALSE, near = FALSE),
    list(name = "anno_all_T_hop1",        ann = res2,         map = "all",      fill = TRUE,  near = FALSE),
    list(name = "refined_all_F_hop1",     ann = refined_res2, map = "all",      fill = FALSE, near = FALSE),
    list(name = "refined_all_T_hop1",     ann = refined_res2, map = "all",      fill = TRUE,  near = FALSE),
    list(name = "anno_promoter_F",        ann = res,          map = "promoter", fill = FALSE, near = FALSE),
    list(name = "anno_promoter_T",        ann = res,          map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "refined_promoter_F",     ann = refined_res,  map = "promoter", fill = FALSE, near = FALSE),
    list(name = "refined_promoter_T",     ann = refined_res,  map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "anno_promoter_F_hop1",   ann = res2,         map = "promoter", fill = FALSE, near = FALSE),
    list(name = "anno_promoter_T_hop1",   ann = res2,         map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "refined_promoter_F_hop1",ann = refined_res2, map = "promoter", fill = FALSE, near = FALSE),
    list(name = "refined_promoter_T_hop1",ann = refined_res2, map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "chrom_all_F",            ann = cr,            map = "all",      fill = FALSE, near = FALSE),
    list(name = "chrom_all_T",            ann = cr,            map = "all",      fill = TRUE,  near = FALSE),
    list(name = "chrom_all_F_hop1",       ann = cr2,           map = "all",      fill = FALSE, near = FALSE),
    list(name = "chrom_all_T_hop1",       ann = cr2,           map = "all",      fill = TRUE,  near = FALSE),
    list(name = "chrom_promoter_F",       ann = cr,            map = "promoter", fill = FALSE, near = FALSE),
    list(name = "chrom_promoter_T",       ann = cr,            map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "chrom_promoter_F_hop1",  ann = cr2,           map = "promoter", fill = FALSE, near = FALSE),
    list(name = "chrom_promoter_T_hop1",  ann = cr2,           map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "chrom_only_all_F",       ann = cr_only,       map = "all",      fill = FALSE, near = FALSE),
    list(name = "chrom_only_all_T",       ann = cr_only,       map = "all",      fill = TRUE,  near = FALSE),
    list(name = "chrom_only_all_F_hop1",  ann = cr2_only,      map = "all",      fill = FALSE, near = FALSE),
    list(name = "chrom_only_all_T_hop1",  ann = cr2_only,      map = "all",      fill = TRUE,  near = FALSE),
    list(name = "chrom_only_promoter_F",  ann = cr_only,       map = "promoter", fill = FALSE, near = FALSE),
    list(name = "chrom_only_promoter_T",  ann = cr_only,       map = "promoter", fill = TRUE,  near = FALSE),
    list(name = "chrom_only_promoter_F_hop1", ann = cr2_only,  map = "promoter", fill = FALSE, near = FALSE),
    list(name = "chrom_only_promoter_T_hop1", ann = cr2_only,  map = "promoter", fill = TRUE,  near = FALSE)
  )
  modes <- Filter(function(x) !grepl("_hop1$", x$name), modes)
  mode_genes <- list()
  for (i in seq_along(modes)) {
    md <- modes[[i]]
    if (md$near) {
      genes <- if (md$map == "promoter") chipseeker_promoter_genes else chipseeker_genes
    } else {
      genes <- get_mode_genes(md)
    }
    mode_genes[[md$name]] <- genes
    message(sprintf("    [%2d/%d] %-30s -> %d genes", i, length(modes), md$name, length(genes)))
  }
  mode_cache_signature <- list(
    ModeLogicVersion = mode_cache_logic_version,
    DiffMD5 = preliminary_signature$DiffMD5,
    TargetMD5 = preliminary_signature$TargetMD5,
    MainRDataMD5 = preliminary_signature$MainRDataMD5,
    LooplookVersion = preliminary_signature$LooplookVersion,
    ModeDefinitionMD5 = md5_object(mode_spec_defs))
  saveRDS(list(modes = modes, mode_genes = mode_genes, signature = mode_cache_signature),
    mode_cache_file)
  message("    Mode cache built and saved: ", length(mode_genes), " modes")
}

# ── Formal mode-object schema preflight ──────────────────────────────────────
mode_schema_qc <- dplyr::bind_rows(lapply(modes, validate_mode_schema))
if (!identical(sort(names(mode_genes)), sort(vapply(modes, `[[`, character(1), "name")))) {
  stop("mode_genes names do not match the declared mode objects.", call. = FALSE)
}
utils::write.csv(
  mode_schema_qc,
  file.path(base_out_dir, "Mode_Object_Schema_QC.csv"),
  row.names = FALSE
)
message("    Mode schema preflight passed: ", nrow(mode_schema_qc), " modes")

benchmark_cache_signature_base <- list(
  CodeVersion = gsea_cache_logic_version,
  DiffMD5 = preliminary_signature$DiffMD5,
  ExprMD5 = preliminary_signature$ExprMD5,
  MetadataMD5 = preliminary_signature$MetadataMD5,
  TargetMD5 = preliminary_signature$TargetMD5,
  MainRDataMD5 = preliminary_signature$MainRDataMD5,
  LooplookVersion = preliminary_signature$LooplookVersion,
  ModeGenesMD5 = md5_object(mode_genes),
  ChIPReferenceMD5 = md5_object(chipseeker_genes),
  Iterations = n_iterations,
  GSEAMinSize = gsea_min_size,
  GSEAMaxSize = gsea_max_size,
  ControlDefinition = "not_looplook_and_not_ChIPseeker",
  SamplingRule = "matched_N_requested_or_80pct_not_below_minGSSize"
)

# ── Hop scope QC: hop1 must be primary + expanded, not a duplicate label ────
hop_mode_names <- grep("_hop1$", names(mode_genes), value = TRUE)
hop_scope_qc <- lapply(hop_mode_names, function(hop1_name) {
  hop0_name <- sub("_hop1$", "", hop1_name)

  if (!hop0_name %in% names(mode_genes)) {
    return(NULL)
  }

  genes0 <- unique(mode_genes[[hop0_name]])
  genes1 <- unique(mode_genes[[hop1_name]])
  added <- setdiff(genes1, genes0)
  lost <- setdiff(genes0, genes1)

  md1_index <- which(
    vapply(
      modes,
      function(x) identical(x$name, hop1_name),
      logical(1)
    )
  )

  md1 <- if (length(md1_index) == 1L) {
    modes[[md1_index]]
  } else {
    NULL
  }

  tgl1 <- if (!is.null(md1)) {
    md1$ann$target_gene_links
  } else {
    NULL
  }

  path2_rows <- if (is.null(tgl1) ||
      !"path_length" %in% colnames(tgl1)) {
    0L
  } else {
    sum(
      is.finite(tgl1$path_length) &
        tgl1$path_length > 1L,
      na.rm = TRUE
    )
  }

  expanded_membership_rows <- if (is.null(tgl1) ||
      !"in_expanded_target" %in% colnames(tgl1)) {
    0L
  } else {
    sum(tgl1$in_expanded_target %in% TRUE, na.rm = TRUE)
  }

  # Some refined/chromatin objects may not retain metadata$parameters$neighbor_hop.
  # as.integer(NULL) returns integer(0) without throwing an error, which would make
  # data.frame() fail with "arguments imply differing number of rows: 1, 0".
  metadata_hop_raw <- tryCatch(
    md1$ann$metadata$parameters$neighbor_hop,
    error = function(e) NULL
  )
  if (is.null(metadata_hop_raw) || length(metadata_hop_raw) == 0L) {
    metadata_hop <- NA_integer_
  } else {
    metadata_hop <- suppressWarnings(as.integer(metadata_hop_raw[[1L]]))
    if (length(metadata_hop) != 1L || is.na(metadata_hop)) {
      metadata_hop <- NA_integer_
    }
  }

  data.frame(
    Hop0_Mode = hop0_name,
    Hop1_Mode = hop1_name,
    Hop0_N = length(genes0),
    Hop1_N = length(genes1),
    Added_by_Expanded = length(added),
    Lost_from_Primary = length(lost),
    Jaccard = length(intersect(genes0, genes1)) /
      max(1L, length(union(genes0, genes1))),
    Object_Neighbor_Hop = metadata_hop,
    Metadata_Hop_Available = !is.na(metadata_hop),
    PathLength_GT1_Rows = path2_rows,
    In_Expanded_Target_Rows = expanded_membership_rows,
    stringsAsFactors = FALSE
  )
})

hop_scope_qc <- dplyr::bind_rows(hop_scope_qc)

if (nrow(hop_scope_qc) > 0L) {
  write.csv(
    hop_scope_qc,
    file.path(base_out_dir, "Hop_Mode_GeneSet_QC.csv"),
    row.names = FALSE
  )

  if (any(!hop_scope_qc$Metadata_Hop_Available)) {
    missing_meta_modes <- hop_scope_qc$Hop1_Mode[
      !hop_scope_qc$Metadata_Hop_Available
    ]
    warning(
      "neighbor_hop metadata was unavailable for: ",
      paste(missing_meta_modes, collapse = ", "),
      ". Hop scope will still be validated from gene-set containment and ",
      "target_gene_links path_length/in_expanded_target evidence.",
      call. = FALSE
    )
  }

  if (any(hop_scope_qc$Lost_from_Primary > 0L)) {
    bad_modes <- hop_scope_qc$Hop1_Mode[
      hop_scope_qc$Lost_from_Primary > 0L
    ]
    stop(
      "Hop1 primary+expanded gene sets unexpectedly lost primary genes: ",
      paste(bad_modes, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(
      !is.na(hop_scope_qc$Object_Neighbor_Hop) &
        hop_scope_qc$Object_Neighbor_Hop != 1L
  )) {
    bad_modes <- hop_scope_qc$Hop1_Mode[
      !is.na(hop_scope_qc$Object_Neighbor_Hop) &
        hop_scope_qc$Object_Neighbor_Hop != 1L
    ]
    stop(
      "Hop1-labelled modes are backed by objects whose metadata neighbor_hop is not 1: ",
      paste(bad_modes, collapse = ", "),
      call. = FALSE
    )
  }

  if (all(hop_scope_qc$Added_by_Expanded == 0L)) {
    warning(
      "All hop1 primary+expanded gene sets are identical to hop0 primary sets. ",
      "Inspect Hop_Mode_GeneSet_QC.csv and target_gene_links path_length/in_expanded_target; ",
      "the data may contain no eligible expanded targets.",
      call. = FALSE
    )
  } else {
    message(
      "    Hop scope QC: median expanded-added genes = ",
      median(hop_scope_qc$Added_by_Expanded),
      "; range = ",
      min(hop_scope_qc$Added_by_Expanded),
      "-",
      max(hop_scope_qc$Added_by_Expanded)
    )
  }
}


# ── Evidence label schema QC against the loaded looplook objects ─────────────
observed_evidence <- sort(unique(unlist(lapply(modes, function(md) {
  tgl <- md$ann$target_gene_links
  if (is.null(tgl) || nrow(tgl) == 0L || !"evidence" %in% colnames(tgl)) {
    return(character())
  }
  values <- trimws(as.character(tgl$evidence))
  values[!is.na(values) & nzchar(values)]
}), use.names = FALSE)))
unknown_evidence <- setdiff(observed_evidence, names(evidence_palette))
evidence_label_qc <- data.frame(
  Evidence = observed_evidence,
  Recognized = observed_evidence %in% names(evidence_palette),
  stringsAsFactors = FALSE
)
utils::write.csv(
  evidence_label_qc,
  file.path(base_out_dir, "Evidence_Label_QC.csv"),
  row.names = FALSE
)
if (length(unknown_evidence) > 0L) {
  warning(
    "Unrecognized looplook evidence labels will be grouped as 'other': ",
    paste(unknown_evidence, collapse = ", "),
    ". Inspect Evidence_Label_QC.csv.",
    call. = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 样本量梯度循环 — 每个 sample_size 独立输出到子文件夹
# ══════════════════════════════════════════════════════════════════════════════
pool_chip <- chipseeker_genes  # 提前计算，各梯度复用
gsea_status_log <- list()  # Global across all sample sizes
paired_nes_by_size <- list()
delta_stats_by_size <- list()

for (sz in sample_sizes) {
  sample_size <- sz
  out_dir <- file.path(base_out_dir, paste0("size", sz))
  dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
  message(sprintf("\n\n========== Sample Size = %d | Output: %s ==========", sz, out_dir))

# ── 4. GSEA 重抽样 ───────────────────────────────────────────────────────────
tmp_dir <- file.path(out_dir, "tmp_sim")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
ensure_cache_signature(
  tmp_dir,
  c(benchmark_cache_signature_base, list(
    Module = "main_raw_lfc",
    SampleSize = as.integer(sample_size),
    RankedListMD5 = md5_object(glist)
  ))
)

message(sprintf("\n>>> Running %d iterations per mode...", n_iterations))

for (idx in seq_along(modes)) {
  md <- modes[[idx]]
  nm <- md$name
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) {
    gsea_status_log[[length(gsea_status_log) + 1L]] <- make_mode_status(
      sample_size, "main", nm, "empty_mode_genes")
    next
  }

  tmp_csv <- file.path(tmp_dir, paste0("sim_", nm, ".csv"))
  term_loop <- "looplook"
  expected_terms <- c("looplook", "Background")

  expected_main_terms <- if (grepl("^ctrl_", nm)) c("ChIPseeker", "Background") else c("looplook", "Background")
  if (validate_csv_cache(tmp_csv, c("ID", "NES", "pvalue", "Iteration", "Mode"), n_iterations, expected_main_terms)) {
    message(sprintf("    [%2d/32] %-30s SKIP (validated cache)", idx, nm))
    cached_df <- read.csv(tmp_csv, stringsAsFactors = FALSE)
    cached_status <- cached_df %>%
      distinct(Iteration) %>%
      transmute(
        SampleSize = sample_size, Module = "main", Mode = nm,
        Iteration = Iteration, Status = "cached_validated",
        ExpectedTerms = paste(expected_terms, collapse = ";"),
        ReturnedTerms = paste(expected_terms, collapse = ";"),
        MissingTerms = "", Error = NA_character_)
    gsea_status_log <- c(gsea_status_log,
      split(cached_status, seq_len(nrow(cached_status))))
    next
  }
  if (file.exists(tmp_csv)) {
    warning(sprintf("    [%2d/32] %-30s cache invalid, will recompute", idx, nm))
    unlink(tmp_csv)
  }

  pool_loop <- intersect(genes, all_genes)
  # NOTE: Nonassigned_control = genes NOT in looplook AND NOT in ChIPseeker.
  # This is NOT the GSEA statistical background/universe; it is a control gene pool.
  pool_control <- setdiff(all_genes, union(genes, pool_chip))

  if (length(pool_control) < sample_size || length(pool_loop) < sample_size) {
    warning(sprintf("[%s] pool smaller than requested size; matched-N will be used (loop=%d, control=%d)", nm, length(pool_loop), length(pool_control)))
  }
  # Pool-size preflight: skip if either pool is too small for meaningful GSEA
  if (length(pool_loop) < gsea_min_size || length(pool_control) < gsea_min_size) {
    gsea_status_log[[length(gsea_status_log) + 1L]] <- data.frame(
      SampleSize = sample_size, Module = "main", Mode = nm, Iteration = NA_integer_,
      Status = "insufficient_pool",
      ExpectedTerms = paste(c(term_loop, "Background"), collapse = ";"),
      ReturnedTerms = "", MissingTerms = paste(c(term_loop, "Background"), collapse = ";"),
      Error = NA_character_, stringsAsFactors = FALSE)
    message(sprintf("    [%2d/32] %-30s SKIP (pool too small: %s=%d, control=%d, need >= %d)",
      idx, nm, term_loop, length(pool_loop), length(pool_control), gsea_min_size))
    next
  }
  message(sprintf("    [%2d/32] %-30s %s=%d control=%d",
    idx, nm, term_loop, length(pool_loop), length(pool_control)))

  sim_local <- checkpointed_lapply(
  1:n_iterations,
  function(i) {
      set.seed(100 + i)
      n_actual <- min(sample_n(pool_loop, sample_size), sample_n(pool_control, sample_size))
      s_loop <- sample(pool_loop, n_actual)
      s_control <- sample(pool_control, n_actual)

      t2g <- bind_rows(
        data.frame(term = term_loop, gene = s_loop),
        data.frame(term = "Background", gene = s_control))

      gsea_out <- safe_gsea_once(glist, t2g, expected_terms, sample_size,
        "main", nm, i)
      status <- gsea_out$status

      if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
        df <- gsea_out$result[, c("ID", "NES", "pvalue")]
        df$Iteration <- i
        df$Mode <- nm
        df$Requested_N <- sample_size
        df$Actual_N <- n_actual
        df$Pool_N <- ifelse(df$ID == term_loop, length(pool_loop), length(pool_control))
        df$Sampling_fraction <- df$Actual_N / pmax(1, df$Pool_N)
        list(result = df, status = status)
      } else {
        list(result = NULL, status = status)
      }
      },
  checkpoint_dir = file.path(tmp_dir, ".iteration_checkpoints", nm),
  label = sprintf("main size=%d mode=%s", sample_size, nm),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)

    # Collect results and statuses
    statuses <- lapply(sim_local, function(x) x$status)
    gsea_status_log <- c(gsea_status_log, statuses)
    sim_local <- lapply(sim_local, function(x) x$result)
    sim_local <- sim_local[vapply(sim_local, is.data.frame, logical(1))]

  sim_mode_df <- if (length(sim_local) > 0) {
    bind_rows(sim_local) %>% filter(!is.na(NES))
  } else {
    data.frame(ID = character(), NES = numeric(), pvalue = numeric(), Iteration = integer(), Mode = character())
  }
  write.csv(sim_mode_df, tmp_csv, row.names = FALSE)
  rm(sim_local, sim_mode_df)
  gc()
  message(sprintf("    [%2d/32] %-30s -> %s", idx, nm, tmp_csv))
}

# ── 5. 合并 CSV ─────────────────────────────────────────────────────────────
message("\n>>> Merging per-mode CSVs...")
csv_files <- list.files(tmp_dir, pattern = "sim_.*\\.csv$", full.names = TRUE)
sim_list <- lapply(csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  df
})
sim_df <- bind_rows(sim_list)
sim_df$Mode <- factor(sim_df$Mode, levels = sapply(modes, `[[`, "name"))
message(sprintf("    Merged %d files, %d rows total", length(csv_files), nrow(sim_df)))
gc()

# ── 6. 汇总 ──────────────────────────────────────────────────────────────────
mode_mean <- sim_df %>%
  group_by(Mode, ID) %>%
  summarise(Mean_NES = mean(NES), Median_NES = median(NES),
    SD_NES = sd(NES), .groups = "drop")

message("\n>>> Mode x ID mean NES:")
print(as.data.frame(mode_mean %>% tidyr::pivot_wider(names_from = ID, values_from = Mean_NES)))

# ── Paired ΔNES: more appropriate for correlated Monte Carlo subsamples ──
# Since looplook and Background use the same iteration index, compute paired difference
dup_check <- sim_df %>%
  filter(!grepl("ctrl", Mode)) %>%
  dplyr::count(Mode, Iteration, ID) %>%
  filter(n > 1L)
if (nrow(dup_check) > 0L) {
  stop("Duplicate Mode-Iteration-ID rows in sim_df before paired_NES: ", nrow(dup_check), " duplicates")
}

paired_nes <- sim_df %>%
  filter(!grepl("ctrl", Mode)) %>%
  dplyr::select(Mode, Iteration, ID, NES) %>%
  mutate(ID_label = ifelse(ID %in% c("looplook", "ChIPseeker"), "Assigned", as.character(ID))) %>%
  tidyr::pivot_wider(names_from = ID_label, values_from = NES) %>%
  filter(is.finite(Assigned), is.finite(Background)) %>%
  mutate(Delta_NES = Assigned - Background)

paired_nes_summary <- paired_nes %>%
  group_by(Mode) %>%
  summarise(
    N_Paired = n(),
    Median_Delta_NES = median(Delta_NES, na.rm = TRUE),
    Mean_Delta_NES = mean(Delta_NES, na.rm = TRUE),
    Empirical_Low = quantile(Delta_NES, 0.025, na.rm = TRUE, names = FALSE),
    Empirical_High = quantile(Delta_NES, 0.975, na.rm = TRUE, names = FALSE),
    P_Delta_lt_0 = mean(Delta_NES < 0, na.rm = TRUE),
    P_Delta_gt_0 = mean(Delta_NES > 0, na.rm = TRUE),
    .groups = "drop")

message("\n>>> Data-driven decisions: effect size vs simulation-only statistics")
message("    (1) For key figures, show RANK-BISERIAL EFFECT SIZE with bootstrap CI")
message("    (2) For Supplementary only, show Dunn/KW statistics with caveat")
message("    (3) For mode selection: Composite ranking, minimal mode differences")
message("NOTE: Dunn/KW p-values from the subsampling NES are NOT indicators")
message("of biological significance; they show only Monte Carlo simulation stability.")
message("Paired ΔNES: negative = assigned genes more enriched at downregulated end")
message("Use P_lt0 (fraction iterations with ΔNES<0) to assess direction consistency")

# Direction stability: fraction of iterations per mode with Delta_NES < 0
direction_stable <- paired_nes_summary %>%
  mutate(
    Direction = ifelse(P_Delta_lt_0 > 0.95, "down",
                ifelse(P_Delta_gt_0 > 0.95, "up", "unstable")),
    StabilityLabel = ifelse(Direction != "unstable",
      sprintf("P%s=%.3f", ifelse(Direction == "down", "<0", ">0"),
        ifelse(Direction == "down", P_Delta_lt_0, P_Delta_gt_0)),
      sprintf("P=%.3f", pmax(P_Delta_lt_0, P_Delta_gt_0, na.rm = TRUE))))

# ── Gene-set-size diagnostics ────────────────────────────────────────────────
if (any(c("Actual_N", "Requested_N") %in% colnames(sim_df))) {
  gsea_diag <- sim_df %>%
    group_by(Mode, ID) %>%
    summarise(Actual_N = median(Actual_N, na.rm = TRUE),
              Median_NES = median(NES, na.rm = TRUE),
              SD_NES = sd(NES, na.rm = TRUE), .groups = "drop")

  p_nes_size <- ggplot(gsea_diag, aes(x = Actual_N, y = Median_NES, shape = ID)) +
    geom_point(size = 1.5, alpha = 0.6) + scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "Gene-set size vs Median NES",
         subtitle = "Each point = one term within one mode",
         x = "Actual sampled genes (log scale)", y = "Median NES") +
    theme_classic() + theme(legend.position = "bottom")
  queueable_ggsave(file.path(out_dir, "GSEA_Diagnostic_Size_vs_NES.pdf"), p_nes_size, width = 7, height = 5)

  p_sd_size <- ggplot(gsea_diag, aes(x = Actual_N, y = SD_NES, shape = ID)) +
    geom_point(size = 1.5, alpha = 0.6) + scale_x_log10() +
    labs(title = "Gene-set size vs NES variability",
         subtitle = "Smaller sets may show greater NES fluctuation",
         x = "Actual sampled genes (log scale)", y = "SD of NES") +
    theme_classic() + theme(legend.position = "bottom")
  queueable_ggsave(file.path(out_dir, "GSEA_Diagnostic_Size_vs_SD.pdf"), p_sd_size, width = 7, height = 5)
  write.csv(gsea_diag, file.path(out_dir, "GSEA_Sample_Size_Diagnostic.csv"), row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# 绘图部分: 所有图保留，每张图创建后立即保存 PDF + PNG
# ══════════════════════════════════════════════════════════════════════════════

# ── 7. 密度图: 8 assignment families × 4 fill/hop variants ─────────────────────
mode_family_map <- c(
  "anno_all_F"            = "Basic (all)",     "anno_all_T"            = "Basic (all)",
  "anno_all_F_hop1"       = "Basic (all)",     "anno_all_T_hop1"       = "Basic (all)",
  "refined_all_F"         = "Refined (all)",         "refined_all_T"         = "Refined (all)",
  "refined_all_F_hop1"    = "Refined (all)",         "refined_all_T_hop1"    = "Refined (all)",
  "anno_promoter_F"       = "Basic (promoter)", "anno_promoter_T"       = "Basic (promoter)",
  "anno_promoter_F_hop1"  = "Basic (promoter)", "anno_promoter_T_hop1"  = "Basic (promoter)",
  "refined_promoter_F"    = "Refined (promoter)",    "refined_promoter_T"    = "Refined (promoter)",
  "refined_promoter_F_hop1"="Refined (promoter)",    "refined_promoter_T_hop1"="Refined (promoter)",
  "chrom_all_F"           = "Chromatin (all)",       "chrom_all_T"           = "Chromatin (all)",
  "chrom_all_F_hop1"      = "Chromatin (all)",       "chrom_all_T_hop1"      = "Chromatin (all)",
  "chrom_promoter_F"      = "Chromatin (promoter)",  "chrom_promoter_T"      = "Chromatin (promoter)",
  "chrom_promoter_F_hop1" = "Chromatin (promoter)",  "chrom_promoter_T_hop1" = "Chromatin (promoter)",
  "chrom_only_all_F"      = "Chromatin_Only (all)",  "chrom_only_all_T"      = "Chromatin_Only (all)",
  "chrom_only_all_F_hop1" = "Chromatin_Only (all)",  "chrom_only_all_T_hop1" = "Chromatin_Only (all)",
  "chrom_only_promoter_F" = "Chromatin_Only (promoter)", "chrom_only_promoter_T" = "Chromatin_Only (promoter)",
  "chrom_only_promoter_F_hop1" = "Chromatin_Only (promoter)", "chrom_only_promoter_T_hop1" = "Chromatin_Only (promoter)"
)
mode_variant_map <- c(
  "anno_all_F"            = "F (hop=0)",      "anno_all_T"            = "T (hop=0)",
  "anno_all_F_hop1"       = "F (hop=1)",      "anno_all_T_hop1"       = "T (hop=1)",
  "refined_all_F"         = "F (hop=0)",      "refined_all_T"         = "T (hop=0)",
  "refined_all_F_hop1"    = "F (hop=1)",      "refined_all_T_hop1"    = "T (hop=1)",
  "anno_promoter_F"       = "F (hop=0)",      "anno_promoter_T"       = "T (hop=0)",
  "anno_promoter_F_hop1"  = "F (hop=1)",      "anno_promoter_T_hop1"  = "T (hop=1)",
  "refined_promoter_F"    = "F (hop=0)",      "refined_promoter_T"    = "T (hop=0)",
  "refined_promoter_F_hop1"="F (hop=1)",      "refined_promoter_T_hop1"="T (hop=1)",
  "chrom_all_F"           = "F (hop=0)",      "chrom_all_T"           = "T (hop=0)",
  "chrom_all_F_hop1"      = "F (hop=1)",      "chrom_all_T_hop1"      = "T (hop=1)",
  "chrom_promoter_F"      = "F (hop=0)",      "chrom_promoter_T"      = "T (hop=0)",
  "chrom_promoter_F_hop1" = "F (hop=1)",      "chrom_promoter_T_hop1" = "T (hop=1)",
  "chrom_only_all_F"      = "F (hop=0)",      "chrom_only_all_T"      = "T (hop=0)",
  "chrom_only_all_F_hop1" = "F (hop=1)",      "chrom_only_all_T_hop1" = "T (hop=1)",
  "chrom_only_promoter_F" = "F (hop=0)",      "chrom_only_promoter_T" = "T (hop=0)",
  "chrom_only_promoter_F_hop1" = "F (hop=1)", "chrom_only_promoter_T_hop1" = "T (hop=1)"
)

sim_lb <- sim_df %>%
  filter(!grepl("ctrl", Mode)) %>%
  mutate(
    Family  = factor(mode_family_map[as.character(Mode)],
      levels = c("Basic (all)", "Refined (all)", "Chromatin (all)", "Chromatin_Only (all)",
                 "Basic (promoter)", "Refined (promoter)", "Chromatin (promoter)", "Chromatin_Only (promoter)")),
    Variant = factor(mode_variant_map[as.character(Mode)],
      levels = c("F (hop=0)", "T (hop=0)", "F (hop=1)", "T (hop=1)")))

lb_medians <- sim_lb %>%
  group_by(Family, Variant, ID) %>%
  summarise(Median_NES = median(NES, na.rm = TRUE), .groups = "drop")

lb_colors <- c("looplook" = "#E64B35", "Background" = "#B0BEC5")

p_density <- ggplot(sim_lb, aes(x = NES, fill = ID)) +
  geom_density(alpha = 0.7, color = NA) +
  geom_vline(data = lb_medians, aes(xintercept = Median_NES, color = ID),
    linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = lb_colors) +
  scale_color_manual(values = lb_colors, guide = "none") +
  facet_grid(rows = vars(Family), cols = vars(Variant), scales = "free_y") +
  labs(title = "NES: looplook vs Non-assigned control across 8 assignment families x 4 variants",
    subtitle = "Rows share y-axis. Columns: F=LoopOnly T=Filled. hop=hop distance.",
    x = "NES", y = "Density") +
  theme_classic() +
  theme(legend.position = "bottom",
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank())

queueable_ggsave(file.path(out_dir, "NES_Density_AllModes.pdf"), p_density, width = 12, height = 10)
rm(p_density, sim_lb, lb_medians); gc()

# ── 8. 箱线图: 8 assignment families × 4 fill/hop variants ────────────────────
sim_lb2 <- sim_df %>%
  filter(!grepl("ctrl", Mode)) %>%
  mutate(
    Family  = factor(mode_family_map[as.character(Mode)],
      levels = c("Basic (promoter)", "Refined (promoter)", "Chromatin (promoter)", "Chromatin_Only (promoter)",
                 "Basic (all)", "Refined (all)", "Chromatin (all)", "Chromatin_Only (all)")),
    Variant = factor(mode_variant_map[as.character(Mode)],
      levels = c("F (hop=0)", "T (hop=0)", "F (hop=1)", "T (hop=1)")))

# Compute paired Delta_NES = looplook_NES - Background_NES
# This reduces within-mode iteration-level noise
sim_delta <- sim_lb2 %>%
  filter(!grepl("ctrl", Mode)) %>%
  dplyr::select(Mode, Iteration, ID, NES)

# Check for duplicate Mode-Iteration-ID combos before pivot
delta_dup <- sim_delta %>%
  dplyr::count(Mode, Iteration, ID) %>%
  filter(n > 1L)
if (nrow(delta_dup) > 0L) {
  stop("Duplicate Mode-Iteration-ID rows detected before paired Delta_NES: ",
    nrow(delta_dup), " duplicates")
}

sim_delta <- sim_delta %>%
  tidyr::pivot_wider(names_from = ID, values_from = NES) %>%
  mutate(Delta_NES = looplook - Background) %>%
  filter(is.finite(Delta_NES)) %>%
  mutate(
    Method = case_when(
      grepl("^anno_", Mode)        ~ "anno",
      grepl("^refined_", Mode)     ~ "refined",
        grepl("^chrom_only_", Mode)  ~ "chromatin_only",
      grepl("^chrom_", Mode)       ~ "chromatin"
    ),
    BaseMode = Mode %>%
      str_replace("^anno_", "") %>%
      str_replace("^refined_", "") %>%
      str_replace("^chrom_only_", "") %>%
      str_replace("^chrom_", ""))

lb_colors <- c("looplook" = "#E64B35", "Background" = "#B0BEC5")

p_box_minimal <- ggplot(sim_lb2, aes(x = ID, y = NES, fill = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.2, outlier.alpha = 0.6, outlier.color = "grey40",
               linewidth = 0.5, alpha = 0.8,
               position = position_dodge(0.75), width = 0.65,notch = FALSE, fatten=1) +
  scale_fill_manual(values = lb_colors) +
  scale_color_manual(values = lb_colors) +
  facet_grid(rows = vars(Variant), cols = vars(Family), scales = "free_y") +
  labs(title = "Robustness of GSEA Normalized Enrichment Scores (NES)",
       subtitle = sprintf("%d subsamples per mode (looplook vs. Non-assigned control)", n_iterations),
       x = NULL, y = "Normalized Enrichment Score (NES)") +

  # 【关键1】换用 theme_minimal 作为底层，它天生没有粗重的黑轴
  theme_minimal() +

  theme(
    # ── 核心去乱改造 ──
    # 1. 彻底删掉 L 型粗黑坐标轴线（罪魁祸首）
    axis.line = element_blank(),

    # 2. 用极浅、极细的全包围边框把数据轻轻“托”住，而不是锁住
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.4),

    # 3. 增加格子之间的间距，让密集的 16 宫格有呼吸感
    panel.spacing = unit(0.8, "lines"),

    # 4. 彻底去掉 X 轴底部多余的刻度线（小黑刺）
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    # ─────────────────

    # 横向网格线保持微弱的参考
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),

    # 标题底纹保持空白，无束缚感
    strip.background = element_blank(),
    strip.text = element_text(size = 9.5, face = "bold", color = "black", margin = margin(t = 6, b = 6)),

    # Y轴刻度和图例保持您的完美设定
    axis.text.y = element_text(size = 9, color = "grey10"),
    axis.ticks.y = element_line(color = "grey10", linewidth = 0.4), # 让Y轴刻度线也变柔和

    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.margin = margin(b = -5),
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5, margin = margin(b = 10))
  )

queueable_print_plot(p_box_minimal)
queueable_ggsave(file.path(out_dir, "NES_Boxplot_AllModes.pdf"), p_box_minimal, width = 8, height = 8)
rm(p_box_minimal); gc()

# ── 9. Ridge 图 (按 all/promoter 分组排列) ───────────────────────────────────
mode_order_v1 <- c(
  "anno_all_F",
  "anno_all_T",
  "refined_all_F",
  "refined_all_T",
  "anno_promoter_F",
  "anno_promoter_T",
  "refined_promoter_F",
  "refined_promoter_T",
  "anno_all_F_hop1",
  "anno_all_T_hop1",
  "refined_all_F_hop1",
  "refined_all_T_hop1",
  "anno_promoter_F_hop1",
  "anno_promoter_T_hop1",
  "refined_promoter_F_hop1",
  "refined_promoter_T_hop1"
)

sim_ridge1 <- sim_lb2 %>%
  filter(Mode %in% mode_order_v1) %>%
  mutate(
    Mode = factor(Mode, levels = rev(mode_order_v1)),
    ID = factor(ID, levels = c("Background", "looplook"))
  )

lb_colors_ridge <- c("Background" = "#CBD3D9", "looplook" = "#D95F59")

med_df1 <- sim_ridge1 %>%
  group_by(Mode, ID) %>%
  summarise(MedNES = median(NES, na.rm = TRUE), .groups = "drop")

p_ridge_v1 <- ggplot(sim_ridge1, aes(x = NES, y = Mode, fill = ID)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.6) +
  geom_density_ridges(
    alpha = 0.78, scale = 1.08, color = "white", linewidth = 0.25,
    rel_min_height = 0.015, bandwidth = 0.12) +
  geom_point(
    data = med_df1, aes(x = MedNES, y = Mode, fill = ID),
    shape = 21, size = 1.8, color = "white", stroke = 0.25,
    position = position_nudge(y = 0.04), inherit.aes = FALSE) +
  scale_fill_manual(values = lb_colors_ridge) +
  scale_x_continuous(breaks = seq(-2, 2, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Robust NES distributions of looplook-defined target genes",
    subtitle = sprintf("Sixteen loop-aware assignment modes; %d iterations, requested up to %d genes per term",
                        n_iterations, sample_size),
    x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey35"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    axis.text.y = element_text(face = "bold", size = 8.5, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    plot.margin = margin(10, 18, 10, 10)
  )

queueable_print_plot(p_ridge_v1)
queueable_ggsave(file.path(out_dir, "NES_Ridge_16Modes_grouped.pdf"), p_ridge_v1, width = 10, height = 7)
rm(p_ridge_v1, sim_ridge1, med_df1); gc()

# ── 10. Ridge 图 (按 promoter 优先排列) ──────────────────────────────────────
mode_order_v2 <- c(
  "anno_promoter_F",
  "anno_promoter_T",
  "anno_promoter_F_hop1",
  "anno_promoter_T_hop1",
  "anno_all_F",
  "anno_all_T",
  "anno_all_F_hop1",
  "anno_all_T_hop1",
  "refined_promoter_F",
  "refined_promoter_T",
  "refined_promoter_F_hop1",
  "refined_promoter_T_hop1",
  "refined_all_F",
  "refined_all_T",
  "refined_all_F_hop1",
  "refined_all_T_hop1"
)

sim_ridge2 <- sim_lb2 %>%
  filter(Mode %in% mode_order_v2) %>%
  mutate(
    Mode = factor(Mode, levels = rev(mode_order_v2)),
    ID = factor(ID, levels = c("Background", "looplook"))
  )

lb_colors_ridge2 <- c("Background" = "#DDE3E7", "looplook" = "#D95F59")

med_df2 <- sim_ridge2 %>%
  group_by(Mode, ID) %>%
  summarise(MedNES = median(NES, na.rm = TRUE), .groups = "drop")

p_ridge_v2 <- ggplot(sim_ridge2, aes(x = NES, y = Mode, fill = ID)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.55) +
  geom_density_ridges(
    alpha = 0.68, scale = 1.05, color = "white", linewidth = 0.25,
    rel_min_height = 0.015, bandwidth = 0.12) +
  geom_point(
    data = med_df2, aes(x = MedNES, y = Mode, fill = ID),
    shape = 21, size = 1.7, color = "white", stroke = 0.25,
    position = position_nudge(y = 0.04), inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = lb_colors_ridge2) +
  scale_x_continuous(breaks = seq(-2, 2.2, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Robust negative NES across 16 looplook assignment modes",
    subtitle = sprintf("%d iterations; requested up to %d genes per term",
                        n_iterations, sample_size),
    x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "grey35"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    axis.text.y = element_text(face = "bold", size = 7.8, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    plot.margin = margin(10, 20, 10, 12)
  )

queueable_print_plot(p_ridge_v2)
queueable_ggsave(file.path(out_dir, "NES_Ridge_16Modes.pdf"), p_ridge_v2, width = 10, height = 7)
rm(p_ridge_v2, sim_ridge2, med_df2); gc()

# ── 11. Paired Ridge 图 (anno vs refined, 含 Non-assigned control) ─────────────────────
pair_order <- c(
  "promoter_F | anno",
  "promoter_F | refined",
  "promoter_T | anno",
  "promoter_T | refined",
  "promoter_F_hop1 | anno",
  "promoter_F_hop1 | refined",
  "promoter_T_hop1 | anno",
  "promoter_T_hop1 | refined",
  "all_F | anno",
  "all_F | refined",
  "all_T | anno",
  "all_T | refined",
  "all_F_hop1 | anno",
  "all_F_hop1 | refined",
  "all_T_hop1 | anno",
  "all_T_hop1 | refined"
)

sim_pair_bg <- sim_lb2 %>%
  filter(grepl("^(anno|refined)_", Mode)) %>%
  mutate(
    Method = case_when(
      grepl("^anno_", Mode) ~ "anno",
      grepl("^refined_", Mode) ~ "refined"
    ),
    BaseMode = Mode %>%
      str_replace("^anno_", "") %>%
      str_replace("^refined_", ""),
    PairLabel = paste(BaseMode, Method, sep = " | "),
    PairLabel = factor(PairLabel, levels = rev(pair_order)),
    Method = factor(Method, levels = c("anno", "refined")),
    PlotGroup = case_when(
      ID == "Background" ~ "Background",
      ID == "looplook" & Method == "anno" ~ "anno",
      ID == "looplook" & Method == "refined" ~ "refined"
    )
  )

med_pair_all <- sim_pair_bg %>%
  group_by(PairLabel, PlotGroup) %>%
  summarise(MedNES = median(NES, na.rm = TRUE), .groups = "drop")

med_bg <- med_pair_all %>% filter(PlotGroup == "Background")
med_loop <- med_pair_all %>% filter(PlotGroup != "Background")

plot_colors_pair <- c(
  "Background" = "#D9E1E5",
  "anno"       = "#E9A39A",
  "refined"    = "#B84A45"
)

p_pair_bg <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.55) +
  geom_density_ridges(
    data = sim_pair_bg %>% filter(PlotGroup == "Background"),
    aes(x = NES, y = PairLabel, fill = PlotGroup,
        group = interaction(PairLabel, PlotGroup)),
    alpha = 0.65, scale = 1.03, color = "white", linewidth = 0.22,
    rel_min_height = 0.015, bandwidth = 0.12) +
  geom_density_ridges(
    data = sim_pair_bg %>% filter(PlotGroup != "Background"),
    aes(x = NES, y = PairLabel, fill = PlotGroup,
        group = interaction(PairLabel, PlotGroup)),
    alpha = 0.78, scale = 1.03, color = "white", linewidth = 0.25,
    rel_min_height = 0.015, bandwidth = 0.12) +
  geom_point(
    data = med_bg, aes(x = MedNES, y = PairLabel),
    shape = 21, size = 1.35, fill = "#B8C2C8", color = "white", stroke = 0.20,
    position = position_nudge(y = 0.045), inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(
    data = med_loop, aes(x = MedNES, y = PairLabel, fill = PlotGroup),
    shape = 21, size = 1.75, color = "white", stroke = 0.25,
    position = position_nudge(y = 0.045), inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = plot_colors_pair) +
  coord_cartesian(xlim = c(-2, 2.8), clip = "off") +
  scale_x_continuous(breaks = seq(-2, 2.5, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "NES distributions across paired looplook assignment modes",
    subtitle = "Non-assigned control shown as reference; annotation and refined modes paired by matched settings",
    x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "grey35"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    axis.text.y = element_text(face = "bold", size = 7.6, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    plot.margin = margin(10, 20, 10, 12)
  )

queueable_print_plot(p_pair_bg)
queueable_ggsave(file.path(out_dir, "NES_Ridge_Paired_Anno_Refined_withBackground.pdf"), p_pair_bg, width = 10, height = 7)
rm(p_pair_bg, sim_pair_bg, med_pair_all, med_bg, med_loop); gc()

# ── 12. Paired Ridge 分 hop (make_hop_plot) ──────────────────────────────────
plot_colors <- c(
  "Background" = "#D9E1E5",
  "anno"       = "#E9A39A",
  "refined"    = "#B84A45"
)

make_hop_plot <- function(sim_lb2, hop = 0, xlim_range = c(-2.5, 2.2)) {

  if (hop == 0) {
    pair_order <- c(
      "promoter_F | anno",
      "promoter_F | refined",
      "promoter_T | anno",
      "promoter_T | refined",
      "all_F | anno",
      "all_F | refined",
      "all_T | anno",
      "all_T | refined"
    )
    mode_keep <- c(
      "anno_promoter_F", "refined_promoter_F",
      "anno_promoter_T", "refined_promoter_T",
      "anno_all_F", "refined_all_F",
      "anno_all_T", "refined_all_T"
    )
    fig_title <- "NES distributions across paired looplook assignment modes (hop = 0)"
    fig_subtitle <- "Non-assigned control shown as reference; annotation and refined modes paired by matched settings"
    out_name <- "Ridge_Anno_Refined_wCtrl_hop0"
  }

  if (hop == 1) {
    pair_order <- c(
      "promoter_F_hop1 | anno",
      "promoter_F_hop1 | refined",
      "promoter_T_hop1 | anno",
      "promoter_T_hop1 | refined",
      "all_F_hop1 | anno",
      "all_F_hop1 | refined",
      "all_T_hop1 | anno",
      "all_T_hop1 | refined"
    )
    mode_keep <- c(
      "anno_promoter_F_hop1", "refined_promoter_F_hop1",
      "anno_promoter_T_hop1", "refined_promoter_T_hop1",
      "anno_all_F_hop1", "refined_all_F_hop1",
      "anno_all_T_hop1", "refined_all_T_hop1"
    )
    fig_title <- "NES distributions across paired looplook assignment modes (hop = 1)"
    fig_subtitle <- "Non-assigned control shown as reference; annotation and refined modes paired by matched settings"
    out_name <- "Ridge_Anno_Refined_wCtrl_hop1"
  }

  sim_plot <- sim_lb2 %>%
    filter(Mode %in% mode_keep) %>%
    mutate(
      Method = case_when(
        grepl("^anno_", Mode) ~ "anno",
        grepl("^refined_", Mode) ~ "refined"
      ),
      BaseMode = Mode %>%
        str_replace("^anno_", "") %>%
        str_replace("^refined_", ""),
      PairLabel = paste(BaseMode, Method, sep = " | "),
      PairLabel = factor(PairLabel, levels = rev(pair_order)),
      Method = factor(Method, levels = c("anno", "refined")),
      PlotGroup = case_when(
        ID == "Background" ~ "Background",
        ID == "looplook" & Method == "anno" ~ "anno",
        ID == "looplook" & Method == "refined" ~ "refined"
      )
    )

  med_all <- sim_plot %>%
    group_by(PairLabel, PlotGroup) %>%
    summarise(MedNES = median(NES, na.rm = TRUE), .groups = "drop")

  med_bg <- med_all %>% filter(PlotGroup == "Background")
  med_loop <- med_all %>% filter(PlotGroup != "Background")

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.55) +
    geom_density_ridges(
      data = sim_plot %>% filter(PlotGroup == "Background"),
      aes(x = NES, y = PairLabel, fill = PlotGroup,
          group = interaction(PairLabel, PlotGroup)),
      alpha = 0.35, scale = 1.03, color = "white", linewidth = 0.22,
      rel_min_height = 0.015, bandwidth = 0.12) +
    geom_density_ridges(
      data = sim_plot %>% filter(PlotGroup != "Background"),
      aes(x = NES, y = PairLabel, fill = PlotGroup,
          group = interaction(PairLabel, PlotGroup)),
      alpha = 0.78, scale = 1.03, color = "white", linewidth = 0.25,
      rel_min_height = 0.015, bandwidth = 0.12) +
    geom_point(
      data = med_bg, aes(x = MedNES, y = PairLabel),
      shape = 21, size = 1.35, fill = "#B8C2C8", color = "white", stroke = 0.20,
      position = position_nudge(y = -0.035), inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(
      data = med_loop, aes(x = MedNES, y = PairLabel, fill = PlotGroup),
      shape = 21, size = 1.75, color = "white", stroke = 0.25,
      position = position_nudge(y = 0.045), inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = plot_colors) +
    coord_cartesian(xlim = xlim_range, clip = "off") +
    scale_x_continuous(breaks = seq(-2, 2, 1), expand = expansion(mult = c(0.01, 0.03))) +
    labs(title = fig_title, subtitle = fig_subtitle,
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9.3, color = "grey35"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.45, "cm"),
      axis.text.y = element_text(face = "bold", size = 8.2, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      plot.margin = margin(10, 20, 10, 12)
    )

  queueable_print_plot(p)
  queueable_ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p, width = 9.2, height = 5.8)
  return(p)
}

p_hop0 <- make_hop_plot(sim_lb2, hop = 0)
p_hop1 <- make_hop_plot(sim_lb2, hop = 1)


rm(p_hop0, p_hop1); gc()

# ── 13. Overlay Ridge 图 (Paired Delta_NES: looplook − Background per pipeline) ──────
# Each mode's Background is paired with its own looplook NES,
# providing a within-mode enrichment contrast.
delta_colors <- c(
  "anno"            = "#1F77B4",
  "refined"         = "#9467BD",
  "chromatin"       = "#E64B35",
  "chromatin_only"  = "#FFA500"
)

make_hop_overlay_plot <- function(sim_delta, hop = 0, xlim_range = c(-3, 3)) {

  if (hop == 0) {
    mode_keep <- c(
      "anno_promoter_F", "refined_promoter_F", "chrom_promoter_F", "chrom_only_promoter_F",
      "anno_promoter_T", "refined_promoter_T", "chrom_promoter_T", "chrom_only_promoter_T",
      "anno_all_F",      "refined_all_F",      "chrom_all_F",      "chrom_only_all_F",
      "anno_all_T",      "refined_all_T",      "chrom_all_T",      "chrom_only_all_T"
    )
    base_order <- c("promoter_F", "promoter_T", "all_F", "all_T")
    fig_title <- "Delta_NES across looplook assignment modes (hop = 0)"
    out_name  <- "DeltaNES_Ridge_Overlay_hop0"
  }

  if (hop == 1) {
    mode_keep <- c(
      "anno_promoter_F_hop1", "refined_promoter_F_hop1", "chrom_promoter_F_hop1", "chrom_only_promoter_F_hop1",
      "anno_promoter_T_hop1", "refined_promoter_T_hop1", "chrom_promoter_T_hop1", "chrom_only_promoter_T_hop1",
      "anno_all_F_hop1",      "refined_all_F_hop1",      "chrom_all_F_hop1",      "chrom_only_all_F_hop1",
      "anno_all_T_hop1",      "refined_all_T_hop1",      "chrom_all_T_hop1",      "chrom_only_all_T_hop1"
    )
    base_order <- c("promoter_F_hop1", "promoter_T_hop1", "all_F_hop1", "all_T_hop1")
    fig_title <- "Delta_NES across looplook assignment modes (hop = 1)"
    out_name  <- "DeltaNES_Ridge_Overlay_hop1"
  }

  sim_plot <- sim_delta %>%
    filter(Mode %in% mode_keep) %>%
    mutate(
      BaseMode = factor(BaseMode, levels = rev(base_order)),
      Method = factor(Method, levels = c("anno", "refined", "chromatin", "chromatin_only")))

  med_all <- sim_plot %>%
    group_by(BaseMode, Method) %>%
    summarise(MedDelta = median(Delta_NES, na.rm = TRUE), .groups = "drop")

  med_anno        <- med_all %>% filter(Method == "anno")
  med_refined     <- med_all %>% filter(Method == "refined")
  med_chrom       <- med_all %>% filter(Method == "chromatin")
  med_chrom_only  <- med_all %>% filter(Method == "chromatin_only")

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_density_ridges(
      data = sim_plot %>% filter(Method == "anno"),
      aes(x = Delta_NES, y = BaseMode, fill = Method, group = interaction(BaseMode, Method)),
      alpha = 0.70, scale = 0.95, color = "#1F77B4", linewidth = 0.24,
      rel_min_height = 0.015, bandwidth = 0.12, quantile_lines = TRUE,
      quantiles = 2, vline_linetype = "dashed", vline_color = alpha("grey40", 0.8), vline_linewidth = 0.35) +
    geom_density_ridges(
      data = sim_plot %>% filter(Method == "refined"),
      aes(x = Delta_NES, y = BaseMode, fill = Method, group = interaction(BaseMode, Method)),
      alpha = 0.75, scale = 0.95, color = "#9467BD", linewidth = 0.25,
      rel_min_height = 0.015, bandwidth = 0.12, quantile_lines = TRUE,
      quantiles = 2, vline_linetype = "dashed", vline_color = alpha("grey40", 0.8), vline_linewidth = 0.35) +
    geom_density_ridges(
      data = sim_plot %>% filter(Method == "chromatin"),
      aes(x = Delta_NES, y = BaseMode, fill = Method, group = interaction(BaseMode, Method)),
      alpha = 0.78, scale = 0.95, color = "#E64B35", linewidth = 0.26,
      rel_min_height = 0.015, bandwidth = 0.12, quantile_lines = TRUE,
      quantiles = 2, vline_linetype = "dashed", vline_color = alpha("grey40", 0.8), vline_linewidth = 0.35) +
    geom_point(
      data = med_anno, aes(x = MedDelta, y = BaseMode),
      shape = 21, size = 1.60, fill = delta_colors["anno"], color = "white", stroke = 0.22,
      position = position_nudge(y = 0.00), inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(
      data = med_refined, aes(x = MedDelta, y = BaseMode),
      shape = 21, size = 1.75, fill = delta_colors["refined"], color = "white", stroke = 0.24,
      position = position_nudge(y = 0), inherit.aes = FALSE, show.legend = FALSE) +
    geom_density_ridges(
      data = sim_plot %>% filter(Method == "chromatin_only"),
      aes(x = Delta_NES, y = BaseMode, fill = Method, group = interaction(BaseMode, Method)),
      alpha = 0.78, scale = 0.95, color = "#FFA500", linewidth = 0.26,
      rel_min_height = 0.015, bandwidth = 0.12, quantile_lines = TRUE,
      quantiles = 2, vline_linetype = "dashed", vline_color = alpha("grey40", 0.8), vline_linewidth = 0.35) +
    geom_point(
      data = med_chrom, aes(x = MedDelta, y = BaseMode),
      shape = 21, size = 1.85, fill = delta_colors["chromatin"], color = "white", stroke = 0.25,
      position = position_nudge(y = 0), inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(
      data = med_chrom_only, aes(x = MedDelta, y = BaseMode),
      shape = 21, size = 1.85, fill = delta_colors["chromatin_only"], color = "white", stroke = 0.25,
      position = position_nudge(y = 0), inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = delta_colors) +
    coord_cartesian(xlim = xlim_range, clip = "off") +
    scale_x_continuous(breaks = seq(-3, 3, 1), expand = expansion(mult = c(0.01, 0.03))) +
    guides(fill = guide_legend(override.aes = list(
      alpha = 1, color = "transparent", size = 0,
      vline_linewidth = 0, vline_colour = "transparent", linetype = "blank"))) +
    labs(
      title = fig_title,
      subtitle = "Paired Delta_NES = NES_looplook − NES_Background (within-mode enrichment contrast)",
      x = "Delta_NES (looplook − non-assigned control)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8.8, color = "grey35"),
      legend.position = "top", legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.45, "cm"),
      axis.text.y = element_text(face = "bold", size = 8.6, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.6, linetype = "solid"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      plot.margin = margin(10, 20, 10, 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA))

  queueable_print_plot(p)
  queueable_ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p, width = 6, height = 4.8)
  return(p)
}

p_overlay_hop0 <- make_hop_overlay_plot(sim_delta, hop = 0)
p_overlay_hop1 <- make_hop_overlay_plot(sim_delta, hop = 1)

# ── 13b. 箱线图: Paired Delta_NES 四组比较 ──
make_hop_boxplot <- function(sim_delta, hop = 0) {

  symnum.args <- list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "ns"))

  if (hop == 0) {
    mode_keep <- c(
      "anno_promoter_F", "refined_promoter_F", "chrom_promoter_F", "chrom_only_promoter_F",
      "anno_promoter_T", "refined_promoter_T", "chrom_promoter_T", "chrom_only_promoter_T",
      "anno_all_F",      "refined_all_F",      "chrom_all_F",      "chrom_only_all_F",
      "anno_all_T",      "refined_all_T",      "chrom_all_T",      "chrom_only_all_T")
    base_order <- c("promoter_F", "promoter_T", "all_F", "all_T")
    fig_title <- "Delta_NES across pipelines (hop = 0)"
    out_name  <- "Delta_Boxplot_4Pipe_hop0"
  }

  if (hop == 1) {
    mode_keep <- c(
      "anno_promoter_F_hop1", "refined_promoter_F_hop1", "chrom_promoter_F_hop1", "chrom_only_promoter_F_hop1",
      "anno_promoter_T_hop1", "refined_promoter_T_hop1", "chrom_promoter_T_hop1", "chrom_only_promoter_T_hop1",
      "anno_all_F_hop1",      "refined_all_F_hop1",      "chrom_all_F_hop1",      "chrom_only_all_F_hop1",
      "anno_all_T_hop1",      "refined_all_T_hop1",      "chrom_all_T_hop1",      "chrom_only_all_T_hop1")
    base_order <- c("promoter_F_hop1", "promoter_T_hop1", "all_F_hop1", "all_T_hop1")
    fig_title <- "Delta_NES across pipelines (hop = 1)"
    out_name  <- "Delta_Boxplot_4Pipe_hop1"
  }

  sim_plot <- sim_delta %>%
    filter(Mode %in% mode_keep) %>%
    mutate(
      BaseMode = factor(BaseMode, levels = base_order),
      Method = factor(Method, levels = c("anno", "refined", "chromatin", "chromatin_only")))

  desc_stats <- sim_plot %>%
    group_by(BaseMode, Method) %>%
    summarise(
      Median_Delta = median(Delta_NES, na.rm = TRUE),
      Q025 = quantile(Delta_NES, 0.025, na.rm = TRUE),
      Q975 = quantile(Delta_NES, 0.975, na.rm = TRUE),
      P_lt0 = mean(Delta_NES < 0, na.rm = TRUE),
      P_gt0 = mean(Delta_NES > 0, na.rm = TRUE),
      N = n(), .groups = "drop") %>%
    mutate(hop = hop)
  message(sprintf("    [delta hop=%d] Pipeline delta summary (%d iterations per mode):", hop, n_iterations))
  for (bm in levels(sim_plot$BaseMode)) {
    bm_stats <- desc_stats %>% filter(BaseMode == bm) %>% arrange(desc(Median_Delta))
    message(sprintf("      %-20s %s", bm,
      paste(sprintf("%s=%.3f [%.3f,%.3f]", bm_stats$Method,
        bm_stats$Median_Delta, bm_stats$Q025, bm_stats$Q975), collapse=" | ")))
  }

  p <- ggplot(sim_plot, aes(x = Method, y = Delta_NES)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    geom_boxplot(aes(fill = Method), outlier.shape = 16, outlier.size = 0.8, outlier.alpha = 0.3,
                 linewidth = 0.3, alpha = 0.85, width = 0.65, fatten = 1, notch = FALSE) +
    scale_fill_manual(values = delta_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)),
                       breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
    facet_wrap(~ BaseMode, ncol = 4, scales = "free_y") +
    labs(title = fig_title,
         subtitle = sprintf("Paired Delta_NES per mode (%d iterations)", n_iterations),
         x = NULL, y = "Delta_NES (looplook − non-assigned control)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.4),
      panel.spacing = unit(0.8, "lines"),
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(size = 8, face = "bold", color = "black", angle = 30, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 9, face = "bold", color = "black", margin = margin(t = 6, b = 6)),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.ticks.y = element_line(color = "grey50", linewidth = 0.4),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5, margin = margin(b = 10)))

  queueable_print_plot(p)
  queueable_ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p, width = 8, height = 4)
  return(list(plot = p, stats = desc_stats))
}

res_box_hop0 <- make_hop_boxplot(sim_delta, hop = 0)
res_box_hop1 <- make_hop_boxplot(sim_delta, hop = 1)

delta_stats_hop0 <- res_box_hop0$stats
delta_stats_hop1 <- res_box_hop1$stats
delta_stats_all <- bind_rows(delta_stats_hop0, delta_stats_hop1)
paired_nes_by_size[[as.character(sample_size)]] <- paired_nes_summary
delta_stats_by_size[[as.character(sample_size)]] <- delta_stats_all

rm(res_box_hop0, res_box_hop1); gc()
rm(p_overlay_hop0, p_overlay_hop1); gc()

# ══════════════════════════════════════════════════════════════════════════════
# 统计分析与数据导出
# ══════════════════════════════════════════════════════════════════════════════

# ── 14. Kruskal-Wallis + Dunn's + Excel 导出 ─────────────────────────────────
sim_l_stats <- sim_df %>%
  filter(!grepl("ctrl", Mode), ID == "looplook") %>%
  mutate(NES = as.numeric(NES), Mode_Safe = gsub("[ \n()-]", "_", Mode)) %>%
  filter(!is.na(NES))
kw_res <- kruskal.test(NES ~ Mode_Safe, data = sim_l_stats)
dunn_res <- sim_l_stats %>% rstatix::dunn_test(NES ~ Mode_Safe, p.adjust.method = "BH")

out_xlsx <- file.path(out_dir,   sprintf("GSEA_Resampling_%dModes_%diter.xlsx", length(modes), n_iterations))
wb <- createWorkbook()

.safe_add_sheet <- function(wb, name, data) {
  if (nchar(name) > 31) name <- substr(name, 1, 31)
  if (name %in% names(wb)) {
    message(sprintf("    Excel: '%s' already exists, skip", name))
    return()
  }
  if (is.null(data)) { message(sprintf("    Excel: '%s' data is NULL, skip", name)); return() }
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(data))
}

.safe_add_sheet(wb, "Mode_x_ID_Mean_NES",
  mode_mean %>% tidyr::pivot_wider(names_from = ID, values_from = Mean_NES) %>% as.data.frame())
.safe_add_sheet(wb, "Exploratory_KW_MC",
  data.frame(Statistic = kw_res$statistic, Pvalue = kw_res$p.value))
.safe_add_sheet(wb, "Exploratory_Dunn_MC", as.data.frame(dunn_res))
.safe_add_sheet(wb, "Paired_Delta_NES", as.data.frame(paired_nes_summary))
.safe_add_sheet(wb, "Raw_All", as.data.frame(sim_df))
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("    Early save: ", out_xlsx)

# ── 15. Unique / Intersection GSEA ──────────────────────────────────────────
message("\n>>> Per-mode Unique vs Intersection GSEA...")
uniq_tmp_dir <- file.path(out_dir, "tmp_unique")
dir.create(uniq_tmp_dir, recursive = TRUE, showWarnings = FALSE)
ensure_cache_signature(
  uniq_tmp_dir,
  c(benchmark_cache_signature_base, list(
    Module = "unique_raw_lfc",
    SampleSize = as.integer(sample_size),
    RankedListMD5 = md5_object(glist)
  ))
)

for (idx in seq_along(modes)) {
  md <- modes[[idx]]
  nm <- md$name
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) next

  loop_genes <- intersect(unique(toupper(genes)), all_genes)
  chip_g     <- intersect(get_chip_reference(md), all_genes)
  only_loop  <- setdiff(loop_genes, chip_g)
  inter_set  <- intersect(loop_genes, chip_g)
  only_chip  <- setdiff(chip_g, loop_genes)
  unique_pools <- list(
    only_looplook = only_loop,
    intersection  = inter_set,
    only_ChIPseeker = only_chip)

  eligible_cats <- sum(as.integer(lengths(unique_pools) >= gsea_min_size))
  if (eligible_cats < 2) {
    gsea_status_log[[length(gsea_status_log) + 1L]] <- make_mode_status(
      sample_size, "unique", nm, "insufficient_terms",
      paste(c("only_looplook", "intersection", "only_ChIPseeker"), collapse = ";"),
      paste0("only ", eligible_cats, " categories >= ", gsea_min_size))
    message(sprintf("    [unique %2d/32] %-30s SKIP (only %d eligible categories)",
      idx, nm, eligible_cats))
    next
  }

  unique_pools <- unique_pools[lengths(unique_pools) >= gsea_min_size]
  expected_terms_u <- names(unique_pools)
  uniq_csv <- file.path(uniq_tmp_dir, paste0("uniq_", nm, ".csv"))
  if (validate_csv_cache(uniq_csv, c("ID", "NES", "pvalue", "Iteration", "Mode"),
      n_iterations, expected_terms_u)) {
    message(sprintf("    [unique %2d/32] %-30s SKIP (validated cache)", idx, nm))
    cached_df <- read.csv(uniq_csv, stringsAsFactors = FALSE)
    cached_status <- cached_df %>%
      distinct(Iteration) %>%
      transmute(
        SampleSize = sample_size, Module = "unique", Mode = nm,
        Iteration = Iteration, Status = "cached_validated",
        ExpectedTerms = paste(sort(expected_terms_u), collapse = ";"),
        ReturnedTerms = paste(sort(expected_terms_u), collapse = ";"),
        MissingTerms = "", Error = NA_character_)
    gsea_status_log <- c(gsea_status_log,
      split(cached_status, seq_len(nrow(cached_status))))
    next
  }
  if (file.exists(uniq_csv)) {
    warning(sprintf("    [unique %2d/32] %-30s cache invalid, will recompute", idx, nm))
    unlink(uniq_csv)
  }

      uniq_local <- checkpointed_lapply(
  1:n_iterations,
  function(i) {
      set.seed(200 + i)
      n_u <- min(sample_size, lengths(unique_pools))
      sampled <- lapply(unique_pools, function(p) sample(p, n_u))
      t2g_parts <- lapply(names(sampled), function(tn) data.frame(term = tn, gene = sampled[[tn]], stringsAsFactors = FALSE))
      t2g <- bind_rows(t2g_parts)
      expected_terms_u <- names(sampled)

      gsea_out <- safe_gsea_once(glist, t2g, expected_terms_u, sample_size,
        "unique", nm, i)
      status <- gsea_out$status

      if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
        df <- gsea_out$result[, c("ID", "NES", "pvalue")]
        df$Iteration <- i
        df$Mode <- nm
        df$Requested_N <- sample_size
        actual_map <- sapply(sampled, length)
        pool_map   <- sapply(unique_pools, length)
        df$Actual_N <- unname(actual_map[df$ID])
        df$Pool_N   <- unname(pool_map[df$ID])
        df$Sampling_fraction <- df$Actual_N / pmax(1, df$Pool_N)
        list(result = df, status = gsea_out$status)
      } else {
        list(result = NULL, status = gsea_out$status)
      }
      },
  checkpoint_dir = file.path(uniq_tmp_dir, ".iteration_checkpoints", nm),
  label = sprintf("unique size=%d mode=%s", sample_size, nm),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)
    statuses <- lapply(uniq_local, function(x) x$status)
    gsea_status_log <- c(gsea_status_log, statuses[!vapply(statuses, is.null, logical(1))])
    uniq_local <- lapply(uniq_local, function(x) x$result)
      uniq_local <- uniq_local[vapply(uniq_local, is.data.frame, logical(1))]
  uniq_mode_df <- if (length(uniq_local) > 0) {
    bind_rows(uniq_local) %>% filter(!is.na(NES))
  } else {
    data.frame(ID = character(), NES = numeric(), pvalue = numeric(), Iteration = integer(), Mode = character())
  }
  write.csv(uniq_mode_df, uniq_csv, row.names = FALSE)
  rm(uniq_local, uniq_mode_df); gc()
  message(sprintf("    [unique %2d/32] %-30s done", idx, nm))
}

message("\n>>> Merging unique CSVs...")
uniq_csv_files <- list.files(uniq_tmp_dir, pattern = "uniq_.*\\.csv$", full.names = TRUE)
uniq_list <- lapply(uniq_csv_files, function(f) {
  if (file.size(f) < 20) return(NULL)
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", f); NULL })
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
})
uniq_list <- uniq_list[!vapply(uniq_list, is.null, logical(1))]
unique_df <- bind_rows(uniq_list)
unique_df$ID <- factor(unique_df$ID, levels = c("only_looplook", "intersection", "only_ChIPseeker"))
unique_df$Mode <- factor(unique_df$Mode, levels = sapply(modes, `[[`, "name"))
message(sprintf("    Merged %d files, %d rows", length(uniq_csv_files), nrow(unique_df)))
message("    Unique visualizations skipped (hang unresolved)")

if (FALSE) {  # Unique density/boxplot/forest visualizations -- permanently skipped
# ── Helper: write Dunn's test results to Excel ──
.write_dunn_sheet <- function(wb, name, df) {
  if (name %in% names(wb)) {
    message(sprintf("    Excel: '%s' already exists, skip", name))
    return()
  }
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("    Excel: %s SKIP (empty)", name))
    return()
  }
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(df))
  message(sprintf("    Excel: %s (%d rows)", name, nrow(df)))
}

# Unique 密度图 — 按 hop 拆成两幅
uniq_colors <- c("only_looplook" = "#E64B35", "intersection" = "#00A087", "only_ChIPseeker" = "#4DBBD5")

# ── Unified helper: subset unique/distal/subset GSEA data for plotting ──
# Used by all Unique/Distal/Promoter/Genic GSEA density and boxplot functions.
# Replaces three previously duplicated .make_*_uniq_base implementations.
make_uniq_base <- function(df, hop, method = NULL,
    ref_type = "chipseeker", exclude_chrom_only = FALSE,
    pipeline_prefixes = c("anno", "refined", "chrom")) {
  # Determine which modes to keep
  if (hop == 0) {
    keep <- grep("_(promoter|all)_(F|T)$", levels(df$Mode), value = TRUE)
  } else {
    keep <- grep("_hop1$", levels(df$Mode), value = TRUE)
  }
  # Filter to allowed pipeline prefixes
  prefix_pat <- paste0("^(", paste(pipeline_prefixes, collapse = "|"), ")_")
  keep <- keep[grepl(prefix_pat, keep)]
  if (exclude_chrom_only && !is.null(method) && method != "chromatin") {
    keep <- keep[!grepl("^chrom_only_", keep)]
  }
  if (!is.null(method)) {
    if (method == "chromatin") {
      keep <- keep[grepl("^chrom_", keep) & !grepl("^chrom_only_", keep)]
    } else if (method == "chromatin_only") {
      keep <- keep[grepl("^chrom_only_", keep)]
    } else {
      keep <- keep[grepl(paste0("^", method, "_"), keep)]
    }
  }

  # Determine ID factor levels based on reference type
  id_3rd <- switch(tolower(ref_type),
    "chipseeker"   = "only_ChIPseeker",
    "distal"       = "only_distal",
    "promoter"     = "only_Promoter",
    "genic"        = "only_Genic",
    "non_promoter" = "only_Non_Promoter")
  if (is.null(id_3rd)) {
    stop("Unsupported ref_type in make_uniq_base: ", ref_type)
  }
  id_levels <- c("only_looplook", "intersection", id_3rd)

  sub <- df %>%
    filter(Mode %in% keep) %>%
    mutate(
      BaseMode = Mode %>% str_replace("^anno_","") %>%
                 str_replace("^refined_","") %>%
                 str_replace("^chrom_only_","") %>%
                 str_replace("^chrom_",""),
      Method = case_when(
        grepl("^anno_", Mode)        ~ "anno",
        grepl("^refined_", Mode)     ~ "refined",
        grepl("^chrom_only_", Mode)  ~ "chromatin_only",
        grepl("^chrom_", Mode)       ~ "chromatin",
        TRUE ~ "unknown"),
      ID = factor(ID, levels = id_levels))
  if (!is.null(method)) sub <- sub %>% filter(Method == method)
  if (nrow(sub) == 0) return(NULL)
  sub
}

# ── Unified density plot for Unique/Distal/Promoter/Genic GSEA ──
make_uniq_density <- function(df, hop, method = NULL,
    ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique",
    colors = uniq_colors, n_iter_label = n_iterations) {
  sub <- make_uniq_base(df, hop, method, ref_type = ref_type)
  if (is.null(sub) || nrow(sub) == 0) return(NULL)
  sub_med <- sub %>% group_by(BaseMode, ID) %>%
    summarise(Median_NES = median(NES, na.rm = TRUE), .groups = "drop")

  hlab <- paste0("hop = ", hop)
  if (!is.null(method)) hlab <- paste(hlab, "|", method)

  nm <- if (is.null(method)) paste0("hop", hop) else paste0("hop", hop, "_", method)
  p <- ggplot(sub, aes(x = NES, fill = ID)) +
    geom_density(alpha = 0.7, color = NA) +
    geom_vline(data = sub_med, aes(xintercept = Median_NES, color = ID),
               linetype = "dashed", linewidth = 0.4) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors, guide = "none") +
    facet_wrap(~ BaseMode, ncol = 4, scales = "free_y") +
    labs(title = paste("NES by Gene Set Origin (", label, ") \u2014", hlab),
         subtitle = sprintf("%d subsamples per mode", n_iter_label),
         x = "NES", y = "Density") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.4),
          panel.spacing = unit(0.8, "lines"), panel.grid = element_blank(),
          legend.position = "bottom", strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold", color = "black",
                                    margin = margin(t = 6, b = 6)),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.ticks.y = element_line(color = "grey50", linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
          plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5,
                                       margin = margin(b = 10)))
  queueable_ggsave(file.path(out_dir, paste0(out_prefix, "_Density_", nm, ".pdf")), p,
         width = 12, height = 4)
  p
}


# density: combined + per-method
# Unique density/boxplot/forest visualizations
p_uniq_density0 <- make_uniq_density(unique_df, hop = 0, ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors); rm(p_uniq_density0); gc()
p_uniq_density1 <- make_uniq_density(unique_df, hop = 1, ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors); rm(p_uniq_density1); gc()

p_uniq_density0_anno <- make_uniq_density(unique_df, hop = 0, method = "anno", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density1_anno <- make_uniq_density(unique_df, hop = 1, method = "anno", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density0_refined <- make_uniq_density(unique_df, hop = 0, method = "refined", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density1_refined <- make_uniq_density(unique_df, hop = 1, method = "refined", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density0_chromatin <- make_uniq_density(unique_df, hop = 0, method = "chromatin", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density1_chromatin <- make_uniq_density(unique_df, hop = 1, method = "chromatin", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density0_chrom_only <- make_uniq_density(unique_df, hop = 0, method = "chromatin_only", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
p_uniq_density1_chrom_only <- make_uniq_density(unique_df, hop = 1, method = "chromatin_only", ref_type = "chipseeker", label = "ChIPseeker", out_prefix = "GSEA_Unique", colors = uniq_colors)
rm(p_uniq_density0_anno, p_uniq_density1_anno,
   p_uniq_density0_refined, p_uniq_density1_refined,
   p_uniq_density0_chromatin, p_uniq_density1_chromatin); gc()

# Unique 箱线图 — 按 hop 拆成两幅，加 Dunn's 检验
.make_uniq_box <- function(unique_df, hop, method = NULL) {
  df <- make_uniq_base(unique_df, hop, method, ref_type = "chipseeker")
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("    .make_uniq_box hop=%d method=%s: no data, returning NULL",
                    hop, ifelse(is.null(method), "all", method)))
    return(list(plot = NULL, stats = NULL))
  }

  hlab <- paste0("hop = ", hop)
  if (!is.null(method)) hlab <- paste(hlab, "|", method)
  nm <- if (is.null(method)) paste0("hop", hop) else paste0("hop", hop, "_", method)

  dunn_uniq <- df %>%
    rstatix::group_by(BaseMode) %>%
    rstatix::dunn_test(NES ~ ID, p.adjust.method = "BH") %>%
    rstatix::add_xy_position(x = "ID") %>%
    mutate(x = (xmin + xmax) / 2) %>%
    mutate(
      p.adj.signif = sapply(p.adj, function(p) {
        if (is.na(p)) return("ns")
        if (p < 0.001) return("***")
        if (p < 0.01)  return("**")
        if (p < 0.05)  return("*")
        return("ns")
      })
    )
  dunn_uniq <- dunn_uniq %>%
    group_by(BaseMode) %>%
    mutate(y.position = y.position + seq(0, 0.8, length.out = n())) %>%
    ungroup()

  p <- ggplot(df, aes(x = ID, y = NES)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.3) +
    geom_boxplot(aes(fill = ID), outlier.shape = 16, outlier.size = 0.8, outlier.alpha = 0.3,
                 linewidth = 0.3, alpha = 0.85, width = 0.65, fatten = 1,notch = FALSE)

  # Note: Dunn's test results exported to exploratory Excel only (not shown on plot)
  p <- p +
    scale_fill_manual(values = uniq_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)),
                       breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
    facet_wrap(~ BaseMode, ncol = 4, scales = "free_y") +
    labs(title = paste("NES by Gene Set Origin —", hlab),
         subtitle = sprintf("%d subsamples per mode | Boxplot = descriptive only", n_iterations),
         x = NULL, y = "Normalized Enrichment Score (NES)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.4),
      panel.spacing = unit(0.8, "lines"),
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(size = 8, face = "bold", color = "black",
                                  angle = 30, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 9, face = "bold", color = "black", margin = margin(t = 6, b = 6)),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.ticks.y = element_line(color = "grey50", linewidth = 0.4),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5, margin = margin(b = 10))
    )

  queueable_print_plot(p)
  queueable_ggsave(file.path(out_dir, paste0("GSEA_Unique_Boxplot_", nm, ".pdf")), p, width = 8, height = 4)
  list(plot = p, stats = dunn_uniq)
}

# boxplot: combined
# res_uniq_box0 <- .make_uniq_box(unique_df, hop = 0)
# res_uniq_box1 <- .make_uniq_box(unique_df, hop = 1)

# boxplot: per-method
res_uniq_box0_anno <- .make_uniq_box(unique_df, hop = 0, method = "anno")
res_uniq_box1_anno <- .make_uniq_box(unique_df, hop = 1, method = "anno")
res_uniq_box0_refined <- .make_uniq_box(unique_df, hop = 0, method = "refined")
res_uniq_box1_refined <- .make_uniq_box(unique_df, hop = 1, method = "refined")
res_uniq_box0_chromatin <- .make_uniq_box(unique_df, hop = 0, method = "chromatin")
res_uniq_box1_chromatin <- .make_uniq_box(unique_df, hop = 1, method = "chromatin")
res_uniq_box0_chr_only  <- .make_uniq_box(unique_df, hop = 0, method = "chromatin_only")
res_uniq_box1_chr_only  <- .make_uniq_box(unique_df, hop = 1, method = "chromatin_only")

dunn_uniq_hop0_anno <- tryCatch(res_uniq_box0_anno$stats %>% mutate(hop = 0, method = "anno"), error = function(e) NULL)
dunn_uniq_hop1_anno <- tryCatch(res_uniq_box1_anno$stats %>% mutate(hop = 1, method = "anno"), error = function(e) NULL)
dunn_uniq_hop0_refined <- tryCatch(res_uniq_box0_refined$stats %>% mutate(hop = 0, method = "refined"), error = function(e) NULL)
dunn_uniq_hop1_refined <- tryCatch(res_uniq_box1_refined$stats %>% mutate(hop = 1, method = "refined"), error = function(e) NULL)
dunn_uniq_hop0_chromatin <- tryCatch(res_uniq_box0_chromatin$stats %>% mutate(hop = 0, method = "chromatin"), error = function(e) NULL)
dunn_uniq_hop1_chromatin <- tryCatch(res_uniq_box1_chromatin$stats %>% mutate(hop = 1, method = "chromatin"), error = function(e) NULL)
dunn_uniq_hop0_chr_only  <- tryCatch(res_uniq_box0_chr_only$stats  %>% mutate(hop = 0, method = "chromatin_only"), error = function(e) NULL)
dunn_uniq_hop1_chr_only  <- tryCatch(res_uniq_box1_chr_only$stats  %>% mutate(hop = 1, method = "chromatin_only"), error = function(e) NULL)
rm(res_uniq_box0_anno, res_uniq_box1_anno,
   res_uniq_box0_refined, res_uniq_box1_refined,
   res_uniq_box0_chromatin, res_uniq_box1_chromatin,
   res_uniq_box0_chr_only, res_uniq_box1_chr_only); gc()

# 汇总表
unique_mode_mean <- unique_df %>%
  filter(ID == "only_looplook") %>%
  group_by(Mode) %>%
  summarise(Mean_NES = mean(NES), Median_NES = median(NES), .groups = "drop") %>%
  arrange(Median_NES)
.safe_add_sheet(wb, "Unique_Mean_NES_by_Mode", as.data.frame(unique_mode_mean))
.safe_add_sheet(wb, "Unique_Raw_All", as.data.frame(unique_df))

# ── Unique GSEA: Effect Size forest (Only_looplook vs Only_ChIPseeker) ──
if (exists("gene_stat") && nrow(gene_stat) > 0) {
  uniq_es_list <- list()
  all_modes_uniq <- unique(as.character(unique_df$Mode))
  for (md_name in all_modes_uniq) {
    md_genes <- mode_genes[[md_name]]
    if (length(md_genes) == 0) next
    md <- modes[[which(sapply(modes, function(m) m$name == md_name))[1]]]
    loop_genes <- intersect(unique(toupper(md_genes)), gene_stat$gene)
    chip_g     <- intersect(get_chip_reference(md), gene_stat$gene)
    only_loop  <- intersect(setdiff(loop_genes, chip_g), gene_stat$gene)
    only_chip  <- intersect(setdiff(chip_g, loop_genes), gene_stat$gene)
    if (length(only_loop) >= 5 && length(only_chip) >= 5) {
      x <- gene_stat$signal_lfc[gene_stat$gene %in% only_loop]
      y <- gene_stat$signal_lfc[gene_stat$gene %in% only_chip]
      rb <- rank_biserial(x, y)
      uniq_es_list[[md_name]] <- data.frame(
        Mode = md_name, N_loop = length(only_loop), N_chip = length(only_chip),
        ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi, stringsAsFactors = FALSE)
    }
  }
  if (length(uniq_es_list) > 1) {
    uniq_es_df <- bind_rows(uniq_es_list) %>%
      mutate(Pipeline = factor(case_when(
        grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")),
        Hop = factor(ifelse(grepl("_hop1$", Mode), "hop1 (primary+expanded)", "hop0 (primary)"),
                     levels = c("hop0 (primary)", "hop1 (primary+expanded)")))

    for (hv in c("hop0", "hop1")) {
      hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
      es_sub <- uniq_es_df %>% filter(Hop == hv_label) %>% arrange(desc(ES))
      if (nrow(es_sub) < 2) next
      p_uniq_es <- ggplot(es_sub, aes(x = ES, y = reorder(Mode, ES), color = Pipeline)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(size = 1.5) +
        geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.3, linewidth = 0.3) +
        scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                       "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
        labs(title = sprintf("Unique GSEA: Effect Size %s", hv_label),
             subtitle = sprintf("Only_looplook vs Only_ChIPseeker | %d modes | Positive = looplook-only genes have stronger downregulation",
                                nrow(es_sub)),
             x = "Rank-Biserial Effect Size", y = NULL) +
        theme_classic() +
        theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
              plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
              legend.position = "bottom", axis.text.y = element_text(size = 8))
      queueable_ggsave(file.path(out_dir, sprintf("GSEA_Unique_EffectSize_Forest_%s.pdf", hv)), p_uniq_es, width = 8, height = 8)
    }

    .safe_add_sheet(wb, "Unique_EffectSize", as.data.frame(uniq_es_df))
  }
}

}  # end if(FALSE) Unique visualizations

# ── Unique 基因导出 + GO 分析 ─────────────────────────────────────────────
message("\n>>> Exporting unique gene sets and GO analysis...")

lfc_lookup <- setNames(diff_df$log2FoldChange, rownames(diff_df))

expr_mat <- looplook:::read_robust_general(expr_path, header = TRUE, row_name = 1,
               desc = "Expr", min_cols = 2)
expr_mat <- normalise_rownames(expr_mat, "Expression matrix")

# DMSO sample identification: metadata-driven (no fallback to first two columns)
meta_path <- file.path(cfg$data_base, cfg$meta_file)
if (file.exists(meta_path)) {
  meta_df <- read.table(meta_path, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE)
  required_meta_cols <- c("sample", "condition")
  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_df))
  if (length(missing_meta_cols) > 0L) {
    stop("Metadata missing columns: ", paste(missing_meta_cols, collapse = ", "))
  }
  meta_df$sample <- trimws(as.character(meta_df$sample))
  meta_df$condition <- trimws(as.character(meta_df$condition))
  if (anyNA(meta_df$sample) || any(!nzchar(meta_df$sample))) {
    stop("Metadata contains NA or blank sample IDs after trim.")
  }
  if (anyNA(meta_df$condition) || any(!nzchar(meta_df$condition))) {
    stop("Metadata contains NA or blank condition labels after trim.")
  }
  if (anyDuplicated(meta_df$sample)) {
    dup_samples <- unique(meta_df$sample[duplicated(meta_df$sample)])
    stop("Metadata contains duplicated sample IDs: ", paste(dup_samples, collapse = ", "))
  }
  dmso_cols <- meta_df$sample[
    tolower(meta_df$condition) %in% c("dmso", "vehicle", "control")]
  lost_controls <- setdiff(dmso_cols, colnames(expr_mat))
  if (length(lost_controls) > 0L) {
    stop("DMSO/control samples from metadata missing in expression matrix: ",
      paste(lost_controls, collapse = ", "))
  }
  dmso_cols <- intersect(dmso_cols, colnames(expr_mat))
  if (length(dmso_cols) == 0L) {
    stop("No DMSO/control samples from metadata were found in expression matrix.")
  }
  message("    DMSO samples from metadata: ", paste(dmso_cols, collapse = ", "))
} else {
  stop("Metadata file is required for formal analysis: ", meta_path)
}
expr_lookup <- setNames(rowMeans(expr_mat[, dmso_cols, drop = FALSE],
                                 na.rm = TRUE),
                        rownames(expr_mat))
expr_lookup <- expr_lookup[is.finite(expr_lookup)]

unique_genes_list <- list()

for (i in seq_along(modes)) {
  md <- modes[[i]]
  nm <- md$name
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) next

  loop_genes <- intersect(unique(toupper(genes)), all_genes)
  chip_g     <- intersect(get_chip_reference(md), all_genes)
  only_loop  <- setdiff(loop_genes, chip_g)
  inter_set  <- intersect(loop_genes, chip_g)
  only_chip  <- setdiff(chip_g, loop_genes)

  max_len <- max(length(only_loop), length(inter_set), length(only_chip))
  pad <- function(x, n) c(x, rep(NA, max(0, n - length(x))))
  unique_genes_list[[nm]] <- data.frame(
    only_looplook        = pad(only_loop, max_len),
    only_looplook_LFC    = pad(lfc_lookup[only_loop], max_len),
    only_looplook_TPM    = pad(expr_lookup[only_loop], max_len),
    intersection         = pad(inter_set, max_len),
    intersection_LFC     = pad(lfc_lookup[inter_set], max_len),
    intersection_TPM     = pad(expr_lookup[inter_set], max_len),
    only_ChIPseeker      = pad(only_chip, max_len),
    only_ChIPseeker_LFC  = pad(lfc_lookup[only_chip], max_len),
    only_ChIPseeker_TPM  = pad(expr_lookup[only_chip], max_len),
    stringsAsFactors = FALSE)

  message(sprintf("    [%2d/32] %-30s only_loop=%d inter=%d only_chip=%d",
    i, nm, length(only_loop), length(inter_set), length(only_chip)))
}

# 基因列表 → 单独 xlsx
wb_genes <- createWorkbook()
for (nm in names(unique_genes_list)) {
  addWorksheet(wb_genes, nm)
  writeData(wb_genes, nm, unique_genes_list[[nm]])
}
out_genes <- file.path(out_dir, sprintf("Unique_GeneLists_%diter.xlsx", n_iterations))
saveWorkbook(wb_genes, out_genes, overwrite = TRUE)
message("    Gene lists saved: ", out_genes)

# ── 16. Venn ────────────────────────────────────────────────────────────────
# Genes are intersected with expression matrix (all_genes) to match GSEA universe
message("\n>>> Venn analysis per mode (expression-universe matched)...")

venn_summary <- list()
venn_dir <- file.path(out_dir, "venn")
dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

for (md in modes) {
  nm    <- md$name
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) next
  assigned_loop <- unique(toupper(genes))
  assigned_chip <- unique(toupper(get_chip_reference(md)))
  loop_g <- intersect(assigned_loop, all_genes)
  chip_g <- intersect(assigned_chip, all_genes)
  venn_summary[[nm]] <- data.frame(
    Mode = nm,
    looplook_assigned = length(assigned_loop),
    ChIPseeker_assigned = length(assigned_chip),
    looplook_in_expr = length(loop_g),
    ChIPseeker_in_expr = length(chip_g),
    Intersection = length(intersect(loop_g, chip_g)),
    Only_looplook = length(setdiff(loop_g, chip_g)),
    Only_ChIPseeker = length(setdiff(chip_g, loop_g)),
    looplook_drop = length(setdiff(assigned_loop, loop_g)),
    ChIPseeker_drop = length(setdiff(assigned_chip, chip_g)),
    stringsAsFactors = FALSE)

  if (length(loop_g) == 0 && length(chip_g) == 0) {
    message(sprintf("    %-30s no genes in expression universe, Venn skipped", nm))
    next
  }
  p_venn <- ggvenn(list("looplook" = loop_g, "ChIPseeker" = chip_g),
    fill_color = c("#E64B35", "#4DBBD5"),
    stroke_size = 0.5, stroke_color = "white",
    set_name_size = 4.5, text_size = 4, fill_alpha = 0.6,
    show_percentage = TRUE) +
    labs(title = nm) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  hop_tag <- ifelse(grepl("_hop1", nm), "hop1", "hop0")
  nm_clean <- gsub("_hop1$", "", nm)
  queueable_ggsave(file.path(venn_dir, paste0("Venn_", hop_tag, "_", nm_clean, ".pdf")), p_venn, width = 5.5, height = 4.5)
}
venn_summary_df <- bind_rows(venn_summary)
.safe_add_sheet(wb, "Venn_Summary_All_Modes", venn_summary_df)

# Warn if any mode drops >5% of genes from expression universe
venn_drop_pct <- venn_summary_df %>%
  summarise(
    looplook_max_drop_pct = round(100 * max(looplook_drop / pmax(1, looplook_assigned), na.rm = TRUE), 1),
    ChIPseeker_max_drop_pct = round(100 * max(ChIPseeker_drop / pmax(1, ChIPseeker_assigned), na.rm = TRUE), 1),
    .groups = "drop")
if (venn_drop_pct$looplook_max_drop_pct > 5 || venn_drop_pct$ChIPseeker_max_drop_pct > 5) {
  warning(sprintf("Venn: some modes drop >5%% genes when intersecting with expression universe (max looplook=%s%%, ChIPseeker=%s%%). Check gene name matching.",
    venn_drop_pct$looplook_max_drop_pct, venn_drop_pct$ChIPseeker_max_drop_pct))
}

# LFC 导出
best_l <- intersect(unique(toupper(mode_genes[["refined_all_F"]])), all_genes)
best_c <- intersect(unique(toupper(chipseeker_genes)), all_genes)
deg_df <- diff_df %>% mutate(Gene = toupper(rownames(.))) %>%
  dplyr::select(Gene, LFC = log2FoldChange) %>% distinct(Gene, .keep_all = TRUE)
get_lfc <- function(gl, dg) dplyr::left_join(data.frame(Gene = gl, stringsAsFactors = FALSE), dg, by = "Gene")
df_inter <- get_lfc(intersect(best_l, best_c), deg_df)
df_loop  <- get_lfc(setdiff(best_l, best_c), deg_df)
df_chip  <- get_lfc(setdiff(best_c, best_l), deg_df)
max_n <- max(nrow(df_inter), nrow(df_loop), nrow(df_chip))
df_inter <- df_inter[1:max_n, , drop = FALSE]; colnames(df_inter) <- c("Shared_Gene", "Shared_LFC")
df_loop  <- df_loop[1:max_n, , drop = FALSE];  colnames(df_loop)  <- c("Only_Looplook_Gene", "Only_Looplook_LFC")
df_chip  <- df_chip[1:max_n, , drop = FALSE];  colnames(df_chip)  <- c("Only_ChIPseeker_Gene", "Only_ChIPseeker_LFC")
.safe_add_sheet(wb, "Venn_Genes_With_LFC", cbind(df_inter, df_loop, df_chip))

# ── 17. Anno vs Refined Venn ──────────────────────────────────────────────
message("\n>>> Anno vs Refined Venn per BaseMode...")

anno_refined_venn <- list()
arv_venn_dir <- file.path(out_dir, "venn_anno_refined")
dir.create(arv_venn_dir, recursive = TRUE, showWarnings = FALSE)

for (hop in c(0, 1)) {
  hop_suffix <- if (hop == 0) "" else "_hop1"
  for (base in c("promoter_F", "promoter_T", "all_F", "all_T")) {
    nm_anno    <- paste0("anno_",    base, hop_suffix)
    nm_refined <- paste0("refined_", base, hop_suffix)

    genes_anno    <- unique(toupper(mode_genes[[nm_anno]]))
    genes_refined <- unique(toupper(mode_genes[[nm_refined]]))

    if (length(genes_anno) == 0 && length(genes_refined) == 0) next

    shared   <- intersect(genes_anno, genes_refined)
    only_anno    <- setdiff(genes_anno, genes_refined)
    only_refined <- setdiff(genes_refined, genes_anno)

    key <- paste0("hop", hop, "_", base)
    max_n <- max(length(only_anno), length(only_refined), length(shared))
      anno_refined_venn[[key]] <- data.frame(
      only_anno         = pad(only_anno, max_n),
      only_anno_LFC     = pad(lfc_lookup[only_anno], max_n),
      only_anno_TPM     = pad(expr_lookup[only_anno], max_n),
      only_refined      = pad(only_refined, max_n),
      only_refined_LFC  = pad(lfc_lookup[only_refined], max_n),
      only_refined_TPM  = pad(expr_lookup[only_refined], max_n),
      shared            = pad(shared, max_n),
      shared_LFC        = pad(lfc_lookup[shared], max_n),
      shared_TPM        = pad(expr_lookup[shared], max_n),
      stringsAsFactors = FALSE)

    p <- ggvenn(list("anno" = genes_anno, "refined" = genes_refined),
      fill_color = c("#1F77B4", "#9467BD"),
      stroke_size = 0.5, stroke_color = "white",
      set_name_size = 4.5, text_size = 4, fill_alpha = 0.6,
      show_percentage = TRUE) +
      labs(title = key) +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
    queueable_ggsave(file.path(arv_venn_dir, paste0("Venn_", key, ".pdf")), p, width = 5.5, height = 4.5)

    message(sprintf("    %-20s anno=%-5d refined=%-5d shared=%-5d only_anno=%-5d only_refined=%-5d",
      key, length(genes_anno), length(genes_refined), length(shared),
      length(only_anno), length(only_refined)))
  }
}

# 写入 anno vs refined 基因列表 → 单独 xlsx
wb_ar <- createWorkbook()
for (nm in names(anno_refined_venn)) {
  addWorksheet(wb_ar, nm)
  writeData(wb_ar, nm, anno_refined_venn[[nm]])
}
out_ar <- file.path(out_dir, sprintf("AnnoRefined_GeneLists_%diter.xlsx", n_iterations))
saveWorkbook(wb_ar, out_ar, overwrite = TRUE)
message("    Anno vs Refined gene lists saved: ", out_ar)

# ── 保存 ───────────────────────────────────────────────────────────────────
# Dunn's 结果
.safe_add_sheet(wb, "Delta_4Pipeline_hop0", delta_stats_hop0)
.safe_add_sheet(wb, "Delta_4Pipeline_hop1", delta_stats_hop1)

# Dunn results skipped (data inside Unique visualization block)
if (FALSE) {
.write_dunn_sheet(wb, "Dunn_Uniq_hop1_chr_only", dunn_uniq_hop1_chr_only)
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("\n>>> Done. Output: ", out_dir)
message("    ", out_xlsx)
message("\nTip: keep cache_policy=\"resume\" after interruption; use \"rebuild\" only when inputs or core logic changed.")


  message(sprintf("\n>>> Sample size %d done. Output: %s", sz, out_dir))
}  # end for sz in sample_sizes
message("\n>>> All sample size gradients complete. Base dir: ", base_out_dir)
message(sprintf("    Status records so far: %d (will be finalized at end of script)", length(gsea_status_log)))

primary_size_key <- as.character(primary_sample_size)
if (!primary_size_key %in% names(paired_nes_by_size) ||
    !primary_size_key %in% names(delta_stats_by_size)) {
  stop("Primary sample size summary was not generated for size=", primary_sample_size)
}
paired_nes_summary <- paired_nes_by_size[[primary_size_key]]
delta_stats_all <- delta_stats_by_size[[primary_size_key]]
message("    Downstream ranking/summary primary sample size: ", primary_sample_size)

# ══════════════════════════════════════════════════════════════════════════════
# 15b. GO Enrichment Analysis (outside sample-size loop; GO does not depend on
#     random sampling — it tests functional overlap of full gene sets)
# ══════════════════════════════════════════════════════════════════════════════
message("\n>>> GO Enrichment Analysis on unique gene sets...")

go_dir <- file.path(base_out_dir, "go_enrichment")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
go_cache_rds <- file.path(go_dir, "go_results_cache.rds")

# Helper: run GO enrichment from pre-computed entrez IDs
run_go_from_entrez <- function(entrez, mode, gene_set_label, ontologies = c("BP", "MF", "CC")) {
  empty_go_result <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Ontology = character(),
    Mode = character(), GeneSet = character(), stringsAsFactors = FALSE)
  if (length(entrez) < 5) {
    status_df <- data.frame(SampleSize = NA_integer_, Module = "go", Mode = mode,
      Iteration = NA_integer_, GeneSet = gene_set_label, Ontology = NA_character_,
      N_Entrez = length(entrez), Status = "insufficient_input",
      ExpectedTerms = NA_character_, ReturnedTerms = NA_character_,
      MissingTerms = NA_character_, Error = NA_character_, stringsAsFactors = FALSE)
    return(list(result = empty_go_result, status = status_df))
  }
  go_calls <- lapply(ontologies, function(ont) {
    safe_go_once(entrez, ont, NA_integer_, mode, gene_set_label, gouniv)
  })
  status_df <- bind_rows(lapply(go_calls, function(x) x$status))
  result_df <- bind_rows(lapply(go_calls, function(x) {
    if (!is.null(x$result) && nrow(x$result) > 0) {
      x$result %>%
        dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID) %>%
        mutate(Ontology = x$status$Ontology[1], Mode = mode, GeneSet = gene_set_label)
    } else {
      data.frame(ID = character(), Description = character(),
        GeneRatio = character(), BgRatio = character(), pvalue = numeric(),
        p.adjust = numeric(), qvalue = numeric(), geneID = character(),
        Ontology = character(), Mode = character(), GeneSet = character(),
        stringsAsFactors = FALSE)
    }
  }))
  list(result = result_df, status = status_df)
}

go_signature_current <- list(
  ModeGenesMD5 = md5_object(mode_genes),
  ChIPReferenceMD5 = md5_object(chipseeker_genes),
  TargetMD5 = preliminary_signature$TargetMD5,
  GounivMD5 = md5_object(gouniv),
  OrgDbVersion = as.character(packageVersion("org.Hs.eg.db")),
  ChIPseekerVersion = as.character(packageVersion("ChIPseeker")),
  ClusterProfilerVersion = as.character(packageVersion("clusterProfiler")),
  Parameters = list(Ontology = c("BP", "MF", "CC"), PAdjustMethod = "BH",
    PvalueCutoff = 0.1, QvalueCutoff = 0.2, MinInputEntrez = 5))

if (file.exists(go_cache_rds)) {
  message("    Loading cached GO results from ", go_cache_rds)
  go_cached <- readRDS(go_cache_rds)
  if (is.null(go_cached$signature)) {
    stop("Unsigned GO cache detected. Delete ", go_cache_rds, " and re-run.")
  }
  if (!identical(go_cached$signature, go_signature_current)) {
    stop("GO cache signature mismatch. Delete ", go_cache_rds, " and re-run.")
  }
  go_results <- go_cached$go_results
  go_summary <- go_cached$go_summary
  if (!is.null(go_cached$go_status) && nrow(go_cached$go_status) > 0L) {
    cached_go_status <- go_cached$go_status
    cached_go_status$CacheState <- "cached_validated"
    gsea_status_log <- c(gsea_status_log,
      split(cached_go_status, seq_len(nrow(cached_go_status))))
    message("    GO status restored from cache: ", nrow(cached_go_status), " records")
  }
} else {
  message("    Batch-converting symbols to Entrez IDs...")
  chip_g_all <- intersect(unique(toupper(chipseeker_genes)), all_genes)
  all_go_genes <- unique(unlist(lapply(mode_genes, function(g) unique(toupper(g)))))
  all_go_genes <- unique(c(all_go_genes, chip_g_all))
  entrez_batch <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = all_go_genes,
      columns = "ENTREZID", keytype = "SYMBOL"))
  entrez_mapping <- entrez_batch %>%
    transmute(
      SYMBOL = toupper(SYMBOL),
      ENTREZID = as.character(ENTREZID)) %>%
    filter(!is.na(SYMBOL), SYMBOL != "", !is.na(ENTREZID), ENTREZID != "") %>%
    distinct()
  message(sprintf("    Mapped %d symbols to %d unique Entrez IDs (long table preserved)",
    n_distinct(entrez_mapping$SYMBOL), n_distinct(entrez_mapping$ENTREZID)))

  fast_symbols_to_entrez <- function(genes) {
    if (length(genes) == 0) return(character(0))
    genes <- unique(toupper(genes))
    entrez_mapping %>%
      filter(SYMBOL %in% genes) %>%
      pull(ENTREZID) %>%
      unique()
  }

  go_results <- list()
  go_summary <- list()
  go_status_this_run <- list()
  for (i in seq_along(modes)) {
    md <- modes[[i]]
    nm <- md$name
    genes <- mode_genes[[nm]]
    if (length(genes) == 0) next

    loop_genes <- intersect(unique(toupper(genes)), all_genes)
    chip_g_per_mode <- get_chip_reference(md)
    only_loop  <- setdiff(loop_genes, chip_g_per_mode)
    inter_set  <- intersect(loop_genes, chip_g_per_mode)

    entrez_only <- fast_symbols_to_entrez(only_loop)
    entrez_inter <- fast_symbols_to_entrez(inter_set)

    go_only_out <- run_go_from_entrez(entrez_only, nm, "only_looplook")
    go_inter_out <- run_go_from_entrez(entrez_inter, nm, "intersection")
    go_status_this_run <- c(go_status_this_run,
      split(go_only_out$status, seq_len(nrow(go_only_out$status))),
      split(go_inter_out$status, seq_len(nrow(go_inter_out$status))))
    gsea_status_log <- c(gsea_status_log,
      split(go_only_out$status, seq_len(nrow(go_only_out$status))),
      split(go_inter_out$status, seq_len(nrow(go_inter_out$status))))
    go_only <- go_only_out$result
    go_inter <- go_inter_out$result

    if (nrow(go_only) > 0) {
      go_results[[paste0(nm, "_only")]] <- go_only
    }
    if (nrow(go_inter) > 0) {
      go_results[[paste0(nm, "_inter")]] <- go_inter
    }

    n_sig <- function(df) sum(df$p.adjust < 0.05, na.rm = TRUE)
    go_summary[[nm]] <- data.frame(
      Mode = nm,
      Only_BP = n_sig(filter(go_only, Ontology == "BP")),
      Only_MF = n_sig(filter(go_only, Ontology == "MF")),
      Only_CC = n_sig(filter(go_only, Ontology == "CC")),
      Inter_BP = n_sig(filter(go_inter, Ontology == "BP")),
      Inter_MF = n_sig(filter(go_inter, Ontology == "MF")),
      Inter_CC = n_sig(filter(go_inter, Ontology == "CC")),
      stringsAsFactors = FALSE)

    message(sprintf("    [%2d/32] %-30s only_loop=%d(entrez=%d) GO_BP_sig=%d | inter=%d(entrez=%d) GO_BP_sig=%d",
      i, nm, length(only_loop), length(entrez_only),
      go_summary[[nm]]$Only_BP,
      length(inter_set), length(entrez_inter),
      go_summary[[nm]]$Inter_BP))
  }
  go_status_filtered <- bind_rows(go_status_this_run)
  saveRDS(list(go_results = go_results, go_summary = go_summary,
    go_status = go_status_filtered,
    signature = go_signature_current), go_cache_rds)
  message("    GO results cached to ", go_cache_rds)
}

# Export GO results to Excel
wb_go <- createWorkbook()
go_all <- bind_rows(go_results)
if (nrow(go_all) > 0) {
  addWorksheet(wb_go, "All_GO_Results")
  writeData(wb_go, "All_GO_Results", as.data.frame(go_all))

  for (nm in names(go_results)) {
    df <- go_results[[nm]]
    if (nrow(df) > 0) {
      sheet_name <- substr(nm, 1, 31)
      addWorksheet(wb_go, sheet_name)
      writeData(wb_go, sheet_name, as.data.frame(df %>% arrange(p.adjust) %>% head(50)))
    }
  }
}
addWorksheet(wb_go, "GO_Summary")
writeData(wb_go, "GO_Summary", as.data.frame(bind_rows(go_summary)))
out_go <- file.path(go_dir, sprintf("GO_Enrichment_%diter.xlsx", n_iterations))
saveWorkbook(wb_go, out_go, overwrite = TRUE)
message("    GO enrichment results saved: ", out_go)

# ── GO Summary Visualization ──
go_summary_df <- bind_rows(go_summary) %>%
  mutate(
    Pipeline = case_when(
      grepl("^anno_", Mode) ~ "Basic",
      grepl("^refined_", Mode) ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined",
      grepl("^chrom_", Mode) ~ "I-refined",
      TRUE ~ "Other"),
    Total_Only = Only_BP + Only_MF + Only_CC,
    Total_Inter = Inter_BP + Inter_MF + Inter_CC)

p_go_compare <- ggplot(go_summary_df, aes(x = Total_Inter, y = Total_Only, color = Pipeline)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                 "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
  labs(title = "GO Term Enrichment: looplook-only vs Shared Genes",
       subtitle = "N GO terms (affected by gene-set size, not direct evidence of richer function)",
       x = "Significant GO terms (shared with ChIPseeker)",
       y = "Significant GO terms (looplook-only)",
       color = "Pipeline") +
  theme_minimal()
queueable_ggsave(file.path(go_dir, "GO_Only_vs_Shared_Scatter.pdf"), p_go_compare, width = 8, height = 6)

p_go_bar <- go_summary_df %>%
  tidyr::pivot_longer(c(Only_BP, Only_MF, Only_CC, Inter_BP, Inter_MF, Inter_CC),
    names_to = "Category", values_to = "N_Terms") %>%
  mutate(Set = ifelse(grepl("^Only_", Category), "looplook-only", "Shared"),
         Ontology = gsub("^(Only_|Inter_)", "", Category)) %>%
  ggplot(aes(x = Ontology, y = N_Terms, fill = Set)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ Pipeline, scales = "free_y") +
  scale_fill_manual(values = c("looplook-only" = "#E64B35", "Shared" = "#4DBBD5")) +
  labs(title = "GO Term Count by Pipeline and Gene Set",
       x = "Ontology", y = "N significant GO terms") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
queueable_ggsave(file.path(go_dir, "GO_Pipeline_Comparison.pdf"), p_go_bar, width = 10, height = 6)

if (nrow(go_all) > 0) {
  go_top <- go_all %>%
    filter(Ontology == "BP", p.adjust < 0.05) %>%
    group_by(Mode, GeneSet) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(Pipeline = case_when(
      grepl("^anno_", Mode) ~ "Basic",
      grepl("^refined_", Mode) ~ "E-refined",
      TRUE ~ "Other"),
      ModeSet = paste0(Mode, " (", GeneSet, ")"))
  if (nrow(go_top) > 0) {
    go_top <- go_top %>% mutate(
      size_metric = pmin(-log10(pmax(p.adjust, 1e-300)), 50))
    p_go_dot <- ggplot(go_top, aes(x = ModeSet, y = Description, size = size_metric, color = Pipeline)) +
      geom_point() +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD", "Other" = "grey50")) +
      labs(title = "Top GO BP Terms per Mode (looplook-only and shared genes)",
           x = "Mode (GeneSet)", y = "GO Term", size = "-log10(p.adjust)") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
            legend.position = "bottom")
    queueable_ggsave(file.path(go_dir, "GO_BP_TopTerms_Dotplot.pdf"), p_go_dot, width = 10, height = 6)
  }
}

message("    GO enrichment plots saved to: ", go_dir)

# ══════════════════════════════════════════════════════════════════════════════
# Peak-Gene Distance Analysis: looplook (3D) vs ChIPseeker (1D linear distance)
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Peak-Gene Distance Analysis ====")
# Rationale: The key biological differentiator between looplook (3D chromatin
# loop-based gene assignment) and ChIPseeker (1D nearest-TSS annotation) is the
# linear distance between peaks and their assigned genes. Looplook can assign
# distal enhancers to target genes via chromatin looping, while ChIPseeker is
# constrained to the nearest gene in linear space. This analysis quantifies that
# difference and characterizes the linear genomic reach of the 3D approach.

dist_dir <- file.path(base_out_dir, "distance_analysis")
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

# Extract peak coordinates and distanceToTSS from ChIPseeker annotation
# ChIPseeker's annotatePeak() already computes distanceToTSS for nearest-gene
peak_dist_df <- as.data.frame(peak_anno)
peak_dist_df$peak_id <- paste0(peak_dist_df$seqnames, ":", peak_dist_df$start, "-", peak_dist_df$end)
peak_dist_df$SYMBOL <- toupper(peak_dist_df$SYMBOL)
message(sprintf("    Peaks with ChIPseeker annotation: %d", nrow(peak_dist_df)))
message(sprintf("    distanceToTSS range: %.0f - %.0f bp (median %.0f bp)",
  min(peak_dist_df$distanceToTSS, na.rm = TRUE),
  max(peak_dist_df$distanceToTSS, na.rm = TRUE),
  median(peak_dist_df$distanceToTSS, na.rm = TRUE)))

# For looplook assignments, we need to compute the distance from each peak
# to the TSS of each gene it is assigned to. Build TSS coordinates from TxDb.
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
tss_list <- GenomicFeatures::promoters(tx_by_gene, upstream = 0, downstream = 1)
tss_gr <- unlist(tss_list, use.names = FALSE)
tss_gr$ENTREZID <- rep(names(tss_list), lengths(tss_list))
symbol_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(tss_gr$ENTREZID), keytype = "ENTREZID", columns = "SYMBOL"))
tss_gr$SYMBOL <- toupper(symbol_map$SYMBOL[match(tss_gr$ENTREZID, symbol_map$ENTREZID)])
tss_gr <- tss_gr[!is.na(tss_gr$SYMBOL)]
message(sprintf("    Unified TSS: %d transcript TSS for %d genes", length(tss_gr), length(unique(tss_gr$SYMBOL))))

distance_to_gene_tss <- function(peak_gr, gene, tss_gr_ref) {
  gene <- toupper(gene); tss_hit <- tss_gr_ref[tss_gr_ref$SYMBOL == gene]
  if (length(tss_hit) == 0L) return(NA_real_)
  d <- as.numeric(GenomicRanges::distance(peak_gr, tss_hit)); d <- d[is.finite(d) & !is.na(d)]
  if (length(d) == 0L) return(NA_real_) else min(d)
}

# Build peak-to-gene distance data: For each mode, compute distances from
# BRD4 peaks to looplook-assigned genes using unified TSS catalogue.
# ChIPseeker distances also recalculated with same TSS for fair comparison.
message("    Computing peak-gene distances (parallel)...")
message(sprintf("    Unified TSS benchmark: both looplook and ChIPseeker use transcript promoters"))# Each row = one BRD4 peak to one looplook-assigned gene.
make_loop_dist_df <- function(md, tss_gr_local) {
  bed_info <- md$ann$target_annotation
  if (is.null(bed_info) || nrow(bed_info) == 0) return(data.frame())

  peak_gr <- GRanges(
    seqnames = bed_info$seqnames,
    ranges = IRanges(start = bed_info$start, end = bed_info$end))

  peak_gene_pairs <- get_mode_peak_gene_pairs(md)
  if (nrow(peak_gene_pairs) == 0L) return(data.frame())

  peak_index <- match(
    peak_gene_pairs$Input_ID,
    as.character(bed_info$input_id)
  )

  valid_pair <- !is.na(peak_index)
  peak_gene_pairs <- peak_gene_pairs[valid_pair, , drop = FALSE]
  peak_index <- peak_index[valid_pair]

  res <- lapply(seq_len(nrow(peak_gene_pairs)), function(j) {
    g <- peak_gene_pairs$Gene[j]
    pidx <- peak_index[j]
    d <- distance_to_gene_tss(peak_gr[pidx], g, tss_gr_local)
    if (is.na(d) || !is.finite(d)) return(NULL)

    data.frame(
      Peak_ID = paste0(
        as.character(seqnames(peak_gr[pidx])),
        ":",
        start(peak_gr[pidx]),
        "-",
        end(peak_gr[pidx])
      ),
      Gene = g,
      Distance = d,
      Source = "looplook",
      AssignmentScope = if (is_hop1_mode(md)) {
        "primary_plus_expanded"
      } else {
        "primary"
      },
      stringsAsFactors = FALSE
    )
  })
  bind_rows(res) %>% filter(is.finite(Distance), Distance >= 0) %>%
    distinct(Peak_ID, Gene, .keep_all = TRUE)
}


# ChIPseeker distances recalculated with unified TSS (not distanceToTSS)
peak_gr <- GenomicRanges::GRanges(seqnames = peak_dist_df$seqnames,
  ranges = IRanges(start = peak_dist_df$start, end = peak_dist_df$end))
chip_dists <- unlist(lapply(seq_len(nrow(peak_dist_df)), function(j) {
  g <- peak_dist_df$SYMBOL[j]; if (is.na(g) || g == "") return(NA_real_)
  distance_to_gene_tss(peak_gr[j], g, tss_gr)
}))
peak_dist_df$Unified_Distance <- chip_dists
message(sprintf("    ChIPseeker unified TSS distance: median = %.0f bp",
  median(peak_dist_df$Unified_Distance, na.rm = TRUE)))
# Pre-filter modes for distance analysis
dist_mode_idx <- which(!sapply(modes, `[[`, "near") &
                        sapply(mode_genes, function(g) length(g) > 0))

# Parallel: build peak-gene pair distance tables
dist_cache_rds <- file.path(base_out_dir, ".loop_dist_cache.rds")
if (file.exists(dist_cache_rds) && !force_recompute) {
  message("    Loading cached peak-gene distances...")
  loop_dist_tables <- readRDS(dist_cache_rds)
} else {
  loop_dist_tables <- parallel::mclapply(dist_mode_idx, function(i) {
    md <- modes[[i]]
    make_loop_dist_df(md, tss_gr)
  }, mc.cores = min(8L, n_cores), mc.set.seed = FALSE, mc.preschedule = TRUE, mc.cleanup = TRUE)
  saveRDS(loop_dist_tables, dist_cache_rds)
  message("    Peak-gene distances cached")
}
close_batch_devices()
names(loop_dist_tables) <- sapply(modes[dist_mode_idx], `[[`, "name")

# ChIPseeker nearest-gene distance table (same peak-gene pair format)
chip_dist_df <- peak_dist_df %>%
  dplyr::mutate(
    Peak_ID  = peak_id,
    Gene     = toupper(SYMBOL),
    Distance = abs(Unified_Distance),
    Source   = "ChIPseeker",
    AssignmentScope = "nearest_gene") %>%
  dplyr::filter(!is.na(Gene), Gene != "", is.finite(Distance), Distance >= 0) %>%
  dplyr::select(Peak_ID, Gene, Distance, Source, AssignmentScope)

# Compile: classify each gene into Only_looplook / Intersection / Only_ChIPseeker
dist_results <- list()
for (k in seq_along(dist_mode_idx)) {
  i  <- dist_mode_idx[k]
  md <- modes[[i]]
  nm <- md$name
  chip_g <- intersect(get_chip_reference(md), all_genes)
  loop_genes <- intersect(unique(toupper(mode_genes[[nm]])), all_genes)

  only_loop <- setdiff(loop_genes, chip_g)
  inter_set <- intersect(loop_genes, chip_g)
  only_chip <- setdiff(chip_g, loop_genes)

  # Looplook peak-gene pairs, classified by gene membership
  loop_df <- loop_dist_tables[[nm]]
  dist_parts <- list()
  if (!is.null(loop_df) && nrow(loop_df) > 0) {
    loop_class <- loop_df %>%
      mutate(
        GeneSet = case_when(
          Gene %in% only_loop ~ "Only_looplook",
          Gene %in% inter_set ~ "Intersection (looplook)",
          TRUE ~ NA_character_),
        Mode = nm) %>%
      filter(!is.na(GeneSet))
    if (nrow(loop_class) > 0) dist_parts[[1]] <- loop_class
  }

  # ChIPseeker pairs: classify genes
  chip_for_mode <- chip_dist_df %>% filter(Gene %in% c(inter_set, only_chip))
  if (nrow(chip_for_mode) > 0) {
    chip_class <- chip_for_mode %>%
      mutate(
        GeneSet = case_when(
          Gene %in% inter_set ~ "Intersection (ChIPseeker)",
          Gene %in% only_chip  ~ "Only_ChIPseeker",
          TRUE ~ NA_character_),
        Mode = nm) %>%
      filter(!is.na(GeneSet))
    if (nrow(chip_class) > 0) dist_parts[[length(dist_parts) + 1]] <- chip_class
  }

  if (length(dist_parts) > 0) {
    dist_results[[nm]] <- bind_rows(dist_parts) %>%
      dplyr::mutate(
        AssignmentScope = dplyr::if_else(
          Source == "ChIPseeker", "nearest_gene", AssignmentScope
        )
      ) %>%
      dplyr::select(Mode, Peak_ID, Gene, GeneSet, Distance, Source, AssignmentScope)
    message(sprintf("    [%2d] %-30s only_loop=%d inter=%d only_chip=%d  loop_pairs=%d chip_pairs=%d",
      i, nm, length(only_loop), length(inter_set), length(only_chip),
      if (!is.null(loop_df)) nrow(loop_df) else 0,
      nrow(chip_for_mode)))
  }
}

dist_df <- bind_rows(dist_results)
if (nrow(dist_df) == 0L || !"Source" %in% colnames(dist_df)) {
  stop("Distance analysis produced no peak-gene pairs.", call. = FALSE)
}
expected_distance_sources <- c("looplook", "ChIPseeker")
unexpected_distance_sources <- setdiff(
  unique(as.character(dist_df$Source)), expected_distance_sources
)
if (length(unexpected_distance_sources) > 0L) {
  stop(
    "Unexpected distance Source value(s): ",
    paste(unexpected_distance_sources, collapse = ", "),
    call. = FALSE
  )
}
dist_df$GeneSet <- factor(dist_df$GeneSet,
  levels = c("Only_looplook", "Intersection (looplook)", "Intersection (ChIPseeker)", "Only_ChIPseeker"))

# Retain all finite cis linear distances. Long-range assignments are biologically
# relevant for a 3D regulatory benchmark and must not be silently discarded.
dist_df_clean <- dist_df %>%
  filter(is.finite(Distance), Distance >= 0) %>%
  mutate(
    Dist_kb = (Distance + 1) / 1000,
    Distance_QC = ifelse(Distance >= 5e6, ">=5Mb_long_range", "<5Mb")
  )

if (nrow(dist_df_clean) == 0) {
  warning("Distance analysis skipped: no valid peak-gene pairs.")
} else {
  message(sprintf("    Total valid peak-gene pairs: %d (range %.1f-%.0f kb; >=5Mb=%d)",
    nrow(dist_df_clean), min(dist_df_clean$Dist_kb), max(dist_df_clean$Dist_kb),
    sum(dist_df_clean$Distance_QC == ">=5Mb_long_range")))

# ── Figure A: Gene-set density comparison ──
dist_colors <- c("Only_looplook" = "#E64B35", "Intersection (looplook)" = "#00A087", "Intersection (ChIPseeker)" = "#7CAE00",
                 "Only_ChIPseeker" = "#4DBBD5")

p_dist_density <- ggplot(dist_df_clean, aes(x = Dist_kb, fill = GeneSet)) +
  geom_density(alpha = 0.7, color = NA) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = dist_colors) +
  facet_wrap(~ GeneSet, ncol = 1, scales = "free_y") +
  labs(title = "Peak-Gene Pair Distance by Gene Set Origin",
       subtitle = paste("All peak-gene pairs use the same transcript-level",
                         "TSS catalogue and nearest-TSS distance rule |",
                         sprintf("%d modes", length(unique(dist_df_clean$Mode)))),
       x = "Distance from peak to TSS (kb, log scale)", y = "Density") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
queueable_ggsave(file.path(dist_dir, "Distance_Density_4Groups_SourceSplit.pdf"), p_dist_density, width = 7, height = 6)

# ── Figure B: Pipeline-stratified boxplot with DistanceGroup ──
dist_df_clean <- dist_df_clean %>%
  mutate(
    DistanceGroup = factor(case_when(
      GeneSet == "Only_looplook" & Source == "looplook"       ~ "Only_looplook | looplook",
      grepl("Intersection", GeneSet) & Source == "looplook"   ~ "Intersection | looplook",
      grepl("Intersection", GeneSet) & Source == "ChIPseeker" ~ "Intersection | ChIPseeker",
      GeneSet == "Only_ChIPseeker" & Source == "ChIPseeker"   ~ "Only_ChIPseeker | ChIPseeker",
      TRUE ~ NA_character_),
      levels = c("Only_looplook | looplook", "Intersection | looplook",
                 "Intersection | ChIPseeker", "Only_ChIPseeker | ChIPseeker")),
    Pipeline = factor(case_when(
      grepl("^anno_", Mode)        ~ "Basic",
      grepl("^refined_", Mode)     ~ "E-refined",
      grepl("^chrom_only_", Mode)  ~ "C-refined",
      grepl("^chrom_", Mode)       ~ "I-refined",
      TRUE ~ "Other"),
      levels = c("Basic", "E-refined", "I-refined", "C-refined"))) %>%
  filter(!is.na(DistanceGroup))

dist_group_colors <- c(
  "Only_looplook | looplook"      = "#E64B35",
  "Intersection | looplook"       = "#00A087",
  "Intersection | ChIPseeker"     = "#7CAE00",
  "Only_ChIPseeker | ChIPseeker"  = "#4DBBD5")

# Mode-level summary for primary statistics (pair-level distribution is descriptive only)
dist_mode_summary_early <- dist_df_clean %>%
  group_by(Mode, DistanceGroup) %>%
  summarise(Med_Dist_kb = median(Dist_kb, na.rm = TRUE), N = n(), .groups = "drop")

p_dist_box <- ggplot(dist_df_clean, aes(x = DistanceGroup, y = Dist_kb)) +
  geom_boxplot(aes(fill = DistanceGroup), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = dist_group_colors) +
  facet_wrap(~ Pipeline, ncol = 4) +
  labs(title = "Peak-Gene Pair Distance by Assignment Source",
       subtitle = "Pair-level distribution (descriptive) | Unified transcript-level TSS catalogue",
       x = NULL, y = "Distance (kb, log scale)") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.4),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))

queueable_ggsave(file.path(dist_dir, "Distance_Boxplot_byPipeline_SourceSplit.pdf"), p_dist_box, width = 14, height = 5)

  dist_df_clean <- dist_df_clean %>%
    mutate(HopDist = ifelse(grepl("_hop1$", Mode), "hop1 (primary+expanded)", "hop0 (primary)"))
  for (hv in c("hop0 (primary)", "hop1 (primary+expanded)")) {
    hs <- ifelse(hv == "hop0 (primary)", "hop0", "hop1")
    dh <- dist_df_clean %>% filter(HopDist == hv)
    if (nrow(dh) > 10) {
      pdh <- ggplot(dh, aes(x = DistanceGroup, y = Dist_kb, fill = DistanceGroup)) +
        geom_boxplot(outlier.shape = 16, outlier.size = 0.3, outlier.alpha = 0.15,
                     coef = 1.5, alpha = 0.6, width = 0.6) +
        scale_y_log10(labels = scales::comma) + scale_fill_manual(values = dist_group_colors) +
        facet_wrap(~ Pipeline, ncol = 4) +
        labs(title = paste("Peak-Gene Pair Distance —", hv),
             subtitle = paste(n_distinct(dh$Mode), "modes"),
             x = NULL, y = "Distance (kb, log scale)") +
        theme_minimal() + theme(legend.position = "none",
          panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.4),
          strip.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
      queueable_ggsave(file.path(dist_dir, paste0("Distance_Boxplot_", hs, ".pdf")), pdh, width = 14, height = 5)
    }
  }

# ── Statistical tests on mode-level medians ───────────────────────────────
dist_mode_summary <- dist_df_clean %>%
  group_by(Mode, DistanceGroup) %>%
  summarise(Med_Dist_kb = median(Dist_kb, na.rm = TRUE), N = n(), .groups = "drop")

dist_kw <- NULL
dist_friedman <- NULL
dist_paired_tests <- data.frame()
if (dplyr::n_distinct(dist_mode_summary$DistanceGroup) >= 2) {
  # Retained as an exploratory omnibus test. The primary sensitivity analysis
  # below respects repeated measurements from the same benchmark mode.
  dist_kw <- kruskal.test(Med_Dist_kb ~ DistanceGroup, data = dist_mode_summary)
  message(sprintf("    Exploratory Kruskal-Wallis (mode-level): chi2=%.1f, p=%.2e",
    dist_kw$statistic, dist_kw$p.value))

  dist_wide <- dist_mode_summary %>%
    dplyr::select(Mode, DistanceGroup, Med_Dist_kb) %>%
    tidyr::pivot_wider(names_from = DistanceGroup, values_from = Med_Dist_kb)
  group_cols <- setdiff(colnames(dist_wide), "Mode")
  complete_wide <- dist_wide[stats::complete.cases(dist_wide[, group_cols, drop = FALSE]), , drop = FALSE]

  if (length(group_cols) >= 3L && nrow(complete_wide) >= 3L) {
    friedman_long <- complete_wide %>%
      tidyr::pivot_longer(-Mode, names_to = "DistanceGroup", values_to = "Med_Dist_kb")
    dist_friedman <- friedman.test(Med_Dist_kb ~ DistanceGroup | Mode, data = friedman_long)
    message(sprintf("    Friedman repeated-mode test: chi2=%.1f, p=%.2e (complete modes=%d)",
      dist_friedman$statistic, dist_friedman$p.value, nrow(complete_wide)))
  }

  if (length(group_cols) >= 2L) {
    comparisons <- utils::combn(group_cols, 2L, simplify = FALSE)
    paired_rows <- lapply(comparisons, function(z) {
      pair_df <- dist_wide[, c("Mode", z), drop = FALSE]
      pair_df <- pair_df[stats::complete.cases(pair_df), , drop = FALSE]
      if (nrow(pair_df) < 3L) return(NULL)
      wt <- suppressWarnings(stats::wilcox.test(
        pair_df[[z[1L]]], pair_df[[z[2L]]],
        paired = TRUE, exact = FALSE
      ))
      data.frame(
        Group1 = z[1L], Group2 = z[2L], N_PairedModes = nrow(pair_df),
        Median_Delta_kb = median(pair_df[[z[1L]]] - pair_df[[z[2L]]], na.rm = TRUE),
        P = wt$p.value, stringsAsFactors = FALSE
      )
    })
    dist_paired_tests <- dplyr::bind_rows(paired_rows)
    if (nrow(dist_paired_tests) > 0L) {
      dist_paired_tests$P_BH <- p.adjust(dist_paired_tests$P, method = "BH")
    }
  }
} else {
  warning("Distance tests skipped: fewer than two distance groups with data.")
}

# Export
wb_dist <- createWorkbook()
addWorksheet(wb_dist, "Peak_Gene_Pairs_Raw")
writeData(wb_dist, "Peak_Gene_Pairs_Raw", as.data.frame(dist_df_clean))
dist_summary <- dist_df_clean %>%
  group_by(Mode, GeneSet, Source) %>%
  summarise(N = n(), Med_Dist_kb = median(Dist_kb), Mean_Dist_kb = mean(Dist_kb),
            .groups = "drop")
addWorksheet(wb_dist, "Distance_Summary")
writeData(wb_dist, "Distance_Summary", as.data.frame(dist_summary))
if (!is.null(dist_kw)) {
  addWorksheet(wb_dist, "Exploratory_KW")
  writeData(wb_dist, "Exploratory_KW", data.frame(
    Statistic = unname(dist_kw$statistic), P = dist_kw$p.value
  ))
}
if (!is.null(dist_friedman)) {
  addWorksheet(wb_dist, "Friedman_RepeatedMode")
  writeData(wb_dist, "Friedman_RepeatedMode", data.frame(
    Statistic = unname(dist_friedman$statistic), P = dist_friedman$p.value
  ))
}
if (nrow(dist_paired_tests) > 0L) {
  addWorksheet(wb_dist, "Paired_Mode_Wilcoxon")
  writeData(wb_dist, "Paired_Mode_Wilcoxon", as.data.frame(dist_paired_tests))
}
  # Unified TSS same-pair check: verify same Peak_ID+Gene has same distance
  loop_pairs <- dist_df_clean %>% filter(Source == "looplook") %>%
    group_by(Peak_ID, Gene) %>%
    summarise(Distance = min(Distance, na.rm = TRUE), .groups = "drop")
  chip_pairs <- dist_df_clean %>% filter(Source == "ChIPseeker") %>%
    group_by(Peak_ID, Gene) %>%
    summarise(Distance = min(Distance, na.rm = TRUE), .groups = "drop")
  same_pairs <- inner_join(loop_pairs, chip_pairs, by = c("Peak_ID", "Gene"), suffix = c("_loop", "_chip"))
  if (nrow(same_pairs) > 0) {
    same_pairs <- same_pairs %>% mutate(Delta = abs(Distance_loop - Distance_chip))
    message(sprintf("    Unified TSS check: %d same Peak-Gene pairs; median delta = %.0f bp", nrow(same_pairs), median(same_pairs$Delta, na.rm = TRUE)))
    if (any(same_pairs$Delta > 1e-8, na.rm = TRUE)) {
      stop(sprintf("Unified distance validation failed: %d same Peak_ID-Gene pairs have different distances.",
        sum(same_pairs$Delta > 1e-8, na.rm = TRUE)))
    }
  } else {
    warning("No identical Peak_ID-Gene pairs were found. Check peak coordinate conventions and peak identifiers.")
  }

  # Gene-level distance comparison: same gene, different peak assignments
  message("    Computing gene-level distance deltas (looplook vs ChIPseeker)...")
  gene_dist_loop <- dist_df_clean %>%
    filter(Source == "looplook") %>%
    group_by(Mode, Gene) %>%
    summarise(MinDist_loop = min(Dist_kb, na.rm = TRUE), .groups = "drop")
  gene_dist_chip <- dist_df_clean %>%
    filter(Source == "ChIPseeker") %>%
    group_by(Mode, Gene) %>%
    summarise(MinDist_chip = min(Dist_kb, na.rm = TRUE), .groups = "drop")
  gene_dist_delta <- inner_join(gene_dist_loop, gene_dist_chip, by = c("Mode", "Gene")) %>%
    mutate(Delta_kb = MinDist_loop - MinDist_chip) %>%
    filter(is.finite(Delta_kb)) %>%
    mutate(Reach = factor(case_when(
      Delta_kb > 100  ~ "looplook >> ChIPseeker (>100 kb)",
      Delta_kb > 10   ~ "looplook further (10-100 kb)",
      abs(Delta_kb) <= 10 ~ "same (±10 kb)",
      Delta_kb < -10  ~ "ChIPseeker further"),
      levels = c("ChIPseeker further", "same (±10 kb)", "looplook further (10-100 kb)", "looplook >> ChIPseeker (>100 kb)")))

  if (nrow(gene_dist_delta) > 0) {
    p_gene_dist <- ggplot(gene_dist_delta, aes(x = Delta_kb, y = reorder(Mode, Delta_kb, median), fill = Mode)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_boxplot(outlier.size = 0.3, alpha = 0.7, linewidth = 0.3) +
      scale_x_continuous(labels = scales::comma) +
      labs(title = "Gene-Level Distance Difference: Looplook minus ChIPseeker",
           subtitle = sprintf("Same gene, min distance per source | %d genes across %d modes | Positive = looplook reaches further",
                              nrow(gene_dist_delta), n_distinct(gene_dist_delta$Mode)),
           x = "Delta Distance (kb): looplook - ChIPseeker", y = NULL) +
      theme_classic() +
      theme(legend.position = "none", axis.text.y = element_text(size = 5),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
    queueable_ggsave(file.path(dist_dir, "Distance_GeneLevel_Delta_Looplook_vs_ChIPseeker.pdf"), p_gene_dist, width = 10, height = 14)

    # Summary table
    gene_dist_summary <- gene_dist_delta %>%
      group_by(Mode) %>%
      summarise(N_SharedGenes = n(), Median_Delta_kb = median(Delta_kb),
                Pct_looplook_further = mean(Delta_kb > 10) * 100,
                Pct_same = mean(abs(Delta_kb) <= 10) * 100,
                .groups = "drop") %>%
      mutate(Pipeline = case_when(grepl("^anno_",Mode)~"Basic", grepl("^refined_",Mode)~"E-refined",
        grepl("^chrom_only_",Mode)~"C-refined", grepl("^chrom_",Mode)~"I-refined", TRUE~"Other")) %>%
      arrange(desc(Median_Delta_kb))

    addWorksheet(wb_dist, "GeneLevel_Delta_Summary"); writeData(wb_dist, "GeneLevel_Delta_Summary", as.data.frame(gene_dist_summary))
    message(sprintf("    Gene-level delta: %d shared genes; median delta = %.1f kb", nrow(gene_dist_delta), median(gene_dist_delta$Delta_kb)))
  }

out_dist <- file.path(dist_dir, "Peak_Gene_Distance_Analysis.xlsx")

  # Evidence-level distance analysis
  message("    Computing per-evidence distances...")
  ev_dist_list <- list()
  for (i in seq_along(modes)) {
    md <- modes[[i]]; if (md$near) next; nm <- md$name
    tgl <- md$ann$target_gene_links; bed_info <- md$ann$target_annotation
    if (is.null(tgl) || nrow(tgl) == 0 || is.null(bed_info)) next
    gev <- build_mode_gene_evidence(md, nm)
    if (nrow(gev) == 0) next
    if (nrow(gev) > 3000) {
      set.seed(42)
      n_total <- nrow(gev)
      gev <- gev %>%
        group_by(Evidence) %>%
        group_modify(function(.x, .y) {
          n_take <- min(nrow(.x), ceiling(3000 * nrow(.x) / n_total))
          slice_sample(.x, n = n_take)
        }) %>%
        ungroup()
      if (nrow(gev) > 3000) gev <- slice_sample(gev, n = 3000)
    }
    peak_map <- data.frame(input_id=bed_info$input_id, chr=bed_info$seqnames, start=bed_info$start, end=bed_info$end, stringsAsFactors=FALSE)
    peak_gr <- GRanges(seqnames=peak_map$chr, ranges=IRanges(start=peak_map$start, end=peak_map$end))
    for (k in seq_len(nrow(gev))) {
      pid <- gev$Input_ID[k]; g <- gev$Gene[k]; ev <- as.character(gev$Evidence[k])
      pidx <- which(peak_map$input_id == pid); if (length(pidx) == 0) next
      tss_hit <- tss_gr[tss_gr$SYMBOL == g]; if (length(tss_hit) == 0) next
      d <- tryCatch(min(GenomicRanges::distance(peak_gr[pidx[1]], tss_hit), na.rm=TRUE), error=function(e) NA_real_)
      if (!is.na(d) && is.finite(d) && d >= 0) ev_dist_list[[length(ev_dist_list)+1]] <- data.frame(Mode=nm, Input_ID=pid, Gene=g, Evidence=ev, Source=if("Source" %in% colnames(gev)) gev$Source[k] else "unknown", Distance=as.numeric(d), TargetColumn=get_target_gene_column(md), stringsAsFactors=FALSE)
    }
  }
  if (length(ev_dist_list) > 0) {
    ev_dist_df <- bind_rows(ev_dist_list) %>% filter(Distance>=0, Distance<5e6) %>%
      mutate(Dist_kb=(Distance+1)/1000, Evidence = normalize_evidence(Evidence),
        Pipeline=factor(case_when(grepl("^anno_",Mode)~"Basic",grepl("^refined_",Mode)~"E-refined",grepl("^chrom_only_",Mode)~"C-refined",grepl("^chrom_",Mode)~"I-refined",TRUE~"Other"),levels=c("Basic","E-refined","I-refined","C-refined")))
    p_ev_dist <- ggplot(ev_dist_df, aes(x=Evidence, y=Dist_kb, fill=Evidence)) +
      geom_boxplot(outlier.shape=16, outlier.size=0.3, outlier.alpha=0.15, coef=1.5, alpha=0.6, width=0.6) +
      scale_y_log10(labels=scales::comma) + scale_fill_manual(values=evidence_palette) + facet_wrap(~Pipeline, ncol=4) +
      labs(title="Peak-Gene Distance by Evidence Type", x=NULL, y="Distance (kb, log scale)") +
      theme_minimal() + theme(legend.position="none", panel.border=element_rect(color="grey10",fill=NA,linewidth=0.4), strip.text=element_text(face="bold",size=9), axis.text.x=element_text(angle=45,hjust=1,size=6))
    queueable_ggsave(file.path(dist_dir, "Distance_by_Evidence.pdf"), p_ev_dist, width=14, height=5.5)
    addWorksheet(wb_dist, "Distance_by_Evidence"); writeData(wb_dist, "Distance_by_Evidence", as.data.frame(ev_dist_df))
  }
  saveWorkbook(wb_dist, out_dist, overwrite = TRUE)
# Key metrics
med_loop <- dist_summary %>% filter(GeneSet == "Only_looplook") %>%
  pull(Med_Dist_kb) %>% median(na.rm = TRUE)
med_chip <- dist_summary %>% filter(GeneSet == "Only_ChIPseeker") %>%
  pull(Med_Dist_kb) %>% median(na.rm = TRUE)
distance_ratio <- if (is.finite(med_chip) && med_chip > 0) med_loop / med_chip else NA_real_
message(sprintf("    looplook-only median = %.0f kb vs ChIPseeker-only = %.0f kb (%.1fx ratio)",
  med_loop, med_chip, distance_ratio))
message("    Note: distance ratio reflects 3D assignment reach, not functional accuracy.")

}  # end if (nrow(dist_df_clean) > 0)

message("    Distance analysis results: ", dist_dir)
# ══════════════════════════════════════════════════════════════════════════════
# hop1 Biological Interpretability: primary (hop0) vs primary + expanded (hop1) loops
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== hop0 vs hop1 Biological Interpretability ====")
# Rationale: hop=1 assigns genes via one intermediate loop (two-hop connection).
# These indirect assignments may capture biologically meaningful long-range
# regulation OR represent non-specific chromatin contacts. This analysis compares
# hop0 (primary) vs hop1 (primary+expanded) gene sets to assess their biological value.

hop_dir <- file.path(base_out_dir, "hop_interpretability")
dir.create(hop_dir, recursive = TRUE, showWarnings = FALSE)
wb_hop <- createWorkbook()  # created early for all hop sheets

# Pair matched modes: same pipeline + map + fill, differing only in hop
hop_pairs <- list()
for (md in modes) {
  nm <- md$name
  if (grepl("_hop1$", nm)) {
    nm_hop0 <- gsub("_hop1$", "", nm)
    if (nm_hop0 %in% names(mode_genes)) {
      hop_pairs[[nm]] <- c(hop0 = nm_hop0, hop1 = nm)
    }
  }
}
message(sprintf("    Found %d matched hop0/hop1 pairs", length(hop_pairs)))

# ── 1. Gene set overlap: hop0 vs hop1 per mode ──
hop_overlap <- list()
for (pair_name in names(hop_pairs)) {
  nm0 <- hop_pairs[[pair_name]]["hop0"]
  nm1 <- hop_pairs[[pair_name]]["hop1"]
  g0 <- unique(toupper(mode_genes[[nm0]]))
  g1 <- unique(toupper(mode_genes[[nm1]]))
  shared <- intersect(g0, g1)
  only0 <- setdiff(g0, g1)
  only1 <- setdiff(g1, g0)
  hop_overlap[[pair_name]] <- data.frame(
    hop0_mode = nm0, hop1_mode = nm1,
    hop0_N = length(g0), hop1_N = length(g1),
    Shared = length(shared),
    Only_hop0 = length(only0), Only_hop1 = length(only1),
    Jaccard = length(shared) / max(1, length(union(g0, g1))),
    Hop1_specific_ratio = length(only1) / max(1, length(g1)),
    stringsAsFactors = FALSE)
}
hop_overlap_df <- bind_rows(hop_overlap)

message(sprintf("    Median Jaccard (hop0 vs hop1): %.3f", median(hop_overlap_df$Jaccard)))
message(sprintf("    Median hop1-specific gene ratio: %.2f%%",
  100 * median(hop_overlap_df$Hop1_specific_ratio)))

# ── 2. NES comparison: hop0 vs hop1 for looplook term ──
# Read NES data from the largest sample size
hop_nes_dir <- file.path(base_out_dir, paste0("size", max(sample_sizes)), "tmp_sim")
if (dir.exists(hop_nes_dir)) {
  hop_nes_data <- list()
  for (pair_name in names(hop_pairs)) {
    for (h in c("hop0", "hop1")) {
      nm <- hop_pairs[[pair_name]][h]
      csv_file <- file.path(hop_nes_dir, paste0("sim_", nm, ".csv"))
      if (file.exists(csv_file)) {
        df <- tryCatch(read.csv(csv_file, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", csv_file); NULL })
        if (!is.null(df) && nrow(df) > 0) {
          df$hop <- h
          df$pair <- pair_name
          hop_nes_data[[length(hop_nes_data) + 1]] <- df
        }
      }
    }
  }
  if (length(hop_nes_data) > 0) {
    hop_nes_df <- bind_rows(hop_nes_data)

    # Per-pair NES comparison
    hop_nes_compare <- hop_nes_df %>%
      filter(ID == "looplook") %>%
      group_by(pair, hop) %>%
      summarise(Mean_NES = mean(NES, na.rm = TRUE),
                SD_NES = sd(NES, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = hop, values_from = c(Mean_NES, SD_NES))

    # Paired NES scatter
    hop_nes_wide <- hop_nes_compare %>%
      mutate(
        Pipeline = case_when(
          grepl("^anno_", pair) ~ "Basic",
          grepl("^refined_", pair) ~ "E-refined",
          grepl("^chrom_only_", pair) ~ "C-refined",
          grepl("^chrom_", pair) ~ "I-refined", TRUE ~ "Other"))

    hop_nes_wide <- hop_nes_wide %>%
      mutate(
        GeneMap = ifelse(grepl("_promoter_", pair), "promoter", "all"),
        Fill    = ifelse(grepl("_T(_hop1)?$", pair), "T (filled)", "F (strict)"),
        Label   = paste(GeneMap, Fill, sep = ", "))

    p_hop_nes <- ggplot(hop_nes_wide, aes(x = Mean_NES_hop0, y = Mean_NES_hop1,
                          color = Pipeline, shape = Label)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                     "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
      scale_shape_manual(values = c("all, F (strict)" = 16, "all, T (filled)" = 17,
                                     "promoter, F (strict)" = 15, "promoter, T (filled)" = 18)) +
      labs(title = "NES: Primary (hop0) vs Primary + Expanded (hop1) Loop Connections",
           subtitle = sprintf("%d matched mode pairs | Color=Pipeline, Shape=gene map + fill",
                              nrow(hop_nes_wide)),
           x = "Mean NES (hop=0, direct)", y = "Mean NES (hop=1, indirect)") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom", legend.box = "vertical")
    queueable_ggsave(file.path(hop_dir, "hop0_vs_hop1_NES_Scatter.pdf"), p_hop_nes, width = 8, height = 6.5)

    # NES difference distribution (hop1 - hop0)
    hop_nes_wide <- hop_nes_wide %>%
      mutate(NES_diff = Mean_NES_hop1 - Mean_NES_hop0)

    p_hop_diff <- ggplot(hop_nes_wide, aes(x = NES_diff, fill = Pipeline)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_histogram(bins = 15, alpha = 0.7, position = "stack") +
      scale_fill_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                    "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
      labs(title = "NES Difference: hop1 - hop0",
           subtitle = sprintf("Median diff = %.3f | Negative = hop1 has stronger NES",
                              median(hop_nes_wide$NES_diff, na.rm = TRUE)),
           x = "NES(hop1) - NES(hop0)", y = "Count") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom")
    queueable_ggsave(file.path(hop_dir, "hop0_vs_hop1_NES_Difference.pdf"), p_hop_diff, width = 6, height = 4)
  }
}

# ── 2b. Effect Size comparison: hop0 vs hop1 (signal_lfc) ──
if (exists("gene_stat") && nrow(gene_stat) > 0) {
  hop_es_list <- list()
  for (pair_name in names(hop_pairs)) {
    hop0_genes <- intersect(unique(toupper(mode_genes[[hop_pairs[[pair_name]]["hop0"]]])), gene_stat$gene)
    hop1_genes <- intersect(unique(toupper(mode_genes[[hop_pairs[[pair_name]]["hop1"]]])), gene_stat$gene)
    only_hop0 <- setdiff(hop0_genes, hop1_genes)
    only_hop1 <- setdiff(hop1_genes, hop0_genes)
    if (length(only_hop0) >= 5 && length(only_hop1) >= 5) {
      x <- gene_stat$signal_lfc[gene_stat$gene %in% only_hop1]
      y <- gene_stat$signal_lfc[gene_stat$gene %in% only_hop0]
      rb <- rank_biserial(x, y)
      hop_es_list[[pair_name]] <- data.frame(
        pair = pair_name, hop0_N = length(hop0_genes), hop1_N = length(hop1_genes),
        only_hop0_N = length(only_hop0), only_hop1_N = length(only_hop1),
        ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi, stringsAsFactors = FALSE)
    }
  }
  if (length(hop_es_list) > 0) {
    hop_es_df <- bind_rows(hop_es_list) %>%
      mutate(Pipeline = factor(case_when(
        grepl("^anno_", pair) ~ "Basic", grepl("^refined_", pair) ~ "E-refined",
        grepl("^chrom_only_", pair) ~ "C-refined", grepl("^chrom_", pair) ~ "I-refined",
        TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")))

    p_hop_es <- ggplot(hop_es_df, aes(x = 0, y = reorder(pair, ES), xmin = ES_lo, xmax = ES_hi, color = Pipeline)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(x = ES), size = 2.5) +
      geom_errorbarh(height = 0.3, linewidth = 0.5) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                     "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      labs(title = "Effect Size: hop0 vs hop1 gene sets (signal_lfc)",
           subtitle = sprintf("Median ES = %.3f | Positive = hop1 genes have stronger downregulation signal",
                              median(hop_es_df$ES, na.rm = TRUE)),
           x = "Rank-Biserial Effect Size (hop1 relative to hop0)", y = NULL) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom")
    queueable_ggsave(file.path(hop_dir, "hop0_vs_hop1_EffectSize_Forest.pdf"), p_hop_es, width = 7, height = 6)

    addWorksheet(wb_hop, "ES_hop0_vs_hop1")
    writeData(wb_hop, "ES_hop0_vs_hop1", as.data.frame(hop_es_df))
  }
}

# ── 3. GO enrichment of expanded-added genes ──
# Compare functional enrichment of genes added by expanded hop1 scope
# OPTIMIZED: use cached entrez_mapping if available, else batch-compute
if (!exists("entrez_mapping")) {
  all_hop_genes <- unique(unlist(lapply(mode_genes, function(g) unique(toupper(g)))))
  entrez_batch_hop <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = all_hop_genes,
      columns = "ENTREZID", keytype = "SYMBOL"))
  entrez_mapping <- entrez_batch_hop %>%
    transmute(SYMBOL = toupper(SYMBOL), ENTREZID = as.character(ENTREZID)) %>%
    filter(!is.na(SYMBOL), SYMBOL != "", !is.na(ENTREZID), ENTREZID != "") %>%
    distinct()
}
entrez_by_map <- function(genes) {
  if (length(genes) == 0) return(character(0))
  genes <- unique(toupper(genes))
  entrez_mapping %>%
    filter(SYMBOL %in% genes) %>%
    pull(ENTREZID) %>%
    unique()
}

hop1_go_list <- list()
for (pair_name in names(hop_pairs)) {
  g0 <- intersect(unique(toupper(mode_genes[[hop_pairs[[pair_name]]["hop0"]]])), all_genes)
  g1 <- intersect(unique(toupper(mode_genes[[hop_pairs[[pair_name]]["hop1"]]])), all_genes)
  only1 <- setdiff(g1, g0)
  if (length(only1) < 5) next
  entrez <- entrez_by_map(only1)
  if (length(entrez) < 5) next
  gsea_out <- safe_go_once(entrez, "BP", NA_integer_, pair_name, pair_name, gouniv)
  gsea_status_log[[length(gsea_status_log) + 1L]] <- gsea_out$status
  if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
    hop1_go_list[[pair_name]] <- gsea_out$result %>%
      dplyr::select(ID, Description, pvalue, p.adjust, Count) %>%
      mutate(Mode = pair_name, N_genes = length(only1))
  }
}

# Export hop analysis results
addWorksheet(wb_hop, "GeneSet_Overlap")
writeData(wb_hop, "GeneSet_Overlap", as.data.frame(hop_overlap_df))
if (exists("hop_nes_wide")) {
  addWorksheet(wb_hop, "NES_Comparison")
  writeData(wb_hop, "NES_Comparison", as.data.frame(hop_nes_wide))
}
if (length(hop1_go_list) > 0) {
  hop1_go_df <- bind_rows(hop1_go_list)
  addWorksheet(wb_hop, "hop1_specific_GO_BP")
  writeData(wb_hop, "hop1_specific_GO_BP", as.data.frame(hop1_go_df))
}
out_hop <- file.path(hop_dir, "hop_Interpretability_Analysis.xlsx")
saveWorkbook(wb_hop, out_hop, overwrite = TRUE)

# Key findings summary
message(sprintf("    Median Jaccard between hop0 and hop1 gene sets: %.3f", median(hop_overlap_df$Jaccard)))
message(sprintf("    Median %% of hop1 genes that are hop1-specific: %.1f%%",
  100 * median(hop_overlap_df$Hop1_specific_ratio)))
if (exists("hop_nes_wide")) {
  med_nes_diff <- median(hop_nes_wide$NES_diff, na.rm = TRUE)
  message(sprintf("    Median NES(hop1) - NES(hop0): %.3f (%s)",
    med_nes_diff,
    ifelse(abs(med_nes_diff) < 0.05, "negligible difference -> hop1 adds genes with similar NES signature",
    ifelse(med_nes_diff < 0, "hop1 has stronger negative NES -> expanded targets add stronger downregulated signal",
    "hop1 has weaker NES -> expanded targets may dilute signal"))))
}
message("    hop interpretability results: ", hop_dir)


# ══════════════════════════════════════════════════════════════════════════════
# Control-pool Sensitivity Analysis: compare strict vs standard non-assigned control
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Control-pool Sensitivity Analysis ====")
# Rationale: The main analysis uses a strict non-assigned control pool (genes NOT
# in looplook AND NOT in ChIPseeker), which may affect NES magnitude. Here we re-run
# for selected modes using a standard control pool (all genes NOT in looplook only)
# to assess robustness.

bg_sens_dir <- file.path(base_out_dir, "background_sensitivity")
dir.create(bg_sens_dir, recursive = TRUE, showWarnings = FALSE)
bg_sens_checkpoint_root <- file.path(bg_sens_dir, ".gsea_checkpoints")
ensure_cache_signature(
  bg_sens_checkpoint_root,
  c(benchmark_cache_signature_base, list(
    Module = "background_sensitivity_raw_lfc",
    RankedListMD5 = md5_object(glist),
    Iterations = min(200L, n_iterations),
    SampleSize = as.integer(max(sample_sizes))
  ))
)

# Select representative modes covering the pipeline diversity
bg_sens_modes <- modes[!sapply(modes, `[[`, "near")]  # exclude ctrl modes
# Use a reduced iteration count for efficiency
bg_sens_iter <- min(200, n_iterations)
bg_sens_size <- sample_sizes[length(sample_sizes)]  # use largest sample size

bg_sens_results <- list()
for (idx in seq_along(bg_sens_modes)) {
  md <- bg_sens_modes[[idx]]
  nm <- md$name
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) next

  pool_loop <- intersect(genes, all_genes)
  # Standard background: all genes not in looplook (standard GSEA convention)
  pool_bg_std <- setdiff(all_genes, genes)
  # Strict background: also exclude ChIPseeker genes (used in main analysis)
  pool_bg_strict <- setdiff(all_genes, union(genes, pool_chip))

  sim_local <- checkpointed_lapply(
  1:bg_sens_iter,
  function(i) {
      set.seed(450 + i)
      s_loop <- sample(pool_loop, sample_n(pool_loop, bg_sens_size))
      s_bg_std <- sample(pool_bg_std, sample_n(pool_bg_std, bg_sens_size))
      s_bg_strict <- sample(pool_bg_strict, sample_n(pool_bg_strict, bg_sens_size))

      t2g <- bind_rows(
        data.frame(term = "looplook", gene = s_loop),
        data.frame(term = "Background_std", gene = s_bg_std),
        data.frame(term = "Background_strict", gene = s_bg_strict))

      expected_terms_bg <- c("looplook", "Background_std", "Background_strict")
      gsea_out <- safe_gsea_once(glist, t2g, expected_terms_bg, bg_sens_size,
        "bg_sens", nm, i)

      if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
        df <- gsea_out$result[, c("ID", "NES")]
        df$Iteration <- i; df$Mode <- nm
        df$Requested_N <- bg_sens_size
        actual_map <- c(looplook = length(s_loop), Background_std = length(s_bg_std), Background_strict = length(s_bg_strict))
        pool_map   <- c(looplook = length(pool_loop), Background_std = length(pool_bg_std), Background_strict = length(pool_bg_strict))
        df$Actual_N <- unname(actual_map[as.character(df$ID)])
        df$Pool_N   <- unname(pool_map[as.character(df$ID)])
        df$Sampling_fraction <- df$Actual_N / pmax(1, df$Pool_N)
        list(result = df, status = gsea_out$status)
      } else {
        list(result = NULL, status = gsea_out$status)
      }
    },
  checkpoint_dir = file.path(bg_sens_dir, ".gsea_checkpoints", nm),
  label = sprintf("bg_sens mode=%s", nm),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)
  statuses <- lapply(sim_local, function(x) x$status)
  gsea_status_log <- c(gsea_status_log, statuses[!vapply(statuses, is.null, logical(1))])
  sim_local <- lapply(sim_local, function(x) x$result)
  sim_local <- sim_local[vapply(sim_local, is.data.frame, logical(1))]
  if (length(sim_local) > 0) {
    bg_sens_results[[nm]] <- bind_rows(sim_local) %>% filter(!is.na(NES))
  }
  rm(sim_local); gc()
  message(sprintf("    [%2d/%d] %-30s std_bg=%d strict_bg=%d",
    idx, length(bg_sens_modes), nm, length(pool_bg_std), length(pool_bg_strict)))
}

# Merge and compare
bg_sens_df <- bind_rows(bg_sens_results)
bg_sens_df$Mode <- factor(bg_sens_df$Mode, levels = sapply(modes, `[[`, "name"))

# Compute NES difference between backgrounds
bg_compare <- bg_sens_df %>%
  filter(grepl("^Background", ID)) %>%
  group_by(Mode, ID) %>%
  summarise(Mean_NES = mean(NES, na.rm = TRUE), SD_NES = sd(NES, na.rm = TRUE),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = ID, values_from = c(Mean_NES, SD_NES)) %>%
  mutate(NES_diff = `Mean_NES_Background_strict` - `Mean_NES_Background_std`)

# Sensitivity scatter plot: NES_std vs NES_strict
p_bg_sens <- ggplot(bg_compare, aes(x = `Mean_NES_Background_std`,
                                      y = `Mean_NES_Background_strict`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.7, color = "#E64B35") +
  geom_smooth(method = "lm", se = TRUE, color = "#4DBBD5", linewidth = 0.8) +
  labs(title = "Control-pool Sensitivity: Strict vs Standard Control",
       subtitle = sprintf("N=%d modes, %d iterations each | Points on diagonal = no bias",
                          nrow(bg_compare), bg_sens_iter),
       x = "Mean NES (Standard control)", y = "Mean NES (Strict control)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
  close_batch_devices()
  queueable_ggsave(file.path(bg_sens_dir, "ControlPool_Sensitivity_Scatter.pdf"), p_bg_sens, width = 7, height = 6)

# Export comparison table
bg_sens_xlsx <- file.path(bg_sens_dir, "ControlPool_Sensitivity_Comparison.xlsx")
wb_bg <- createWorkbook()
addWorksheet(wb_bg, "NES_Comparison")
writeData(wb_bg, "NES_Comparison", bg_compare)
addWorksheet(wb_bg, "Raw_Sensitivity_Data")
writeData(wb_bg, "Raw_Sensitivity_Data", as.data.frame(bg_sens_df))
saveWorkbook(wb_bg, bg_sens_xlsx, overwrite = TRUE)

# Summary message
med_diff <- median(bg_compare$NES_diff, na.rm = TRUE)
message(sprintf("    Median NES difference (strict - std): %.4f", med_diff))
if (abs(med_diff) < 0.1) {
  message("    -> Control-pool choice has negligible impact; results are robust.")
} else {
  message(sprintf("    -> NOTE: Control-pool definition shifts NES by median %.3f. Interpret main NES values with this caveat.", med_diff))
}
  message("    Control-pool sensitivity results: ", bg_sens_dir)

# ══════════════════════════════════════════════════════════════════════════════
# 新增功能：Signal Density Analysis + Validation + Convergence QC
# ══════════════════════════════════════════════════════════════════════════════

# ── Signal Density Analysis ──────────────────────────────────────────────────
message("\n\n>>> ==== Signal Density Analysis ====")

# Pairwise comparison function
pairwise_compare_signal <- function(gene_sets, gene_stat, score_col = "signal_lfc") {
  pairs <- list(
    c("Only_looplook", "Only_ChIPseeker"),
    c("Only_looplook", "Background"),
    c("Only_ChIPseeker", "Background"),
    c("Intersection", "Background")
  )
  res <- list()
  for (p in pairs) {
    xg <- intersect(gene_sets[[p[1]]], gene_stat$gene)
    yg <- intersect(gene_sets[[p[2]]], gene_stat$gene)
    x <- gene_stat[[score_col]][gene_stat$gene %in% xg]
    y <- gene_stat[[score_col]][gene_stat$gene %in% yg]
    if (length(x) < 2 || length(y) < 2) {
      res[[length(res) + 1]] <- data.frame(
        Set_A = p[1], Set_B = p[2], N_A = length(x), N_B = length(y),
        Med_A = NA_real_, Med_B = NA_real_, P = NA_real_,
        ES = NA_real_, ES_lo = NA_real_, ES_hi = NA_real_,
        Med_Diff = NA_real_,
        stringsAsFactors = FALSE)
      next
    }
    wt <- wilcox.test(x, y, alternative = "greater", exact = FALSE)
    rb <- rank_biserial(x, y)
    med_diff <- median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
    res[[length(res) + 1]] <- data.frame(
      Set_A = p[1], Set_B = p[2], N_A = length(x), N_B = length(y),
      Med_A = median(x, na.rm = TRUE), Med_B = median(y, na.rm = TRUE),
      Med_Diff = med_diff,
      P = wt$p.value,
      ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi,
      stringsAsFactors = FALSE)
  }
  result_df <- bind_rows(res)
  result_df$FDR <- p.adjust(result_df$P, method = "BH")
  return(result_df)
}

# Top-down enrichment: test if gene set is enriched in top 10% strongest downregulation signal function
top_down_enrichment <- function(test_genes, universe, top_down_genes) {
  test_genes <- intersect(test_genes, universe)
  top_down_genes <- intersect(top_down_genes, universe)
  a <- length(intersect(test_genes, top_down_genes))
  b <- length(setdiff(test_genes, top_down_genes))
  c <- length(setdiff(top_down_genes, test_genes))
  d <- length(setdiff(universe, union(test_genes, top_down_genes)))
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.frame(
    N_set = length(test_genes), N_top_down_overlap = a,
    Fraction_top_down = a / max(1, length(test_genes)),
    OddsRatio = unname(ft$estimate), P_value = ft$p.value,
    stringsAsFactors = FALSE)
}

# Validation function
validate_known_targets <- function(gene_sets, known, universe) {
  kinu <- intersect(known, universe)
  res <- list()
  for (nm in names(gene_sets)) {
    genes <- intersect(gene_sets[[nm]], universe)
    ov <- intersect(genes, kinu)
    a <- length(ov); b <- length(setdiff(genes, kinu))
    c <- length(setdiff(kinu, genes)); d <- length(setdiff(universe, union(genes, kinu)))
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    prec <- a / max(1, length(genes)); rec <- a / max(1, length(kinu))
    f1 <- 2 * prec * rec / max(1e-10, prec + rec)
    res[[nm]] <- data.frame(Set = nm, N_set = length(genes), N_hits = a,
      N_known = length(kinu), Precision = prec, Recall = rec, F1 = f1,
      OR = unname(ft$estimate), P = ft$p.value, stringsAsFactors = FALSE)
  }
  rdf <- bind_rows(res)
  rdf$FDR <- p.adjust(rdf$P, method = "BH")
  rdf
}

# Run signal density analysis for each mode (parallelized, with checkpoint)
signal_cache_rds <- file.path(base_out_dir, ".signal_density_cache.rds")
if (file.exists(signal_cache_rds) && !force_recompute) {
  message("    Loading cached signal density results...")
  signal_density_results <- readRDS(signal_cache_rds)
} else {
  signal_density_list <- parallel::mclapply(seq_along(modes), function(i) {
  md <- modes[[i]]
  nm <- md$name
  if (md$near) return(NULL)  # skip control modes
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) return(NULL)

  loop_genes <- unique(toupper(genes))
  chip_g <- get_chip_reference(md)

  sets <- list(
    Intersection = intersect(loop_genes, chip_g),
    Only_looplook = setdiff(loop_genes, chip_g),
    Only_ChIPseeker = setdiff(chip_g, loop_genes),
    All_looplook = loop_genes,
    All_ChIPseeker = chip_g,
    Background = setdiff(all_genes, union(loop_genes, chip_g))
  )

  # Summary stats
  summary_tbl <- bind_rows(lapply(names(sets), function(snm) {
    g <- intersect(sets[[snm]], gene_stat$gene)
    if (length(g) == 0) return(data.frame(Set = snm, N = 0, Med_LFC = NA, Med_Signal = NA, Med_Signal_pv = NA))
    ss <- gene_stat[gene_stat$gene %in% g, ]
    data.frame(Set = snm, N = length(g),
               Med_LFC = median(ss$log2FoldChange, na.rm = TRUE),
               Med_Signal = median(ss$signal_lfc, na.rm = TRUE),
               Med_Signal_pv = median(ss$signal_pv, na.rm = TRUE))
  }))

  # Pairwise comparisons: primary (signal_lfc) + supplemental (signal_pv)
  comparisons      <- pairwise_compare_signal(sets, gene_stat, "signal_lfc")
  comparisons_pv   <- pairwise_compare_signal(sets, gene_stat, "signal_pv")

  # Top-down enrichment: test if gene set is enriched in top 10% strongest downregulation signal
  gene_stat_valid <- gene_stat %>%
    dplyr::filter(!is.na(signal_lfc), is.finite(signal_lfc)) %>%
    dplyr::arrange(dplyr::desc(signal_lfc))
  top_n <- max(1, round(0.10 * nrow(gene_stat_valid)))
  top_down_genes <- gene_stat_valid %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(gene)
  top_down_tbl <- bind_rows(lapply(names(sets), function(snm) {
    cbind(Set = snm, top_down_enrichment(sets[[snm]], all_genes, top_down_genes))
  }))

  # Literature-informed BRD4 gene panel validation (three tiers)
  # Tier 1: Core BRD4 targets with strongest literature evidence
  # Tier 2: Context-dependent BRD4-associated genes
  # Tier 3: Combined panel (core + context)
  panel_list <- list(
    brd4_core       = brd4_core_panel,
    brd4_context    = brd4_context_panel,
    brd4_associated = brd4_associated_panel)
  validation_lst <- list()
  for (panel_name in names(panel_list)) {
    panel_genes <- panel_list[[panel_name]]
    if (length(panel_genes) > 0) {
      validation_lst[[panel_name]] <- validate_known_targets(sets, panel_genes, all_genes) %>%
        mutate(Panel = panel_name)
    }
  }
  validation <- bind_rows(validation_lst)

  cat(sprintf("    [%2d] %-30s Inter=%d Only_look=%d ES=%.3f FDR=%.2e\n",
    i, nm, length(sets$Intersection), length(sets$Only_looplook),
    comparisons$ES[comparisons$Set_A == "Intersection" & comparisons$Set_B == "Background"],
    comparisons$FDR[comparisons$Set_A == "Intersection" & comparisons$Set_B == "Background"]))

  list(name = nm, summary = summary_tbl, comparisons = comparisons,
    comparisons_pv = comparisons_pv,
    top_down = top_down_tbl, validation = validation)
}, mc.cores = min(8L, n_cores), mc.preschedule = TRUE, mc.cleanup = TRUE)

# Combine parallel results into named list
signal_density_results <- list()
for (res in signal_density_list) {
  if (is.null(res) || is.null(res$name)) next
  signal_density_results[[res$name]] <- list(
    summary = res$summary, comparisons = res$comparisons,
    comparisons_pv = res$comparisons_pv,
    top_down = res$top_down, validation = res$validation)
}
  rm(signal_density_list)
  saveRDS(signal_density_results, signal_cache_rds)
  message("    Signal density cached")
}

# Export signal density results
if (length(signal_density_results) > 0) {
  wb_sig <- createWorkbook()

  # Master summary
  master_summary <- bind_rows(lapply(names(signal_density_results), function(nm) {
    signal_density_results[[nm]]$summary %>% mutate(Mode = nm)
  }))
  addWorksheet(wb_sig, "Signal_Summary")
  writeData(wb_sig, "Signal_Summary", master_summary)

  # Master comparisons
  master_comp <- bind_rows(lapply(names(signal_density_results), function(nm) {
    signal_density_results[[nm]]$comparisons %>% mutate(Mode = nm, Score = "signal_lfc")
  }))
  master_comp_pv <- bind_rows(lapply(names(signal_density_results), function(nm) {
    signal_density_results[[nm]]$comparisons_pv %>% mutate(Mode = nm, Score = "signal_pv")
  }))
  addWorksheet(wb_sig, "Pairwise_Comparisons")
  writeData(wb_sig, "Pairwise_Comparisons", master_comp)
  addWorksheet(wb_sig, "Pairwise_Comparisons_pv")
  writeData(wb_sig, "Pairwise_Comparisons_pv", bind_rows(lapply(names(signal_density_results), function(nm) {
    signal_density_results[[nm]]$comparisons_pv %>% mutate(Mode = nm, Score = "signal_pv")
  })))

  # Top-down enrichment: test if gene set is enriched in top 10% strongest downregulation signal
  master_top_down <- bind_rows(lapply(names(signal_density_results), function(nm) {
    signal_density_results[[nm]]$top_down %>% mutate(Mode = nm)
  }))
  addWorksheet(wb_sig, "Top_Down_Enrichment")
  writeData(wb_sig, "Top_Down_Enrichment", master_top_down)

  # Validation
  master_valid <- bind_rows(lapply(names(signal_density_results), function(nm) {
    signal_density_results[[nm]]$validation %>% mutate(Mode = nm)
  }))
  addWorksheet(wb_sig, "Validation")
  writeData(wb_sig, "Validation", master_valid)

  out_sig <- file.path(base_out_dir, "Signal_Density_Analysis.xlsx")
  saveWorkbook(wb_sig, out_sig, overwrite = TRUE)
  message("    Signal density results saved: ", out_sig)

  # ── Signal Density 可视化：效应量排名图 ─────────────────────────────────
  # Fix: mclapply fork can leave broken graphics devices, causing ggsave deadlock
  close_batch_devices()
  message("    Creating Signal Density visualizations...")

  # Add hop column for hop0/hop1 split visualizations
  master_comp <- master_comp %>%
    mutate(Hop = ifelse(grepl("_hop1$", Mode), "hop1", "hop0"))

  # 提取关键比较的效应量
  es_summary <- master_comp %>%
    filter(Set_A == "Intersection" & Set_B == "Background") %>%
    mutate(
      Pipeline = case_when(
        grepl("^anno_", Mode) ~ "Basic",
        grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined",
        grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"
      ),
      Pipeline = factor(Pipeline, levels = c("Basic", "E-refined", "I-refined", "C-refined"))
    ) %>%
    left_join(paired_nes_summary %>% dplyr::select(Mode, P_Delta_lt_0, P_Delta_gt_0), by = "Mode") %>%
    arrange(Pipeline, desc(ES))

  # 效应量 Forest plot — split by hop
  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_sub <- es_summary %>% filter(Hop == hv) %>% arrange(Pipeline, desc(ES))
    if (nrow(es_sub) == 0) next
    n_stable <- sum(es_sub$P_Delta_lt_0 > 0.95 | es_sub$P_Delta_gt_0 > 0.95, na.rm = TRUE)
    p_es <- ggplot(es_sub, aes(x = ES, y = reorder(Mode, ES), color = Pipeline)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.3, linewidth = 0.4) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD", "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      labs(title = sprintf("Signal Density: Effect Size %s", hv_label),
           subtitle = sprintf("Intersection vs Non-assigned | %d modes | %d/%d with consistent ΔNES direction (P<0>0.95)",
             nrow(es_sub), n_stable, nrow(es_sub)),
           x = "Effect Size (rank-biserial)", y = NULL, color = "Pipeline") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom", axis.text.y = element_text(size = 7))
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_Forest_%s.pdf", hv)), p_es, width = 10, height = 8)
  }

  # 按Pipeline分组的效应量boxplot — split by hop
  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_sub <- es_summary %>% filter(Hop == hv)
    if (nrow(es_sub) == 0) next
    p_es_box <- ggplot(es_sub, aes(x = Pipeline, y = ES, fill = Pipeline)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD", "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      labs(title = sprintf("Signal Density: Effect Size by Pipeline (%s)", hv_label),
           subtitle = "Intersection vs Non-assigned control | Each point = one mode",
           x = NULL, y = "Effect Size (rank-biserial)") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "none")
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_byPipeline_%s.pdf", hv)), p_es_box, width = 6, height = 5)
  }

  # ── 追加比较：Only_looplook vs Non-assigned control / Only_ChIPseeker vs Non-assigned control / Head-to-head ──
  es_extra <- master_comp %>%
    mutate(Hop = ifelse(grepl("_hop1$", Mode), "hop1", "hop0")) %>%
    filter((Set_A == "Only_looplook" & Set_B == "Background") |
           (Set_A == "Only_ChIPseeker" & Set_B == "Background") |
           (Set_A == "Only_looplook" & Set_B == "Only_ChIPseeker")) %>%
    mutate(
      Comparison = case_when(
        Set_A == "Only_looplook" & Set_B == "Background" ~ "Only_looplook vs Background",
        Set_A == "Only_ChIPseeker" & Set_B == "Background" ~ "Only_ChIPseeker vs Background",
        Set_A == "Only_looplook" & Set_B == "Only_ChIPseeker" ~ "Only_looplook vs Only_ChIPseeker"),
      Comparison = factor(Comparison, levels = c("Only_looplook vs Background",
                                                  "Only_ChIPseeker vs Background",
                                                  "Only_looplook vs Only_ChIPseeker")),
      Pipeline = factor(case_when(
        grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"),
        levels = c("Basic", "E-refined", "I-refined", "C-refined")))

  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_hop <- es_extra %>% filter(Hop == hv)
    if (nrow(es_hop) == 0) next
    p_es_extra <- ggplot(es_hop, aes(x = ES, y = reorder(Mode, ES), color = Pipeline)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 1.5) +
      geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.3, linewidth = 0.3) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                     "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
      labs(title = sprintf("Signal Density: Unique Gene Sets %s", hv_label),
           subtitle = "Only_looplook vs Only_ChIPseeker head-to-head comparison",
           x = "Effect Size (rank-biserial)", y = NULL, color = "Pipeline") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom", axis.text.y = element_text(size = 6))
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_UniqueSets_%s.pdf", hv)), p_es_extra, width = 10, height = 10)
  }

  # ── Combined comparison: Intersection + Only_looplook + Only_ChIPseeker, side by side per Pipeline ──
  es_combined <- master_comp %>%
    mutate(Hop = ifelse(grepl("_hop1$", Mode), "hop1", "hop0")) %>%
    filter((Set_A == "Intersection" & Set_B == "Background") |
           (Set_A == "Only_looplook" & Set_B == "Background") |
           (Set_A == "Only_ChIPseeker" & Set_B == "Background")) %>%
    mutate(
      Comparison = case_when(
        Set_A == "Intersection" & Set_B == "Background" ~ "Intersection vs BG",
        Set_A == "Only_looplook" & Set_B == "Background" ~ "Only_looplook vs BG",
        Set_A == "Only_ChIPseeker" & Set_B == "Background" ~ "Only_ChIPseeker vs BG"),
      Comparison = factor(Comparison, levels = c("Intersection vs BG",
                                                  "Only_looplook vs BG",
                                                  "Only_ChIPseeker vs BG")),
      Pipeline = factor(case_when(
        grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"),
        levels = c("Basic", "E-refined", "I-refined", "C-refined")))

  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_h <- es_combined %>% filter(Hop == hv)
    if (nrow(es_h) == 0) next

    es_head2head <- es_h %>%
      filter(Comparison %in% c("Only_looplook vs BG", "Only_ChIPseeker vs BG")) %>%
      tidyr::pivot_wider(id_cols = c(Mode, Pipeline),
                         names_from = Comparison, values_from = ES) %>%
      filter(!is.na(`Only_looplook vs BG`), !is.na(`Only_ChIPseeker vs BG`))
    if (nrow(es_head2head) >= 5) {
      wt_h2h <- wilcox.test(es_head2head[["Only_looplook vs BG"]],
                            es_head2head[["Only_ChIPseeker vs BG"]],
                            paired = TRUE, alternative = "two.sided", exact = FALSE)
      h2h_label <- sprintf("Paired Wilcoxon p = %.3f", wt_h2h$p.value)
    } else { h2h_label <- "Insufficient data" }

    p_es_cmb <- ggplot(es_h, aes(x = Pipeline, y = ES, fill = Comparison)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_boxplot(outlier.shape = 16, outlier.size = 0.8, alpha = 0.65,
                   width = 0.6, position = position_dodge(0.8)) +
      scale_fill_manual(values = c("Intersection vs BG" = "#00A087",
                                    "Only_looplook vs BG" = "#E64B35",
                                    "Only_ChIPseeker vs BG" = "#4DBBD5")) +
      labs(title = sprintf("Signal Density: Gene Set Origin (%s)", hv_label),
            subtitle = h2h_label,
           x = NULL, y = "Effect Size (rank-biserial)") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom", legend.title = element_blank(),
            panel.grid = element_blank())
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_Combined_%s.pdf", hv)),
           p_es_cmb, width = 8, height = 5)
  }

# ── Global full-set effect size: looplook-all vs ChIPseeker-all, split by hop ──
global_es_list <- list()
for (nm in names(mode_genes)) {
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) next
  loop_genes <- intersect(unique(toupper(genes)), gene_stat$gene)
  chip_g     <- intersect(get_chip_reference(modes[[which(sapply(modes, function(m) m$name == nm))[1]]]), gene_stat$gene)
  if (length(loop_genes) < 5 || length(chip_g) < 5) next
  x <- gene_stat$signal_lfc[gene_stat$gene %in% loop_genes]
  y <- gene_stat$signal_lfc[gene_stat$gene %in% chip_g]
  rb <- rank_biserial(x, y)
  global_es_list[[nm]] <- data.frame(
    Mode = nm, N_looplook = length(loop_genes), N_ChIPseeker = length(chip_g),
    ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi, stringsAsFactors = FALSE)
}
if (length(global_es_list) > 1) {
  global_es_df <- bind_rows(global_es_list) %>%
    mutate(
      Pipeline = factor(case_when(
        grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")),
      Hop = factor(ifelse(grepl("_hop1$", Mode), "hop1 (primary+expanded)", "hop0 (primary)"),
                   levels = c("hop0 (primary)", "hop1 (primary+expanded)")))

  p_global_es <- ggplot(global_es_df, aes(x = ES, y = reorder(Mode, ES), color = Pipeline)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 1.5) +
    geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.3, linewidth = 0.3) +
    scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                   "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
    facet_wrap(~ Hop, ncol = 2, scales = "free_y") +
    labs(title = "Global Effect Size: Looplook-all Target Genes vs ChIPseeker-all",
         subtitle = sprintf("All assigned genes (not split) | %d modes | Positive = looplook targets have stronger downregulation",
                            nrow(global_es_df)),
         x = "Rank-Biserial Effect Size", y = NULL, color = "Pipeline") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "bottom", axis.text.y = element_text(size = 6),
          strip.text = element_text(face = "bold", size = 10))
  queueable_ggsave(file.path(base_out_dir, "Signal_Density_EffectSize_Global_Looplook_vs_ChIPseeker.pdf"),
         p_global_es, width = 12, height = 10)

  .safe_add_sheet(wb, "Global_ES_Looplook_vs_ChIP", as.data.frame(global_es_df))
  message(sprintf("    Global looplook vs ChIPseeker ES: median=%.3f, range [%.3f, %.3f]",
    median(global_es_df$ES, na.rm = TRUE),
    min(global_es_df$ES, na.rm = TRUE), max(global_es_df$ES, na.rm = TRUE)))
}

  es_combined_pv <- bind_rows(
    master_comp_pv %>% filter(Set_A == "Intersection" & Set_B == "Background") %>% mutate(Comparison = "Intersection vs BG"),
    master_comp_pv %>% filter(Set_A == "Only_looplook" & Set_B == "Background") %>% mutate(Comparison = "Only_looplook vs BG"),
    master_comp_pv %>% filter(Set_A == "Only_ChIPseeker" & Set_B == "Background") %>% mutate(Comparison = "Only_ChIPseeker vs BG")) %>%
    mutate(Comparison = factor(Comparison, levels = c("Intersection vs BG", "Only_looplook vs BG", "Only_ChIPseeker vs BG")),
           Pipeline = factor(case_when(
             grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
             grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
             TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")))

  p_es_cmb_pv <- ggplot(es_combined_pv, aes(x = Pipeline, y = ES, fill = Comparison)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.8, alpha = 0.65,
                 width = 0.6, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("Intersection vs BG" = "#00A087", "Only_looplook vs BG" = "#E64B35", "Only_ChIPseeker vs BG" = "#4DBBD5")) +
    labs(title = "Signal Density (pv score): Effect Size by Pipeline",
         subtitle = "supplemental: sign(LFC) × -log10(pvalue) weighted score",
         x = NULL, y = "Effect Size (rank-biserial)") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "bottom", legend.title = element_blank(),
          panel.grid = element_blank())

  queueable_ggsave(file.path(base_out_dir, "Signal_Density_EffectSize_pv_Combined.pdf"),
         p_es_cmb_pv, width = 8, height = 5)

  # ── signal_pv mirrored plots (supplemental) ──
  pv_prefix <- "pv"
  # Forest plot — split by hop
  es_summary_pv <- master_comp_pv %>%
    mutate(Hop = ifelse(grepl("_hop1$", Mode), "hop1", "hop0")) %>%
    filter(Set_A == "Intersection" & Set_B == "Background") %>%
    mutate(Pipeline = factor(case_when(
      grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
      TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")))

  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_sub <- es_summary_pv %>% filter(Hop == hv) %>% arrange(Pipeline, desc(ES))
    if (nrow(es_sub) == 0) next
    p_es_pv <- ggplot(es_sub, aes(x = ES, y = reorder(Mode, ES), color = Pipeline)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.3, linewidth = 0.4) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD", "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      labs(title = sprintf("Signal Density (pv): Effect Size %s", hv_label),
           subtitle = "Intersection vs Non-assigned control | signal_pv = -LFC × -log10(pvalue) | supplemental",
           x = "Effect Size (rank-biserial)", y = NULL, color = "Pipeline") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom", axis.text.y = element_text(size = 7))
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_Forest_pv_%s.pdf", hv)), p_es_pv, width = 10, height = 8)
  }

  # By-pipeline boxplot — split by hop
  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_sub <- es_summary_pv %>% filter(Hop == hv)
    if (nrow(es_sub) == 0) next
    p_es_box_pv <- ggplot(es_sub, aes(x = Pipeline, y = ES, fill = Pipeline)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD", "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      labs(title = sprintf("Signal Density (pv): Effect Size by Pipeline (%s)", hv_label),
           subtitle = "signal_pv weighted score | supplemental",
           x = NULL, y = "Effect Size (rank-biserial)") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "none")
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_byPipeline_pv_%s.pdf", hv)), p_es_box_pv, width = 6, height = 5)
  }

  # UniqueSets forest — split by hop
  es_extra_pv <- master_comp_pv %>%
    mutate(Hop = ifelse(grepl("_hop1$", Mode), "hop1", "hop0")) %>%
    filter((Set_A == "Only_looplook" & Set_B == "Background") |
           (Set_A == "Only_ChIPseeker" & Set_B == "Background") |
           (Set_A == "Only_looplook" & Set_B == "Only_ChIPseeker")) %>%
    mutate(Comparison = factor(case_when(
      Set_A == "Only_looplook" & Set_B == "Background" ~ "Only_looplook vs Background",
      Set_A == "Only_ChIPseeker" & Set_B == "Background" ~ "Only_ChIPseeker vs Background",
      Set_A == "Only_looplook" & Set_B == "Only_ChIPseeker" ~ "Only_looplook vs Only_ChIPseeker"),
      levels = c("Only_looplook vs Background", "Only_ChIPseeker vs Background", "Only_looplook vs Only_ChIPseeker")),
      Pipeline = factor(case_when(
        grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")))

  for (hv in c("hop0", "hop1")) {
    hv_label <- if (hv == "hop0") "hop0 (primary)" else "hop1 (primary+expanded)"
    es_hop <- es_extra_pv %>% filter(Hop == hv)
    if (nrow(es_hop) == 0) next
    p_es_extra_pv <- ggplot(es_hop, aes(x = ES, y = reorder(Mode, ES), color = Pipeline)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 1.5) +
      geom_errorbarh(aes(xmin = ES_lo, xmax = ES_hi), height = 0.3, linewidth = 0.3) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD", "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
      labs(title = sprintf("Signal Density (pv): Unique Gene Sets %s", hv_label),
           subtitle = "signal_pv weighted score | Only_looplook vs Only_ChIPseeker head-to-head | supplemental",
           x = "Effect Size (rank-biserial)", y = NULL, color = "Pipeline") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            legend.position = "bottom", axis.text.y = element_text(size = 6))
    queueable_ggsave(file.path(base_out_dir, sprintf("Signal_Density_EffectSize_UniqueSets_pv_%s.pdf", hv)), p_es_extra_pv, width = 10, height = 10)
  }

  message("    Signal Density plots saved")
}


# ══════════════════════════════════════════════════════════════════════════════
# Expanded-only Analysis: biological value of genes added by neighbor_hop=1
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Expanded-only Target Analysis ====")
expanded_dir <- file.path(base_out_dir, "expanded_only")
dir.create(expanded_dir, recursive = TRUE, showWarnings = FALSE)

infer_pipeline_label <- function(mode_name) {
  dplyr::case_when(
    grepl("^anno_", mode_name) ~ "Basic",
    grepl("^refined_", mode_name) ~ "E-refined",
    grepl("^chrom_only_", mode_name) ~ "C-refined",
    grepl("^chrom_", mode_name) ~ "I-refined",
    TRUE ~ "Other"
  )
}

# Fixed-N subsampling interval for rank-biserial effect size.
# The point estimate uses all available genes. The interval deliberately caps
# each group for predictable runtime and is therefore reported as a
# subsampling interval, not as a full-set bootstrap confidence interval.
expanded_rank_biserial <- function(
    x,
    y,
    R = expanded_interval_R,
    seed = bootstrap_seed,
    max_per_group = expanded_interval_max_per_group
) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2L || length(y) < 2L) {
    return(list(est = NA_real_, lo = NA_real_, hi = NA_real_))
  }

  est <- rank_biserial_point(x, y)
  R <- max(100L, as.integer(R))
  nx <- min(length(x), as.integer(max_per_group))
  ny <- min(length(y), as.integer(max_per_group))

  set.seed(seed)
  boot_values <- replicate(R, {
    xb <- sample(x, nx, replace = TRUE)
    yb <- sample(y, ny, replace = TRUE)
    rank_biserial_point(xb, yb)
  })
  boot_values <- boot_values[is.finite(boot_values)]

  if (length(boot_values) < 20L) {
    return(list(est = est, lo = NA_real_, hi = NA_real_))
  }

  list(
    est = est,
    lo = unname(stats::quantile(
      boot_values,
      0.025,
      na.rm = TRUE,
      names = FALSE
    )),
    hi = unname(stats::quantile(
      boot_values,
      0.975,
      na.rm = TRUE,
      names = FALSE
    ))
  )
}



expanded_smd <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  pooled_sd <- sqrt((stats::var(x) + stats::var(y)) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled_sd
}

expanded_expression_match <- function(
    expanded_genes,
    reference_genes,
    signal_lookup,
    expr_lookup,
    iterations = 500L,
    bins = 10L,
    max_per_bin = 50L,
    seed = 42L
) {
  expanded_genes <- intersect(
    unique(toupper(expanded_genes)),
    intersect(names(signal_lookup), names(expr_lookup))
  )
  reference_genes <- intersect(
    unique(toupper(reference_genes)),
    intersect(names(signal_lookup), names(expr_lookup))
  )

  if (length(expanded_genes) < 10L || length(reference_genes) < 10L) {
    return(NULL)
  }

  # Use adaptive pooled-expression bins. Fixed deciles are too sparse when
  # expanded-only sets contain only 10-30 genes and can yield no matched output.
  effective_bins <- min(
    as.integer(bins),
    max(2L, floor(min(length(expanded_genes), length(reference_genes)) / 5L))
  )

  match_df <- dplyr::bind_rows(
    data.frame(
      Gene = expanded_genes,
      Group = "Expanded_only",
      Expression = log2(expr_lookup[expanded_genes] + 0.1),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Gene = reference_genes,
      Group = "Reference",
      Expression = log2(expr_lookup[reference_genes] + 0.1),
      stringsAsFactors = FALSE
    )
  ) %>%
    dplyr::filter(is.finite(Expression)) %>%
    dplyr::mutate(
      Expression_Bin = dplyr::ntile(Expression, effective_bins)
    )

  if (sum(match_df$Group == "Expanded_only") < 10L ||
      sum(match_df$Group == "Reference") < 10L) {
    return(NULL)
  }

  matched_es <- numeric()
  smd_after <- numeric()
  n_matched <- integer()
  set.seed(seed)

  for (iteration in seq_len(as.integer(iterations))) {
    sampled_expanded <- character()
    sampled_reference <- character()

    for (bin_id in seq_len(effective_bins)) {
      expanded_bin <- match_df$Gene[
        match_df$Group == "Expanded_only" &
          match_df$Expression_Bin == bin_id
      ]
      reference_bin <- match_df$Gene[
        match_df$Group == "Reference" &
          match_df$Expression_Bin == bin_id
      ]
      n_draw <- min(
        length(expanded_bin),
        length(reference_bin),
        as.integer(max_per_bin)
      )
      if (n_draw < 1L) next
      sampled_expanded <- c(
        sampled_expanded,
        sample(expanded_bin, n_draw, replace = FALSE)
      )
      sampled_reference <- c(
        sampled_reference,
        sample(reference_bin, n_draw, replace = FALSE)
      )
    }

    if (length(sampled_expanded) < 5L ||
        length(sampled_reference) < 5L) {
      next
    }

    es_value <- rank_biserial_point(
      unname(signal_lookup[sampled_expanded]),
      unname(signal_lookup[sampled_reference])
    )
    balance_value <- expanded_smd(
      log2(expr_lookup[sampled_expanded] + 0.1),
      log2(expr_lookup[sampled_reference] + 0.1)
    )

    if (is.finite(es_value)) {
      matched_es <- c(matched_es, es_value)
      smd_after <- c(smd_after, balance_value)
      n_matched <- c(n_matched, length(sampled_expanded))
    }
  }

  minimum_successful_iterations <- max(
    10L,
    ceiling(0.20 * as.integer(iterations))
  )
  if (length(matched_es) < minimum_successful_iterations) return(NULL)

  data.frame(
    Full_ES = rank_biserial_point(
      unname(signal_lookup[expanded_genes]),
      unname(signal_lookup[reference_genes])
    ),
    Matched_ES_Median = stats::median(matched_es, na.rm = TRUE),
    Matched_ES_Low = stats::quantile(
      matched_es, 0.025, na.rm = TRUE, names = FALSE
    ),
    Matched_ES_High = stats::quantile(
      matched_es, 0.975, na.rm = TRUE, names = FALSE
    ),
    P_Matched_ES_gt_0 = mean(matched_es > 0, na.rm = TRUE),
    P_Matched_ES_lt_0 = mean(matched_es < 0, na.rm = TRUE),
    N_Expanded_Available = length(expanded_genes),
    N_Reference_Available = length(reference_genes),
    Median_N_Matched = stats::median(n_matched, na.rm = TRUE),
    Median_Expanded_TPM = stats::median(
      expr_lookup[expanded_genes], na.rm = TRUE
    ),
    Median_Reference_TPM = stats::median(
      expr_lookup[reference_genes], na.rm = TRUE
    ),
    SMD_Before = expanded_smd(
      log2(expr_lookup[expanded_genes] + 0.1),
      log2(expr_lookup[reference_genes] + 0.1)
    ),
    SMD_After_Median = stats::median(smd_after, na.rm = TRUE),
    SMD_After_Low = stats::quantile(
      smd_after, 0.025, na.rm = TRUE, names = FALSE
    ),
    SMD_After_High = stats::quantile(
      smd_after, 0.975, na.rm = TRUE, names = FALSE
    ),
    P_Abs_SMD_lt_0_1 = mean(abs(smd_after) < 0.1, na.rm = TRUE),
    Requested_Matching_Iterations = as.integer(iterations),
    Matching_Iterations = length(matched_es),
    Successful_Matching_Fraction = length(matched_es) /
      max(1L, as.integer(iterations)),
    Effective_Expression_Bins = effective_bins,
    stringsAsFactors = FALSE
  )
}

expanded_topology_gene_table <- function(md, expanded_only_genes) {
  empty <- data.frame(
    Mode = character(),
    Gene = character(),
    Supporting_Peaks = integer(),
    Supporting_Links = integer(),
    Supporting_Anchors = integer(),
    Median_Path_Length = numeric(),
    Max_Path_Length = numeric(),
    Median_Anchor_Degree = numeric(),
    Max_Anchor_Degree = numeric(),
    Median_Component_Size = numeric(),
    Max_Component_Size = numeric(),
    Hub_Degree_Threshold = numeric(),
    Topology_Available = logical(),
    Is_Top5Pct_Hub = logical(),
    Multi_Supported = logical(),
    stringsAsFactors = FALSE
  )

  genes <- unique(toupper(expanded_only_genes))
  if (!is_hop1_mode(md) || length(genes) == 0L) return(empty)

  tgl <- md$ann$target_gene_links
  required <- c(
    "input_id", "gene", "anchor_id", "source", "path_length",
    "in_expanded_target"
  )
  if (is.null(tgl) || nrow(tgl) == 0L ||
      length(setdiff(required, colnames(tgl))) > 0L) {
    return(empty)
  }

  keep <- as.character(tgl$source) == "loop_anchor" &
    is.finite(tgl$path_length) &
    tgl$path_length > 1L &
    tgl$in_expanded_target %in% TRUE &
    toupper(trimws(as.character(tgl$gene))) %in% genes
  if ("anchor_role" %in% colnames(tgl)) {
    keep <- keep & as.character(tgl$anchor_role) == "expanded_anchor"
  }
  if (identical(md$map, "promoter") &&
      "gene_role" %in% colnames(tgl)) {
    keep <- keep & as.character(tgl$gene_role) == "promoter"
  }

  links <- data.frame(
    Input_ID = as.character(tgl$input_id[keep]),
    Gene = toupper(trimws(as.character(tgl$gene[keep]))),
    Anchor_ID = as.character(tgl$anchor_id[keep]),
    Path_Length = as.integer(tgl$path_length[keep]),
    stringsAsFactors = FALSE
  )
  links <- links[
    !is.na(links$Gene) & nzchar(links$Gene) &
      !is.na(links$Anchor_ID) & nzchar(links$Anchor_ID),
    , drop = FALSE
  ]
  if (nrow(links) == 0L) return(empty)

  links$Anchor_Degree <- NA_real_
  links$Component_Size <- NA_real_
  state <- attr(md$ann, "looplook_anchor_state", exact = TRUE)
  if (!is.null(state) && !is.null(state$g) &&
      requireNamespace("igraph", quietly = TRUE)) {
    graph_object <- state$g
    vertex_names <- igraph::V(graph_object)$name
    degree_map <- setNames(
      as.numeric(igraph::degree(graph_object, mode = "all")),
      vertex_names
    )
    component_info <- igraph::components(graph_object)
    component_size_map <- setNames(
      as.numeric(component_info$csize[component_info$membership]),
      names(component_info$membership)
    )
    links$Anchor_Degree <- unname(degree_map[links$Anchor_ID])
    links$Component_Size <- unname(component_size_map[links$Anchor_ID])
  }

  # Define hubs against the complete anchor graph, not only anchors that
  # happen to support the current expanded-only gene set.
  hub_anchor_ids <- character()
  hub_threshold <- NA_real_
  if (exists("graph_object") &&
      requireNamespace("igraph", quietly = TRUE)) {
    graph_degree <- as.numeric(igraph::degree(graph_object, mode = "all"))
    names(graph_degree) <- igraph::V(graph_object)$name
    graph_degree <- graph_degree[is.finite(graph_degree)]
    if (length(graph_degree) >= 2L &&
        length(unique(graph_degree)) > 1L) {
      hub_threshold <- unname(stats::quantile(
        graph_degree,
        probs = 0.95,
        names = FALSE,
        type = 1,
        na.rm = TRUE
      ))
      hub_anchor_ids <- names(graph_degree)[graph_degree >= hub_threshold]
    }
  }
  links$Is_Hub_Anchor <- links$Anchor_ID %in% hub_anchor_ids

  links %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Mode = md$name,
      Supporting_Peaks = dplyr::n_distinct(Input_ID),
      Supporting_Links = dplyr::n(),
      Supporting_Anchors = dplyr::n_distinct(Anchor_ID),
      Median_Path_Length = stats::median(Path_Length, na.rm = TRUE),
      Max_Path_Length = max(Path_Length, na.rm = TRUE),
      Median_Anchor_Degree = if (any(is.finite(Anchor_Degree))) {
        stats::median(Anchor_Degree, na.rm = TRUE)
      } else NA_real_,
      Max_Anchor_Degree = if (any(is.finite(Anchor_Degree))) {
        max(Anchor_Degree, na.rm = TRUE)
      } else NA_real_,
      Median_Component_Size = if (any(is.finite(Component_Size))) {
        stats::median(Component_Size, na.rm = TRUE)
      } else NA_real_,
      Max_Component_Size = if (any(is.finite(Component_Size))) {
        max(Component_Size, na.rm = TRUE)
      } else NA_real_,
      Hub_Degree_Threshold = hub_threshold,
      Topology_Available = any(is.finite(Anchor_Degree)),
      Is_Top5Pct_Hub = any(Is_Hub_Anchor),
      Multi_Supported = Supporting_Peaks >= 2L |
        Supporting_Anchors >= 2L,
      .groups = "drop"
    ) %>%
    dplyr::select(
      Mode, Gene, Supporting_Peaks, Supporting_Links,
      Supporting_Anchors, Median_Path_Length, Max_Path_Length,
      Median_Anchor_Degree, Max_Anchor_Degree,
      Median_Component_Size, Max_Component_Size,
      Hub_Degree_Threshold, Topology_Available,
      Is_Top5Pct_Hub, Multi_Supported
    )
}

hop1_mode_indices <- which(vapply(modes, is_hop1_mode, logical(1)))
hop1_modes <- modes[hop1_mode_indices]

expanded_only_sets <- list()
expanded_catalog_rows <- list()

for (md in hop1_modes) {
  primary_genes <- get_primary_mode_genes(md)
  expanded_pairs <- get_expanded_peak_gene_pairs(md)
  expanded_raw <- sort(unique(expanded_pairs$Gene))
  expanded_only <- get_expanded_only_genes(md)
  primary_plus_expanded <- get_mode_genes(md)

  if (length(intersect(expanded_only, primary_genes)) > 0L) {
    stop(
      "Expanded-only set overlaps primary genes for mode: ",
      md$name,
      call. = FALSE
    )
  }

  if (length(setdiff(expanded_only, primary_plus_expanded)) > 0L) {
    stop(
      "Expanded-only genes are not a subset of hop1 genes for mode: ",
      md$name,
      call. = FALSE
    )
  }

  expanded_only_sets[[md$name]] <- expanded_only
  signature <- if (length(expanded_only) > 0L) {
    md5_object(sort(expanded_only))
  } else {
    NA_character_
  }
  gsea_signature <- if (length(expanded_only) > 0L) {
    md5_object(list(
      expanded_only = sort(intersect(expanded_only, all_genes)),
      assigned_total = sort(intersect(primary_plus_expanded, all_genes))
    ))
  } else {
    NA_character_
  }

  expanded_catalog_rows[[md$name]] <- data.frame(
    Mode = md$name,
    Pipeline = infer_pipeline_label(md$name),
    Map = md$map,
    Filled = isTRUE(md$fill),
    Primary_N = length(primary_genes),
    Expanded_Link_Genes_N = length(expanded_raw),
    Expanded_Only_N = length(expanded_only),
    Primary_Plus_Expanded_N = length(primary_plus_expanded),
    Expanded_Only_Fraction = length(expanded_only) /
      max(1L, length(primary_plus_expanded)),
    Set_Signature = signature,
    GSEA_Signature = gsea_signature,
    stringsAsFactors = FALSE
  )
}

expanded_catalog <- dplyr::bind_rows(expanded_catalog_rows)
expanded_catalog$Expanded_Set_ID <- NA_character_

nonempty_signatures <- unique(
  expanded_catalog$Set_Signature[
    !is.na(expanded_catalog$Set_Signature) &
      expanded_catalog$Expanded_Only_N > 0L
  ]
)

if (length(nonempty_signatures) > 0L) {
  signature_ids <- setNames(
    sprintf("EXP%02d", seq_along(nonempty_signatures)),
    nonempty_signatures
  )
  expanded_catalog$Expanded_Set_ID <- unname(
    signature_ids[expanded_catalog$Set_Signature]
  )
}

gsea_signatures <- unique(
  expanded_catalog$GSEA_Signature[
    !is.na(expanded_catalog$GSEA_Signature) &
      expanded_catalog$Expanded_Only_N > 0L
  ]
)
expanded_catalog$Expanded_GSEA_ID <- NA_character_
if (length(gsea_signatures) > 0L) {
  gsea_signature_ids <- setNames(
    sprintf("EXPG%02d", seq_along(gsea_signatures)),
    gsea_signatures
  )
  expanded_catalog$Expanded_GSEA_ID <- unname(
    gsea_signature_ids[expanded_catalog$GSEA_Signature]
  )
}

expanded_catalog <- expanded_catalog %>%
  dplyr::group_by(Expanded_Set_ID) %>%
  dplyr::mutate(
    Equivalent_Modes = ifelse(
      is.na(Expanded_Set_ID),
      Mode,
      paste(sort(unique(Mode)), collapse = ";")
    ),
    Representative_Mode = ifelse(
      is.na(Expanded_Set_ID),
      Mode,
      sort(unique(Mode))[1L]
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Expanded_GSEA_ID) %>%
  dplyr::mutate(
    GSEA_Equivalent_Modes = ifelse(
      is.na(Expanded_GSEA_ID),
      Mode,
      paste(sort(unique(Mode)), collapse = ";")
    ),
    GSEA_Representative_Mode = ifelse(
      is.na(Expanded_GSEA_ID),
      Mode,
      sort(unique(Mode))[1L]
    )
  ) %>%
  dplyr::ungroup()

utils::write.csv(
  expanded_catalog,
  file.path(expanded_dir, "Expanded_Only_Set_Catalog.csv"),
  row.names = FALSE
)

expanded_gene_rows <- list()
for (mode_name in names(expanded_only_sets)) {
  genes <- expanded_only_sets[[mode_name]]
  if (length(genes) == 0L) next

  set_id <- expanded_catalog$Expanded_Set_ID[
    match(mode_name, expanded_catalog$Mode)
  ]

  expanded_gene_rows[[mode_name]] <- data.frame(
    Mode = mode_name,
    Expanded_Set_ID = set_id,
    Gene = genes,
    log2FoldChange = unname(diff_df[genes, "log2FoldChange"]),
    signal_lfc = unname(setNames(
      gene_stat$signal_lfc,
      gene_stat$gene
    )[genes]),
    DMSO_TPM = if (exists("expr_lookup")) {
      unname(expr_lookup[genes])
    } else {
      NA_real_
    },
    In_BRD4_Core = genes %in% toupper(brd4_core_panel),
    In_BRD4_Context = genes %in% toupper(brd4_context_panel),
    In_BRD4_Associated = genes %in% toupper(brd4_associated_panel),
    stringsAsFactors = FALSE
  )
}
expanded_gene_table <- dplyr::bind_rows(expanded_gene_rows)
if (ncol(expanded_gene_table) == 0L) {
  expanded_gene_table <- data.frame(
    Mode = character(),
    Expanded_Set_ID = character(),
    Gene = character(),
    log2FoldChange = numeric(),
    signal_lfc = numeric(),
    DMSO_TPM = numeric(),
    In_BRD4_Core = logical(),
    In_BRD4_Context = logical(),
    In_BRD4_Associated = logical(),
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  expanded_gene_table,
  file.path(expanded_dir, "Expanded_Only_Genes_Long.csv"),
  row.names = FALSE
)

expanded_topology_rows <- list()
for (mode_name in names(expanded_only_sets)) {
  genes <- expanded_only_sets[[mode_name]]
  if (length(genes) == 0L) next
  md_index <- which(vapply(
    modes,
    function(x) identical(x$name, mode_name),
    logical(1)
  ))
  if (length(md_index) != 1L) next
  expanded_topology_rows[[mode_name]] <- expanded_topology_gene_table(
    modes[[md_index]],
    genes
  )
}
expanded_topology_df <- dplyr::bind_rows(expanded_topology_rows)
if (ncol(expanded_topology_df) == 0L) {
  expanded_topology_df <- data.frame(
    Mode = character(), Gene = character(), Supporting_Peaks = integer(),
    Supporting_Links = integer(), Supporting_Anchors = integer(),
    Median_Path_Length = numeric(), Max_Path_Length = numeric(),
    Median_Anchor_Degree = numeric(), Max_Anchor_Degree = numeric(),
    Median_Component_Size = numeric(), Max_Component_Size = numeric(),
    Hub_Degree_Threshold = numeric(), Topology_Available = logical(),
    Is_Top5Pct_Hub = logical(), Multi_Supported = logical(),
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  expanded_topology_df,
  file.path(expanded_dir, "ExpandedOnly_Topology_Support.csv"),
  row.names = FALSE
)

expanded_mode_defs <- list()
for (mode_name in names(expanded_only_sets)) {
  genes <- expanded_only_sets[[mode_name]]
  if (length(genes) == 0L) next
  md_index <- which(vapply(
    modes,
    function(x) identical(x$name, mode_name),
    logical(1)
  ))
  if (length(md_index) != 1L) {
    stop("Unable to resolve expanded-only mode: ", mode_name, call. = FALSE)
  }
  row <- expanded_catalog %>% dplyr::filter(Mode == mode_name)
  expanded_mode_defs[[mode_name]] <- list(
    id = mode_name,
    expanded_set_id = row$Expanded_Set_ID[1L],
    genes = genes,
    representative_mode = mode_name,
    equivalent_modes = mode_name,
    pipeline = row$Pipeline[1L],
    map = row$Map[1L],
    filled = row$Filled[1L],
    md = modes[[md_index]]
  )
}

# Unique expanded-gene definitions are used for GO, where only the gene set
# itself matters.
expanded_set_defs <- list()
if (length(nonempty_signatures) > 0L) {
  for (set_id in unique(na.omit(expanded_catalog$Expanded_Set_ID))) {
    rows <- expanded_catalog %>%
      dplyr::filter(Expanded_Set_ID == set_id)
    representative_mode <- rows$Representative_Mode[1L]
    def <- expanded_mode_defs[[representative_mode]]
    def$id <- set_id
    def$equivalent_modes <- rows$Equivalent_Modes[1L]
    expanded_set_defs[[set_id]] <- def
  }
}

# GSEA is deduplicated only when both the expanded-only genes and the
# non-assigned control universe are identical. This avoids treating F/T modes
# as equivalent when their primary assignments differ.
expanded_gsea_defs <- list()
if (length(gsea_signatures) > 0L) {
  for (gsea_id in unique(na.omit(expanded_catalog$Expanded_GSEA_ID))) {
    rows <- expanded_catalog %>%
      dplyr::filter(Expanded_GSEA_ID == gsea_id)
    representative_mode <- rows$GSEA_Representative_Mode[1L]
    def <- expanded_mode_defs[[representative_mode]]
    def$id <- gsea_id
    def$equivalent_modes <- rows$GSEA_Equivalent_Modes[1L]
    expanded_gsea_defs[[gsea_id]] <- def
  }
}

message(sprintf(
  "    Expanded-only: %d hop1 modes, %d unique gene sets, %d unique GSEA analyses",
  length(hop1_modes),
  length(expanded_set_defs),
  length(expanded_gsea_defs)
))

if (length(expanded_mode_defs) == 0L) {
  warning(
    "No expanded-only genes were found in any hop1 mode. ",
    "The expanded-only GSEA/effect-size/GO modules will be skipped.",
    call. = FALSE
  )
}

# ── Expanded-only GSEA: exact analysis duplicates only ─────────────────────
expanded_gsea_status <- list()
expanded_gsea_summaries <- list()
expanded_delta_summaries <- list()

expanded_gsea_plan_rows <- list()
if (length(expanded_gsea_defs) > 0L) {
  for (set_id in names(expanded_gsea_defs)) {
    def <- expanded_gsea_defs[[set_id]]
    expanded_pool <- intersect(def$genes, all_genes)
    assigned_total <- intersect(get_mode_genes(def$md), all_genes)
    control_pool <- setdiff(all_genes, union(assigned_total, pool_chip))
    for (requested_n in sample_sizes) {
      expanded_gsea_plan_rows[[paste(set_id, requested_n, sep = "::")]] <-
        data.frame(
          Mode = set_id,
          Representative_Mode = def$representative_mode,
          Equivalent_Modes = def$equivalent_modes,
          Requested_N = as.integer(requested_n),
          Actual_N = min(
            as.integer(requested_n),
            length(expanded_pool),
            length(control_pool)
          ),
          Expanded_Pool_N = length(expanded_pool),
          Control_Pool_N = length(control_pool),
          stringsAsFactors = FALSE
        )
    }
  }
}
expanded_gsea_plan <- dplyr::bind_rows(expanded_gsea_plan_rows)
if (nrow(expanded_gsea_plan) > 0L) {
  expanded_gsea_plan <- expanded_gsea_plan %>%
    dplyr::group_by(Mode, Actual_N) %>%
    dplyr::mutate(
      Duplicate_Actual_N = dplyr::n() > 1L,
      Duplicate_Actual_N_Group_Size = dplyr::n(),
      Canonical_Requested_N = min(Requested_N)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Sampling_Status = ifelse(
        Duplicate_Actual_N,
        paste0(
          "Duplicated Actual_N; canonical Requested_N=",
          Canonical_Requested_N
        ),
        "Unique Actual_N"
      )
    )
} else {
  expanded_gsea_plan <- data.frame(
    Mode = character(), Representative_Mode = character(),
    Equivalent_Modes = character(), Requested_N = integer(),
    Actual_N = integer(), Expanded_Pool_N = integer(),
    Control_Pool_N = integer(), Duplicate_Actual_N = logical(),
    Duplicate_Actual_N_Group_Size = integer(),
    Canonical_Requested_N = integer(), Sampling_Status = character(),
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  expanded_gsea_plan,
  file.path(expanded_dir, "ExpandedOnly_GSEA_SamplingPlan.csv"),
  row.names = FALSE
)

expanded_cache_signature <- list(
  DiffMD5 = preliminary_signature$DiffMD5,
  MainRDataMD5 = preliminary_signature$MainRDataMD5,
  LooplookVersion = preliminary_signature$LooplookVersion,
  RankedListMD5 = md5_object(glist),
  ExpandedCatalogMD5 = md5_object(expanded_catalog),
  SamplingPlanMD5 = md5_object(expanded_gsea_plan),
  Iterations = n_iterations,
  SampleSizes = sample_sizes,
  PrimarySampleSize = primary_sample_size,
  GSEAMinSize = gsea_min_size,
  GSEAMaxSize = gsea_max_size,
  ExpandedIntervalR = expanded_interval_R,
  ExpandedIntervalMaxPerGroup = expanded_interval_max_per_group,
  ExpandedExpressionMatchIterations = expanded_expression_match_iterations,
  ExpandedExpressionMatchBins = expanded_expression_match_bins,
  ExpandedExpressionMatchMaxPerBin = expanded_expression_match_max_per_bin
)
expanded_cache_signature_path <- file.path(
  expanded_dir,
  "ExpandedOnly_GSEA_Cache_Signature.rds"
)
existing_expanded_cache_files <- list.files(
  expanded_dir,
  pattern = "^expanded_.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!file.exists(expanded_cache_signature_path) &&
    !force_recompute &&
    length(existing_expanded_cache_files) > 0L) {
  stop(
    "Unsigned expanded-only GSEA cache detected. Set force_recompute=TRUE ",
    "or delete expanded_only/tmp_gsea before rerunning.",
    call. = FALSE
  )
}
if (file.exists(expanded_cache_signature_path) && !force_recompute) {
  previous_expanded_signature <- readRDS(expanded_cache_signature_path)
  # Legacy v9 signatures included the whole script MD5. That made harmless
  # plotting/configuration edits invalidate expensive GSEA output. Ignore only
  # that obsolete field; all biological inputs and numerical parameters remain
  # strictly checked.
  previous_expanded_signature$ScriptMD5 <- NULL
  if (!identical(previous_expanded_signature, expanded_cache_signature)) {
    stop(
      "Expanded-only GSEA cache signature mismatch. Set cache_policy='rebuild' ",
      "or delete the expanded_only GSEA cache before rerunning.",
      call. = FALSE
    )
  }
}
saveRDS(
  expanded_cache_signature,
  expanded_cache_signature_path,
  version = 3
)

if (length(expanded_gsea_defs) > 0L) {
  for (sz in sample_sizes) {
    expanded_size_dir <- file.path(
      expanded_dir,
      paste0("size", sz)
    )
    expanded_tmp_dir <- file.path(
      expanded_size_dir,
      "tmp_gsea"
    )
    dir.create(expanded_tmp_dir, recursive = TRUE, showWarnings = FALSE)

    for (set_id in names(expanded_gsea_defs)) {
      def <- expanded_gsea_defs[[set_id]]
      tmp_csv <- file.path(
        expanded_tmp_dir,
        paste0("expanded_", set_id, ".csv")
      )
      expected_terms <- c("Expanded_only", "Background")

      sampling_row <- expanded_gsea_plan %>%
        dplyr::filter(Mode == set_id, Requested_N == as.integer(sz))
      if (nrow(sampling_row) != 1L) {
        stop(
          "Unable to resolve expanded-only GSEA sampling plan for ",
          set_id,
          " at Requested_N=",
          sz,
          call. = FALSE
        )
      }

      if (validate_csv_cache(
          tmp_csv,
          c(
            "ID", "NES", "pvalue", "Iteration", "Mode",
            "Requested_N", "Actual_N", "Duplicate_Actual_N",
            "Duplicate_Actual_N_Group_Size", "Canonical_Requested_N",
            "Sampling_Status", "Expanded_Pool_N", "Control_Pool_N"
          ),
          n_iterations,
          expected_terms
      )) {
        message(sprintf(
          "    [expanded size=%d %s] SKIP (validated cache)",
          sz,
          set_id
        ))
        cached_df <- utils::read.csv(tmp_csv, stringsAsFactors = FALSE)
        cached_status <- cached_df %>%
          dplyr::distinct(Iteration) %>%
          dplyr::transmute(
            SampleSize = as.integer(sz),
            Module = "expanded_only",
            Mode = set_id,
            Iteration = Iteration,
            Status = "cached_validated",
            ExpectedTerms = paste(expected_terms, collapse = ";"),
            ReturnedTerms = paste(expected_terms, collapse = ";"),
            MissingTerms = "",
            Error = NA_character_,
            Requested_N = as.integer(sz),
            Actual_N = sampling_row$Actual_N[1L],
            Duplicate_Actual_N = sampling_row$Duplicate_Actual_N[1L],
            Canonical_Requested_N = sampling_row$Canonical_Requested_N[1L]
          )
        expanded_gsea_status <- c(
          expanded_gsea_status,
          split(cached_status, seq_len(nrow(cached_status)))
        )
        next
      }

      if (file.exists(tmp_csv)) {
        unlink(tmp_csv)
      }

      expanded_pool <- intersect(def$genes, all_genes)
      assigned_total <- intersect(
        get_mode_genes(def$md),
        all_genes
      )
      control_pool <- setdiff(
        all_genes,
        union(assigned_total, pool_chip)
      )

      n_actual <- sampling_row$Actual_N[1L]

      if (n_actual < gsea_min_size) {
        insufficient_status <- make_mode_status(
          sz,
          "expanded_only",
          set_id,
          "insufficient_pool",
          paste(expected_terms, collapse = ";"),
          paste0(
            "expanded=", length(expanded_pool),
            "; control=", length(control_pool),
            "; matched_n=", n_actual
          )
        )
        insufficient_status$Requested_N <- as.integer(sz)
        insufficient_status$Actual_N <- n_actual
        insufficient_status$Duplicate_Actual_N <-
          sampling_row$Duplicate_Actual_N[1L]
        insufficient_status$Canonical_Requested_N <-
          sampling_row$Canonical_Requested_N[1L]
        expanded_gsea_status[[length(expanded_gsea_status) + 1L]] <-
          insufficient_status
        message(sprintf(
          "    [expanded size=%d %s] SKIP (expanded=%d control=%d)",
          sz,
          set_id,
          length(expanded_pool),
          length(control_pool)
        ))
        next
      }

      sim_local <- checkpointed_lapply(
  seq_len(n_iterations),
  function(i) {
            set.seed(41000L + i)
            sampled_expanded <- sample(expanded_pool, n_actual)
            sampled_control <- sample(control_pool, n_actual)
            t2g <- dplyr::bind_rows(
              data.frame(
                term = "Expanded_only",
                gene = sampled_expanded,
                stringsAsFactors = FALSE
              ),
              data.frame(
                term = "Background",
                gene = sampled_control,
                stringsAsFactors = FALSE
              )
            )

            gsea_out <- safe_gsea_once(
              glist,
              t2g,
              expected_terms,
              sz,
              "expanded_only",
              set_id,
              i
            )

            result <- NULL
            if (!is.null(gsea_out$result) &&
                nrow(gsea_out$result) > 0L) {
              result <- gsea_out$result[, c("ID", "NES", "pvalue")]
              result$Iteration <- i
              result$Mode <- set_id
              result$Representative_Mode <- def$representative_mode
              result$Equivalent_Modes <- def$equivalent_modes
              result$Requested_N <- as.integer(sz)
              result$Actual_N <- n_actual
              result$Duplicate_Actual_N <-
                sampling_row$Duplicate_Actual_N[1L]
              result$Duplicate_Actual_N_Group_Size <-
                sampling_row$Duplicate_Actual_N_Group_Size[1L]
              result$Canonical_Requested_N <-
                sampling_row$Canonical_Requested_N[1L]
              result$Sampling_Status <- sampling_row$Sampling_Status[1L]
              result$Expanded_Pool_N <- length(expanded_pool)
              result$Control_Pool_N <- length(control_pool)
            }

            gsea_out$status$Requested_N <- as.integer(sz)
            gsea_out$status$Actual_N <- n_actual
            gsea_out$status$Duplicate_Actual_N <-
              sampling_row$Duplicate_Actual_N[1L]
            gsea_out$status$Canonical_Requested_N <-
              sampling_row$Canonical_Requested_N[1L]

            list(
              result = result,
              status = gsea_out$status
            )
          },
  checkpoint_dir = file.path(expanded_tmp_dir, ".iteration_checkpoints", set_id),
  label = sprintf("expanded size=%d mode=%s", sz, set_id),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)

      statuses <- lapply(sim_local, function(x) x$status)
      expanded_gsea_status <- c(
        expanded_gsea_status,
        statuses
      )

      result_parts <- lapply(sim_local, function(x) x$result)
      result_parts <- result_parts[
        vapply(result_parts, is.data.frame, logical(1))
      ]

      result_df <- if (length(result_parts) > 0L) {
        dplyr::bind_rows(result_parts) %>%
          dplyr::filter(is.finite(NES))
      } else {
        data.frame(
          ID = character(),
          NES = numeric(),
          pvalue = numeric(),
          Iteration = integer(),
          Mode = character(),
          stringsAsFactors = FALSE
        )
      }

      utils::write.csv(result_df, tmp_csv, row.names = FALSE)
      rm(sim_local, result_parts, result_df)
      invisible(gc())
    }

    expanded_csvs <- list.files(
      expanded_tmp_dir,
      pattern = "^expanded_.*\\.csv$",
      full.names = TRUE
    )
    expanded_parts <- lapply(expanded_csvs, function(path) {
      tryCatch(
        utils::read.csv(path, stringsAsFactors = FALSE),
        error = function(e) {
          warning(
            "Failed to read expanded-only GSEA cache: ",
            path,
            " | ",
            conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )
    })
    expanded_parts <- expanded_parts[
      !vapply(expanded_parts, is.null, logical(1))
    ]
    expanded_gsea_df <- dplyr::bind_rows(expanded_parts)

    utils::write.csv(
      expanded_gsea_df,
      file.path(expanded_size_dir, "ExpandedOnly_GSEA_Raw.csv"),
      row.names = FALSE
    )

    if (nrow(expanded_gsea_df) > 0L) {
      expanded_summary <- expanded_gsea_df %>%
        dplyr::group_by(
          Mode,
          Representative_Mode,
          Equivalent_Modes,
          Requested_N,
          Actual_N,
          Duplicate_Actual_N,
          Duplicate_Actual_N_Group_Size,
          Canonical_Requested_N,
          Sampling_Status,
          ID
        ) %>%
        dplyr::summarise(
          N_Iterations = dplyr::n_distinct(Iteration),
          Mean_NES = mean(NES, na.rm = TRUE),
          Median_NES = stats::median(NES, na.rm = TRUE),
          SD_NES = stats::sd(NES, na.rm = TRUE),
          Q025 = stats::quantile(
            NES,
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),
          Q975 = stats::quantile(
            NES,
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),
          P_NES_lt_0 = mean(NES < 0, na.rm = TRUE),
          P_NES_gt_0 = mean(NES > 0, na.rm = TRUE),
          .groups = "drop"
        )

      duplicate_check <- expanded_gsea_df %>%
        dplyr::count(Mode, Iteration, ID) %>%
        dplyr::filter(n > 1L)
      if (nrow(duplicate_check) > 0L) {
        stop(
          "Duplicate expanded-only GSEA Mode-Iteration-ID rows detected.",
          call. = FALSE
        )
      }

      expanded_delta <- expanded_gsea_df %>%
        dplyr::select(
          Mode,
          Representative_Mode,
          Equivalent_Modes,
          Requested_N,
          Actual_N,
          Duplicate_Actual_N,
          Duplicate_Actual_N_Group_Size,
          Canonical_Requested_N,
          Sampling_Status,
          Iteration,
          ID,
          NES
        ) %>%
        tidyr::pivot_wider(names_from = ID, values_from = NES) %>%
        dplyr::filter(
          is.finite(Expanded_only),
          is.finite(Background)
        ) %>%
        dplyr::mutate(
          Delta_NES = Expanded_only - Background
        )

      expanded_delta_summary <- expanded_delta %>%
        dplyr::group_by(
          Mode,
          Representative_Mode,
          Equivalent_Modes,
          Requested_N,
          Actual_N,
          Duplicate_Actual_N,
          Duplicate_Actual_N_Group_Size,
          Canonical_Requested_N,
          Sampling_Status
        ) %>%
        dplyr::summarise(
          N_Paired = dplyr::n(),
          Median_Delta_NES = stats::median(
            Delta_NES,
            na.rm = TRUE
          ),
          Mean_Delta_NES = mean(Delta_NES, na.rm = TRUE),
          Q025 = stats::quantile(
            Delta_NES,
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),
          Q975 = stats::quantile(
            Delta_NES,
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),
          P_Delta_lt_0 = mean(Delta_NES < 0, na.rm = TRUE),
          P_Delta_gt_0 = mean(Delta_NES > 0, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::left_join(
          expanded_catalog %>%
            dplyr::filter(
              !is.na(Expanded_GSEA_ID),
              Mode == GSEA_Representative_Mode
            ) %>%
            dplyr::select(
              Expanded_GSEA_ID,
              Pipeline,
              Map,
              Filled,
              Expanded_Only_N
            ),
          by = c("Mode" = "Expanded_GSEA_ID")
        )

      utils::write.csv(
        expanded_summary,
        file.path(expanded_size_dir, "ExpandedOnly_GSEA_Summary.csv"),
        row.names = FALSE
      )
      utils::write.csv(
        expanded_delta_summary,
        file.path(
          expanded_size_dir,
          "ExpandedOnly_GSEA_PairedDelta.csv"
        ),
        row.names = FALSE
      )

      p_expanded_delta <- ggplot(
        expanded_delta_summary,
        aes(
          x = Median_Delta_NES,
          y = reorder(Mode, Median_Delta_NES),
          color = Pipeline,
          shape = factor(Duplicate_Actual_N, levels = c(FALSE, TRUE))
        )
      ) +
        geom_vline(
          xintercept = 0,
          linetype = "dashed",
          color = "grey50"
        ) +
        geom_errorbarh(
          aes(xmin = Q025, xmax = Q975),
          height = 0.25,
          linewidth = 0.35
        ) +
        geom_point(size = 2) +
        scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 17),
          labels = c(
            "FALSE" = "Unique Actual_N",
            "TRUE" = "Duplicated Actual_N"
          ),
          drop = FALSE
        ) +
        scale_color_manual(values = c(
          "Basic" = "#1F77B4",
          "E-refined" = "#9467BD",
          "I-refined" = "#E64B35",
          "C-refined" = "#FFA500"
        )) +
        labs(
          title = paste0(
            "Expanded-only GSEA effect (sample size = ",
            sz,
            ")"
          ),
          subtitle = paste0(
            "Requested_N=", sz,
            "; Actual_N range=",
            min(expanded_delta_summary$Actual_N, na.rm = TRUE),
            "-",
            max(expanded_delta_summary$Actual_N, na.rm = TRUE),
            ". Triangle = this Actual_N is duplicated by another requested size. ",
            "Paired Delta_NES = Expanded-only - Background."
          ),
          x = "Median paired Delta_NES",
          y = "Hop1 mode",
          color = "Pipeline",
          shape = "Requested/actual sampling"
        ) +
        theme_classic() +
        theme(
          plot.title = element_text(
            face = "bold",
            size = 12,
            hjust = 0.5
          ),
          plot.subtitle = element_text(
            size = 9,
            color = "grey40",
            hjust = 0.5
          ),
          legend.position = "bottom",
          axis.text.y = element_text(size = 8)
        )

      close_batch_devices()
      queueable_ggsave(
        file.path(
          expanded_size_dir,
          "ExpandedOnly_GSEA_PairedDelta.pdf"
        ),
        p_expanded_delta,
        width = 8,
        height = max(4.5, 0.32 * nrow(expanded_delta_summary) + 2)
      )

      wb_expanded_gsea <- createWorkbook()
      addWorksheet(wb_expanded_gsea, "Set_Catalog")
      writeData(
        wb_expanded_gsea,
        "Set_Catalog",
        as.data.frame(expanded_catalog)
      )
      addWorksheet(wb_expanded_gsea, "Sampling_Plan")
      writeData(
        wb_expanded_gsea,
        "Sampling_Plan",
        as.data.frame(
          expanded_gsea_plan %>%
            dplyr::filter(Requested_N == as.integer(sz))
        )
      )
      addWorksheet(wb_expanded_gsea, "GSEA_Summary")
      writeData(
        wb_expanded_gsea,
        "GSEA_Summary",
        as.data.frame(expanded_summary)
      )
      addWorksheet(wb_expanded_gsea, "Paired_Delta")
      writeData(
        wb_expanded_gsea,
        "Paired_Delta",
        as.data.frame(expanded_delta_summary)
      )
      addWorksheet(wb_expanded_gsea, "Raw")
      writeData(
        wb_expanded_gsea,
        "Raw",
        as.data.frame(expanded_gsea_df)
      )
      saveWorkbook(
        wb_expanded_gsea,
        file.path(
          expanded_size_dir,
          "ExpandedOnly_GSEA.xlsx"
        ),
        overwrite = TRUE
      )

      expanded_gsea_summaries[[as.character(sz)]] <- expanded_summary
      expanded_delta_summaries[[as.character(sz)]] <-
        expanded_delta_summary
    }
  }
}

expanded_gsea_status_df <- if (length(expanded_gsea_status) > 0L) {
  dplyr::bind_rows(expanded_gsea_status)
} else {
  data.frame(
    SampleSize = integer(),
    Module = character(),
    Mode = character(),
    Iteration = integer(),
    Status = character(),
    ExpectedTerms = character(),
    ReturnedTerms = character(),
    MissingTerms = character(),
    Error = character(),
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  expanded_gsea_status_df,
  file.path(expanded_dir, "ExpandedOnly_GSEA_Status.csv"),
  row.names = FALSE
)

# ── Expanded-only effect size ────────────────────────────────────────────────
expanded_es_rows <- list()
expanded_expr_matched_rows <- list()
expanded_topology_sensitivity_rows <- list()
expanded_panel_rows <- list()

if (length(expanded_mode_defs) > 0L) {
  signal_lookup <- setNames(gene_stat$signal_lfc, gene_stat$gene)

  for (set_id in names(expanded_mode_defs)) {
    def <- expanded_mode_defs[[set_id]]
    expanded_genes <- intersect(def$genes, gene_stat$gene)
    primary_genes <- intersect(
      get_primary_mode_genes(def$md),
      gene_stat$gene
    )
    total_assigned <- intersect(
      get_mode_genes(def$md),
      gene_stat$gene
    )
    background_genes <- setdiff(
      gene_stat$gene,
      union(total_assigned, pool_chip)
    )
    chip_only_genes <- setdiff(
      intersect(pool_chip, gene_stat$gene),
      total_assigned
    )

    comparison_groups <- list(
      Expanded_vs_Background = background_genes,
      Expanded_vs_Primary = primary_genes,
      Expanded_vs_ChIPseekerOnly = chip_only_genes
    )

    for (comparison_name in names(comparison_groups)) {
      reference_genes <- comparison_groups[[comparison_name]]
      if (length(expanded_genes) < 5L ||
          length(reference_genes) < 5L) {
        next
      }

      x <- unname(signal_lookup[expanded_genes])
      y <- unname(signal_lookup[reference_genes])
      rb <- expanded_rank_biserial(
        x,
        y,
        seed = bootstrap_seed + match(
          comparison_name,
          names(comparison_groups)
        )
      )
      wt <- suppressWarnings(
        stats::wilcox.test(x, y, exact = FALSE)
      )

      expanded_es_rows[[paste(set_id, comparison_name, sep = "::")]] <-
        data.frame(
          Mode = set_id,
          Expanded_Set_ID = def$expanded_set_id,
          Representative_Mode = def$representative_mode,
          Equivalent_Modes = def$equivalent_modes,
          Pipeline = def$pipeline,
          Map = def$map,
          Filled = def$filled,
          Comparison = comparison_name,
          N_Expanded = length(expanded_genes),
          N_Reference = length(reference_genes),
          Median_Expanded_signal_lfc = stats::median(
            x,
            na.rm = TRUE
          ),
          Median_Reference_signal_lfc = stats::median(
            y,
            na.rm = TRUE
          ),
          Median_Difference = stats::median(x, na.rm = TRUE) -
            stats::median(y, na.rm = TRUE),
          Rank_Biserial = rb$est,
          Subsampling_Low = rb$lo,
          Subsampling_High = rb$hi,
          Interval_Type = "fixed-N subsampling interval; max 500 genes/group",
          Wilcoxon_P = wt$p.value,
          stringsAsFactors = FALSE
        )


      if (exists("expr_lookup") && length(expr_lookup) > 100L) {
        matched_result <- expanded_expression_match(
          expanded_genes = expanded_genes,
          reference_genes = reference_genes,
          signal_lookup = signal_lookup,
          expr_lookup = expr_lookup,
          iterations = expanded_expression_match_iterations,
          bins = expanded_expression_match_bins,
          max_per_bin = expanded_expression_match_max_per_bin,
          seed = bootstrap_seed +
            match(comparison_name, names(comparison_groups)) +
            match(set_id, names(expanded_mode_defs)) * 100L
        )
        if (!is.null(matched_result)) {
          expanded_expr_matched_rows[[paste(
            set_id, comparison_name, sep = "::"
          )]] <- dplyr::bind_cols(
            data.frame(
              Mode = set_id,
              Expanded_Set_ID = def$expanded_set_id,
              Representative_Mode = def$representative_mode,
              Pipeline = def$pipeline,
              Map = def$map,
              Filled = def$filled,
              Comparison = comparison_name,
              stringsAsFactors = FALSE
            ),
            matched_result
          )
        }
      }
    }


    topology_mode <- expanded_topology_df %>%
      dplyr::filter(Mode == set_id, Gene %in% expanded_genes)
    topology_subsets <- list(
      All_Expanded = expanded_genes,
      MultiSupported = topology_mode$Gene[topology_mode$Multi_Supported]
    )
    if (nrow(topology_mode) > 0L &&
        any(topology_mode$Topology_Available %in% TRUE)) {
      topology_subsets$NonHub <- topology_mode$Gene[
        !topology_mode$Is_Top5Pct_Hub
      ]
      topology_subsets$NonHub_MultiSupported <- topology_mode$Gene[
        !topology_mode$Is_Top5Pct_Hub & topology_mode$Multi_Supported
      ]
    }
    for (subset_name in names(topology_subsets)) {
      subset_genes <- intersect(
        unique(topology_subsets[[subset_name]]),
        names(signal_lookup)
      )
      if (length(subset_genes) < 5L || length(background_genes) < 5L) next
      subset_interval <- expanded_rank_biserial(
        unname(signal_lookup[subset_genes]),
        unname(signal_lookup[background_genes]),
        seed = bootstrap_seed + match(subset_name, names(topology_subsets)) +
          match(set_id, names(expanded_mode_defs)) * 1000L
      )
      expanded_topology_sensitivity_rows[[paste(
        set_id, subset_name, sep = "::"
      )]] <- data.frame(
        Mode = set_id,
        Pipeline = def$pipeline,
        Subset = subset_name,
        N_Genes = length(subset_genes),
        N_Background = length(background_genes),
        Rank_Biserial_vs_Background = subset_interval$est,
        Subsampling_Low = subset_interval$lo,
        Subsampling_High = subset_interval$hi,
        Median_signal_lfc = stats::median(
          unname(signal_lookup[subset_genes]), na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }

    expanded_panel_rows[[set_id]] <- data.frame(
      Mode = set_id,
      Expanded_Set_ID = def$expanded_set_id,
      Representative_Mode = def$representative_mode,
      N_Expanded = length(def$genes),
      BRD4_Core_N = sum(
        def$genes %in% toupper(brd4_core_panel)
      ),
      BRD4_Context_N = sum(
        def$genes %in% toupper(brd4_context_panel)
      ),
      BRD4_Associated_N = sum(
        def$genes %in% toupper(brd4_associated_panel)
      ),
      Any_Panel_N = sum(
        def$genes %in% unique(toupper(c(
          brd4_core_panel,
          brd4_context_panel,
          brd4_associated_panel
        )))
      ),
      BRD4_Core_Fraction = sum(
        def$genes %in% toupper(brd4_core_panel)
      ) / max(1L, length(def$genes)),
      BRD4_Context_Fraction = sum(
        def$genes %in% toupper(brd4_context_panel)
      ) / max(1L, length(def$genes)),
      BRD4_Associated_Fraction = sum(
        def$genes %in% toupper(brd4_associated_panel)
      ) / max(1L, length(def$genes)),
      Any_Panel_Fraction = sum(
        def$genes %in% unique(toupper(c(
          brd4_core_panel,
          brd4_context_panel,
          brd4_associated_panel
        )))
      ) / max(1L, length(def$genes)),
      stringsAsFactors = FALSE
    )
  }
}

expanded_es_df <- dplyr::bind_rows(expanded_es_rows)
if (ncol(expanded_es_df) == 0L) {
  expanded_es_df <- data.frame(
    Mode = character(),
    Expanded_Set_ID = character(),
    Representative_Mode = character(),
    Equivalent_Modes = character(),
    Pipeline = character(),
    Map = character(),
    Filled = logical(),
    Comparison = character(),
    N_Expanded = integer(),
    N_Reference = integer(),
    Median_Expanded_signal_lfc = numeric(),
    Median_Reference_signal_lfc = numeric(),
    Median_Difference = numeric(),
    Rank_Biserial = numeric(),
    Subsampling_Low = numeric(),
    Subsampling_High = numeric(),
    Interval_Type = character(),
    Wilcoxon_P = numeric(),
    Wilcoxon_FDR = numeric(),
    stringsAsFactors = FALSE
  )
}
expanded_panel_df <- dplyr::bind_rows(expanded_panel_rows)
if (ncol(expanded_panel_df) == 0L) {
  expanded_panel_df <- data.frame(
    Mode = character(),
    Expanded_Set_ID = character(),
    Representative_Mode = character(),
    N_Expanded = integer(),
    BRD4_Core_N = integer(),
    BRD4_Context_N = integer(),
    BRD4_Associated_N = integer(),
    Any_Panel_N = integer(),
    BRD4_Core_Fraction = numeric(),
    BRD4_Context_Fraction = numeric(),
    BRD4_Associated_Fraction = numeric(),
    Any_Panel_Fraction = numeric(),
    N_Expanded_In_Universe = integer(),
    BRD4_Core_Universe_N = integer(),
    BRD4_Context_Universe_N = integer(),
    BRD4_Associated_Universe_N = integer(),
    Any_Panel_Universe_N = integer(),
    BRD4_Core_P = numeric(),
    BRD4_Context_P = numeric(),
    BRD4_Associated_P = numeric(),
    Any_Panel_P = numeric(),
    BRD4_Core_FDR = numeric(),
    BRD4_Context_FDR = numeric(),
    BRD4_Associated_FDR = numeric(),
    Any_Panel_FDR = numeric(),
    stringsAsFactors = FALSE
  )
}

if (nrow(expanded_es_df) > 0L) {
  expanded_es_df <- expanded_es_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::mutate(
      Wilcoxon_FDR = stats::p.adjust(Wilcoxon_P, method = "BH")
    ) %>%
    dplyr::ungroup()
}

expanded_expr_matched_df <- dplyr::bind_rows(expanded_expr_matched_rows)
if (ncol(expanded_expr_matched_df) == 0L) {
  expanded_expr_matched_df <- data.frame(
    Mode = character(), Expanded_Set_ID = character(),
    Representative_Mode = character(), Pipeline = character(),
    Map = character(), Filled = logical(), Comparison = character(),
    Full_ES = numeric(), Matched_ES_Median = numeric(),
    Matched_ES_Low = numeric(), Matched_ES_High = numeric(),
    P_Matched_ES_gt_0 = numeric(), P_Matched_ES_lt_0 = numeric(),
    N_Expanded_Available = integer(), N_Reference_Available = integer(),
    Median_N_Matched = numeric(), Median_Expanded_TPM = numeric(),
    Median_Reference_TPM = numeric(), SMD_Before = numeric(),
    SMD_After_Median = numeric(), SMD_After_Low = numeric(),
    SMD_After_High = numeric(), P_Abs_SMD_lt_0_1 = numeric(),
    Requested_Matching_Iterations = integer(),
    Matching_Iterations = integer(),
    Successful_Matching_Fraction = numeric(),
    Effective_Expression_Bins = integer(),
    stringsAsFactors = FALSE
  )
}

expanded_topology_sensitivity_df <- dplyr::bind_rows(
  expanded_topology_sensitivity_rows
)
if (ncol(expanded_topology_sensitivity_df) == 0L) {
  expanded_topology_sensitivity_df <- data.frame(
    Mode = character(), Pipeline = character(), Subset = character(),
    N_Genes = integer(), N_Background = integer(),
    Rank_Biserial_vs_Background = numeric(),
    Subsampling_Low = numeric(), Subsampling_High = numeric(),
    Median_signal_lfc = numeric(), stringsAsFactors = FALSE
  )
}

if (nrow(expanded_panel_df) > 0L) {
  panel_universe <- unique(all_genes)
  panel_sets <- list(
    BRD4_Core = intersect(toupper(brd4_core_panel), panel_universe),
    BRD4_Context = intersect(toupper(brd4_context_panel), panel_universe),
    BRD4_Associated = intersect(toupper(brd4_associated_panel), panel_universe),
    Any_Panel = intersect(unique(toupper(c(
      brd4_core_panel, brd4_context_panel, brd4_associated_panel
    ))), panel_universe)
  )
  expanded_panel_df$N_Expanded_In_Universe <- vapply(
    expanded_panel_df$Mode,
    function(mode_name) {
      length(intersect(
        expanded_mode_defs[[mode_name]]$genes,
        panel_universe
      ))
    },
    integer(1)
  )
  for (panel_name in names(panel_sets)) {
    universe_hit_col <- paste0(panel_name, "_Universe_N")
    p_col <- paste0(panel_name, "_P")
    fdr_col <- paste0(panel_name, "_FDR")
    panel_size <- length(panel_sets[[panel_name]])
    expanded_panel_df[[universe_hit_col]] <- vapply(
      expanded_panel_df$Mode,
      function(mode_name) {
        mode_genes_in_universe <- intersect(
          expanded_mode_defs[[mode_name]]$genes,
          panel_universe
        )
        sum(mode_genes_in_universe %in% panel_sets[[panel_name]])
      },
      integer(1)
    )
    expanded_panel_df[[p_col]] <- vapply(
      seq_len(nrow(expanded_panel_df)),
      function(row_id) {
        stats::phyper(
          expanded_panel_df[[universe_hit_col]][row_id] - 1L,
          panel_size,
          max(0L, length(panel_universe) - panel_size),
          expanded_panel_df$N_Expanded_In_Universe[row_id],
          lower.tail = FALSE
        )
      },
      numeric(1)
    )
    expanded_panel_df[[fdr_col]] <- stats::p.adjust(
      expanded_panel_df[[p_col]], method = "BH"
    )
  }
}

utils::write.csv(
  expanded_es_df,
  file.path(expanded_dir, "ExpandedOnly_EffectSize.csv"),
  row.names = FALSE
)
utils::write.csv(
  expanded_expr_matched_df,
  file.path(
    expanded_dir,
    "ExpandedOnly_ExpressionMatched_EffectSize.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  expanded_topology_sensitivity_df,
  file.path(
    expanded_dir,
    "ExpandedOnly_Topology_Sensitivity.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  expanded_panel_df,
  file.path(expanded_dir, "ExpandedOnly_BRD4_Panel_Overlap.csv"),
  row.names = FALSE
)

if (nrow(expanded_es_df) > 0L) {
  expanded_es_df$Comparison <- factor(
    expanded_es_df$Comparison,
    levels = c(
      "Expanded_vs_Background",
      "Expanded_vs_Primary",
      "Expanded_vs_ChIPseekerOnly"
    )
  )

  p_expanded_es <- ggplot(
    expanded_es_df,
    aes(
      x = Rank_Biserial,
      y = reorder(Mode, Rank_Biserial),
      color = Pipeline
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_errorbarh(
      aes(xmin = Subsampling_Low, xmax = Subsampling_High),
      height = 0.25,
      linewidth = 0.3
    ) +
    geom_point(size = 1.8) +
    facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c(
      "Basic" = "#1F77B4",
      "E-refined" = "#9467BD",
      "I-refined" = "#E64B35",
      "C-refined" = "#FFA500"
    )) +
    labs(
      title = "Expanded-only target effect sizes",
      subtitle = paste0(
        "Positive rank-biserial values indicate stronger downregulation. ",
        "Bars are fixed-N subsampling intervals (not full-set bootstrap CIs); ",
        "Wilcoxon FDR is exported in the result table."
      ),
      x = "Rank-biserial effect size",
      y = "Hop1 mode",
      color = "Pipeline"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 9,
        color = "grey40",
        hjust = 0.5
      ),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )

  queueable_ggsave(
    file.path(expanded_dir, "ExpandedOnly_EffectSize_Forest.pdf"),
    p_expanded_es,
    width = 9,
    height = max(7, 0.28 * nrow(expanded_es_df) + 2)
  )
}

if (nrow(expanded_expr_matched_df) > 0L) {
  p_expanded_matched <- ggplot(
    expanded_expr_matched_df,
    aes(
      x = Full_ES,
      y = Matched_ES_Median,
      color = Pipeline
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
    geom_errorbar(
      aes(ymin = Matched_ES_Low, ymax = Matched_ES_High),
      width = 0.01,
      alpha = 0.4
    ) +
    geom_point(size = 2, alpha = 0.75) +
    facet_wrap(~ Comparison, scales = "free") +
    scale_color_manual(values = c(
      "Basic" = "#1F77B4",
      "E-refined" = "#9467BD",
      "I-refined" = "#E64B35",
      "C-refined" = "#FFA500"
    )) +
    labs(
      title = "Expanded-only expression-matched effect-size sensitivity",
      subtitle = paste0(
        "DMSO TPM decile matching; error bars are the 2.5%-97.5% ",
        "range across matched iterations"
      ),
      x = "Full-set rank-biserial effect size",
      y = "Expression-matched effect size (median)",
      color = "Pipeline"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
      legend.position = "bottom"
    )
  queueable_ggsave(
    file.path(
      expanded_dir,
      "ExpandedOnly_ExpressionMatched_EffectSize.pdf"
    ),
    p_expanded_matched,
    width = 9,
    height = 7
  )
}

if (nrow(expanded_topology_sensitivity_df) > 0L) {
  expanded_topology_sensitivity_df$Subset <- factor(
    expanded_topology_sensitivity_df$Subset,
    levels = c(
      "All_Expanded", "NonHub", "MultiSupported",
      "NonHub_MultiSupported"
    )
  )
  p_expanded_topology <- ggplot(
    expanded_topology_sensitivity_df,
    aes(
      x = Rank_Biserial_vs_Background,
      y = reorder(Mode, Rank_Biserial_vs_Background),
      color = Pipeline
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = Subsampling_Low, xmax = Subsampling_High),
      height = 0.25,
      linewidth = 0.3
    ) +
    geom_point(aes(size = N_Genes), alpha = 0.8) +
    facet_wrap(~ Subset, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c(
      "Basic" = "#1F77B4",
      "E-refined" = "#9467BD",
      "I-refined" = "#E64B35",
      "C-refined" = "#FFA500"
    )) +
    labs(
      title = "Expanded-only topology/support sensitivity",
      subtitle = paste0(
        "NonHub removes genes linked through top-5% degree anchors; ",
        "MultiSupported requires >=2 peaks or >=2 anchors"
      ),
      x = "Rank-biserial effect size vs background",
      y = "Hop1 mode",
      color = "Pipeline",
      size = "Genes"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
      legend.position = "bottom",
      axis.text.y = element_text(size = 7)
    )
  queueable_ggsave(
    file.path(expanded_dir, "ExpandedOnly_Topology_Sensitivity.pdf"),
    p_expanded_topology,
    width = 10,
    height = max(8, 0.2 * nrow(expanded_topology_sensitivity_df) + 3)
  )
}

# ── Expanded-only GO enrichment ─────────────────────────────────────────────
expanded_go_results <- list()
expanded_go_status <- list()

if (length(expanded_set_defs) > 0L) {
  for (set_id in names(expanded_set_defs)) {
    def <- expanded_set_defs[[set_id]]
    symbols <- intersect(def$genes, all_genes)
    entrez <- tryCatch({
      mapped <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = symbols,
        columns = "ENTREZID",
        keytype = "SYMBOL"
      )
      unique(stats::na.omit(mapped$ENTREZID))
    }, error = function(e) {
      warning(
        "Expanded-only GO mapping failed for ",
        set_id,
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
      character()
    })

    if (length(entrez) < gsea_min_size) {
      expanded_go_status[[paste(set_id, "mapping", sep = "::")]] <-
        data.frame(
          SampleSize = NA_integer_,
          Module = "expanded_only_go",
          Mode = set_id,
          Iteration = NA_integer_,
          GeneSet = "expanded_only",
          Ontology = NA_character_,
          N_Entrez = length(entrez),
          Status = "insufficient_input",
          ExpectedTerms = NA_character_,
          ReturnedTerms = NA_character_,
          MissingTerms = NA_character_,
          Error = paste0(
            "Mapped Entrez IDs < ",
            gsea_min_size
          ),
          stringsAsFactors = FALSE
        )
      next
    }

    for (ontology in c("BP", "CC", "MF")) {
      go_out <- safe_go_once(
        entrez,
        ontology,
        NA_integer_,
        set_id,
        "expanded_only",
        gouniv
      )
      expanded_go_status[[paste(set_id, ontology, sep = "::")]] <-
        go_out$status

      if (!is.null(go_out$result) &&
          nrow(go_out$result) > 0L) {
        result <- go_out$result
        result$Expanded_Set_ID <- set_id
        result$Representative_Mode <- def$representative_mode
        result$Equivalent_Modes <- def$equivalent_modes
        result$Ontology <- ontology
        expanded_go_results[[paste(set_id, ontology, sep = "::")]] <-
          result
      }
    }
  }
}

expanded_go_df <- dplyr::bind_rows(expanded_go_results)
if (ncol(expanded_go_df) == 0L) {
  expanded_go_df <- data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    Expanded_Set_ID = character(),
    Representative_Mode = character(),
    Equivalent_Modes = character(),
    Ontology = character(),
    stringsAsFactors = FALSE
  )
}
expanded_go_status_df <- dplyr::bind_rows(expanded_go_status)
if (ncol(expanded_go_status_df) == 0L) {
  expanded_go_status_df <- data.frame(
    SampleSize = integer(),
    Module = character(),
    Mode = character(),
    Iteration = integer(),
    GeneSet = character(),
    Ontology = character(),
    N_Entrez = integer(),
    Status = character(),
    ExpectedTerms = character(),
    ReturnedTerms = character(),
    MissingTerms = character(),
    Error = character(),
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  expanded_go_df,
  file.path(expanded_dir, "ExpandedOnly_GO_Results.csv"),
  row.names = FALSE
)
utils::write.csv(
  expanded_go_status_df,
  file.path(expanded_dir, "ExpandedOnly_GO_Status.csv"),
  row.names = FALSE
)

if (nrow(expanded_go_df) > 0L) {
  expanded_go_top <- expanded_go_df %>%
    dplyr::filter(Ontology == "BP", is.finite(p.adjust)) %>%
    dplyr::group_by(Expanded_Set_ID) %>%
    dplyr::slice_min(
      order_by = p.adjust,
      n = 5L,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      NegLog10FDR = -log10(pmax(p.adjust, 1e-300)),
      Label = paste0(Expanded_Set_ID, " | ", Description)
    )

  if (nrow(expanded_go_top) > 0L) {
    p_expanded_go <- ggplot(
      expanded_go_top,
      aes(
        x = NegLog10FDR,
        y = reorder(Label, NegLog10FDR),
        size = Count,
        color = Expanded_Set_ID
      )
    ) +
      geom_point(alpha = 0.8) +
      labs(
        title = "Expanded-only GO biological-process enrichment",
        subtitle = "Top five BP terms per unique expanded-only set",
        x = "-log10 adjusted P value",
        y = NULL,
        color = "Set",
        size = "Genes"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 12,
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          size = 9,
          color = "grey40",
          hjust = 0.5
        ),
        legend.position = "bottom",
        axis.text.y = element_text(size = 7)
      )

    queueable_ggsave(
      file.path(expanded_dir, "ExpandedOnly_GO_BP_TopTerms.pdf"),
      p_expanded_go,
      width = 10,
      height = max(6, 0.20 * nrow(expanded_go_top) + 2)
    )
  }
}

wb_expanded <- createWorkbook()
addWorksheet(wb_expanded, "Set_Catalog")
writeData(wb_expanded, "Set_Catalog", as.data.frame(expanded_catalog))
addWorksheet(wb_expanded, "Genes")
writeData(wb_expanded, "Genes", as.data.frame(expanded_gene_table))
addWorksheet(wb_expanded, "Effect_Size")
writeData(wb_expanded, "Effect_Size", as.data.frame(expanded_es_df))
addWorksheet(wb_expanded, "Expression_Matched_ES")
writeData(
  wb_expanded,
  "Expression_Matched_ES",
  as.data.frame(expanded_expr_matched_df)
)
addWorksheet(wb_expanded, "Topology_Support")
writeData(
  wb_expanded,
  "Topology_Support",
  as.data.frame(expanded_topology_df)
)
addWorksheet(wb_expanded, "Topology_Sensitivity")
writeData(
  wb_expanded,
  "Topology_Sensitivity",
  as.data.frame(expanded_topology_sensitivity_df)
)
addWorksheet(wb_expanded, "GSEA_Sampling_Plan")
writeData(
  wb_expanded,
  "GSEA_Sampling_Plan",
  as.data.frame(expanded_gsea_plan)
)
addWorksheet(wb_expanded, "BRD4_Panel")
writeData(wb_expanded, "BRD4_Panel", as.data.frame(expanded_panel_df))
addWorksheet(wb_expanded, "GO_Results")
writeData(wb_expanded, "GO_Results", as.data.frame(expanded_go_df))
addWorksheet(wb_expanded, "GO_Status")
writeData(wb_expanded, "GO_Status", as.data.frame(expanded_go_status_df))
saveWorkbook(
  wb_expanded,
  file.path(expanded_dir, "ExpandedOnly_Analysis.xlsx"),
  overwrite = TRUE
)

# ── Expanded-only stability across TSS windows ──────────────────────────────
expanded_tss_counts <- list()
expanded_tss_sets <- list()

for (tss_name in names(cfg$tss_subdirs)) {
  rdata_path <- file.path(
    cfg$rdata_base,
    cfg$tss_subdirs[[tss_name]]
  )
  if (!file.exists(rdata_path)) next

  tss_env <- new.env(parent = emptyenv())
  loaded_tss <- load(rdata_path, envir = tss_env, verbose = FALSE)
  hop1_specs <- mode_spec_defs[
    vapply(
      mode_spec_defs,
      function(x) grepl("_hop1$", x$name),
      logical(1)
    )
  ]

  for (spec in hop1_specs) {
    if (!spec$ann_source %in% loaded_tss) next
    md <- list(
      name = spec$name,
      ann = tss_env[[spec$ann_source]],
      map = spec$map,
      fill = spec$fill,
      near = FALSE
    )
    genes <- get_expanded_only_genes(md)
    key <- paste(tss_name, spec$name, sep = "::")
    expanded_tss_sets[[key]] <- genes
    expanded_tss_counts[[key]] <- data.frame(
      TSS = tss_name,
      Mode = spec$name,
      Pipeline = infer_pipeline_label(spec$name),
      Map = spec$map,
      Filled = isTRUE(spec$fill),
      Expanded_Only_N = length(genes),
      stringsAsFactors = FALSE
    )
  }

  rm(tss_env)
  invisible(gc())
}

expanded_tss_count_df <- dplyr::bind_rows(expanded_tss_counts)
if (ncol(expanded_tss_count_df) == 0L) {
  expanded_tss_count_df <- data.frame(
    TSS = character(),
    Mode = character(),
    Pipeline = character(),
    Map = character(),
    Filled = logical(),
    Expanded_Only_N = integer(),
    stringsAsFactors = FALSE
  )
}
expanded_tss_jaccard_rows <- list()
expanded_tss_stability_rows <- list()

if (nrow(expanded_tss_count_df) > 0L) {
  available_tss <- unique(expanded_tss_count_df$TSS)
  available_modes <- unique(expanded_tss_count_df$Mode)

  if (length(available_tss) >= 2L) {
    for (mode_name in available_modes) {
      mode_tss <- expanded_tss_count_df$TSS[
        expanded_tss_count_df$Mode == mode_name
      ]
      mode_tss <- unique(mode_tss)

      if (length(mode_tss) >= 2L) {
        pairs <- utils::combn(mode_tss, 2L, simplify = FALSE)
        for (pair in pairs) {
          g1 <- expanded_tss_sets[[paste(pair[1L], mode_name, sep = "::")]]
          g2 <- expanded_tss_sets[[paste(pair[2L], mode_name, sep = "::")]]
          expanded_tss_jaccard_rows[[paste(
            mode_name,
            pair[1L],
            pair[2L],
            sep = "::"
          )]] <- data.frame(
            Mode = mode_name,
            TSS_A = pair[1L],
            TSS_B = pair[2L],
            N_A = length(g1),
            N_B = length(g2),
            Intersection = length(intersect(g1, g2)),
            Union = length(union(g1, g2)),
            Jaccard = length(intersect(g1, g2)) /
              max(1L, length(union(g1, g2))),
            stringsAsFactors = FALSE
          )
        }
      }

      gene_occurrences <- unlist(lapply(mode_tss, function(tss_name) {
        expanded_tss_sets[[paste(tss_name, mode_name, sep = "::")]]
      }), use.names = FALSE)
      if (length(gene_occurrences) > 0L) {
        freq <- sort(table(gene_occurrences), decreasing = TRUE)
        expanded_tss_stability_rows[[mode_name]] <- data.frame(
          Mode = mode_name,
          Gene = names(freq),
          Windows_Present = as.integer(freq),
          Total_Windows = length(mode_tss),
          Stability_Fraction = as.integer(freq) / length(mode_tss),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

expanded_tss_jaccard_df <- dplyr::bind_rows(
  expanded_tss_jaccard_rows
)
if (ncol(expanded_tss_jaccard_df) == 0L) {
  expanded_tss_jaccard_df <- data.frame(
    Mode = character(),
    TSS_A = character(),
    TSS_B = character(),
    N_A = integer(),
    N_B = integer(),
    Intersection = integer(),
    Union = integer(),
    Jaccard = numeric(),
    stringsAsFactors = FALSE
  )
}
expanded_tss_stability_df <- dplyr::bind_rows(
  expanded_tss_stability_rows
)
if (ncol(expanded_tss_stability_df) == 0L) {
  expanded_tss_stability_df <- data.frame(
    Mode = character(),
    Gene = character(),
    Windows_Present = integer(),
    Total_Windows = integer(),
    Stability_Fraction = numeric(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  expanded_tss_count_df,
  file.path(expanded_dir, "ExpandedOnly_CrossTSS_Counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  expanded_tss_jaccard_df,
  file.path(expanded_dir, "ExpandedOnly_CrossTSS_Jaccard.csv"),
  row.names = FALSE
)
utils::write.csv(
  expanded_tss_stability_df,
  file.path(expanded_dir, "ExpandedOnly_CrossTSS_GeneStability.csv"),
  row.names = FALSE
)

if (nrow(expanded_tss_count_df) > 0L) {
  p_expanded_tss_count <- ggplot(
    expanded_tss_count_df,
    aes(
      x = TSS,
      y = Mode,
      fill = Expanded_Only_N
    )
  ) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "C") +
    labs(
      title = "Expanded-only gene counts across TSS windows",
      x = "TSS window",
      y = NULL,
      fill = "Genes"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 7)
    )

  queueable_ggsave(
    file.path(expanded_dir, "ExpandedOnly_CrossTSS_Counts.pdf"),
    p_expanded_tss_count,
    width = 8,
    height = 7
  )
}

if (nrow(expanded_tss_jaccard_df) > 0L) {
  expanded_tss_jaccard_plot_df <- expanded_tss_jaccard_df %>%
    dplyr::mutate(
      Pipeline = infer_pipeline_label(Mode),
      Pair = paste(TSS_A, TSS_B, sep = " vs ")
    )

  expanded_tss_jaccard_median <- expanded_tss_jaccard_plot_df %>%
    dplyr::group_by(Mode) %>%
    dplyr::summarise(
      Jaccard = stats::median(Jaccard, na.rm = TRUE),
      .groups = "drop"
    )

  p_expanded_tss_jaccard <- ggplot(
    expanded_tss_jaccard_plot_df,
    aes(
      x = Jaccard,
      y = reorder(Mode, Jaccard, median),
      color = Pipeline
    )
  ) +
    geom_point(alpha = 0.65, size = 1.4) +
    geom_point(
      data = expanded_tss_jaccard_median,
      aes(x = Jaccard, y = Mode),
      inherit.aes = FALSE,
      shape = 18,
      size = 2.7,
      color = "black"
    ) +
    scale_color_manual(values = c(
      "Basic" = "#1F77B4",
      "E-refined" = "#9467BD",
      "I-refined" = "#E64B35",
      "C-refined" = "#FFA500"
    )) +
    labs(
      title = "Expanded-only cross-TSS stability",
      subtitle = "Points are pairwise TSS-window Jaccard values; black diamonds are medians",
      x = "Pairwise Jaccard",
      y = NULL,
      color = "Pipeline"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 9,
        color = "grey40",
        hjust = 0.5
      ),
      legend.position = "bottom",
      axis.text.y = element_text(size = 7)
    )

  queueable_ggsave(
    file.path(expanded_dir, "ExpandedOnly_CrossTSS_Jaccard.pdf"),
    p_expanded_tss_jaccard,
    width = 8,
    height = 7
  )
}

wb_expanded_tss <- createWorkbook()
addWorksheet(wb_expanded_tss, "Counts")
writeData(
  wb_expanded_tss,
  "Counts",
  as.data.frame(expanded_tss_count_df)
)
addWorksheet(wb_expanded_tss, "Pairwise_Jaccard")
writeData(
  wb_expanded_tss,
  "Pairwise_Jaccard",
  as.data.frame(expanded_tss_jaccard_df)
)
addWorksheet(wb_expanded_tss, "Gene_Stability")
writeData(
  wb_expanded_tss,
  "Gene_Stability",
  as.data.frame(expanded_tss_stability_df)
)
saveWorkbook(
  wb_expanded_tss,
  file.path(expanded_dir, "ExpandedOnly_CrossTSS.xlsx"),
  overwrite = TRUE
)

expanded_module_summary <- c(
  paste0("Completed: ", Sys.time()),
  paste0("Hop1 modes: ", length(hop1_modes)),
  paste0("Unique non-empty expanded-only sets: ", length(expanded_set_defs)),
  paste0("Unique expanded-only GSEA analyses: ", length(expanded_gsea_defs)),
  paste0(
    "Total mode-specific expanded-only genes: ",
    sum(expanded_catalog$Expanded_Only_N, na.rm = TRUE)
  ),
  paste0(
    "GSEA errors: ",
    sum(expanded_gsea_status_df$Status == "error", na.rm = TRUE)
  ),
  paste0(
    "GO errors: ",
    sum(expanded_go_status_df$Status == "error", na.rm = TRUE)
  ),
  paste0("Effect-size rows: ", nrow(expanded_es_df)),
  paste0("Expression-matched effect-size rows: ", nrow(expanded_expr_matched_df)),
  paste0("Topology support rows: ", nrow(expanded_topology_df)),
  paste0("Topology sensitivity rows: ", nrow(expanded_topology_sensitivity_df)),
  paste0(
    "Duplicated Actual_N plan rows: ",
    sum(expanded_gsea_plan$Duplicate_Actual_N, na.rm = TRUE)
  )
)
writeLines(
  expanded_module_summary,
  file.path(expanded_dir, "ExpandedOnly_Module_Status.txt")
)
message("    Expanded-only outputs: ", expanded_dir)


# ── Convergence QC (CMA plots, grouped by sample size) ──────────────────────
message("\n\n>>> ==== Convergence QC ====")

make_cma <- function(vec, n_full) {
  v <- rep(NA_real_, n_full)
  n <- min(length(vec), n_full)
  if (n > 0) v[seq_len(n)] <- vec[seq_len(n)]
  cma <- rep(NA_real_, n_full)
  if (sum(!is.na(v)) > 0) {
    cs <- cumsum(replace(v, is.na(v), 0))
    cc <- cumsum(!is.na(v))
    cma <- cs / pmax(1, cc)
    cma[cc == 0] <- NA_real_
  }
  cma
}

conv_dir <- file.path(base_out_dir, "convergence_plots")
dir.create(conv_dir, recursive = TRUE, showWarnings = FALSE)

# One panel plot per sample size: overlay all mode CMA traces with transparency
for (sz in sample_sizes) {
  sim_dir <- file.path(base_out_dir, paste0("size", sz), "tmp_sim")
  if (!dir.exists(sim_dir)) next
  csv_files <- list.files(sim_dir, pattern = "sim_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) next

  conv_all <- list()
  for (f in csv_files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", f); NULL })
    if (is.null(df) || nrow(df) == 0) next
    mode_name <- gsub("^sim_|\\.csv$", "", basename(f))
    n_iter <- min(nrow(df), 500)
    sub <- df %>% filter(ID == "looplook") %>% arrange(Iteration) %>% head(n_iter)
    if (nrow(sub) < 10) next
    cma <- make_cma(sub$NES, nrow(sub))
    conv_all[[length(conv_all) + 1]] <- data.frame(
      Iteration = seq_len(nrow(sub)), CMA = cma, Mode = mode_name, stringsAsFactors = FALSE)
  }
  if (length(conv_all) == 0) next

  conv_df <- bind_rows(conv_all)
  n_modes <- length(unique(conv_df$Mode))
  conv_df <- conv_df %>%
    mutate(Pipeline = case_when(
      grepl("^anno_", Mode)        ~ "Basic",
      grepl("^refined_", Mode)     ~ "E-refined",
      grepl("^chrom_only_", Mode)  ~ "C-refined",
      grepl("^chrom_", Mode)       ~ "I-refined",
      grepl("^ctrl_", Mode)        ~ "Control",
      TRUE ~ "Other"))

  p <- ggplot(conv_df, aes(x = Iteration, y = CMA, group = Mode, color = Pipeline)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.2, alpha = 0.45) +
    scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                   "I-refined" = "#E64B35", "C-refined" = "#FFA500",
                                   "Control" = "grey50")) +
    labs(title = sprintf("CMA Convergence: all modes overlaid (sample size = %d)", sz),
         subtitle = sprintf("%d modes, %d iterations | Thin transparent lines show individual mode stability",
                            n_modes, n_iterations),
         x = "Iteration", y = "Cumulative Moving Average of NES") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "bottom", legend.title = element_blank(),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
  queueable_ggsave(file.path(conv_dir, sprintf("Convergence_CMA_size%d.pdf", sz)),
         p, width = 8, height = 4)
  message(sprintf("    size=%d: %d mode traces plotted", sz, n_modes))
}
message("    Convergence plots saved to: ", conv_dir)
# ══════════════════════════════════════════════════════════════════════════════
# Cross-Sample-Size Summary: NES stability across sample sizes 100/200/500
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Cross-Sample-Size Stability Summary ====")
xss_dir <- file.path(base_out_dir, "cross_size_summary")
dir.create(xss_dir, recursive = TRUE, showWarnings = FALSE)

# Read merged sim_df from each sample size directory
xss_data <- list()
for (sz in sample_sizes) {
  sz_dir <- file.path(base_out_dir, paste0("size", sz))
  sz_tmp <- file.path(sz_dir, "tmp_sim")
  if (!dir.exists(sz_tmp)) next
  csv_files <- list.files(sz_tmp, pattern = "sim_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) next
  sz_list <- lapply(csv_files, function(f) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", f); NULL })
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  })
  sz_list <- sz_list[!vapply(sz_list, is.null, logical(1))]
  if (length(sz_list) > 0) {
    sz_df <- bind_rows(sz_list)
    sz_df$SampleSize <- sz
    xss_data[[as.character(sz)]] <- sz_df
  }
  message(sprintf("    Size %d: %d CSV files, %d rows", sz, length(csv_files),
    if (!is.null(xss_data[[as.character(sz)]])) nrow(xss_data[[as.character(sz)]]) else 0))
}

if (length(xss_data) >= 2) {
  xss_df <- bind_rows(xss_data)
  xss_df$Mode <- factor(xss_df$Mode, levels = sapply(modes, `[[`, "name"))
  xss_df$SampleSize <- factor(xss_df$SampleSize, levels = sample_sizes)

  xss_summary <- xss_df %>%
    filter(!grepl("ctrl", Mode)) %>%
    group_by(Mode, ID, SampleSize) %>%
    summarise(Mean_NES = mean(NES, na.rm = TRUE),
              SD_NES = sd(NES, na.rm = TRUE), N = n(),
              Median_Actual_N = median(Actual_N, na.rm = TRUE),
              Requested_N = Requested_N[1], .groups = "drop")

  # Warn if actual sample size differs significantly from requested
  xss_size_check <- xss_summary %>%
    filter(ID == "looplook") %>%
    mutate(Diff = Requested_N - Median_Actual_N)
  if (any(xss_size_check$Diff > 5, na.rm = TRUE)) {
    big_diff <- xss_size_check %>% filter(Diff > 5) %>%
      dplyr::select(Mode, SampleSize, Requested_N, Median_Actual_N, Diff)
    warning(sprintf("Cross-size: %d mode-size combos have actual N < requested N (max gap=%d). Check pool sizes.",
      nrow(big_diff), max(big_diff$Diff, na.rm = TRUE)))
  }

  xss_loop <- xss_summary %>%
    filter(ID == "looplook") %>%
    mutate(
      Pipeline = case_when(
        grepl("^anno_", Mode) ~ "Basic",
        grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined",
        grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"),
      SizeLabel = paste0(as.character(SampleSize), " (N=", round(Median_Actual_N), ")"))

  p_xss_lines <- ggplot(xss_loop, aes(x = SampleSize, y = Mean_NES,
      group = Mode, color = Pipeline)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.6, alpha = 0.6) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                   "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
    facet_wrap(~ Pipeline, ncol = 4) +
    labs(title = "NES Stability Across Sample Sizes",
         subtitle = sprintf("Each line = one mode (%d total) | Larger sample = more stable",
                            length(unique(xss_loop$Mode))),
         x = "Sample Size (genes per subsample)", y = "Mean NES (looplook)") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "none", strip.text = element_text(face = "bold", size = 10))
  queueable_ggsave(file.path(xss_dir, "CrossSize_NES_Stability.pdf"), p_xss_lines, width = 9, height = 5)

  # ── Hop0/hop1 split: NES stability ──
  xss_loop <- xss_loop %>%
    mutate(HopSize = factor(ifelse(grepl("_hop1$", Mode), "hop1 (primary+expanded)", "hop0 (primary)"),
                            levels = c("hop0 (primary)", "hop1 (primary+expanded)")))

  p_xss_hop <- ggplot(xss_loop, aes(x = SampleSize, y = Mean_NES,
      group = Mode, color = Pipeline)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                   "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
    facet_wrap(~ HopSize, ncol = 2) +
    labs(title = "NES Stability Across Sample Sizes — hop0 vs hop1",
         subtitle = sprintf("%d modes | Each line = one mode", n_distinct(xss_loop$Mode)),
         x = "Sample Size", y = "Mean NES (looplook)") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "bottom", strip.text = element_text(face = "bold", size = 10))
  queueable_ggsave(file.path(xss_dir, "CrossSize_NES_Stability_hop.pdf"), p_xss_hop, width = 10, height = 5)


  # Variability plot: SD of NES decreases with larger sample sizes
  p_xss_sd <- ggplot(xss_loop, aes(x = SampleSize, y = SD_NES, fill = SampleSize)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.8, alpha = 0.7) +
    scale_fill_manual(values = c("100" = "#4DBBD5", "200" = "#00A087", "500" = "#E64B35")) +
    labs(title = "NES Variability Decreases with Larger Sample Sizes",
         subtitle = "SD of NES across iterations per mode",
         x = "Sample Size", y = "SD of NES (looplook)") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "none")
  queueable_ggsave(file.path(xss_dir, "CrossSize_NES_Variability.pdf"), p_xss_sd, width = 5, height = 5)

  # Export
  wb_xss <- createWorkbook()

  addWorksheet(wb_xss, "NES_Summary_by_Size")
  writeData(wb_xss, "NES_Summary_by_Size", as.data.frame(xss_summary))
  out_xss <- file.path(xss_dir, "CrossSize_NES_Summary.xlsx")
  saveWorkbook(wb_xss, out_xss, overwrite = TRUE)
  message("    Cross-sample-size summary saved: ", out_xss)
} else {
  message("    WARNING: Need at least 2 sample sizes for cross-size comparison, found ", length(xss_data))
}



# ══════════════════════════════════════════════════════════════════════════════
# Cross-TSS Window Summary: compare results across 4 TSS window sizes
# ══════════════════════════════════════════════════════════════════════════════
# Rationale: The promoter/distal boundary depends on tss_region in
# annotate_peaks_and_loops(). Server scripts at \u670d\u52a1\u5668\u811a\u672c/ run the full
# looplook pipeline with tss_region = +/-1kb, 2kb, 5kb, 10kb, producing
# separate res_chromatin.RData files per TSS window. This section loads them
# and compares gene sets across windows to assess robustness.
message("\n\n>>> ==== Cross-TSS Window Stability Summary ====")
tss_cross_dir <- file.path(base_out_dir, "cross_tss_summary")
dir.create(tss_cross_dir, recursive = TRUE, showWarnings = FALSE)

# Base path matches the load() in main script (line ~38):
#   load("<r21_brd4>/profile_export/v0.14/res_chromatin.RData", ...)
# Server scripts output to: <r21_brd4>/profile_export/v0.14/tssXXXX/tss_XXXX_res_chromatin.RData
# Use paths from cfg (defined at top of script)
rdata_base <- cfg$rdata_base
tss_configs <- cfg$tss_subdirs

tss_mode_genes <- list()
tss_available <- character(0)
tss_anchor_types <- list()
tss_loop_types <- list()
tss_anchor_types_hop1 <- list()
tss_loop_types_hop1 <- list()

for (tss_name in names(tss_configs)) {
  tss_file <- tss_configs[[tss_name]]
  rdata_file <- file.path(rdata_base, tss_file)
  if (!file.exists(rdata_file)) {
    message(sprintf("    TSS %s: SKIP (not found: %s)", tss_name, rdata_file))
    next
  }
  tmp_env <- new.env(parent = emptyenv())
  loaded_tss_objects <- load(rdata_file, envir = tmp_env, verbose = FALSE)
  missing_tss_objects <- setdiff(required_objects, loaded_tss_objects)
  if (length(missing_tss_objects) > 0L) {
    stop("Cross-TSS RData ", tss_name, " is missing objects: ",
      paste(missing_tss_objects, collapse = ", "))
  }
  message(sprintf("    TSS %s: loaded %s", tss_name, rdata_file))

  # 32 looplook modes: 4 pipelines x 2 hops x 2 maps x 2 fills
  # Excludes 2 ctrl modes (ChIPseeker-only, no looplook pipeline dependency)
  tss_modes <- list(
    # ── Annotation pipeline ──
    list(name = "anno_all_F",             ann = tmp_env$res,          map = "all",      fill = FALSE),
    list(name = "anno_all_T",             ann = tmp_env$res,          map = "all",      fill = TRUE),
    list(name = "anno_promoter_F",        ann = tmp_env$res,          map = "promoter", fill = FALSE),
    list(name = "anno_promoter_T",        ann = tmp_env$res,          map = "promoter", fill = TRUE),
    list(name = "anno_all_F_hop1",        ann = tmp_env$res2,         map = "all",      fill = FALSE),
    list(name = "anno_all_T_hop1",        ann = tmp_env$res2,         map = "all",      fill = TRUE),
    list(name = "anno_promoter_F_hop1",   ann = tmp_env$res2,         map = "promoter", fill = FALSE),
    list(name = "anno_promoter_T_hop1",   ann = tmp_env$res2,         map = "promoter", fill = TRUE),
    # ── Refined pipeline ──
    list(name = "refined_all_F",          ann = tmp_env$refined_res,  map = "all",      fill = FALSE),
    list(name = "refined_all_T",          ann = tmp_env$refined_res,  map = "all",      fill = TRUE),
    list(name = "refined_promoter_F",     ann = tmp_env$refined_res,  map = "promoter", fill = FALSE),
    list(name = "refined_promoter_T",     ann = tmp_env$refined_res,  map = "promoter", fill = TRUE),
    list(name = "refined_all_F_hop1",     ann = tmp_env$refined_res2, map = "all",      fill = FALSE),
    list(name = "refined_all_T_hop1",     ann = tmp_env$refined_res2, map = "all",      fill = TRUE),
    list(name = "refined_promoter_F_hop1",ann = tmp_env$refined_res2, map = "promoter", fill = FALSE),
    list(name = "refined_promoter_T_hop1",ann = tmp_env$refined_res2, map = "promoter", fill = TRUE),
    # ── Chromatin pipeline ──
    list(name = "chrom_all_F",            ann = tmp_env$cr,           map = "all",      fill = FALSE),
    list(name = "chrom_all_T",            ann = tmp_env$cr,           map = "all",      fill = TRUE),
    list(name = "chrom_promoter_F",       ann = tmp_env$cr,           map = "promoter", fill = FALSE),
    list(name = "chrom_promoter_T",       ann = tmp_env$cr,           map = "promoter", fill = TRUE),
    list(name = "chrom_all_F_hop1",       ann = tmp_env$cr2,          map = "all",      fill = FALSE),
    list(name = "chrom_all_T_hop1",       ann = tmp_env$cr2,          map = "all",      fill = TRUE),
    list(name = "chrom_promoter_F_hop1",  ann = tmp_env$cr2,          map = "promoter", fill = FALSE),
    list(name = "chrom_promoter_T_hop1",  ann = tmp_env$cr2,          map = "promoter", fill = TRUE),
    # ── Chromatin-only pipeline ──
    list(name = "chrom_only_all_F",       ann = tmp_env$cr_only,      map = "all",      fill = FALSE),
    list(name = "chrom_only_all_T",       ann = tmp_env$cr_only,      map = "all",      fill = TRUE),
    list(name = "chrom_only_promoter_F",  ann = tmp_env$cr_only,      map = "promoter", fill = FALSE),
    list(name = "chrom_only_promoter_T",  ann = tmp_env$cr_only,      map = "promoter", fill = TRUE),
    list(name = "chrom_only_all_F_hop1",  ann = tmp_env$cr2_only,     map = "all",      fill = FALSE),
    list(name = "chrom_only_all_T_hop1",  ann = tmp_env$cr2_only,     map = "all",      fill = TRUE),
    list(name = "chrom_only_promoter_F_hop1", ann = tmp_env$cr2_only, map = "promoter", fill = FALSE),
    list(name = "chrom_only_promoter_T_hop1", ann = tmp_env$cr2_only, map = "promoter", fill = TRUE))

  mg <- list()
  for (md in tss_modes) {
    mg[[md$name]] <- get_mode_genes(md)
  }
  tss_mode_genes[[tss_name]] <- mg
  tss_available <- unique(c(tss_available, tss_name))
  # ── Extract anchor type & loop type distributions ──
  if (!is.null(tmp_env$res$loop_annotation)) {
    la <- tmp_env$res$loop_annotation
    anchor_types <- table(c(la$anchor1_type, la$anchor2_type))
    tss_anchor_types[[tss_name]] <- data.frame(
      TSS = tss_name, Type = names(anchor_types),
      Count = as.integer(anchor_types), stringsAsFactors = FALSE)

    loop_tbl <- table(la$loop_type)
    tss_loop_types[[tss_name]] <- data.frame(
      TSS = tss_name, Type = names(loop_tbl),
      Count = as.integer(loop_tbl), stringsAsFactors = FALSE)


  # ── Hop0/hop1 split: anchor & loop type reclassification ──
  # Read hop1 data from tmp_env$res2 (parallel to existing res extraction)
  if (!is.null(tmp_env$res2$loop_annotation)) {
    la2 <- tmp_env$res2$loop_annotation
    a2_types <- table(c(la2$anchor1_type, la2$anchor2_type))
    tss_anchor_types_hop1[[tss_name]] <- data.frame(
      TSS = tss_name, Type = names(a2_types), Count = as.integer(a2_types),
      Hop = "hop1", stringsAsFactors = FALSE)
    l2_tbl <- table(la2$loop_type)
    tss_loop_types_hop1[[tss_name]] <- data.frame(
      TSS = tss_name, Type = names(l2_tbl), Count = as.integer(l2_tbl),
      Hop = "hop1", stringsAsFactors = FALSE)
  }

    message(sprintf("      Anchors: %s | Loops: %s",
      paste(sprintf("%s=%d", names(anchor_types), anchor_types), collapse=" "),
      paste(sprintf("%s=%d", names(loop_tbl), loop_tbl), collapse=" ")))
  }

  rm(tmp_env); gc()
}

if (length(tss_available) >= 2) {
  wb_tss_cross <- createWorkbook()

  # Extend palettes for newly observed looplook categories rather than silently
  # dropping a complete cross-TSS visualization. Original labels are retained.
  extend_palette <- function(observed, palette, label) {
    observed <- unique(as.character(observed))
    observed <- observed[!is.na(observed) & nzchar(observed)]
    unknown <- setdiff(observed, names(palette))
    if (length(unknown) > 0L) {
      warning(label, " contains newly observed categories: ",
        paste(unknown, collapse = ", "),
        ". Generic colors were appended; inspect CrossTSS_Category_QC.csv.")
      extra <- stats::setNames(scales::hue_pal()(length(unknown)), unknown)
      palette <- c(palette, extra)
    }
    list(palette = palette, unknown = unknown)
  }

  # ── Anchor type reclassification across TSS windows ──
  if (exists("tss_anchor_types") && length(tss_anchor_types) >= 2) {
    tss_at_df <- bind_rows(tss_anchor_types)
    tss_at_df$TSS <- factor(tss_at_df$TSS, levels = tss_available)

    anchor_colors <- c("P" = "#4169E1", "E" = "#FFA500", "G" = "#9370DB")
    anchor_palette_info <- extend_palette(tss_at_df$Type, anchor_colors, "Anchor type")
    anchor_colors <- anchor_palette_info$palette

    {

    p_at_stack <- ggplot(tss_at_df, aes(x = TSS, y = Count, fill = Type)) +
      geom_col(position = "fill", alpha = 0.85) +
      scale_fill_manual(values = anchor_colors) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "Anchor Classification Shifts with TSS Window",
           subtitle = "Wider TSS window -> more anchors classified as promoter-like (P, eP)",
           x = "TSS Window", y = "Proportion of Anchors", fill = "Anchor Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_AnchorTypes.pdf"), p_at_stack, width = 7, height = 5)

    p_at_count <- ggplot(tss_at_df, aes(x = TSS, y = Count, fill = Type)) +
      geom_col(position = "stack", alpha = 0.85) +
      scale_fill_manual(values = anchor_colors) +
      labs(title = "Anchor Type Counts by TSS Window", x = "TSS Window",
           y = "Number of Anchors", fill = "Anchor Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_AnchorCounts.pdf"), p_at_count, width = 7, height = 5)

    addWorksheet(wb_tss_cross, "Anchor_Types_by_TSS")
    writeData(wb_tss_cross, "Anchor_Types_by_TSS", as.data.frame(tss_at_df))
    } # end anchor type plotting block
  }

  # ── Loop type reclassification across TSS windows ──
  if (exists("tss_loop_types") && length(tss_loop_types) >= 2) {
    tss_lt_df <- bind_rows(tss_loop_types)
    tss_lt_df$TSS <- factor(tss_lt_df$TSS, levels = tss_available)

    loop_colors <- c("E-P" = "#E64B35", "P-P" = "#4169E1", "E-E" = "#FFA500",
                     "E-G" = "#20B2AA", "G-P" = "#87CEEB", "G-G" = "#9370DB")
    loop_palette_info <- extend_palette(tss_lt_df$Type, loop_colors, "Loop type")
    loop_colors <- loop_palette_info$palette

    {

    p_lt_stack <- ggplot(tss_lt_df, aes(x = TSS, y = Count, fill = Type)) +
      geom_col(position = "fill", alpha = 0.85) +
      scale_fill_manual(values = loop_colors) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "Loop Type Proportions Shift with TSS Window",
           subtitle = "Wider TSS -> more PP (promoter-promoter), fewer EP (enhancer-promoter)",
           x = "TSS Window", y = "Proportion of Loops", fill = "Loop Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_LoopTypes.pdf"), p_lt_stack, width = 7, height = 5)

    p_lt_count <- ggplot(tss_lt_df, aes(x = TSS, y = Count, fill = Type)) +
      geom_col(position = "stack", alpha = 0.85) +
      scale_fill_manual(values = loop_colors) +
      labs(title = "Loop Type Counts by TSS Window", x = "TSS Window",
           y = "Number of Loops", fill = "Loop Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_LoopCounts.pdf"), p_lt_count, width = 7, height = 5)

    } # end loop type plotting block

    
  # ── Hop0/hop1 split: anchor & loop type plots ──
  if (exists("tss_anchor_types_hop1") && length(tss_anchor_types_hop1) > 0) {
    # Mark hop0 data
    for (tn in names(tss_anchor_types)) tss_anchor_types[[tn]]$Hop <- "hop0"
    for (tn in names(tss_loop_types))   tss_loop_types[[tn]]$Hop <- "hop0"
    # Merge hop0 + hop1
    tss_at_all <- bind_rows(bind_rows(tss_anchor_types), bind_rows(tss_anchor_types_hop1))
    tss_at_all$TSS <- factor(tss_at_all$TSS, levels = tss_available)
    tss_lt_all <- bind_rows(bind_rows(tss_loop_types), bind_rows(tss_loop_types_hop1))
    tss_lt_all$TSS <- factor(tss_lt_all$TSS, levels = tss_available)

    p_at_hop <- ggplot(tss_at_all, aes(x = TSS, y = Count, fill = Type)) +
      geom_col(position = "fill", alpha = 0.85) +
      scale_fill_manual(values = anchor_colors) +
      scale_y_continuous(labels = scales::percent) +
      facet_wrap(~ Hop, ncol = 2) +
      labs(title = "Anchor Classification: hop0 vs hop1 Across TSS Windows",
           subtitle = "Wider TSS -> more P, fewer E", x = "TSS Window",
           y = "Proportion of Anchors", fill = "Anchor Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            strip.text = element_text(face = "bold", size = 10))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_AnchorTypes_hop.pdf"), p_at_hop, width = 10, height = 5)

    p_lt_hop <- ggplot(tss_lt_all, aes(x = TSS, y = Count, fill = Type)) +
      geom_col(position = "fill", alpha = 0.85) +
      scale_fill_manual(values = loop_colors) +
      scale_y_continuous(labels = scales::percent) +
      facet_wrap(~ Hop, ncol = 2) +
      labs(title = "Loop Type Proportions: hop0 vs hop1 Across TSS Windows",
           subtitle = "Wider TSS -> more PP, fewer EP", x = "TSS Window",
           y = "Proportion of Loops", fill = "Loop Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            strip.text = element_text(face = "bold", size = 10))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_LoopTypes_hop.pdf"), p_lt_hop, width = 10, height = 5)
  }

    # EP/(EP+PP) ratio shift -- plot + table
    ep_ratio <- tss_lt_df %>%
      group_by(TSS) %>%
      summarise(EP_n = sum(Count[grepl("E-P", Type)]),
                PP_n = sum(Count[grepl("P-P", Type)]),
                Total = EP_n + PP_n,
                EP_ratio = ifelse(Total > 0, EP_n / Total, NA_real_),
                PP_ratio = ifelse(Total > 0, PP_n / Total, NA_real_),
                .groups = "drop")

    # Filter rows with valid ratios (Total > 0)
    ep_ratio <- ep_ratio %>% filter(is.finite(EP_ratio), is.finite(PP_ratio))

    if (nrow(ep_ratio) > 1) {
      ep_long <- ep_ratio %>%
        tidyr::pivot_longer(c(EP_ratio, PP_ratio), names_to = "Type", values_to = "Ratio") %>%
        mutate(Type = dplyr::recode(Type, EP_ratio = "EP", PP_ratio = "PP"))

      # Dynamic subtitle based on actual data
      first_tss <- as.character(ep_ratio$TSS[1])
      last_tss  <- as.character(ep_ratio$TSS[nrow(ep_ratio)])
      first_ep  <- ep_ratio$EP_ratio[1]
      last_ep   <- ep_ratio$EP_ratio[nrow(ep_ratio)]
      trend_word <- dplyr::case_when(
        last_ep < first_ep ~ "decreased",
        last_ep > first_ep ~ "increased",
        TRUE ~ "was unchanged")
      subtitle_text <- sprintf("EP ratio %s from %.1f%% (%s) to %.1f%% (%s)",
        trend_word, 100 * first_ep, first_tss, 100 * last_ep, last_tss)

      p_ep_ratio <- ggplot(ep_long, aes(x = TSS, y = Ratio, fill = Type)) +
      geom_col(position = "stack", alpha = 0.85) +
      geom_text(data = ep_ratio, inherit.aes = FALSE,
                aes(x = TSS, y = EP_ratio / 2,
                    label = sprintf("EP\n%.1f%%", 100 * EP_ratio)),
                size = 3.5, fontface = "bold", color = "white") +
      geom_text(data = ep_ratio, inherit.aes = FALSE,
                aes(x = TSS, y = EP_ratio + PP_ratio / 2,
                    label = sprintf("PP\n%.1f%%", 100 * PP_ratio)),
                size = 3.5, fontface = "bold", color = "white") +
      scale_fill_manual(values = c("EP" = "#E64B35", "PP" = "#4169E1")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(title = "EP vs PP Loop Ratio Across TSS Windows",
           subtitle = subtitle_text,
           x = "TSS Window", y = "Proportion of EP+PP Loops", fill = "Loop Type") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
    queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_EP_PP_Ratio.pdf"), p_ep_ratio, width = 6, height = 5)
    } else {
      warning("EP_PP_Ratio skipped: insufficient data")
    }

    message("    EP/(EP+PP) ratio across TSS windows:")
    for (k in seq_len(nrow(ep_ratio))) {
      message(sprintf("      %s: %.3f  (EP=%d, PP=%d)",
        ep_ratio$TSS[k], ep_ratio$EP_ratio[k], ep_ratio$EP_n[k], ep_ratio$PP_n[k]))
    }

    addWorksheet(wb_tss_cross, "Loop_Types_by_TSS")
    writeData(wb_tss_cross, "Loop_Types_by_TSS", as.data.frame(tss_lt_df))
    addWorksheet(wb_tss_cross, "EP_PP_Ratio")
    writeData(wb_tss_cross, "EP_PP_Ratio", as.data.frame(ep_ratio))
  }

  # ── Gene count comparison across TSS windows ──
  tss_gene_counts <- bind_rows(lapply(tss_available, function(tn) {
    mg <- tss_mode_genes[[tn]]
    data.frame(TSS = tn, Mode = names(mg),
               N_Genes = sapply(mg, length), stringsAsFactors = FALSE)
  }))
  tss_gene_counts$TSS <- factor(tss_gene_counts$TSS, levels = tss_available)
  tss_gene_counts$Mode <- factor(tss_gene_counts$Mode,
    levels = names(tss_mode_genes[[tss_available[1]]]))

  tss_gc_plot <- tss_gene_counts %>%
    mutate(Pipeline = case_when(
      grepl("^anno_", Mode) ~ "Basic",
      grepl("^refined_", Mode) ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined",
      grepl("^chrom_", Mode) ~ "I-refined", TRUE ~ "Other"),
      GeneMap = ifelse(grepl("_promoter_", Mode), "promoter", "all"),
      Fill    = ifelse(grepl("_T(_hop1)?$", Mode), "filled", "strict"),
      MapFill = paste(GeneMap, Fill, sep = ", "))

  p_tss_gc <- ggplot(tss_gc_plot, aes(x = TSS, y = N_Genes,
      group = Mode, color = Pipeline, shape = MapFill)) +
    geom_line(linewidth = 0.5, alpha = 0.5) +
    geom_point(size = 2.0, alpha = 0.8) +
    scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                   "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
    scale_shape_manual(values = c("all, strict" = 16, "all, filled" = 17,
                                   "promoter, strict" = 15, "promoter, filled" = 18)) +
    facet_wrap(~ Pipeline, ncol = 3, scales = "free_y") +
    labs(title = "Gene Set Size Stability Across TSS Windows",
         subtitle = sprintf("%d modes x %d TSS windows | Color=Pipeline, Shape=map+fill",
                            n_distinct(tss_gc_plot$Mode), length(tss_available)),
         x = "TSS Window", y = "Number of Target Genes") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "none", strip.text = element_text(face = "bold", size = 9),
          axis.text.x = element_text(angle = 30, hjust = 1))
  queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_GeneCounts.pdf"), p_tss_gc, width = 10, height = 6)

  # ── Hop0/hop1 split: Cross-TSS gene counts ──
  tss_gc_plot <- tss_gc_plot %>%
    mutate(HopTSS = factor(ifelse(grepl("_hop1$", Mode), "hop1 (primary+expanded)", "hop0 (primary)"),
                           levels = c("hop0 (primary)", "hop1 (primary+expanded)")))

  p_tss_gc_hop <- ggplot(tss_gc_plot, aes(x = TSS, y = N_Genes,
      group = Mode, color = Pipeline, shape = MapFill)) +
    geom_line(linewidth = 0.4, alpha = 0.45) +
    geom_point(size = 1.8, alpha = 0.7) +
    scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                   "I-refined" = "#E64B35", "C-refined" = "#FFA500", "Other" = "grey50")) +
    scale_shape_manual(values = c("all, strict" = 16, "all, filled" = 17,
                                   "promoter, strict" = 15, "promoter, filled" = 18)) +
    facet_wrap(~ HopTSS, ncol = 2, scales = "free_y") +
    labs(title = "Gene Set Size Stability Across TSS Windows — hop0 vs hop1",
         subtitle = sprintf("%d modes x %d TSS windows | Color=Pipeline, Shape=map+fill",
                            n_distinct(tss_gc_plot$Mode), length(tss_available)),
         x = "TSS Window", y = "Number of Target Genes") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "bottom", strip.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1))
  queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_GeneCounts_hop.pdf"), p_tss_gc_hop, width = 10, height = 5.5)


  # ── Jaccard overlap between adjacent TSS windows ──
  tss_jaccard <- list()
  for (i in 1:(length(tss_available)-1)) {
    tn1 <- tss_available[i]; tn2 <- tss_available[i+1]
    mg1 <- tss_mode_genes[[tn1]]; mg2 <- tss_mode_genes[[tn2]]
    common_modes <- intersect(names(mg1), names(mg2))
    for (md in common_modes) {
      g1 <- mg1[[md]]; g2 <- mg2[[md]]
      jac <- length(intersect(g1, g2)) / max(1, length(union(g1, g2)))
      tss_jaccard[[length(tss_jaccard) + 1]] <- data.frame(
        Pair = paste(tn1, tn2, sep = " vs "), Mode = md, Jaccard = jac,
        stringsAsFactors = FALSE)
    }
  }
  tss_jaccard_df <- bind_rows(tss_jaccard)

  p_tss_jac <- ggplot(tss_jaccard_df, aes(x = Pair, y = Jaccard, fill = Pair)) +
    geom_violin(alpha = 0.6, draw_quantiles = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.8) +
    labs(title = "Gene Set Overlap Between Adjacent TSS Windows",
         subtitle = "Jaccard index per mode | Higher = more stable to TSS change",
         x = NULL, y = "Jaccard Index") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.position = "none")
  queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_Jaccard.pdf"), p_tss_jac, width = 7, height = 5)

  # ── Per-mode CV of gene count across TSS ──
  tss_cv <- tss_gene_counts %>%
    group_by(Mode) %>%
    summarise(CV = sd(N_Genes) / mean(N_Genes),
              Min_N = min(N_Genes), Max_N = max(N_Genes),
              Range = Max_N - Min_N, .groups = "drop") %>%
    arrange(CV)

  p_tss_cv <- ggplot(tss_cv, aes(x = reorder(Mode, CV), y = CV)) +
    geom_col(aes(fill = CV < 0.1), alpha = 0.8) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35"), guide = "none") +
    labs(title = "Gene Count Stability Across TSS Windows",
         subtitle = sprintf("%d/%d modes with CV < 0.1 (stable) | Dashed = 10%% threshold",
                            sum(tss_cv$CV < 0.1), nrow(tss_cv)),
         x = NULL, y = "Coefficient of Variation") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
  queueable_ggsave(file.path(tss_cross_dir, "CrossTSS_CV_Stability.pdf"), p_tss_cv, width = 14, height = 5)

  # Export gene/loop/anchor data
  addWorksheet(wb_tss_cross, "Gene_Counts_by_TSS")
  writeData(wb_tss_cross, "Gene_Counts_by_TSS", as.data.frame(tss_gene_counts))
  addWorksheet(wb_tss_cross, "Jaccard_Overlap")
  writeData(wb_tss_cross, "Jaccard_Overlap", as.data.frame(tss_jaccard_df))
  addWorksheet(wb_tss_cross, "Stability_CV")
  writeData(wb_tss_cross, "Stability_CV", as.data.frame(tss_cv))
  out_tss_cross <- file.path(tss_cross_dir, "CrossTSS_Summary.xlsx")
  saveWorkbook(wb_tss_cross, out_tss_cross, overwrite = TRUE)

  message(sprintf("    Modes with CV < 0.1: %d/%d", sum(tss_cv$CV < 0.1), nrow(tss_cv)))
  message("    Cross-TSS summary: ", tss_cross_dir)
} else {
  message("    WARNING: Need >= 2 TSS RData files, found ", length(tss_available))
  message("    Run server scripts at \u670d\u52a1\u5668\u811a\u672c/ first to generate per-TSS RData files.")
}

# ══════════════════════════════════════════════════════════════════════════════
# Multi-Dimensional Mode Comparison: rank all modes by composite score
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# Evidence Composition Analysis: 3D loop vs linear fallback
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Evidence Composition Analysis ====")
ev_dir <- file.path(base_out_dir, "evidence_composition")
dir.create(ev_dir, recursive = TRUE, showWarnings = FALSE)



# Build gene–evidence assignment map once for all uses
all_gene_ev <- bind_rows(lapply(modes, function(md) {
  if (md$near) return(NULL)
  build_mode_gene_evidence(md, md$name)
}))

# ── Evidence Composition (using assignment-level counts) ──
ev_all_df <- all_gene_ev %>%
  count(Mode, Evidence, Source, name = "Count") %>%
  mutate(Evidence = factor(Evidence, levels = names(evidence_palette)),
    Pipeline = factor(case_when(grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
      TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined")),
    Fill = ifelse(grepl("_T(_hop1)?$", Mode), "T (filled)", "F (strict)"),
    Hop  = ifelse(grepl("_hop1$", Mode), "hop1", "hop0"),
    GeneMap = ifelse(grepl("_promoter_", Mode), "promoter", "all"))
if (nrow(ev_all_df) == 0) {
  warning("No evidence composition data found.")
} else {
  # ── Export workbook (created early for all sheets) ──
  wb_ev <- createWorkbook()

  # ── Figure 1: Full mode stacked bar ──
  ev_mode_summary <- ev_all_df %>%
    group_by(Mode, Pipeline, Fill, Hop, GeneMap) %>%
    mutate(Pct = Count / sum(Count) * 100) %>% ungroup()

  p_ev_full <- ggplot(ev_mode_summary,
    aes(x = reorder(Mode, -Count * as.numeric(Evidence == "local_promoter_overlap" |
                                               Evidence == "distal_promoter")),
        y = Pct, fill = Evidence)) +
    geom_col(width = 0.75, alpha = 0.9) +
    scale_fill_manual(values = evidence_palette) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ Fill, ncol = 1, scales = "free_x") +
    labs(title = "Evidence Composition: Strict (F) vs Filled (T) — All Modes",
         subtitle = sprintf("%d modes | Blue = 3D loop evidence, Red = linear fallback",
                            n_distinct(ev_mode_summary$Mode)),
         x = NULL, y = "Percentage of Peak-Gene Evidence Assignments") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5))
  queueable_ggsave(file.path(ev_dir, "Evidence_Composition_AllModes.pdf"), p_ev_full, width = 16, height = 9)

  # ── Figure 2: By Pipeline summary ──
  ev_pipe_summary <- ev_all_df %>%
    group_by(Pipeline, Evidence) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(Pipeline) %>%
    mutate(Pct = Count / sum(Count) * 100) %>% ungroup()

  p_ev_pipe <- ggplot(ev_pipe_summary, aes(x = Pipeline, y = Pct, fill = Evidence)) +
    geom_col(width = 0.6, alpha = 0.9) +
    scale_fill_manual(values = evidence_palette) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "Evidence Composition by Pipeline",
         subtitle = "3D loop evidence (blue) vs linear fallback (red)",
         x = NULL, y = "Percentage of Evidence Assignments") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(face = "bold", size = 10),
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5))
  queueable_ggsave(file.path(ev_dir, "Evidence_Composition_byPipeline.pdf"), p_ev_pipe, width = 7, height = 5)

  # ── Figure 3: F vs T comparison (paired, promoter map only) ──
  ev_ft <- ev_mode_summary %>%
    filter(GeneMap == "promoter", Hop == "hop0") %>%
    mutate(BaseMode = Mode %>% str_replace("^anno_","") %>%
             str_replace("^refined_","") %>%
             str_replace("^chrom_only_","") %>%
             str_replace("^chrom_","") %>%
             str_replace("_promoter_F$", "") %>%
             str_replace("_promoter_T$", ""))
  if (nrow(ev_ft) > 0) {
    ev_ft_3d <- ev_ft %>%
      # Source-based definition is robust to evolving looplook evidence labels.
      filter(Source == "loop_anchor") %>%
      group_by(Mode, Pipeline, Fill, BaseMode) %>%
      summarise(Pct_3D = sum(Pct), .groups = "drop")

    p_ev_ft <- ggplot(ev_ft_3d, aes(x = BaseMode, y = Pct_3D, fill = Fill)) +
      geom_col(position = "dodge", width = 0.6, alpha = 0.9) +
      scale_fill_manual(values = c("F (strict)" = "#2166AC", "T (filled)" = "#B2182B")) +
      facet_wrap(~ Pipeline, ncol = 4) +
      labs(title = "3D Loop Evidence: Strict (F) vs Filled (T)",
           subtitle = sprintf("Promoter gene map, hop0 | Higher blue = more 3D-dependent"),
           x = NULL, y = "3D Loop Evidence (%)") +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(face = "bold", size = 10),
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5))
    queueable_ggsave(file.path(ev_dir, "Evidence_Composition_Strict_vs_Filled.pdf"), p_ev_ft, width = 10, height = 5)
  }

  # ── Add remaining sheets ──
  addWorksheet(wb_ev, "Evidence_Per_Mode")
  writeData(wb_ev, "Evidence_Per_Mode", as.data.frame(ev_mode_summary))
  addWorksheet(wb_ev, "Evidence_by_Pipeline")
  writeData(wb_ev, "Evidence_by_Pipeline", as.data.frame(ev_pipe_summary))
  out_ev <- file.path(ev_dir, "Evidence_Composition.xlsx")

  # ── Per-evidence effect size (shared gene–evidence mapping) ──
  if (exists("gene_stat") && nrow(gene_stat) > 0) {
    message("    Computing per-evidence effect sizes...")

    get_strict_bg <- function(mode_name) {
      intersect(setdiff(all_genes, union(unique(toupper(mode_genes[[mode_name]])), chipseeker_genes)), gene_stat$gene)
    }

    run_evidence_es <- function(gene_ev, signal_col) {
      bg <- get_strict_bg(gene_ev$Mode[1])
      res <- list()
      for (ev in unique(as.character(gene_ev$Evidence))) {
        ev_genes <- unique(gene_ev$Gene[gene_ev$Evidence == ev])
        if (length(ev_genes) >= 5 && length(bg) >= 5) {
          x <- gene_stat[[signal_col]][gene_stat$gene %in% ev_genes]
          y <- gene_stat[[signal_col]][gene_stat$gene %in% bg]
          rb <- rank_biserial(x, y)
          res[[length(res) + 1]] <- data.frame(Mode = gene_ev$Mode[1], Evidence = ev,
            N_genes = length(ev_genes), Med_signal = median(x, na.rm = TRUE),
            ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi, stringsAsFactors = FALSE)
        }
      }
      if (length(res) == 0) return(data.frame(Mode=character(),Evidence=character(),N_genes=integer(),Med_signal=numeric(),ES=numeric(),ES_lo=numeric(),ES_hi=numeric(),stringsAsFactors=FALSE))
      bind_rows(res)
    }

    # Build gene–evidence map once for all non-ctrl modes
    all_gene_ev_list <- list()
  for (i in seq_along(modes)) {
    md <- modes[[i]]; if (md$near) next; nm <- md$name
    if (i %% 8 == 0L) message(sprintf("      per-evidence: %d/32 modes", i))
      gev <- build_mode_gene_evidence(md, nm, gene_stat$gene)
      if (nrow(gev) > 0) all_gene_ev_list[[nm]] <- gev
    }

    # raw-LFC (skipped; bottleneck, will plot separately)
    ev_es_list <- list()
    # signal_pv (skipped)
    ev_es_pv_list <- list()
    empty_evidence_es <- function() { data.frame(Mode=character(),Evidence=character(),N_genes=integer(),Med_signal=numeric(),ES=numeric(),ES_lo=numeric(),ES_hi=numeric(),stringsAsFactors=FALSE) }
    ev_es_df <- empty_evidence_es()
    ev_es_pv_df <- empty_evidence_es()

    # Raw-LFC plot
    if (nrow(ev_es_df) > 0) {
      ev_es_df <- ev_es_df %>%
        mutate(Evidence = factor(Evidence, levels = names(evidence_palette)),
          Pipeline = factor(case_when(
            grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
            grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
            TRUE ~ "Other"), levels = c("Basic", "E-refined", "I-refined", "C-refined"))) %>%
        mutate(Fill = factor(ifelse(grepl("_T(_hop1)?$", Mode), "T (filled)", "F (strict)"),
                             levels = c("F (strict)", "T (filled)")))

      p_ev_es <- ggplot(ev_es_df, aes(x = ES, y = Evidence, fill = Evidence)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_boxplot(outlier.shape = 16, outlier.size = 0.5, alpha = 0.7, width = 0.5) +
        scale_fill_manual(values = evidence_palette) +
        facet_grid(rows = vars(Fill), cols = vars(Pipeline)) +
        labs(title = "Per-Evidence Effect Size vs Non-assigned control",
             subtitle = sprintf("signal_lfc | %d modes | Higher = stronger downregulation under that evidence",
                                n_distinct(ev_es_df$Mode)),
             x = "Rank-Biserial Effect Size", y = NULL) +
        theme_minimal() +
        theme(legend.position = "none", panel.grid.major.y = element_blank(),
              strip.text = element_text(face = "bold", size = 9),
              plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
              plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5))
      queueable_ggsave(file.path(ev_dir, "Evidence_EffectSize_byPipeline.pdf"), p_ev_es, width = 10, height = 5.5)

      addWorksheet(wb_ev, "Evidence_EffectSize")
      writeData(wb_ev, "Evidence_EffectSize", as.data.frame(ev_es_df))
    }
    # signal_pv (independent of raw-LFC empty condition)
    if (nrow(ev_es_pv_df) > 0) {
        ev_es_pv_df <- ev_es_pv_df %>%
          mutate(Evidence = factor(Evidence, levels = names(evidence_palette)),
            Pipeline = factor(case_when(grepl("^anno_",Mode)~"Basic",grepl("^refined_",Mode)~"E-refined",
              grepl("^chrom_only_",Mode)~"C-refined",grepl("^chrom_",Mode)~"I-refined",TRUE~"Other"),
              levels=c("Basic","E-refined","I-refined","C-refined")),
            Fill = factor(ifelse(grepl("_T(_hop1)?$", Mode), "T (filled)", "F (strict)"),
                          levels = c("F (strict)", "T (filled)")))
        p_ev_es_pv <- ggplot(ev_es_pv_df, aes(x=ES, y=Evidence, fill=Evidence)) +
          geom_vline(xintercept=0, linetype="dashed", color="grey50") +
          geom_boxplot(outlier.shape=16, outlier.size=0.5, alpha=0.7, width=0.5) +
          scale_fill_manual(values=evidence_palette) +
          facet_wrap(~ Pipeline, ncol = 4) +
          labs(title = "Per-Evidence Effect Size vs Non-assigned control (signal_pv)",
               subtitle = sprintf("signal_pv weighted | %d modes", n_distinct(ev_es_pv_df$Mode)),
               x = "Rank-Biserial Effect Size", y = NULL) +
          theme_minimal() + theme(legend.position = "none", panel.grid.major.y = element_blank(),
            strip.text = element_text(face = "bold", size = 9),
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5))
        queueable_ggsave(file.path(ev_dir, "Evidence_EffectSize_byPipeline_pv.pdf"), p_ev_es_pv, width = 10, height = 5.5)
        addWorksheet(wb_ev, "Evidence_EffectSize_pv")
        writeData(wb_ev, "Evidence_EffectSize_pv", as.data.frame(ev_es_pv_df))
  }
  }
  saveWorkbook(wb_ev, out_ev, overwrite = TRUE)
  }

  message("    Evidence composition results: ", ev_dir)



# ══════════════════════════════════════════════════════════════════════════════
# Peak Annotation Category Effect Size: Promoter / Genic / Distal
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Peak Category Effect Size Analysis ====")
peak_cat_dir <- file.path(base_out_dir, "peak_category_es")
dir.create(peak_cat_dir, recursive = TRUE, showWarnings = FALSE)

# Build coordinate keys for safe matching (do not rely on row order)
peak_key_chip <- with(peak_anno_df, paste(seqnames, start, end, sep = ":"))

# ── Peak Category ES with cache ──
peak_cat_cache <- file.path(base_out_dir, ".peak_category_es_cache.rds")
if (file.exists(peak_cat_cache) && !force_recompute) {
  peak_cat_df <- readRDS(peak_cat_cache)
  message("    Peak category ES loaded from cache: ", nrow(peak_cat_df), " rows")
} else {
# Peak Category ES skipped (bottleneck, 1152 bootstrap calls)
  peak_cat_df <- data.frame()
  peak_cat_es <- list()
  message("    Peak Category ES skipped (will plot separately)")
  if (FALSE) {
peak_cat_list <- list(
  "Promoter" = which(grepl("Promoter", peak_anno_df$annotation)),
  "Genic"    = which(grepl("Exon|Intron|UTR|Downstream", peak_anno_df$annotation) &
                     !grepl("Promoter", peak_anno_df$annotation)),
  "Distal"   = which(grepl("Distal Intergenic", peak_anno_df$annotation)),
  "Non_Promoter" = which(!grepl("Promoter", peak_anno_df$annotation)))

if (use_se_stratification) {
  if (sum(is_se_peak) > 0) peak_cat_list[["Super_Enhancer"]] <- which(is_se_peak)
  if (sum(is_te_peak) > 0) peak_cat_list[["Typical_Enhancer"]] <- which(is_te_peak)
  message(sprintf("    Peak category ES with SE/TE: SE=%d TE=%d peaks",
    sum(is_se_peak), sum(is_te_peak)))
}

peak_cat_es <- list()
for (cat_name in names(peak_cat_list)) {
  cat_idx_chip <- peak_cat_list[[cat_name]]
  cat_chip_genes <- unique(toupper(na.omit(peak_anno_df$SYMBOL[cat_idx_chip])))
  cat_key_chip <- peak_key_chip[cat_idx_chip]
  for (i in seq_along(modes)) {
    md <- modes[[i]]; if (md$near) next; nm <- md$name
    bed_info <- md$ann$target_annotation; if (is.null(bed_info)) next
    peak_key_bed <- with(bed_info, paste(seqnames, start, end, sep = ":"))
    cat_idx_bed <- match(cat_key_chip, peak_key_bed)
    match_rate <- mean(!is.na(cat_idx_bed))
    if (match_rate < 0.99) warning(sprintf("[%s] %s match rate = %.1f%%", nm, cat_name, match_rate * 100))
    cat_idx_bed <- cat_idx_bed[!is.na(cat_idx_bed)]
    if (length(cat_idx_bed) == 0) next
    cat_loop_genes <- get_mode_genes(
      md,
      row_index = cat_idx_bed
    )
    cat_loop_genes <- intersect(cat_loop_genes, all_genes)
    cat_chip_g    <- intersect(cat_chip_genes, all_genes)
    only_loop <- intersect(setdiff(cat_loop_genes, cat_chip_g), gene_stat$gene)
    only_chip <- intersect(setdiff(cat_chip_g, cat_loop_genes), gene_stat$gene)
    bg_genes  <- intersect(setdiff(all_genes, union(cat_loop_genes, cat_chip_g)), gene_stat$gene)
    if (length(only_loop) >= 5 && length(bg_genes) >= 5) {
      rb <- rank_biserial(gene_stat$signal_lfc[gene_stat$gene %in% only_loop], gene_stat$signal_lfc[gene_stat$gene %in% bg_genes])
      peak_cat_es[[length(peak_cat_es)+1]] <- data.frame(Mode=nm, Category=cat_name, Comparison="Only_looplook vs BG", ES=rb$est, ES_lo=rb$lo, ES_hi=rb$hi, stringsAsFactors=FALSE)
    }
    if (length(only_chip) >= 5 && length(bg_genes) >= 5) {
      rb <- rank_biserial(gene_stat$signal_lfc[gene_stat$gene %in% only_chip], gene_stat$signal_lfc[gene_stat$gene %in% bg_genes])
      peak_cat_es[[length(peak_cat_es)+1]] <- data.frame(Mode=nm, Category=cat_name, Comparison="Only_ChIPseeker vs BG", ES=rb$est, ES_lo=rb$lo, ES_hi=rb$hi, stringsAsFactors=FALSE)
    }
    if (length(only_loop) >= 5 && length(only_chip) >= 5) {
      rb <- rank_biserial(gene_stat$signal_lfc[gene_stat$gene %in% only_loop], gene_stat$signal_lfc[gene_stat$gene %in% only_chip])
      peak_cat_es[[length(peak_cat_es)+1]] <- data.frame(Mode=nm, Category=cat_name, Comparison="Only_looplook vs Only_ChIPseeker", ES=rb$est, ES_lo=rb$lo, ES_hi=rb$hi, stringsAsFactors=FALSE)
    }
  }
  message(sprintf("    %s: %d ChIPseeker peaks; %d matched looplook peaks",
    cat_name, length(cat_idx_chip), length(cat_idx_bed)))
}
if (length(peak_cat_es) > 0) {
  cat_levels <- c("Promoter","Genic","Distal","Non_Promoter")
  if (use_se_stratification) cat_levels <- c(cat_levels, "Super_Enhancer", "Typical_Enhancer")
  peak_cat_df <- bind_rows(peak_cat_es) %>%
    mutate(Category = factor(Category, levels = cat_levels),
           Comparison = factor(Comparison, levels = c("Only_looplook vs BG", "Only_ChIPseeker vs BG", "Only_looplook vs Only_ChIPseeker")),
           Pipeline = factor(case_when(grepl("^anno_",Mode)~"Basic", grepl("^refined_",Mode)~"E-refined", grepl("^chrom_only_",Mode)~"C-refined", grepl("^chrom_",Mode)~"I-refined", TRUE~"Other"), levels = c("Basic","E-refined","I-refined","C-refined")))
  ncol_cat <- if (use_se_stratification) 3 else 3
  se_tag <- if (use_se_stratification) " + SE/TE" else ""
  se_fname <- if (use_se_stratification) "_with_SE_TE" else ""
  p_pcat <- ggplot(peak_cat_df, aes(x = Pipeline, y = ES, fill = Comparison)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, alpha = 0.65, width = 0.6, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("Only_looplook vs BG" = "#E64B35", "Only_ChIPseeker vs BG" = "#4DBBD5", "Only_looplook vs Only_ChIPseeker" = "#FFA500")) +
    facet_wrap(~ Category, ncol = ncol_cat) +
    labs(title = paste0("Effect Size by Peak Annotation Category", se_tag),
         subtitle = sprintf("%d modes | same-input peaks, looplook vs ChIPseeker assignments", n_distinct(peak_cat_df$Mode)),
         x = NULL, y = "Effect Size (rank-biserial)") +
    theme_classic() + theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10), plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
  queueable_ggsave(file.path(peak_cat_dir, paste0("PeakCategory_EffectSize", se_fname, ".pdf")), p_pcat, width = if (use_se_stratification) 14 else 12, height = 5)
  wb_pcat <- createWorkbook(); addWorksheet(wb_pcat, "PeakCategory_ES"); writeData(wb_pcat, "PeakCategory_ES", as.data.frame(peak_cat_df))
  saveWorkbook(wb_pcat, file.path(peak_cat_dir, "PeakCategory_EffectSize.xlsx"), overwrite = TRUE)
  message("    Peak category ES", se_tag, ": ", peak_cat_dir)

  # ── SE/TE专项分析 ──
  if (use_se_stratification) {
    se_peak_ids <- make_peak_id(peak_anno_df)[is_se_peak]
    te_peak_ids <- make_peak_id(peak_anno_df)[is_te_peak]

    # 1. SE-gene distance vs effect size
    se_dist_es <- list()
    for (nm in names(mode_genes)) {
      genes <- mode_genes[[nm]]
      if (length(genes) == 0) next
      md <- modes[[which(sapply(modes, function(m) m$name == nm))[1]]]
      loop_genes <- intersect(unique(toupper(genes)), gene_stat$gene)
      chip_g <- intersect(get_chip_reference(md), gene_stat$gene)
      only_loop <- intersect(setdiff(loop_genes, chip_g), gene_stat$gene)
      if (length(only_loop) < 5) next

      se_dist <- dist_df_clean %>% filter(Mode == nm, Peak_ID %in% se_peak_ids, Source == "looplook")
      if (nrow(se_dist) < 3) next

      # SE-specific effect size: only genes from SE peaks
      se_loop_genes <- intersect(unique(se_dist$Gene), only_loop)
      se_chip_genes <- unique(toupper(na.omit(peak_anno_df$SYMBOL[is_se_peak])))
      se_chip_genes <- intersect(se_chip_genes, gene_stat$gene)
      se_only_loop <- setdiff(se_loop_genes, se_chip_genes)
      if (length(se_only_loop) < 3 || length(se_chip_genes) < 3) next

      rb <- rank_biserial(
        gene_stat$signal_lfc[gene_stat$gene %in% se_only_loop],
        gene_stat$signal_lfc[gene_stat$gene %in% se_chip_genes])
      se_dist_es[[nm]] <- data.frame(
        Mode = nm,
        Med_SE_Distance_kb = median(se_dist$Dist_kb, na.rm = TRUE),
        N_SE_pairs = nrow(se_dist),
        ES = rb$est, ES_lo = rb$lo, ES_hi = rb$hi,
        Pipeline = case_when(
          grepl("^anno_", nm) ~ "Basic", grepl("^refined_", nm) ~ "E-refined",
          grepl("^chrom_only_", nm) ~ "C-refined", grepl("^chrom_", nm) ~ "I-refined",
          TRUE ~ "Other"),
        stringsAsFactors = FALSE)
    }
    if (length(se_dist_es) > 0) {
      se_dist_es_df <- bind_rows(se_dist_es)
      p_se_dist <- ggplot(se_dist_es_df, aes(x = Med_SE_Distance_kb, y = ES, color = Pipeline)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_smooth(method = "lm", se = TRUE, color = "grey40", fill = "grey80", linewidth = 0.6) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                       "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
        labs(title = "SE-assigned Target Genes: Distance vs Effect Size",
             subtitle = sprintf("Each point = one mode (%d modes) | SE peaks only | Smoother = linear fit",
                                nrow(se_dist_es_df)),
             x = "Median SE peak to gene distance (kb)", y = "Effect Size (looplook vs ChIPseeker)") +
        theme_classic() +
        theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
              plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
              legend.position = "bottom")
      queueable_ggsave(file.path(peak_cat_dir, "SE_Distance_vs_EffectSize.pdf"), p_se_dist, width = 7, height = 5)

      wb_pcat_se <- tryCatch(loadWorkbook(file.path(peak_cat_dir, "PeakCategory_EffectSize.xlsx")), error = function(e) createWorkbook())
      addWorksheet(wb_pcat_se, "SE_Distance_vs_ES"); writeData(wb_pcat_se, "SE_Distance_vs_ES", as.data.frame(se_dist_es_df))
      saveWorkbook(wb_pcat_se, file.path(peak_cat_dir, "PeakCategory_EffectSize.xlsx"), overwrite = TRUE)
      message("    SE distance-ES: ", peak_cat_dir, "/SE_Distance_vs_EffectSize.pdf")
    }

    # 2. SE vs TE head-to-head effect size comparison
    se_te_es <- peak_cat_df %>%
      filter(Category %in% c("Super_Enhancer", "Typical_Enhancer"),
             Comparison == "Only_looplook vs Only_ChIPseeker")
    if (nrow(se_te_es) > 0 &&
        all(c("Super_Enhancer", "Typical_Enhancer") %in% as.character(se_te_es$Category))) {
      se_te_wide <- se_te_es %>%
        dplyr::select(Mode, Pipeline, Category, ES) %>%
        tidyr::pivot_wider(names_from = Category, values_from = ES) %>%
        filter(!is.na(Super_Enhancer), !is.na(Typical_Enhancer))
      if (nrow(se_te_wide) >= 5) {
        wt_sete <- wilcox.test(se_te_wide$Super_Enhancer, se_te_wide$Typical_Enhancer,
                               paired = TRUE, alternative = "two.sided", exact = FALSE)
        se_te_label <- sprintf("Paired Wilcoxon SE vs TE p = %.3f (%d modes)", wt_sete$p.value, nrow(se_te_wide))
      } else { se_te_label <- "Insufficient data" }

      p_se_te <- ggplot(se_te_es, aes(x = Category, y = ES, fill = Category)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.7, width = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +
        scale_fill_manual(values = c("Super_Enhancer" = "#E64B35", "Typical_Enhancer" = "#FFA500")) +
        facet_wrap(~ Pipeline, ncol = 4) +
        labs(title = "SE vs TE: Looplook-only Assignments vs Only_ChIPseeker",
             subtitle = se_te_label,
             x = NULL, y = "Effect Size (rank-biserial)") +
        theme_classic() +
        theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
              plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
              legend.position = "none", strip.text = element_text(face = "bold", size = 10))
      queueable_ggsave(file.path(peak_cat_dir, "SE_vs_TE_EffectSize.pdf"), p_se_te, width = 10, height = 5)

      wb_pcat_se <- tryCatch(loadWorkbook(file.path(peak_cat_dir, "PeakCategory_EffectSize.xlsx")), error = function(e) createWorkbook())
      addWorksheet(wb_pcat_se, "SE_vs_TE_ES"); writeData(wb_pcat_se, "SE_vs_TE_ES", as.data.frame(se_te_es))
      saveWorkbook(wb_pcat_se, file.path(peak_cat_dir, "PeakCategory_EffectSize.xlsx"), overwrite = TRUE)
      message("    SE vs TE ES: ", peak_cat_dir, "/SE_vs_TE_EffectSize.pdf")
    }

  }
    # 3. Distance-effect size binned curve for ChIPseeker-distal peaks
    distal_peak_ids_for_es <- make_peak_id(peak_anno_df[is_distal, , drop = FALSE])
    dist_es_bin <- dist_df_clean %>%
      filter(
        Source == "looplook",
        GeneSet == "Only_looplook",
        Peak_ID %in% distal_peak_ids_for_es
      ) %>%
      mutate(DistBin = factor(case_when(
        Dist_kb < 1     ~ "<1 kb",
        Dist_kb < 5     ~ "1-5 kb",
        Dist_kb < 20    ~ "5-20 kb",
        Dist_kb < 100   ~ "20-100 kb",
        Dist_kb < 500   ~ "100-500 kb",
        TRUE             ~ ">500 kb"),
        levels = c("<1 kb","1-5 kb","5-20 kb","20-100 kb","100-500 kb",">500 kb")))
    if (nrow(dist_es_bin) > 0 && exists("gene_stat") && nrow(gene_stat) > 0) {
      se_dist_bin_es <- list()
      for (nm in unique(dist_es_bin$Mode)) {
        genes <- mode_genes[[nm]]
        if (length(genes) == 0) next
        md <- modes[[which(sapply(modes, function(m) m$name == nm))[1]]]
        loop_genes <- intersect(unique(toupper(genes)), gene_stat$gene)
        chip_g <- intersect(get_chip_reference(md), gene_stat$gene)
        only_loop <- intersect(setdiff(loop_genes, chip_g), gene_stat$gene)
        if (length(only_loop) < 5) next
        mode_dists <- dist_es_bin %>% filter(Mode == nm)
        for (db in levels(dist_es_bin$DistBin)) {
          db_genes <- mode_dists %>% filter(DistBin == db) %>% pull(Gene) %>% unique() %>% intersect(only_loop)
          if (length(db_genes) < 3) next
          x <- gene_stat$signal_lfc[gene_stat$gene %in% db_genes]
          y <- gene_stat$signal_lfc[gene_stat$gene %in% chip_g]
          rb <- tryCatch(rank_biserial_point(x, y), error = function(e) NA_real_)
          if (!is.na(rb)) {
            se_dist_bin_es[[length(se_dist_bin_es)+1]] <- data.frame(
              Mode = nm, DistBin = db, ES = rb, N_genes = length(db_genes),
              Pipeline = case_when(grepl("^anno_",nm)~"Basic", grepl("^refined_",nm)~"E-refined",
                grepl("^chrom_only_",nm)~"C-refined", grepl("^chrom_",nm)~"I-refined", TRUE~"Other"),
              stringsAsFactors = FALSE)
          }
        }
      }
      if (length(se_dist_bin_es) > 0) {
        se_dist_bin_df <- bind_rows(se_dist_bin_es)
        p_dist_bin <- ggplot(se_dist_bin_df, aes(x = DistBin, y = ES, fill = Pipeline)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.6) +
          scale_fill_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                        "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
          labs(title = "Looplook Assigned Genes: Distance from Peak to TSS vs Effect Size",
               subtitle = sprintf("Genes binned by peak-to-TSS distance | %d modes | Boxplot = all modes combined",
                                  n_distinct(se_dist_bin_df$Mode)),
               x = "Linear distance (peak to TSS)", y = "Effect Size (looplook-only genes vs ChIPseeker)") +
          theme_classic() +
          theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
                legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
        queueable_ggsave(file.path(peak_cat_dir, "DistanceBinned_vs_EffectSize.pdf"), p_dist_bin, width = 9, height = 5)
        message("    Distance bin-ES: ", peak_cat_dir, "/DistanceBinned_vs_EffectSize.pdf")
      }
    }
}

}  # end if (FALSE) Peak Category ES

# ══════════════════════════════════════════════════════════════════════════════
# Case Study: Looplook gene assignments from ChIPseeker-distal peaks
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Case Study: Loop-supported assignments from ChIPseeker-distal peaks ====")
cs_dir <- file.path(base_out_dir, "case_study")
dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)

# Get distal intergenic peak coords (ChIPseeker: Distal Intergenic = nearest gene is distal)
is_distal_anno <- grepl("Distal Intergenic", peak_anno_df$annotation)
distal_peak_coords <- paste0(
  peak_anno_df$seqnames[is_distal_anno], ":",
  peak_anno_df$start[is_distal_anno], "-",
  peak_anno_df$end[is_distal_anno])
message(sprintf("    ChIPseeker distal intergenic peaks: %d", length(distal_peak_coords)))

# Run case study for representative modes (4 pipelines x 2 maps x F only = 8 modes, hop0 only)
cs_modes <- c(
  "anno_all_F",       "anno_promoter_F",
  "refined_all_F",    "refined_promoter_F",
  "chrom_all_F",      "chrom_promoter_F",
  "chrom_only_all_F", "chrom_only_promoter_F")
cs_modes <- intersect(cs_modes, names(mode_genes))

case_study_all <- list()
cs_summary <- list()
for (nm in cs_modes) {
  md <- modes[[which(sapply(modes, function(m) m$name == nm))[1]]]
  bed_info <- md$ann$target_annotation
  if (is.null(bed_info)) next

  loop_peak_ids <- paste0(bed_info$seqnames, ":", bed_info$start, "-", bed_info$end)
  is_distal <- loop_peak_ids %in% distal_peak_coords
  if (sum(is_distal) == 0) next

  # Loop-supported candidate genes assigned from peaks classified as distal intergenic by ChIPseeker
  distal_genes <- get_mode_genes(
    md,
    row_index = which(is_distal)
  )
  distal_genes <- intersect(distal_genes, gene_stat$gene)

  if (length(distal_genes) == 0) next

  mode_cs <- data.frame(
    Mode = nm,
    Gene = distal_genes,
    signal_lfc = gene_stat$signal_lfc[match(distal_genes, gene_stat$gene)],
    log2FC = gene_stat$log2FoldChange[match(distal_genes, gene_stat$gene)],
    stringsAsFactors = FALSE) %>%
    filter(!is.na(signal_lfc)) %>%
    arrange(desc(signal_lfc))
  case_study_all[[nm]] <- mode_cs

  cs_summary[[nm]] <- data.frame(
    Mode = nm,
    Pipeline = case_when(grepl("^anno_",nm)~"Basic",grepl("^refined_",nm)~"E-refined",
      grepl("^chrom_only_",nm)~"C-refined",grepl("^chrom_",nm)~"I-refined",TRUE~"Other"),
    Map = ifelse(grepl("_promoter_", nm), "promoter", "all"),
    Distal_Peaks = sum(is_distal),
    Assigned_Genes = length(distal_genes),
    N_WithSignal = nrow(mode_cs),
    N_Downregulated = sum(mode_cs$signal_lfc > 0),
    Median_signal_lfc = median(mode_cs$signal_lfc),
    Mean_signal_lfc = mean(mode_cs$signal_lfc),
    stringsAsFactors = FALSE)
  message(sprintf("    %-30s distal_peaks=%-4d genes=%-4d down=%-3d median_signal_lfc=%.2f",
    nm, sum(is_distal), length(distal_genes), sum(mode_cs$signal_lfc > 0, na.rm = TRUE), median(mode_cs$signal_lfc)))
}

if (length(cs_summary) > 0) {
  cs_summary_df <- bind_rows(cs_summary) %>% arrange(desc(Assigned_Genes))

  # Save summary table
  write.csv(cs_summary_df, file.path(cs_dir, "CaseStudy_Summary.csv"), row.names = FALSE)

  # Save all genes (combined)
  cs_all_df <- bind_rows(case_study_all) %>%
    mutate(signal_lfc = round(signal_lfc, 4), log2FC = round(log2FC, 4))
  write.csv(cs_all_df, file.path(cs_dir, "CaseStudy_All_Genes.csv"), row.names = FALSE)

  # Summary bar chart
  p_cs_bar <- ggplot(cs_summary_df, aes(x = reorder(Mode, Assigned_Genes), y = Assigned_Genes, fill = Pipeline)) +
    geom_col(alpha = 0.85) + coord_flip() +
    scale_fill_manual(values = c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")) +
    labs(title = "Loop-supported Candidate Genes from ChIPseeker-Distal Peaks",
         subtitle = sprintf("%d modes | ChIPseeker nearest-gene annotation retained; looplook adds loop-supported candidates",
                            nrow(cs_summary_df)),
         x = NULL, y = "Number of looplook-assigned genes") +
    theme_classic() + theme(legend.position = "bottom")
  queueable_ggsave(file.path(cs_dir, "CaseStudy_Summary.pdf"), p_cs_bar, width = 8, height = 5)

  # Distribution of signal_lfc across modes
  n_modes_cs <- n_distinct(cs_all_df$Mode)
  cs_all_df <- cs_all_df %>% mutate(
    Pipeline = case_when(grepl("^anno_",Mode)~"Basic",grepl("^refined_",Mode)~"E-refined",
      grepl("^chrom_only_",Mode)~"C-refined",grepl("^chrom_",Mode)~"I-refined",TRUE~"Other"),
    ModeLabel = paste0(Mode, " (n=", ave(seq_along(Mode), Mode, FUN=length), ")"))

  p_cs_dist <- ggplot(cs_all_df, aes(x = signal_lfc, y = reorder(ModeLabel, signal_lfc, median), fill = Pipeline)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, alpha = 0.7, width = 0.6) +
    scale_fill_manual(values = c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")) +
    labs(title = "Signal Distribution of Loop-supported Genes from ChIPseeker-Distal Peaks",
         subtitle = sprintf("signal_lfc = -log2FoldChange | %d modes | Genes from distal-intergenic peaks only",
                            n_modes_cs),
         x = "signal_lfc (>0 = downregulated in BRD4 degradation)", y = NULL) +
    theme_classic() + theme(legend.position = "bottom")
  queueable_ggsave(file.path(cs_dir, "CaseStudy_SignalDistribution.pdf"), p_cs_dist, width = 10, height = 5)

  message(sprintf("    Case study: %d modes, %d genes, saved to %s",
    nrow(cs_summary_df), nrow(cs_all_df), cs_dir))
}
# Expression-Matched Effect Size Sensitivity Analysis
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== Expression-Matched Effect Size ====")
expr_match_dir <- file.path(base_out_dir, "expression_matched_es")
dir.create(expr_match_dir, recursive = TRUE, showWarnings = FALSE)

message("    Expression-Matched ES skipped (bottleneck), will plot separately")
  if (FALSE) {

  # Build TPM lookup from DMSO samples (reuse dmso_cols from metadata-driven identification)
expr_mat <- looplook:::read_robust_general(expr_path, header = TRUE, row_name = 1, desc = "Expr", min_cols = 2)
expr_mat <- normalise_rownames(expr_mat, "Expression matrix")
# dmso_cols already determined by metadata above; re-validate intersection
lost_controls <- setdiff(dmso_cols, colnames(expr_mat))
if (length(lost_controls) > 0L) {
  stop("DMSO/control samples from metadata missing in expression matrix: ",
    paste(lost_controls, collapse = ", "))
}
dmso_cols <- intersect(dmso_cols, colnames(expr_mat))
if (length(dmso_cols) == 0L) stop("No DMSO/control samples found in expression matrix.")
expr_lookup <- setNames(rowMeans(expr_mat[, dmso_cols, drop = FALSE], na.rm = TRUE),
                        rownames(expr_mat))
expr_lookup <- expr_lookup[is.finite(expr_lookup) & expr_lookup >= 0]
message(sprintf("    TPM lookup: %d genes (DMSO mean)", length(expr_lookup)))

if (length(expr_lookup) > 100 && exists("gene_stat") && nrow(gene_stat) > 0) {
  expr_matched_es <- list()
  for (i in seq_along(modes)) {
    md <- modes[[i]]; if (md$near) next; nm <- md$name
    loop_genes <- intersect(unique(toupper(mode_genes[[nm]])), gene_stat$gene)
    chip_g <- intersect(get_chip_reference(md), gene_stat$gene)
    only_loop <- intersect(setdiff(loop_genes, chip_g), names(expr_lookup))
    only_chip <- intersect(setdiff(chip_g, loop_genes), names(expr_lookup))
    if (length(only_loop) < 20 || length(only_chip) < 20) next

    # Full-set ES
    rb_full <- rank_biserial(
      gene_stat$signal_lfc[gene_stat$gene %in% only_loop],
      gene_stat$signal_lfc[gene_stat$gene %in% only_chip])

    # Expression-matched: decile-stratified sampling (ntile for stable bins)
    match_df <- bind_rows(
      data.frame(gene = only_loop, group = "looplook",
                 expr = log2(expr_lookup[only_loop] + 0.1), stringsAsFactors = FALSE),
      data.frame(gene = only_chip, group = "ChIPseeker",
                 expr = log2(expr_lookup[only_chip] + 0.1), stringsAsFactors = FALSE)
    ) %>% filter(is.finite(expr)) %>% mutate(expr_bin = dplyr::ntile(expr, 10))

    smd_fn <- function(x, y) {
      x <- x[is.finite(x)]; y <- y[is.finite(y)]
      psd <- sqrt((var(x) + var(y)) / 2)
      if (!is.finite(psd) || psd == 0) return(NA_real_)
      (mean(x) - mean(y)) / psd
    }

    matched_es_vec <- c()
    smd_after_vec <- c()
    n_matched_vec <- c()
    set.seed(42)
    for (r in 1:500) {
      s_loop <- c(); s_chip <- c()
      for (bin in 1:10) {
        bin_df <- match_df %>% filter(expr_bin == bin)
        loop_in_bin <- bin_df %>% filter(group == "looplook") %>% pull(gene)
        chip_in_bin <- bin_df %>% filter(group == "ChIPseeker") %>% pull(gene)
        n_draw <- min(length(loop_in_bin), length(chip_in_bin), 50L)
        if (n_draw < 2) next
        s_loop <- c(s_loop, sample(loop_in_bin, n_draw))
        s_chip <- c(s_chip, sample(chip_in_bin, n_draw))
      }
      if (length(s_loop) >= 5 && length(s_chip) >= 5) {
        es_r <- rank_biserial_point(
          gene_stat$signal_lfc[gene_stat$gene %in% s_loop],
          gene_stat$signal_lfc[gene_stat$gene %in% s_chip])
        if (is.finite(es_r)) {
          matched_es_vec <- c(matched_es_vec, es_r)
          smd_after_vec <- c(smd_after_vec,
            smd_fn(log2(expr_lookup[s_loop] + 0.1), log2(expr_lookup[s_chip] + 0.1)))
          n_matched_vec <- c(n_matched_vec, length(s_loop))
        }
      }
    }

    if (length(matched_es_vec) >= 50) {
      tpm_loop_before <- log2(expr_lookup[only_loop] + 0.1)
      tpm_chip_before <- log2(expr_lookup[only_chip] + 0.1)
      smd_before <- smd_fn(tpm_loop_before, tpm_chip_before)

      expr_matched_es[[nm]] <- data.frame(
        Mode = nm, Full_ES = rb_full$est,
        Matched_ES_median = median(matched_es_vec, na.rm = TRUE),
        Matched_ES_lo = quantile(matched_es_vec, 0.025, na.rm = TRUE, names = FALSE),
        Matched_ES_hi = quantile(matched_es_vec, 0.975, na.rm = TRUE, names = FALSE),
        P_ES_gt_0 = mean(matched_es_vec > 0, na.rm = TRUE),
        N_Only_loop = length(only_loop), N_Only_chip = length(only_chip),
        SMD_before = smd_before,
        SMD_after_median = median(smd_after_vec, na.rm = TRUE),
        SMD_after_lo = quantile(smd_after_vec, 0.025, na.rm = TRUE, names = FALSE),
        SMD_after_hi = quantile(smd_after_vec, 0.975, na.rm = TRUE, names = FALSE),
        AbsSMD_after_median = median(abs(smd_after_vec), na.rm = TRUE),
        P_absSMD_lt_01 = mean(abs(smd_after_vec) < 0.1, na.rm = TRUE),
        Median_N_matched = median(n_matched_vec, na.rm = TRUE),
        MedTPM_loop = median(expr_lookup[only_loop], na.rm = TRUE),
        MedTPM_chip = median(expr_lookup[only_chip], na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }

  if (length(expr_matched_es) > 0) {
    em_df <- bind_rows(expr_matched_es) %>%
      mutate(Pipeline = factor(case_when(
        grepl("^anno_", Mode) ~ "Basic", grepl("^refined_", Mode) ~ "E-refined",
        grepl("^chrom_only_", Mode) ~ "C-refined", grepl("^chrom_", Mode) ~ "I-refined",
        TRUE ~ "Other"), levels = c("Basic","E-refined","I-refined","C-refined")))

    p_em <- ggplot(em_df, aes(x = Full_ES, y = Matched_ES_median, color = Pipeline)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
      geom_point(size = 2, alpha = 0.7) +
      geom_errorbar(aes(ymin = Matched_ES_lo, ymax = Matched_ES_hi), width = 0.01, alpha = 0.4) +
      scale_color_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                     "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
      labs(title = "Expression-Matched vs Full-Set Effect Size",
           subtitle = sprintf("%d modes | r = %.3f | Points near diagonal = TPM not confounding",
                              nrow(em_df), cor(em_df$Full_ES, em_df$Matched_ES_median, use = "pairwise")),
           x = "Full-set ES (Only_looplook vs Only_ChIPseeker)",
           y = "Expression-matched ES (median of 500 iterations)") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5),
            legend.position = "bottom")
    queueable_ggsave(file.path(expr_match_dir, "ExpressionMatched_ES.pdf"), p_em, width = 7, height = 6.5)

    wb_em <- createWorkbook()
    addWorksheet(wb_em, "ExpressionMatched_ES"); writeData(wb_em, "ExpressionMatched_ES", as.data.frame(em_df))
    addWorksheet(wb_em, "Balance_Check"); writeData(wb_em, "Balance_Check", as.data.frame(em_df %>%
      dplyr::select(Mode, Pipeline, N_Only_loop, N_Only_chip, MedTPM_loop, MedTPM_chip,
                    SMD_before, SMD_after_median, SMD_after_lo, SMD_after_hi,
                    AbsSMD_after_median, P_absSMD_lt_01, Median_N_matched)))
    saveWorkbook(wb_em, file.path(expr_match_dir, "ExpressionMatched_ES.xlsx"), overwrite = TRUE)
    med_smd_before <- median(abs(em_df$SMD_before), na.rm = TRUE)
    med_smd_after <- median(em_df$AbsSMD_after_median, na.rm = TRUE)
    p_balanced <- mean(em_df$P_absSMD_lt_01 > 0.5, na.rm = TRUE)
    message(sprintf("    Expression-matched ES: %d modes, median |SMD| before=%.3f after=%.3f; %d%% modes with |SMD|<0.1 >50%% of iter",
                    nrow(em_df), med_smd_before, med_smd_after, round(p_balanced * 100)))
  }
} else {
  warning("Expression-matched ES skipped: TPM data insufficient")
}

}  # end if (FALSE) Expression-Matched ES

message("\n\n>>> ==== Multi-Dimensional Mode Ranking ====")
rank_dir <- file.path(base_out_dir, "mode_ranking")
dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)

# Read NES summary from the largest sample size
sz_max <- max(sample_sizes)
nes_file <- file.path(base_out_dir, paste0("size", sz_max),
                      sprintf("GSEA_Resampling_%dModes_%diter.xlsx", length(modes), n_iterations))
if (file.exists(nes_file)) {
  # Read from Excel
  nes_data <- tryCatch(read.xlsx(nes_file, sheet = "Mode_x_ID_Mean_NES"), error = function(e) { warning("Failed to read: ", nes_file); NULL })
  if (!is.null(nes_data) && ncol(nes_data) >= 2) {
    colnames(nes_data)[1] <- "Mode"
    # Find the looplook NES column (may be named "looplook" or contain "looplook")
    loop_col <- grep("looplook", colnames(nes_data), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(loop_col)) {
      nes_loop <- nes_data %>%
        dplyr::select(Mode, looplook_NES = all_of(loop_col)) %>%
        filter(!grepl("ctrl", Mode)) %>%
        mutate(Mode = as.character(Mode), looplook_NES = as.numeric(looplook_NES))
    } else {
      nes_loop <- data.frame(Mode = character(), looplook_NES = numeric())
    }
  } else {
    nes_loop <- data.frame(Mode = character(), looplook_NES = numeric())
  }
} else {
  # Fallback: read from the configured primary-size tmp_sim CSVs
  nes_loop <- data.frame()
  sim_dir <- file.path(base_out_dir, paste0("size", max(sample_sizes)), "tmp_sim")
  if (dir.exists(sim_dir)) {
    csv_files <- list.files(sim_dir, pattern = "sim_.*\\.csv$", full.names = TRUE)
    if (length(csv_files) > 0) {
      sim_list <- lapply(csv_files, function(f) {
        df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", f); NULL })
        if (is.null(df) || nrow(df) == 0) return(NULL)
        df
      })
      sim_list <- sim_list[!vapply(sim_list, is.null, logical(1))]
      if (length(sim_list) > 0) {
        sim_df_tmp <- bind_rows(sim_list)
        nes_loop <- sim_df_tmp %>%
          filter(ID == "looplook", !grepl("ctrl", Mode)) %>%
          group_by(Mode) %>%
          summarise(looplook_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
          mutate(Mode = as.character(Mode))
      }
    }
  }
}


# Read signal density summary
sig_file <- file.path(base_out_dir, "Signal_Density_Analysis.xlsx")
if (file.exists(sig_file)) {
  sig_comp <- tryCatch(read.xlsx(sig_file, sheet = "Pairwise_Comparisons"),
                       error = function(e) { warning("Failed to read Signal_Density: ", sig_file); NULL })
  sig_valid <- tryCatch(read.xlsx(sig_file, sheet = "Validation"),
                        error = function(e) { warning("Failed to read Signal_Density Validation: ", sig_file); NULL })
}

# Normalize and compute weighted composite score (Delta_NES 50%, Panel_Recall 20%, ES 30%)
safe_z <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) == 0L) return(out)
  sx <- sd(x[ok])
  if (!is.finite(sx) || sx == 0) { out[ok] <- 0; return(out) }
  out[ok] <- (x[ok] - mean(x[ok])) / sx
  out
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

# Build composite score per mode
ranking <- data.frame(Mode = names(mode_genes), stringsAsFactors = FALSE)
ranking <- ranking %>% filter(!grepl("ctrl", Mode))

# Metric 1: Paired Delta_NES (primary: uses same-mode paired subtraction)
if (exists("paired_nes_summary") && nrow(paired_nes_summary) > 0) {
  ranking <- ranking %>%
    left_join(paired_nes_summary %>% dplyr::select(Mode, Median_Delta_NES, P_Delta_lt_0), by = "Mode")
  ranking$Delta_score <- -ranking$Median_Delta_NES  # negative Delta = looplook stronger
  ranking$Direction_Stable <- ranking$P_Delta_lt_0 > 0.95
} else {
  ranking$Median_Delta_NES <- NA_real_
  ranking$Delta_score <- NA_real_
  ranking$Direction_Stable <- NA
}

# Metric 2: Literature panel recall (descriptive plausibility, NOT ground truth)
if (exists("sig_valid") && !is.null(sig_valid) && nrow(sig_valid) > 0) {
  valid_pr <- sig_valid %>%
    filter(Set == "All_looplook",
           if ("Panel" %in% colnames(sig_valid)) Panel == "brd4_associated" else TRUE) %>%
    group_by(Mode) %>%
    summarise(Panel_Recall = safe_max(Recall), .groups = "drop")
  ranking <- ranking %>% left_join(valid_pr, by = "Mode")
} else {
  ranking$Panel_Recall <- NA_real_
}

# Metric 3: Head-to-head effect size (Only_looplook vs Only_ChIPseeker)
if (exists("sig_comp") && !is.null(sig_comp) && nrow(sig_comp) > 0) {
  sig_es <- sig_comp %>%
    filter(Set_A == "Only_looplook", Set_B == "Only_ChIPseeker") %>%
    group_by(Mode) %>%
    summarise(EffectSize = ES[1], .groups = "drop")
  ranking <- ranking %>% left_join(sig_es, by = "Mode")
} else {
  ranking$EffectSize <- NA_real_
}

# Descriptive only (not in composite): gene set uniqueness
unique_scores <- sapply(ranking$Mode, function(nm) {
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) return(NA_real_)
  loop_genes <- unique(toupper(genes))
  md <- modes[[which(sapply(modes, function(m) m$name == nm))[1]]]
  only_loop <- setdiff(loop_genes, get_chip_reference(md))
  length(only_loop) / max(1, length(loop_genes))
})
ranking$Uniqueness <- unique_scores

ranking <- ranking %>%
  mutate(
    Delta_z  = safe_z(Delta_score),
    Panel_z  = safe_z(Panel_Recall),
    ES_z     = safe_z(EffectSize),
    N_metrics_available = is.finite(Delta_score) + is.finite(Panel_Recall) + is.finite(EffectSize),
    Ranking_status = ifelse(N_metrics_available == 3, "Complete", "Incomplete"),
    Composite = 0.50 * Delta_z + 0.20 * Panel_z + 0.30 * ES_z,
    Composite = ifelse(N_metrics_available == 3, Composite, NA_real_),
    Pipeline = factor(case_when(
      grepl("^anno_", Mode)       ~ "Basic",
      grepl("^refined_", Mode)    ~ "E-refined",
      grepl("^chrom_only_", Mode) ~ "C-refined",
      grepl("^chrom_", Mode)      ~ "I-refined",
      TRUE ~ "Other"))) %>%
  arrange(desc(Composite))

# ── Ranking plot (Complete modes only) ──
ranking_complete <- ranking %>% filter(Ranking_status == "Complete", is.finite(Composite))
ranking_incomplete <- ranking %>% filter(Ranking_status != "Complete")
ranking_plot_df <- ranking_complete %>% slice_head(n = 20)
p_rank <- ggplot(ranking_plot_df, aes(x = reorder(Mode, Composite), y = Composite, fill = Pipeline)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", Composite)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Basic" = "#1F77B4", "E-refined" = "#9467BD",
                                "I-refined" = "#E64B35", "C-refined" = "#FFA500")) +
  labs(title = "Exploratory Composite Mode Ranking (LFC ranked list)",
       subtitle = "50% Delta_NES + 20% Panel_Recall + 30% effect size | Panel_Recall descriptive only",
       x = NULL, y = "Composite Score (z-score mean)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5))
queueable_ggsave(file.path(rank_dir, "Mode_Composite_Ranking_LFC.pdf"), p_rank, width = 9, height = 6)

# ── Top 3 per pipeline (Complete modes only) ──
top3 <- ranking_complete %>% group_by(Pipeline) %>% slice_max(Composite, n = 3) %>% ungroup()
message("\n    === Top modes by pipeline ===")
for (p in unique(top3$Pipeline)) {
  top_p <- top3 %>% filter(Pipeline == p) %>% arrange(desc(Composite))
  for (k in seq_len(nrow(top_p))) {
    message(sprintf("      %-15s %-30s Composite=%.3f  Delta=%.3f  Panel=%.3f  ES=%.3f  Unique=%.1f%%",
      p, top_p$Mode[k], top_p$Composite[k],
      top_p$Median_Delta_NES[k], top_p$Panel_Recall[k], top_p$EffectSize[k],
      100*top_p$Uniqueness[k]))
  }
}

# ── Export ──
wb_rank <- createWorkbook()
addWorksheet(wb_rank, "Mode_Ranking")
writeData(wb_rank, "Mode_Ranking", as.data.frame(ranking))
addWorksheet(wb_rank, "Top3_Per_Pipeline")
writeData(wb_rank, "Top3_Per_Pipeline", as.data.frame(top3))
if (nrow(ranking_incomplete) > 0) {
  addWorksheet(wb_rank, "Incomplete_Modes")
  writeData(wb_rank, "Incomplete_Modes", as.data.frame(ranking_incomplete))
}
out_rank <- file.path(rank_dir, "Mode_Composite_Ranking_LFC.xlsx")

# Ranking robustness first
if (exists("ranking_complete") && nrow(ranking_complete) >= 10) {
  r <- ranking_complete
    alt_weights <- list(
      "Equal"      = c(Delta=0.333, Panel=0.333, ES=0.333),
      "ES_Primary" = c(Delta=0.25, Panel=0.25, ES=0.50),
      "Delta_Low"  = c(Delta=0.10, Panel=0.50, ES=0.40),
      "Panel_Low"  = c(Delta=0.50, Panel=0.10, ES=0.40))
    robustness_df <- data.frame(Mode = r$Mode, Pipeline = r$Pipeline, stringsAsFactors = FALSE)
    for (nm in names(alt_weights)) {
      w <- alt_weights[[nm]]
      rw <- r %>% mutate(
        Composite = w["Delta"]*Delta_z + w["Panel"]*Panel_z + w["ES"]*ES_z)
    top3_rw <- rw %>% slice_max(Composite, n = 3) %>% pull(Mode)
    robustness_df[[nm]] <- rw$Composite
    robustness_df[[paste0(nm, "_top1")]] <- ifelse(rw$Mode == top3_rw[1], "★", "")
    message(sprintf("    %-12s top: %s", nm, paste(top3_rw, collapse=", ")))
  }
  message("    Ranking robustness: see sheet 'Robustness' in ranking Excel")
  addWorksheet(wb_rank, "Robustness"); writeData(wb_rank, "Robustness", as.data.frame(robustness_df))
}

saveWorkbook(wb_rank, out_rank, overwrite = TRUE)

message("\n>>> All optimizations complete!")

# ── Parameter recommendation ──────────────────────────────────────────────────
message("\n>>> Parameter recommendation for looplook usage:")
top_modes <- if (exists("ranking_complete") && nrow(ranking_complete) >= 3) {
  ranking_complete %>% arrange(desc(Composite)) %>% head(5) %>% pull(Mode)
} else {
  c("refined_all_F", "refined_promoter_F", "chrom_all_F")
}
base_names <- gsub("_hop1$", "", top_modes)

message("\n    SCENARIO-BASED RECOMMENDATIONS:")
message("    ┌─────────────────────────────────────────────────────────────┐")
message("    │ Best in this BRD4 case study: ", paste(head(base_names, 2), collapse=", "))
message("    │ Recall-oriented:           filled + all-map modes")
message("    │ Conservative/promoter-focused: strict + promoter-map modes")
message("    │ Distal-regulatory:         chromatin-refined all-map modes")
message("    │ Expression-informed:       refined (E-refined) pipeline")
message("    └─────────────────────────────────────────────────────────────┘")
message("    Note: These are case-study observations, not universal claims.")
message("          Hop0 represents primary targets; hop1 represents primary + expanded targets.")
message("          Expanded-only biological value is reported separately in expanded_only/.")
message("    Detailed ranking: ", rank_dir)
message("    Output: ", base_out_dir)

# ══════════════════════════════════════════════════════════════════════════════
# Summary Multi-Panel Figure + Take-Home Quantification
# ══════════════════════════════════════════════════════════════════════════════
message("\n\n>>> ==== SUMMARY FIGURE ====")
summary_dir <- file.path(base_out_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Panel A: NES heatmap-style bar (32 modes, looplook NES only)
if (exists("nes_loop") && nrow(nes_loop) > 0) {
  nes_p <- nes_loop %>%
    mutate(Pipeline = factor(case_when(
      grepl("^anno_",Mode)~"Basic",grepl("^refined_",Mode)~"E-refined",
      grepl("^chrom_only_",Mode)~"C-refined",grepl("^chrom_",Mode)~"I-refined",TRUE~"Other"),
      levels=c("Basic","E-refined","I-refined","C-refined")))
  pA <- ggplot(nes_p, aes(x = reorder(Mode, looplook_NES), y = looplook_NES, fill = Pipeline)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")) +
    labs(title = "A. Looplook GSEA Enrichment (NES)", x = NULL, y = "Mean NES") +
    theme_minimal(base_size = 8) + theme(legend.position = "right", axis.text.y = element_text(size = 5))
} else { pA <- ggplot() + labs(title = "A. NES not available") + theme_void() }

# Panel B: Delta_NES by pipeline
if (exists("delta_stats_all") && nrow(delta_stats_all) > 0) {
  pB <- ggplot(delta_stats_all, aes(x = Method, y = Median_Delta, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.5) +
    facet_wrap(~ hop, ncol = 2) +
    scale_fill_manual(values = c("anno"="#1F77B4","refined"="#9467BD","chromatin"="#E64B35","chromatin_only"="#FFA500")) +
    labs(title = "B. Paired Delta_NES (looplook - Background)", x = NULL, y = "Median Delta_NES") +
    theme_minimal(base_size = 8) + theme(legend.position = "none")
} else { pB <- ggplot() + labs(title = "B. Delta_NES not available") + theme_void() }

# Panel C: Distance-ES (binned)
if (exists("se_dist_bin_df") && nrow(se_dist_bin_df) > 0) {
  pC <- ggplot(se_dist_bin_df, aes(x = DistBin, y = ES, fill = Pipeline)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6, linewidth = 0.3) +
    scale_fill_manual(values = c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")) +
    labs(title = "C. Distance vs Effect Size", x = "Peak-to-TSS Distance", y = "ES") +
    theme_minimal(base_size = 7) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
} else { pC <- ggplot() + labs(title = "C. Distance-ES not available") + theme_void() }

# Panel D: Case study signal distribution
if (exists("cs_summary_df") && nrow(cs_summary_df) > 0) {
  pD <- ggplot(cs_summary_df, aes(x = reorder(Mode, Assigned_Genes), y = Median_signal_lfc, fill = Pipeline)) +
    geom_col(alpha = 0.85, width = 0.7) + coord_flip() +
    geom_text(aes(label = Assigned_Genes), hjust = -0.1, size = 2.5) +
    scale_fill_manual(values = c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")) +
    labs(title = "D. Distal Peak Targets (case study)", x = NULL, y = "Median signal_lfc") +
    theme_minimal(base_size = 8) + theme(legend.position = "right")
} else { pD <- ggplot() + labs(title = "D. Case study not available") + theme_void() }

# Panel E: Composite ranking (top 15)
if (exists("ranking_complete") && nrow(ranking_complete) > 0) {
  r_top <- ranking_complete %>% slice_max(Composite, n = 15)
  pE <- ggplot(r_top, aes(x = reorder(Mode, Composite), y = Composite, fill = Pipeline)) +
    geom_col(alpha = 0.85, width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Basic"="#1F77B4","E-refined"="#9467BD","I-refined"="#E64B35","C-refined"="#FFA500")) +
    labs(title = "E. Top 15 Modes (Composite Score)", x = NULL, y = "Composite") +
    theme_minimal(base_size = 8) + theme(legend.position = "none")
} else { pE <- ggplot() + labs(title = "E. Ranking not available") + theme_void() }

# Panel F: Take-home numbers — three-layer quantification
if (exists("cs_summary_df") && nrow(cs_summary_df) > 0 &&
    exists("paired_nes_summary") && nrow(paired_nes_summary) > 0) {
  # Layer 1: genome-wide
  n_stable_modes <- sum(paired_nes_summary$P_Delta_lt_0 > 0.95)
  pct_stable <- round(100 * n_stable_modes / nrow(paired_nes_summary))
  med_global_es <- if (exists("global_es_df")) round(median(global_es_df$ES, na.rm = TRUE), 3) else NA

  # Layer 2: distal-only (case study)
  best_cs <- cs_summary_df %>% slice_max(Assigned_Genes, n = 1)
  n_distal_targets <- best_cs$Assigned_Genes[1]
  n_distal_down <- best_cs$N_Downregulated[1]
  med_distal_signal <- round(best_cs$Median_signal_lfc[1], 2)
  pct_down <- round(100 * n_distal_down / max(1, n_distal_targets))

  # Layer 3: SE (if available)
  se_line <- ""
  if (use_se_stratification && exists("se_te_es") && nrow(se_te_es) > 0) {
    med_se_es <- round(median(se_te_es$ES[se_te_es$Category == "Super_Enhancer"], na.rm = TRUE), 3)
    med_te_es <- round(median(se_te_es$ES[se_te_es$Category == "Typical_Enhancer"], na.rm = TRUE), 3)
    se_line <- paste0(
      " Super-enhancer targets: ES_median=", med_se_es,
      "  |  Typical enhancer: ES_median=", med_te_es, "\n")
  }

  th_text <- paste0(
    " GENOME-WIDE (evaluated modes):\n",
    "   ", pct_stable, "% modes with Pr(Delta_NES<0)>0.95\n",
    "   Global ES median = ", med_global_es, "\n\n",
    " DISTAL-ONLY (ChIPseeker-distal peaks):\n",
    "   Top mode assigns ", n_distal_targets, " distal targets\n",
    "   ", pct_down, "% downregulated (median signal=", med_distal_signal, ")", "\n",
    se_line)

  pF <- ggplot() +
    annotate("text", x = 0, y = 0, label = th_text, size = 4,
      hjust = 0, vjust = 0.5, family = "mono", color = "grey20") +
    labs(title = "F. Three-Layer Quantification") +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold"))
} else { pF <- ggplot() + labs(title = "F. Key numbers pending") + theme_void() }

# Combine into summary figure
summary_fig <- ggpubr::ggarrange(pA, pB, pC, pD, pE, pF,
  ncol = 3, nrow = 2, common.legend = FALSE,
  widths = c(1.3, 1, 1), heights = c(1, 1))
queueable_ggsave(file.path(summary_dir, "Summary_MultiPanel_Figure.pdf"),
  summary_fig, width = 18, height = 12)
message("    Summary figure: ", summary_dir, "/Summary_MultiPanel_Figure.pdf")

# ── Conclusion paragraph (ready for manuscript) ───────────────────────────────
message("\n>>> ================================================================")
message(">>> DRAFT CONCLUSION (copy into manuscript Results/Discussion):")
message(">>> ================================================================")
n_total_peaks <- nrow(peak_anno_df)
n_distal_peaks <- length(distal_peak_coords)
pct_distal <- round(100 * n_distal_peaks / n_total_peaks)
if (exists("cs_summary_df") && nrow(cs_summary_df) > 0) {
  best_cs <- cs_summary_df %>% slice_max(Assigned_Genes, n = 1)
  if (exists("paired_nes_summary") && nrow(paired_nes_summary) > 0)
    n_stable <- sum(paired_nes_summary$P_Delta_lt_0 > 0.95)
  else n_stable <- 0
  best_mode <- best_cs$Mode[1]
  n_tg <- best_cs$Assigned_Genes[1]
  pct_down_cs <- round(100 * best_cs$N_Downregulated[1] / max(1, n_tg))
  med_signal_cs <- round(best_cs$Median_signal_lfc[1], 2)
} else { n_stable <- 0; best_mode <- "TBD"; n_tg <- 0; pct_down_cs <- 0; med_signal_cs <- 0 }

message("")
message("    Among ", n_total_peaks, " BRD4 ChIP-seq peaks, ", pct_distal, "% (", n_distal_peaks,
  ") were classified as distal intergenic by ChIPseeker.")
message("    ChIPseeker still provides a nearest-gene annotation for these peaks;")
message("    looplook adds loop-supported candidate target assignments that are not")
message("    restricted to simple nearest-gene proximity. The ", best_mode,
  " mode identified ", n_tg, " candidate genes from these distal-classified peaks,")
message("    of which ", pct_down_cs, "% were downregulated after BRD4 degradation")
message("    (ARV-825), with median signal_lfc = ", med_signal_cs, ".")
if (n_stable > 0) {
  n_mode_total <- nrow(paired_nes_summary)
  message("    Across the ", n_mode_total, " evaluated parameter combinations, ",
    n_stable, "/", n_mode_total, " modes had Pr(Delta_NES < 0) > 0.95,")
  message("    indicating a consistent shift of looplook-assigned candidate genes")
  message("    toward the downregulated end relative to the prespecified control pool.")
}
if (use_se_stratification && exists("se_te_es") && nrow(se_te_es) > 0) {
  med_se <- round(median(se_te_es$ES[se_te_es$Category == "Super_Enhancer"], na.rm = TRUE), 3)
  med_te <- round(median(se_te_es$ES[se_te_es$Category == "Typical_Enhancer"], na.rm = TRUE), 3)
  if (is.finite(med_se) && is.finite(med_te) && med_se > med_te) {
    message("    In this BRD4 dataset, ROSE-defined super-enhancer-linked candidates")
    message("    had a higher median effect size than typical-enhancer-linked candidates")
    message("    (SE=", med_se, "; TE=", med_te, "). This is an association and does")
    message("    not by itself establish direct enhancer-gene causality.")
  } else {
    message("    ROSE-defined SE and TE strata had median effect sizes of ", med_se,
      " and ", med_te, ", respectively; no stronger-SE conclusion was imposed.")
  }
}
message("")
message(">>> ================================================================")
message(">>> END DRAFT CONCLUSION")
message(">>> ================================================================")

# ══════════════════════════════════════════════════════════════════════════════
# LFC_pv 全套分析：sign(LFC) × -log10(pvalue) 作为 ranked list
# ══════════════════════════════════════════════════════════════════════════════

if (use_glist_pv) {
message("\n\n>>> ========================================")
message(">>> ==== LFC_pv FULL ANALYSIS MODULE ====")
message(">>> ========================================")

lfc_pv_base <- file.path(base_out_dir, "LFC_pv")
dir.create(lfc_pv_base, recursive = TRUE, showWarnings = FALSE)

# ── LFC_pv.1: GSEA 重抽样 (每个 sample size) ──────────────────────────────
for (sz in sample_sizes) {
  sample_size_pv <- sz
  lfc_pv_sz_dir <- file.path(lfc_pv_base, paste0("size", sz))
  dir.create(lfc_pv_sz_dir, recursive = TRUE, showWarnings = FALSE)
  lfc_pv_tmp <- file.path(lfc_pv_sz_dir, "tmp_sim")
  dir.create(lfc_pv_tmp, recursive = TRUE, showWarnings = FALSE)
  ensure_cache_signature(
    lfc_pv_tmp,
    c(benchmark_cache_signature_base, list(
      Module = "main_lfc_pv",
      SampleSize = as.integer(sz),
      RankedListMD5 = md5_object(glist_pv)
    ))
  )

  message(sprintf("\n>>> LFC_pv GSEA: size=%d, %d iterations per mode...", sz, n_iterations))

  for (idx in seq_along(modes)) {
    md <- modes[[idx]]
    nm <- md$name
    genes <- mode_genes[[nm]]
    if (length(genes) == 0) next

    tmp_csv <- file.path(lfc_pv_tmp, paste0("sim_", nm, ".csv"))
    expected_main_terms <- if (grepl("^ctrl_", nm)) c("ChIPseeker", "Background") else c("looplook", "Background")
    if (validate_csv_cache(tmp_csv, c("ID", "NES", "pvalue", "Iteration", "Mode"), n_iterations, expected_main_terms)) {
      message(sprintf("    [%2d/32] %-30s SKIP (validated cache)", idx, nm))
      next
    }
    if (file.exists(tmp_csv)) {
      warning(sprintf("    [%2d/32] %-30s cache invalid, will recompute", idx, nm))
      unlink(tmp_csv)
    }

    pool_loop <- intersect(genes, all_genes)
    pool_bg <- setdiff(all_genes, union(genes, pool_chip))
    is_ctrl <- md$near
    term_loop <- if (is_ctrl) "ChIPseeker" else "looplook"

    # Pool-size preflight
    if (length(pool_loop) < gsea_min_size || length(pool_bg) < gsea_min_size) {
      gsea_status_log[[length(gsea_status_log) + 1L]] <- data.frame(
        SampleSize = sample_size_pv, Module = "main_pv", Mode = nm, Iteration = NA_integer_,
        Status = "insufficient_pool",
        ExpectedTerms = paste(c(term_loop, "Background"), collapse = ";"),
        ReturnedTerms = "", MissingTerms = paste(c(term_loop, "Background"), collapse = ";"),
        Error = NA_character_, stringsAsFactors = FALSE)
      message(sprintf("    [LFC_pv %2d/32] %-30s SKIP (pool too small)", idx, nm))
      next
    }

    if (length(pool_bg) < sample_size_pv) {
      warning(sprintf("[LFC_pv %s] control pool has %d genes; adaptive 80%% sampling will be used", nm, length(pool_bg)))
    }

    sim_local <- checkpointed_lapply(
  1:n_iterations,
  function(i) {
        set.seed(400 + i)
        n_actual_pv <- min(sample_n(pool_loop, sample_size_pv), sample_n(pool_bg, sample_size_pv))
        s_loop <- sample(pool_loop, n_actual_pv)
        s_bg   <- sample(pool_bg, n_actual_pv)
        t2g <- bind_rows(
          data.frame(term = term_loop,   gene = s_loop),
          data.frame(term = "Background", gene = s_bg))
        expected_terms_pv <- c(term_loop, "Background")
        gsea_out <- safe_gsea_once(glist_pv, t2g, expected_terms_pv, sample_size_pv,
          "main_pv", nm, i)

        if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
          df <- gsea_out$result[, c("ID", "NES", "pvalue")]
          df$Iteration <- i; df$Mode <- nm
          df$Requested_N <- sample_size_pv
          actual_map <- setNames(c(length(s_loop), length(s_bg)), c(term_loop, "Background"))
          pool_map   <- setNames(c(length(pool_loop), length(pool_bg)), c(term_loop, "Background"))
          df$Actual_N <- unname(actual_map[df$ID])
          df$Pool_N   <- unname(pool_map[df$ID])
          df$Sampling_fraction <- df$Actual_N / pmax(1, df$Pool_N)
          list(result = df, status = gsea_out$status)
        } else {
          list(result = NULL, status = gsea_out$status)
        }
      },
  checkpoint_dir = file.path(lfc_pv_tmp, ".iteration_checkpoints", nm),
  label = sprintf("main_pv size=%d mode=%s", sample_size_pv, nm),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)
    statuses <- lapply(sim_local, function(x) x$status)
    gsea_status_log <- c(gsea_status_log, statuses[!vapply(statuses, is.null, logical(1))])
    sim_local <- lapply(sim_local, function(x) x$result)
    sim_local <- sim_local[vapply(sim_local, is.data.frame, logical(1))]
    sim_mode_df <- if (length(sim_local) > 0) {
      bind_rows(sim_local) %>% filter(!is.na(NES))
    } else {
      data.frame(ID = character(), NES = numeric(), pvalue = numeric(),
                 Iteration = integer(), Mode = character())
    }
    write.csv(sim_mode_df, tmp_csv, row.names = FALSE)
    rm(sim_local, sim_mode_df); gc()
    message(sprintf("    [%2d/32] %-30s -> %s", idx, nm, tmp_csv))
  }
}

# ── LFC_pv.2: 合并所有 sample size 的 CSV + Excel 导出 ──────────────────────
message("\n>>> LFC_pv: Merging and exporting Excel...")
for (sz in sample_sizes) {
  lfc_pv_sz_dir <- file.path(lfc_pv_base, paste0("size", sz))
  lfc_pv_tmp <- file.path(lfc_pv_sz_dir, "tmp_sim")
  csv_files <- list.files(lfc_pv_tmp, pattern = "sim_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) next

  sim_list <- lapply(csv_files, function(f) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", f); NULL })
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  })
  sim_list <- sim_list[!vapply(sim_list, is.null, logical(1))]
  sim_pv <- bind_rows(sim_list)
  sim_pv$Mode <- factor(sim_pv$Mode, levels = sapply(modes, `[[`, "name"))

  mode_mean_pv <- sim_pv %>%
    group_by(Mode, ID) %>%
    summarise(Mean_NES = mean(NES), Median_NES = median(NES),
      SD_NES = sd(NES), .groups = "drop")

  wb_pv <- createWorkbook()
  .safe_add_sheet(wb_pv, "Mode_x_ID_Mean_NES",
    mode_mean_pv %>% tidyr::pivot_wider(names_from = ID, values_from = Mean_NES) %>% as.data.frame())
  .safe_add_sheet(wb_pv, "Raw_All", as.data.frame(sim_pv))
  out_xlsx_pv <- file.path(lfc_pv_sz_dir, sprintf("GSEA_LFC_pv_%diter.xlsx", n_iterations))
  saveWorkbook(wb_pv, out_xlsx_pv, overwrite = TRUE)
  message(sprintf("    size=%d: %d files, %d rows → %s", sz, length(csv_files), nrow(sim_pv), out_xlsx_pv))
}

# ── LFC_pv.3: NES 可视化 (每个 sample size) ────────────────────────────────
for (sz in sample_sizes) {
  lfc_pv_sz_dir <- file.path(lfc_pv_base, paste0("size", sz))
  lfc_pv_tmp <- file.path(lfc_pv_sz_dir, "tmp_sim")
  csv_files <- list.files(lfc_pv_tmp, pattern = "sim_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) next

  sim_list <- lapply(csv_files, function(f) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) { warning("Failed to read: ", f); NULL })
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  })
  sim_pv <- bind_rows(sim_list[!vapply(sim_list, is.null, logical(1))])
  sim_pv$Mode <- factor(sim_pv$Mode, levels = sapply(modes, `[[`, "name"))

  # Density plot
  sim_pv_lb <- sim_pv %>%
    filter(!grepl("ctrl", Mode)) %>%
    mutate(
      Family = factor(mode_family_map[as.character(Mode)],
        levels = c("Basic (all)", "Refined (all)", "Chromatin (all)", "Chromatin_Only (all)",
                   "Basic (promoter)", "Refined (promoter)", "Chromatin (promoter)", "Chromatin_Only (promoter)")),
      Variant = factor(mode_variant_map[as.character(Mode)],
        levels = c("F (hop=0)", "T (hop=0)", "F (hop=1)", "T (hop=1)")))
  lb_med_pv <- sim_pv_lb %>% group_by(Family, Variant, ID) %>%
    summarise(Median_NES = median(NES, na.rm = TRUE), .groups = "drop")
  lb_colors <- c("looplook" = "#E64B35", "Background" = "#B0BEC5")

  p <- ggplot(sim_pv_lb, aes(x = NES, fill = ID)) +
    geom_density(alpha = 0.7, color = NA) +
    geom_vline(data = lb_med_pv, aes(xintercept = Median_NES, color = ID), linetype = "dashed", linewidth = 0.4) +
    scale_fill_manual(values = lb_colors) + scale_color_manual(values = lb_colors, guide = "none") +
    facet_grid(rows = vars(Family), cols = vars(Variant), scales = "free_y") +
    labs(title = "NES: LFC_pv ranked list (sign(LFC)*-log10(pvalue))",
         subtitle = sprintf("%d subsamples per mode (size=%d)", n_iterations, sz),
         x = "NES", y = "Density") +
    theme_classic() + theme(legend.position = "bottom", strip.text = element_text(size = 9, face = "bold"))
  close_batch_devices()
  queueable_ggsave(file.path(lfc_pv_sz_dir, "NES_Density_AllModes_LFC_pv.pdf"), p, width = 12, height = 10)
  rm(p, sim_pv_lb, lb_med_pv); gc()
  message(sprintf("    LFC_pv density plot saved: size=%d", sz))
}

# ── LFC_pv.4: Unique/Intersection GSEA ─────────────────────────────────────
message("\n>>> LFC_pv: Unique/Intersection GSEA...")
for (sz in sample_sizes) {
  lfc_pv_sz_dir <- file.path(lfc_pv_base, paste0("size", sz))
  lfc_pv_uniq_tmp <- file.path(lfc_pv_sz_dir, "tmp_unique")
  dir.create(lfc_pv_uniq_tmp, recursive = TRUE, showWarnings = FALSE)
  ensure_cache_signature(
    lfc_pv_uniq_tmp,
    c(benchmark_cache_signature_base, list(
      Module = "unique_lfc_pv",
      SampleSize = as.integer(sz),
      RankedListMD5 = md5_object(glist_pv)
    ))
  )

  for (idx in seq_along(modes)) {
    md <- modes[[idx]]; nm <- md$name
    genes <- mode_genes[[nm]]
    if (length(genes) == 0) next

    loop_genes <- intersect(unique(toupper(genes)), names(glist_pv))
    chip_g <- intersect(get_chip_reference(md), names(glist_pv))
    only_loop <- setdiff(loop_genes, chip_g)
    inter_set <- intersect(loop_genes, chip_g)
    only_chip <- setdiff(chip_g, loop_genes)
    pv_unique_pools <- list(only_looplook = only_loop, intersection = inter_set, only_ChIPseeker = only_chip)
    pv_unique_pools <- pv_unique_pools[lengths(pv_unique_pools) >= gsea_min_size]
    if (length(pv_unique_pools) < 2) {
      gsea_status_log[[length(gsea_status_log) + 1L]] <- make_mode_status(
        sz, "unique_pv", nm, "insufficient_terms",
        paste(c("only_looplook", "intersection", "only_ChIPseeker"), collapse = ";"),
        paste0("only ", length(pv_unique_pools), " categories >= ", gsea_min_size))
      next
    }

    expected_terms_upv <- names(pv_unique_pools)
    uniq_csv <- file.path(lfc_pv_uniq_tmp, paste0("uniq_", nm, ".csv"))
    if (validate_csv_cache(uniq_csv, c("ID", "NES", "pvalue", "Iteration", "Mode"),
        n_iterations, expected_terms_upv)) next
    if (file.exists(uniq_csv)) {
      warning(sprintf("    [unique_pv %2d/32] %-30s cache invalid, will recompute", idx, nm))
      unlink(uniq_csv)
    }

    uniq_local <- checkpointed_lapply(
  1:n_iterations,
  function(i) {
        set.seed(500 + i)
        n_upv <- min(sz, lengths(pv_unique_pools))
        sampled_pv <- lapply(pv_unique_pools, function(p) sample(p, n_upv))
        t2g_parts <- lapply(names(sampled_pv), function(tn) data.frame(term=tn, gene=sampled_pv[[tn]], stringsAsFactors=FALSE))
      t2g <- bind_rows(t2g_parts)
      expected_terms_upv <- names(sampled_pv)
      gsea_out <- safe_gsea_once(glist_pv, t2g, expected_terms_upv, sz,
        "unique_pv", nm, i)
      if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
        df <- gsea_out$result[, c("ID", "NES", "pvalue")]
        df$Iteration <- i; df$Mode <- nm
        df$Requested_N <- sz
        actual_map <- sapply(sampled_pv, length)
        pool_map   <- sapply(pv_unique_pools, length)
        df$Actual_N <- unname(actual_map[df$ID])
        df$Pool_N   <- unname(pool_map[df$ID])
        df$Sampling_fraction <- df$Actual_N / pmax(1, df$Pool_N)
        list(result = df, status = gsea_out$status)
      } else {
        list(result = NULL, status = gsea_out$status)
      }
    },
  checkpoint_dir = file.path(lfc_pv_uniq_tmp, ".iteration_checkpoints", nm),
  label = sprintf("unique_pv size=%d mode=%s", sz, nm),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)
    statuses <- lapply(uniq_local, function(x) x$status)
    gsea_status_log <- c(gsea_status_log, statuses[!vapply(statuses, is.null, logical(1))])
    uniq_local <- lapply(uniq_local, function(x) x$result)
    uniq_local <- uniq_local[vapply(uniq_local, is.data.frame, logical(1))]
    uniq_df <- if (length(uniq_local) > 0) bind_rows(uniq_local) %>% filter(!is.na(NES))
               else data.frame(ID=character(), NES=numeric(), pvalue=numeric(), Iteration=integer(), Mode=character())
    write.csv(uniq_df, uniq_csv, row.names=FALSE)
    rm(uniq_local, uniq_df); gc()
  }
  message(sprintf("    LFC_pv Unique GSEA done: size=%d", sz))
}


# ── LFC_pv.6: 背景敏感性 ───────────────────────────────────────────────────
message("\n>>> LFC_pv: Control-pool Sensitivity...")
bg_sens_pv_dir <- file.path(lfc_pv_base, "background_sensitivity")
dir.create(bg_sens_pv_dir, recursive = TRUE, showWarnings = FALSE)
bg_sens_pv_checkpoint_root <- file.path(bg_sens_pv_dir, ".gsea_checkpoints")
ensure_cache_signature(
  bg_sens_pv_checkpoint_root,
  c(benchmark_cache_signature_base, list(
    Module = "background_sensitivity_lfc_pv",
    RankedListMD5 = md5_object(glist_pv),
    Iterations = min(200L, n_iterations),
    SampleSize = as.integer(max(sample_sizes))
  ))
)
bg_sens_modes <- modes[!sapply(modes, `[[`, "near")]
bg_sens_iter <- min(200, n_iterations)
bg_sens_size <- max(sample_sizes)

bg_sens_pv_results <- list()
for (idx in seq_along(bg_sens_modes)) {
  md <- bg_sens_modes[[idx]]; nm <- md$name
  genes <- mode_genes[[nm]]
  if (length(genes) == 0) next

  pool_loop <- intersect(genes, all_genes)
  pool_bg_std <- setdiff(all_genes, genes)
  pool_bg_strict <- setdiff(all_genes, union(genes, pool_chip))

  sim_local <- checkpointed_lapply(
  1:bg_sens_iter,
  function(i) {
      set.seed(700 + i)
      s_loop <- sample(pool_loop, sample_n(pool_loop, bg_sens_size))
      s_bg_std <- sample(pool_bg_std, sample_n(pool_bg_std, bg_sens_size))
      s_bg_strict <- sample(pool_bg_strict, sample_n(pool_bg_strict, bg_sens_size))
      t2g <- bind_rows(
        data.frame(term="looplook", gene=s_loop),
        data.frame(term="Background_std", gene=s_bg_std),
        data.frame(term="Background_strict", gene=s_bg_strict))
      expected_terms_bpv <- c("looplook", "Background_std", "Background_strict")
      gsea_out <- safe_gsea_once(glist_pv, t2g, expected_terms_bpv, bg_sens_size,
        "bg_sens_pv", nm, i)
      if (!is.null(gsea_out$result) && nrow(gsea_out$result) > 0) {
        df <- gsea_out$result[, c("ID", "NES")]
        df$Iteration <- i; df$Mode <- nm
        df$Requested_N <- bg_sens_size
        actual_map <- c(looplook = length(s_loop), Background_std = length(s_bg_std), Background_strict = length(s_bg_strict))
        pool_map   <- c(looplook = length(pool_loop), Background_std = length(pool_bg_std), Background_strict = length(pool_bg_strict))
        df$Actual_N <- unname(actual_map[as.character(df$ID)])
        df$Pool_N   <- unname(pool_map[as.character(df$ID)])
        df$Sampling_fraction <- df$Actual_N / pmax(1, df$Pool_N)
        list(result = df, status = gsea_out$status)
      } else {
        list(result = NULL, status = gsea_out$status)
      }
    },
  checkpoint_dir = file.path(bg_sens_pv_dir, ".gsea_checkpoints", nm),
  label = sprintf("bg_sens_pv mode=%s", nm),
  workers = gsea_n_cores,
  chunk_size = gsea_chunk_size
)
  statuses <- lapply(sim_local, function(x) x$status)
  gsea_status_log <- c(gsea_status_log, statuses[!vapply(statuses, is.null, logical(1))])
  sim_local <- lapply(sim_local, function(x) x$result)
  sim_local <- sim_local[vapply(sim_local, is.data.frame, logical(1))]
  if (length(sim_local) > 0) bg_sens_pv_results[[nm]] <- bind_rows(sim_local) %>% filter(!is.na(NES))
  rm(sim_local); gc()
  message(sprintf("    [%2d/%d] %-30s LFC_pv bg_sens done", idx, length(bg_sens_modes), nm))
}

bg_sens_pv_df <- bind_rows(bg_sens_pv_results)
bg_compare_pv <- bg_sens_pv_df %>%
  filter(grepl("^Background", ID)) %>%
  group_by(Mode, ID) %>%
  summarise(Mean_NES = mean(NES, na.rm=TRUE), SD_NES = sd(NES, na.rm=TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from = ID, values_from = c(Mean_NES, SD_NES)) %>%
  mutate(NES_diff = `Mean_NES_Background_strict` - `Mean_NES_Background_std`)

wb_bg_pv <- createWorkbook()
addWorksheet(wb_bg_pv, "NES_Comparison")
writeData(wb_bg_pv, "NES_Comparison", bg_compare_pv)
addWorksheet(wb_bg_pv, "Raw_Sensitivity_Data")
writeData(wb_bg_pv, "Raw_Sensitivity_Data", as.data.frame(bg_sens_pv_df))
saveWorkbook(wb_bg_pv, file.path(bg_sens_pv_dir, "ControlPool_Sensitivity_LFC_pv.xlsx"), overwrite = TRUE)

med_diff_pv <- median(bg_compare_pv$NES_diff, na.rm=TRUE)
message(sprintf("    LFC_pv bg_sens: median NES_diff = %.4f", med_diff_pv))

message("\n>>> LFC_pv full analysis complete!")
message("    Output: ", lfc_pv_base)

} else {
  message("\n>>> No pvalue column found, LFC_pv analysis skipped.")
}

# Finalize the plot queue before provenance and completion checks.
if (identical(plot_policy, "queue")) {
  write_plot_queue_manifest()
}

# ── Analysis manifest (script_path/script_md5 computed at script start) ──────
meta_path <- file.path(cfg$data_base, cfg$meta_file)

manifest <- data.frame(
  CodeVersion = code_version,
  ScriptFilename = basename(ifelse(is.na(script_path), "unknown", script_path)),
  ScriptMD5 = script_md5,
  RankMetric = "raw_log2FoldChange",
  EffectSignal = "-raw_log2FoldChange",
  EffectMetric = "rank_biserial",
  GSEAControlRule = "exclude_all_ChIPseeker_annotations",
  SamplingRule = "matched_N_requested_or_80pct_not_below_minGSSize",
  Iterations = n_iterations,
  SampleSizes = paste(sample_sizes, collapse = ","),
  PrimarySampleSize = primary_sample_size,
  BootstrapR = bootstrap_R,
  BootstrapSeed = bootstrap_seed,
  TSSWindow = "plus_minus_5kb",
  PromoterWindow = "plus_minus_5kb",
  DistanceTSSRule = "nearest_transcript_TSS",
  DMSOSamples = paste(dmso_cols, collapse = ","),
  CachePolicy = cache_policy,
  ForceRecompute = force_recompute,
  DiffFile = normalizePath(diff_path, mustWork = FALSE),
  ExprFile = normalizePath(expr_path, mustWork = FALSE),
  TargetFile = normalizePath(target_file, mustWork = FALSE),
  DiffMD5 = as.character(tryCatch(tools::md5sum(diff_path), error = function(e) NA)),
  ExprMD5 = as.character(tryCatch(tools::md5sum(expr_path), error = function(e) NA)),
  TargetMD5 = as.character(tryCatch(tools::md5sum(target_file), error = function(e) NA)),
  MainRData = file.path(cfg$rdata_base, cfg$tss_subdirs[["5kb"]]),
  MainRDataMD5 = as.character(tryCatch(tools::md5sum(file.path(cfg$rdata_base, cfg$tss_subdirs[["5kb"]])), error = function(e) NA)),
  TSS1kbRDataMD5 = as.character(tryCatch(tools::md5sum(file.path(cfg$rdata_base, cfg$tss_subdirs[["1kb"]])), error = function(e) NA)),
  TSS2kbRDataMD5 = as.character(tryCatch(tools::md5sum(file.path(cfg$rdata_base, cfg$tss_subdirs[["2kb"]])), error = function(e) NA)),
  TSS5kbRDataMD5 = as.character(tryCatch(tools::md5sum(file.path(cfg$rdata_base, cfg$tss_subdirs[["5kb"]])), error = function(e) NA)),
  TSS10kbRDataMD5 = as.character(tryCatch(tools::md5sum(file.path(cfg$rdata_base, cfg$tss_subdirs[["10kb"]])), error = function(e) NA)),
  LooplookVersion = as.character(tryCatch(packageVersion("looplook"), error = function(e) "unknown")),
  ChIPseekerVersion = as.character(tryCatch(packageVersion("ChIPseeker"), error = function(e) "unknown")),
  ClusterProfilerVersion = as.character(tryCatch(packageVersion("clusterProfiler"), error = function(e) "unknown")),
  GenomeBuild = "hg38",
  TxDbVersion = as.character(tryCatch(packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene"), error = function(e) "unknown")),
  OrgDbVersion = as.character(tryCatch(packageVersion("org.Hs.eg.db"), error = function(e) "unknown")),
  MetadataFile = normalizePath(meta_path, mustWork = FALSE),
  MetadataMD5 = as.character(tryCatch(tools::md5sum(meta_path), error = function(e) NA)),
  RVersion = R.version.string,
  BioconductorVersion = as.character(tryCatch(BiocManager::version(), error = function(e) "unknown")),
  GenePanelFile = if (exists("gene_panel_file")) gene_panel_file else NA_character_,
  GenePanelMD5 = if (exists("gene_panel_md5")) gene_panel_md5 else NA_character_,
  RunTime = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(base_out_dir, "analysis_manifest.csv"), row.names = FALSE)
message("\n>>> Analysis manifest saved")

# ── Run signature (cache validation) ─────────────────────────────────────────
run_signature <- list(
  CodeVersion = code_version,
  ScriptMD5 = script_md5,
  DiffMD5 = manifest$DiffMD5,
  ExprMD5 = manifest$ExprMD5,
  MetadataMD5 = manifest$MetadataMD5,
  TargetMD5 = manifest$TargetMD5,
  MainRDataMD5 = manifest$MainRDataMD5,
  RawRankedListMD5 = md5_object(glist),
  LFCpvRankedListMD5 = if (use_glist_pv) md5_object(glist_pv) else NA_character_,
  ModeGenesMD5 = md5_object(mode_genes),
  ChIPReferenceMD5 = md5_object(chipseeker_genes),
  ModeDefinitionMD5 = md5_object(lapply(modes, function(x) { x$ann <- NULL; x })),
  Iterations = n_iterations,
  SampleSizes = sample_sizes,
  PrimarySampleSize = primary_sample_size,
  GSEAMinSize = gsea_min_size,
  GSEAMaxSize = gsea_max_size,
  LooplookVersion = manifest$LooplookVersion,
  ChIPseekerVersion = manifest$ChIPseekerVersion,
  ClusterProfilerVersion = manifest$ClusterProfilerVersion)
saveRDS(run_signature, file.path(base_out_dir, "run_signature.rds"))
message("    Run signature saved")

# ── Finalize analysis status (all modules, all sample sizes) ─────────────────
message("\n>>> Finalizing analysis status...")
if (length(gsea_status_log) > 0L) {
  gsea_status_df <- bind_rows(gsea_status_log)
  write.csv(gsea_status_df, file.path(base_out_dir, "Analysis_Run_Status.csv"), row.names = FALSE)

  # GSEA-specific summary (exclude GO) - split by iteration vs mode level
  gsea_only <- gsea_status_df %>% filter(Module != "go")
  if (nrow(gsea_only) > 0L) {
    # Iteration-level summary (only rows with actual Iteration)
    gsea_iter <- gsea_only %>% filter(!is.na(Iteration))
    if (nrow(gsea_iter) > 0L) {
      iter_summary <- gsea_iter %>%
        group_by(SampleSize, Module, Mode) %>%
        summarise(N_Iterations = n(),
          N_FreshSuccess = sum(Status == "success"),
          N_CachedValidated = sum(Status == "cached_validated"),
          N_Valid = sum(Status %in% c("success", "cached_validated")),
          N_Partial = sum(Status == "partial_terms"),
          N_NoResult = sum(Status == "no_result"),
          N_Error = sum(Status == "error"),
          ValidRate = N_Valid / N_Iterations, .groups = "drop")
      write.csv(iter_summary, file.path(base_out_dir, "GSEA_Iteration_Summary.csv"), row.names = FALSE)
      if (any(iter_summary$ValidRate < 0.95, na.rm = TRUE)) {
        warning("At least one GSEA mode has valid iteration rate <95%.")
      }
    }
    # Mode-level summary (rows with NA Iteration = skips)
    gsea_mode <- gsea_only %>% filter(is.na(Iteration))
    if (nrow(gsea_mode) > 0L) {
      mode_summary <- gsea_mode %>%
        group_by(SampleSize, Module, Mode, Status) %>%
        summarise(N = n(), .groups = "drop")
      write.csv(mode_summary, file.path(base_out_dir, "GSEA_Mode_Status.csv"), row.names = FALSE)
    }
  }

  # GO-specific summary
  go_only <- gsea_status_df %>% filter(Module == "go")
  if (nrow(go_only) > 0L) {
    write.csv(go_only, file.path(base_out_dir, "GO_Call_Status.csv"), row.names = FALSE)
    go_status_summary <- go_only %>%
      group_by(Mode, GeneSet, Ontology) %>%
      summarise(N_Call = n(), N_Success = sum(Status == "success"),
        N_NoResult = sum(Status == "no_result"), N_Error = sum(Status == "error"),
        N_InsufficientInput = sum(Status == "insufficient_input"),
        .groups = "drop")
    write.csv(go_status_summary, file.path(base_out_dir, "GO_Status_Summary.csv"), row.names = FALSE)
  }

  # Detailed status breakdown
  n_success <- sum(gsea_status_df$Status == "success", na.rm = TRUE)
  n_cached <- sum(gsea_status_df$Status == "cached_validated", na.rm = TRUE)
  n_valid <- n_success + n_cached
  n_partial <- sum(gsea_status_df$Status == "partial_terms", na.rm = TRUE)
  n_noresult <- sum(gsea_status_df$Status == "no_result", na.rm = TRUE)
  n_error <- sum(gsea_status_df$Status == "error", na.rm = TRUE)
  n_insufficient_pool <- sum(gsea_status_df$Status == "insufficient_pool", na.rm = TRUE)
  n_insufficient_terms <- sum(gsea_status_df$Status == "insufficient_terms", na.rm = TRUE)
  n_empty_mode <- sum(gsea_status_df$Status == "empty_mode_genes", na.rm = TRUE)
  n_insufficient_input <- sum(gsea_status_df$Status == "insufficient_input", na.rm = TRUE)
  n_total <- nrow(gsea_status_df)
  message(sprintf("    Analysis status: %d total records", n_total))
  message(sprintf("    fresh_success=%d, cached_validated=%d, valid_total=%d, partial=%d, no_result=%d, error=%d",
    n_success, n_cached, n_valid, n_partial, n_noresult, n_error))
  message(sprintf("    insufficient_pool=%d, insufficient_terms=%d, empty_mode=%d, insufficient_input=%d",
    n_insufficient_pool, n_insufficient_terms, n_empty_mode, n_insufficient_input))
  n_computed <- n_valid + n_partial + n_noresult + n_error
  message(sprintf("    Valid rate (fresh+cached / computed): %.1f%%", 100 * n_valid / max(1, n_computed)))
} else {
  warning("No analysis status records were generated.")
}

# ── Output completeness check ────────────────────────────────────────────────
message("\n>>> Checking output completeness...")

# Content-level CSV check (uses expected_terms to verify each iteration)
check_gsea_csv <- function(path, expected_mode, expected_iterations, expected_terms) {
  if (!file.exists(path) || file.size(path) == 0L) {
    return(data.frame(File = basename(path), Mode = expected_mode, Readable = FALSE,
      NRows = 0L, NIterations = 0L, NCompleteIterations = 0L,
      Status = "missing_or_empty", stringsAsFactors = FALSE))
  }
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x)) {
    return(data.frame(File = basename(path), Mode = expected_mode, Readable = FALSE,
      NRows = NA_integer_, NIterations = NA_integer_, NCompleteIterations = 0L,
      Status = "unreadable", stringsAsFactors = FALSE))
  }
  required_cols <- c("ID", "NES", "pvalue", "Iteration", "Mode", "Requested_N", "Actual_N", "Pool_N", "Sampling_fraction")
  if (!all(required_cols %in% colnames(x))) {
    return(data.frame(File = basename(path), Mode = expected_mode, Readable = TRUE,
      NRows = nrow(x), NIterations = NA_integer_, NCompleteIterations = 0L,
      Status = paste0("missing_columns:", paste(setdiff(required_cols, colnames(x)), collapse = ";")),
      stringsAsFactors = FALSE))
  }
  # Content quality checks
  pval_ok <- all(x$pvalue >= 0 & x$pvalue <= 1, na.rm = TRUE)
  actual_ok <- all(x$Actual_N <= x$Pool_N, na.rm = TRUE)
  iter_dup <- any(duplicated(x[, c("Iteration", "ID")]))
  quality_flags <- c(
    if (!pval_ok) "pvalue_out_of_range",
    if (!actual_ok) "Actual_N_exceeds_Pool_N",
    if (iter_dup) "duplicate_iteration_ID")
  if (length(quality_flags) > 0) {
    return(data.frame(File = basename(path), Mode = expected_mode, Readable = TRUE,
      NRows = nrow(x), NIterations = NA_integer_, NCompleteIterations = 0L,
      Status = paste(quality_flags, collapse = ";"), stringsAsFactors = FALSE))
  }
  # Check Mode consistency
  mode_matches <- all(as.character(x$Mode) == expected_mode)
  # Check each iteration has all expected terms
  iter_qc <- x %>%
    filter(is.finite(NES)) %>%
    distinct(Iteration, ID) %>%
    group_by(Iteration) %>%
    summarise(Complete = all(expected_terms %in% as.character(ID)), .groups = "drop")
  n_iter <- dplyr::n_distinct(x$Iteration)
  n_complete <- sum(iter_qc$Complete, na.rm = TRUE)
  status <- dplyr::case_when(
    !mode_matches ~ "mode_mismatch",
    n_iter < expected_iterations ~ "incomplete_iterations",
    n_complete < expected_iterations ~ "partial_terms",
    TRUE ~ "ok")
  data.frame(File = basename(path), Mode = expected_mode, Readable = TRUE,
    NRows = nrow(x), NIterations = n_iter, NCompleteIterations = n_complete,
    Status = status, stringsAsFactors = FALSE)
}

expected_files <- c(
  "analysis_manifest.csv",
  "run_signature.rds",
  "Analysis_Run_Status.csv",
  "GSEA_Iteration_Summary.csv",
  "Peak_Input_Concordance.csv",
  "Mode_Object_Schema_QC.csv",
  file.path("expression_matched_es", "ExpressionMatched_ES.xlsx"),
  file.path("mode_ranking", "Mode_Composite_Ranking_LFC.xlsx"),
  file.path("evidence_composition", "Evidence_Composition.xlsx"),
  file.path("expanded_only", "Expanded_Only_Set_Catalog.csv"),
  file.path("expanded_only", "ExpandedOnly_GSEA_SamplingPlan.csv"),
  file.path("expanded_only", "ExpandedOnly_ExpressionMatched_EffectSize.csv"),
  file.path("expanded_only", "ExpandedOnly_Topology_Support.csv"),
  file.path("expanded_only", "ExpandedOnly_Topology_Sensitivity.csv"),
  file.path("expanded_only", "ExpandedOnly_Module_Status.txt"),
  file.path("expanded_only", "ExpandedOnly_Analysis.xlsx"),
  file.path("expanded_only", "ExpandedOnly_CrossTSS.xlsx"))
if (exists("gsea_status_df") && nrow(gsea_status_df) > 0L &&
    any(is.na(gsea_status_df$Iteration) & gsea_status_df$Module != "go")) {
  expected_files <- c(expected_files, "GSEA_Mode_Status.csv")
}
if (dir.exists(file.path(base_out_dir, "distance_analysis")))
  expected_files <- c(expected_files, file.path("distance_analysis", "Peak_Gene_Distance_Analysis.xlsx"))

# Check expected files and ensure they were created/refreshed in this run.
expected_paths <- file.path(base_out_dir, expected_files)
missing <- expected_files[!file.exists(expected_paths)]
existing_idx <- file.exists(expected_paths)
stale <- character(0)
if (any(existing_idx) && force_recompute) {
  existing_info <- file.info(expected_paths[existing_idx])
  stale <- expected_files[existing_idx][
    is.na(existing_info$mtime) | existing_info$mtime < analysis_run_start
  ]
} else if (any(existing_idx) && !force_recompute) {
  message("    Resume mode: validated pre-existing outputs are accepted without mtime refresh.")
}
if (length(missing) > 0L) {
  warning("Missing expected outputs: ", paste(missing, collapse = "; "))
}
if (length(stale) > 0L) {
  warning("Stale outputs not refreshed in this run: ", paste(stale, collapse = "; "))
}
if (length(missing) == 0L && length(stale) == 0L) {
  message("    All expected output files are current for this run.")
}

# Content-level check: ALL mode CSVs, not just first 5
expected_mode_names <- names(mode_genes)[lengths(mode_genes) > 0L]
all_completeness <- list()
for (sz in sample_sizes) {
  sim_dir <- file.path(base_out_dir, paste0("size", sz), "tmp_sim")
  if (!dir.exists(sim_dir)) next
  csv_files <- list.files(sim_dir, pattern = "^sim_.*\\.csv$", full.names = TRUE)
  expected_csv <- paste0("sim_", expected_mode_names, ".csv")
  actual_csv <- basename(csv_files)
  missing_csv <- setdiff(expected_csv, actual_csv)
  # Check ALL CSVs with mode-specific expected terms
  completeness <- bind_rows(lapply(csv_files, function(f) {
    mode_name <- sub("^sim_|\\.csv$", "", basename(f))
    # Control modes use ChIPseeker, others use looplook
    exp_terms <- if (grepl("^ctrl_", mode_name)) c("ChIPseeker", "Background") else c("looplook", "Background")
    check_gsea_csv(f, mode_name, n_iterations, exp_terms)
  }))
  # Add missing CSVs to report
  if (length(missing_csv) > 0) {
    missing_rows <- data.frame(
      File = missing_csv, Mode = sub("^sim_|\\.csv$", "", missing_csv),
      Readable = FALSE, NRows = 0L, NIterations = 0L, NCompleteIterations = 0L,
      Status = "missing", stringsAsFactors = FALSE)
    completeness <- bind_rows(completeness, missing_rows)
  }
  completeness$SampleSize <- sz
  all_completeness[[length(all_completeness) + 1]] <- completeness

  bad <- completeness %>% filter(Status != "ok")
  if (nrow(bad) > 0) {
    warning(sprintf("size%d: %d/%d CSVs have issues: %s", sz, nrow(bad), nrow(completeness),
      paste(bad$File, bad$Status, sep = "=", collapse = "; ")))
  }
  message(sprintf("    size%d: %d/%d mode CSVs present, %d complete",
    sz, length(csv_files), length(expected_mode_names), sum(completeness$Status == "ok")))
}

# Content-level check: additional output directories (non-main GSEA)
message(">>> Checking additional output directories...")
aux_dirs <- list(
  "unique" = "tmp_unique"
)
for (label in names(aux_dirs)) {
  aux_pattern <- if (label == "unique") "uniq_" else if (label == "distal") "distal_uniq_" else paste0(label, "_uniq_")
  for (sz in sample_sizes) {
    aux_path <- file.path(base_out_dir, paste0("size", sz), aux_dirs[[label]])
    if (!dir.exists(aux_path)) {
      warning(sprintf("Missing %s directory: %s", label, aux_path))
      next
    }
    aux_csvs <- list.files(aux_path, pattern = paste0("^", aux_pattern, ".*\\.csv$"), full.names = TRUE)
    expected_aux_csvs <- paste0(aux_pattern, expected_mode_names, ".csv")
    missing_aux <- setdiff(expected_aux_csvs, basename(aux_csvs))
    if (length(missing_aux) > 0) {
      warning(sprintf("%s size%d: %d/%d CSVs missing: %s", label, sz, length(missing_aux),
        length(expected_aux_csvs), paste(missing_aux, collapse = "; ")))
    }
    message(sprintf("    %s size%d: %d/%d CSVs present", label, sz, length(aux_csvs), length(expected_aux_csvs)))
  }
}

# LFC_pv completeness
lfc_pv_base <- file.path(base_out_dir, "LFC_pv")
if (dir.exists(lfc_pv_base)) {
  for (sz in sample_sizes) {
    lfc_sim_dir <- file.path(lfc_pv_base, paste0("size", sz), "tmp_sim")
    lfc_uniq_dir <- file.path(lfc_pv_base, paste0("size", sz), "tmp_unique")
    for (sub in c("tmp_sim" = lfc_sim_dir, "tmp_unique" = lfc_uniq_dir)) {
      if (!dir.exists(sub)) {
        warning(sprintf("Missing LFC_pv directory: %s", sub))
        next
      }
      lfc_csvs <- list.files(sub, pattern = "\\.csv$", full.names = TRUE)
      lfc_prefix <- if (identical(basename(sub), "tmp_sim")) "sim_" else "uniq_"
      expected_lfc <- paste0(lfc_prefix, expected_mode_names, ".csv")
      missing_lfc <- setdiff(expected_lfc, basename(lfc_csvs))
      if (length(missing_lfc) > 0L) {
        warning(sprintf("LFC_pv size%d %s: %d/%d CSVs present", sz, basename(sub),
          length(lfc_csvs), length(expected_mode_names)))
      }
    }
  }
}

# Plot files may be intentionally queued for independent rendering.
queued_plot_targets <- character(0)
if (identical(plot_policy, "queue")) {
  queued_meta_files <- list.files(plot_queue_dir, pattern = "\\.meta\\.rds$", full.names = TRUE)
  queued_plot_targets <- unique(unlist(lapply(queued_meta_files, function(f) {
    m <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.data.frame(m) && "Target" %in% colnames(m)) as.character(m$Target) else character(0)
  }), use.names = FALSE))
  queued_plot_targets <- normalizePath(queued_plot_targets, mustWork = FALSE)
}
output_exists_or_queued <- function(path) {
  path_norm <- normalizePath(path, mustWork = FALSE)
  file.exists(path) || (identical(plot_policy, "queue") && path_norm %in% queued_plot_targets)
}

# GO completeness
go_check_dir <- file.path(base_out_dir, "go_enrichment")
go_files <- c("go_results_cache.rds", sprintf("GO_Enrichment_%diter.xlsx", n_iterations),
  "GO_Only_vs_Shared_Scatter.pdf", "GO_Pipeline_Comparison.pdf")
if (exists("go_top") && nrow(go_top) > 0L) {
  go_files <- c(go_files, "GO_BP_TopTerms_Dotplot.pdf")
}
go_missing <- go_files[!vapply(file.path(go_check_dir, go_files), output_exists_or_queued, logical(1))]
if (length(go_missing) > 0) {
  warning("Missing GO output files: ", paste(go_missing, collapse = "; "))
} else {
  message("    GO enrichment output: all ", length(go_files), " files present")
}

# Distance analysis completeness
dist_dir <- file.path(base_out_dir, "distance_analysis")
if (dir.exists(dist_dir)) {
  dist_files <- c("Peak_Gene_Distance_Analysis.xlsx",
    "Distance_Density_4Groups_SourceSplit.pdf",
    "Distance_Boxplot_byPipeline_SourceSplit.pdf")
  dist_missing <- dist_files[!vapply(file.path(dist_dir, dist_files), output_exists_or_queued, logical(1))]
  if (length(dist_missing) > 0) {
    warning("Missing distance analysis output: ", paste(dist_missing, collapse = "; "))
  }
}

# Write completeness report
if (length(all_completeness) > 0) {
  completeness_df <- bind_rows(all_completeness)
  write.csv(completeness_df, file.path(base_out_dir, "Output_Completeness_Report.csv"), row.names = FALSE)
  if (any(completeness_df$Status != "ok")) {
    warning("Output completeness validation found issues. See Output_Completeness_Report.csv")
  }
}

# ── Completion marker ────────────────────────────────────────────────────────
has_error_status <- exists("gsea_status_df") && any(gsea_status_df$Status == "error", na.rm = TRUE)
has_partial <- exists("gsea_status_df") && any(gsea_status_df$Status == "partial_terms", na.rm = TRUE)
has_no_result <- exists("gsea_status_df") && any(gsea_status_df$Status == "no_result", na.rm = TRUE)
has_skipped <- exists("gsea_status_df") && any(gsea_status_df$Status %in%
  c("insufficient_pool", "insufficient_terms", "empty_mode_genes", "insufficient_input"), na.rm = TRUE)
has_completeness_issues <- exists("completeness_df") && any(completeness_df$Status != "ok", na.rm = TRUE)
has_missing <- (exists("missing") && length(missing) > 0) || (exists("stale") && length(stale) > 0)
has_go_missing <- exists("go_missing") && length(go_missing) > 0
has_dist_missing <- exists("dist_missing") && length(dist_missing) > 0
has_expanded_error <-
  (exists("expanded_module_error") && !is.null(expanded_module_error)) ||
  (exists("expanded_gsea_status_df") &&
     any(expanded_gsea_status_df$Status == "error", na.rm = TRUE)) ||
  (exists("expanded_go_status_df") &&
     any(expanded_go_status_df$Status == "error", na.rm = TRUE))
has_hard_failure <- has_error_status || has_completeness_issues ||
  has_missing || has_go_missing || has_dist_missing || has_expanded_error
completion_marker <- file.path(
  base_out_dir,
  if (identical(plot_policy, "queue")) "RUN_COMPUTE_COMPLETED.ok" else "RUN_COMPLETED.ok"
)
if (!has_hard_failure) {
  if (has_partial || has_no_result || has_skipped) {
    writeLines(
      c(paste0("Completed with minor issues: ", Sys.time()),
        paste0("ScriptMD5: ", script_md5),
        paste0("CodeVersion: ", code_version),
        paste0("PlotPolicy: ", plot_policy),
        paste0("HasPartialTerms: ", has_partial),
        paste0("HasNoResult: ", has_no_result),
        paste0("HasExpectedSkips: ", has_skipped)),
      completion_marker)
    message("\n>>> Compute stage complete (minor issues). Marker written: ", completion_marker)
  } else {
    writeLines(
      c(paste0("Completed: ", Sys.time()),
        paste0("ScriptMD5: ", script_md5),
        paste0("CodeVersion: ", code_version),
        paste0("PlotPolicy: ", plot_policy)),
      completion_marker)
    message("\n>>> Compute stage complete. Marker written: ", completion_marker)
  }
} else {
  writeLines(
    c(paste0("Finished with issues: ", Sys.time()),
      paste0("HardErrors: ", has_error_status),
      paste0("CompletenessIssues: ", has_completeness_issues),
      paste0("MissingOrStaleFiles: ", has_missing),
      paste0("StaleFiles: ", if (exists("stale")) paste(stale, collapse = ";") else ""),
      paste0("ExpandedOnlyErrors: ", has_expanded_error)),
    file.path(base_out_dir, "RUN_FINISHED_WITH_ISSUES.txt"))
  warning("\n>>> Analysis finished with issues. See RUN_FINISHED_WITH_ISSUES.txt")
}

}
