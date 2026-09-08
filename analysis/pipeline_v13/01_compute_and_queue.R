#!/usr/bin/env Rscript
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1L]) else "01_compute_and_queue.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "00_pipeline_config.R"))
source(pipeline_cfg$plot_style_file)
source(file.path(script_dir, "plot_functions", "plot_style_helpers.R"))

Sys.setenv(
  LOOPLOOK_DATA_BASE = pipeline_cfg$data_base,
  LOOPLOOK_RDATA_BASE = pipeline_cfg$rdata_base,
  LOOPLOOK_OUT_BASE = pipeline_cfg$out_base,
  LOOPLOOK_CACHE_POLICY = pipeline_cfg$cache_policy,
  LOOPLOOK_PLOT_POLICY = pipeline_cfg$plot_policy
)

if (!file.exists(pipeline_cfg$master_script)) {
  stop("Master script not found: ", pipeline_cfg$master_script)
}

message("=== looplook stage 01: compute + queue plots ===")
message("Output: ", pipeline_cfg$out_base)
message("Cache policy: ", pipeline_cfg$cache_policy)
message("Plot policy: ", pipeline_cfg$plot_policy)
source(pipeline_cfg$master_script, chdir = TRUE)
