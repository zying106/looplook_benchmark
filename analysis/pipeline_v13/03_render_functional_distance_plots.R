#!/usr/bin/env Rscript
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1L]) else "03_render_functional_distance_plots.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "00_pipeline_config.R"))
source(file.path(script_dir, "plot_queue_runner.R"))
message("=== looplook plot stage: functional_distance ===")
render_plot_group("functional_distance")
