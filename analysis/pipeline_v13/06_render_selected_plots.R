#!/usr/bin/env Rscript
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1L]) else "06_render_selected_plots.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "00_pipeline_config.R"))
source(file.path(script_dir, "plot_queue_runner.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(paste(
    "Usage:",
    "Rscript 06_render_selected_plots.R <group|all> [filename_regex]",
    "Groups: core_gsea, functional_distance, expanded_evidence, summary_lfcpv"
  ))
}
group <- args[1L]
if (identical(group, "all")) {
  group <- c("core_gsea", "functional_distance", "expanded_evidence", "summary_lfcpv")
}
regex <- if (length(args) >= 2L) args[2L] else NULL
render_plot_group(group, filename_regex = regex)
