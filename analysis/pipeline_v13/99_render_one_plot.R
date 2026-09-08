#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 99_render_one_plot.R <task.plot.rds>")
}

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1L]) else "99_render_one_plot.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "00_pipeline_config.R"))
source(pipeline_cfg$plot_style_file)
source(file.path(script_dir, "plot_functions", "plot_style_helpers.R"))

task_path <- normalizePath(args[1L], mustWork = TRUE)
task <- readRDS(task_path)
required_fields <- c("Target", "Plot", "GgsaveArgs")
missing_fields <- setdiff(required_fields, names(task))
if (length(missing_fields) > 0L) {
  stop("Invalid plot task; missing: ", paste(missing_fields, collapse = ", "))
}

required_packages <- unique(c("ggplot2", as.character(task$RequiredPackages)))
required_packages <- setdiff(required_packages, c("grid"))
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required plotting package is not installed: ", pkg)
  }
}

if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")

custom <- customize_queued_plot(task, style = plot_style)
p <- custom$plot
save_args_extra <- custom$ggsave_args

target <- normalizePath(as.character(task$Target), mustWork = FALSE)
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
ext <- tools::file_ext(target)
base_no_ext <- if (nzchar(ext)) tools::file_path_sans_ext(basename(target)) else basename(target)
tmp_name <- if (nzchar(ext)) {
  paste0(".", base_no_ext, ".rendering.", Sys.getpid(), ".", ext)
} else {
  paste0(".", base_no_ext, ".rendering.", Sys.getpid())
}
tmp_target <- file.path(dirname(target), tmp_name)
on.exit(unlink(tmp_target), add = TRUE)

save_args <- c(list(filename = tmp_target, plot = p), save_args_extra)
do.call(ggplot2::ggsave, save_args)

if (!file.exists(tmp_target) || is.na(file.info(tmp_target)$size) || file.info(tmp_target)$size < 100L) {
  stop("Rendered output is missing or too small: ", tmp_target)
}

if (file.exists(target)) unlink(target)
if (!file.rename(tmp_target, target)) {
  stop("Unable to atomically move rendered file to: ", target)
}

message("Rendered: ", target)
message("Style family: ", task$StyleFamily %||% infer_plot_family(target))
