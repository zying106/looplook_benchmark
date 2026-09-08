#!/usr/bin/env Rscript
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1L]) else "08_export_plot_data.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "00_pipeline_config.R"))
`%||%` <- function(x, y) if (is.null(x)) y else x

dir.create(pipeline_cfg$plot_data_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(pipeline_cfg$plot_queue_dir, pattern = "\\.plot\\.rds$", full.names = TRUE)
if (!length(files)) stop("No plot tasks found in: ", pipeline_cfg$plot_queue_dir)

manifest <- list()
for (f in files) {
  task <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(task) || is.null(task$TaskID)) next
  dat <- NULL
  if (!is.null(task$PlotDataRDS) && file.exists(task$PlotDataRDS)) {
    dat <- tryCatch(readRDS(task$PlotDataRDS), error = function(e) NULL)
  }
  if (is.null(dat)) {
    dat <- list(plot_data = task$Plot$data, layer_data = lapply(task$Plot$layers, function(z) z$data))
  }
  out_rds <- file.path(pipeline_cfg$plot_data_dir, paste0(task$TaskID, ".plot_data.rds"))
  saveRDS(dat, out_rds, version = 3)

  csv_files <- character()
  pieces <- c(list(plot_data = dat$plot_data), dat$layer_data)
  for (nm in names(pieces)) {
    x <- pieces[[nm]]
    if (is.data.frame(x) && nrow(x) > 0L) {
      out_csv <- file.path(pipeline_cfg$plot_data_dir, paste0(task$TaskID, "__", nm, ".csv"))
      try(utils::write.csv(x, out_csv, row.names = FALSE), silent = TRUE)
      if (file.exists(out_csv)) csv_files <- c(csv_files, out_csv)
    }
  }
  manifest[[length(manifest) + 1L]] <- data.frame(
    TaskID = task$TaskID,
    Group = task$Group %||% NA_character_,
    StyleFamily = task$StyleFamily %||% NA_character_,
    Target = task$Target,
    PlotDataRDS = out_rds,
    CSVCount = length(csv_files),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, manifest)
utils::write.csv(out, file.path(pipeline_cfg$plot_data_dir, "plot_data_manifest.csv"), row.names = FALSE)
message("Exported plot-ready data for ", nrow(out), " task(s) to: ", pipeline_cfg$plot_data_dir)
