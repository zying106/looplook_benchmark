# looplook modular benchmark configuration
# Edit only this file when moving the pipeline to a new server/project.

.required_env_path <- function(name) {
  value <- trimws(Sys.getenv(name, ""))
  if (!nzchar(value)) stop("Set environment variable ", name, " to an absolute path.", call. = FALSE)
  path.expand(value)
}
.r_lib <- trimws(Sys.getenv("LOOPLOOK_R_LIB", ""))
if (nzchar(.r_lib)) .libPaths(c(path.expand(.r_lib), .libPaths()))

.pipeline_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
    return(dirname(normalizePath(ofile, mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

pipeline_dir <- .pipeline_script_dir()

pipeline_cfg <- list(
  data_base = .required_env_path("LOOPLOOK_DATA_BASE"),
  rdata_base = .required_env_path("LOOPLOOK_RDATA_BASE"),
  out_base = .required_env_path("LOOPLOOK_OUT_BASE"),

  # resume: reuse validated caches/checkpoints; rebuild: true clean recomputation
  cache_policy = "resume",

  # queue is strongly recommended. immediate reproduces legacy in-process plotting.
  plot_policy = "queue",

  # Maximum time allowed for one isolated plot-render process.
  plot_timeout_sec = 900L,

  # Existing up-to-date targets are skipped unless TRUE.
  rerender_existing = FALSE,

  # Centralized figure-style control. Editing this does not invalidate analysis caches.
  plot_style_file = file.path(pipeline_dir, "00_plot_style.R"),

  # Export source data stored with each queued plot task for custom replotting.
  export_plot_data = TRUE
)

stopifnot(pipeline_cfg$cache_policy %in% c("resume", "rebuild"))
stopifnot(pipeline_cfg$plot_policy %in% c("queue", "immediate", "off"))

pipeline_cfg$master_script <- file.path(
  dirname(pipeline_dir),
  "looplook_benchmark_v12_plot_data_master.R"
)
pipeline_cfg$plot_queue_dir <- file.path(pipeline_cfg$out_base, "plot_queue")
pipeline_cfg$render_status_dir <- file.path(pipeline_cfg$plot_queue_dir, "render_status")
pipeline_cfg$plot_data_dir <- file.path(pipeline_cfg$out_base, "plot_data")

invisible(pipeline_cfg)
