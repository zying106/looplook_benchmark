# Shared, portable path configuration for manuscript figure scripts.

.looplook_script_dir <- function() {
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

looplook_env_path <- function(name, required = TRUE, default = NULL) {
  value <- trimws(Sys.getenv(name, ""))
  if (!nzchar(value) && !is.null(default)) value <- default
  if (!nzchar(value) && isTRUE(required)) {
    stop("Set environment variable ", name, " to an absolute path.", call. = FALSE)
  }
  if (!nzchar(value)) return(NULL)
  path.expand(value)
}

configure_looplook_library <- function() {
  lib <- looplook_env_path("LOOPLOOK_R_LIB", required = FALSE)
  if (!is.null(lib)) .libPaths(c(lib, .libPaths()))
  invisible(.libPaths())
}
