# Common plot-queue runner. Source this file from a group entry-point script.

.list_plot_meta <- function(queue_dir) {
  if (!dir.exists(queue_dir)) {
    stop("Plot queue directory does not exist: ", queue_dir)
  }
  files <- list.files(queue_dir, pattern = "\\.meta\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.data.frame(x) || nrow(x) != 1L) return(NULL)
    x
  })
  rows <- rows[vapply(rows, is.data.frame, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    miss <- setdiff(all_names, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$TaskID, fromLast = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}


.render_dependency_mtime <- function(task_path) {
  deps <- c(
    task_path,
    pipeline_cfg$plot_style_file,
    file.path(pipeline_dir, "plot_functions", "plot_style_helpers.R"),
    file.path(pipeline_dir, "99_render_one_plot.R")
  )
  deps <- deps[file.exists(deps)]
  if (!length(deps)) return(as.POSIXct(NA))
  max(file.info(deps)$mtime, na.rm = TRUE)
}

.write_marker <- function(path, text) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, path)
}

render_plot_group <- function(
    group,
    filename_regex = NULL,
    timeout_sec = pipeline_cfg$plot_timeout_sec,
    rerender_existing = pipeline_cfg$rerender_existing) {
  queue_dir <- pipeline_cfg$plot_queue_dir
  status_dir <- pipeline_cfg$render_status_dir
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- .list_plot_meta(queue_dir)
  if (nrow(meta) == 0L) {
    stop("No queued plot tasks found in: ", queue_dir)
  }
  meta <- meta[meta$Group %in% group, , drop = FALSE]
  if (!is.null(filename_regex) && nzchar(filename_regex)) {
    meta <- meta[grepl(filename_regex, basename(meta$Target), ignore.case = TRUE), , drop = FALSE]
  }
  if (nrow(meta) == 0L) {
    message("No tasks matched group/regex.")
    return(invisible(data.frame()))
  }
  meta <- meta[order(meta$Target), , drop = FALSE]

  renderer <- file.path(pipeline_dir, "99_render_one_plot.R")
  if (!file.exists(renderer)) stop("Renderer not found: ", renderer)
  rscript <- file.path(R.home("bin"), "Rscript")
  timeout_bin <- Sys.which("timeout")
  if (!nzchar(timeout_bin)) {
    warning("GNU timeout was not found. Plot processes remain isolated, but cannot be force-killed by time limit.")
  }

  results <- vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    task_id <- as.character(row$TaskID)
    task_path <- as.character(row$TaskRDS)
    target <- as.character(row$Target)
    log_path <- file.path(status_dir, paste0(task_id, ".log"))
    ok_path <- file.path(status_dir, paste0(task_id, ".ok"))
    fail_path <- file.path(status_dir, paste0(task_id, ".failed"))
    timeout_path <- file.path(status_dir, paste0(task_id, ".timeout"))
    unlink(c(fail_path, timeout_path))

    dep_mtime <- tryCatch(.render_dependency_mtime(task_path), error = function(e) as.POSIXct(NA))
    target_mtime <- tryCatch(file.info(target)$mtime, error = function(e) as.POSIXct(NA))
    target_current <- file.exists(target) && !is.na(target_mtime) && !is.na(dep_mtime) && target_mtime >= dep_mtime

    if (target_current && !isTRUE(rerender_existing)) {
      status <- "skipped_current"
      message(sprintf("[%d/%d] SKIP current: %s", i, nrow(meta), basename(target)))
      .write_marker(ok_path, paste("Current target reused:", target))
      results[[i]] <- data.frame(TaskID = task_id, Group = row$Group,
        Target = target, Status = status, ExitCode = 0L,
        Log = log_path, stringsAsFactors = FALSE)
      next
    }

    message(sprintf("[%d/%d] Rendering: %s", i, nrow(meta), basename(target)))
    start <- Sys.time()
    if (nzchar(timeout_bin)) {
      cmd_args <- c(
        "--signal=TERM",
        paste0(as.integer(timeout_sec), "s"),
        shQuote(rscript),
        shQuote(renderer),
        shQuote(task_path)
      )
      exit_code <- system2(timeout_bin, cmd_args, stdout = log_path, stderr = log_path)
    } else {
      exit_code <- system2(rscript, c(shQuote(renderer), shQuote(task_path)),
        stdout = log_path, stderr = log_path)
    }
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    if (identical(as.integer(exit_code), 0L) && file.exists(target)) {
      status <- "success"
      .write_marker(ok_path, c(
        paste("Rendered:", target),
        paste("ElapsedSeconds:", round(elapsed, 2))
      ))
      unlink(c(fail_path, timeout_path))
    } else if (identical(as.integer(exit_code), 124L)) {
      status <- "timeout"
      .write_marker(timeout_path, c(
        paste("Timed out:", target),
        paste("TimeoutSeconds:", timeout_sec),
        paste("Log:", log_path)
      ))
      warning("Plot timed out; continuing to next task: ", basename(target), call. = FALSE)
    } else {
      status <- "failed"
      .write_marker(fail_path, c(
        paste("Failed:", target),
        paste("ExitCode:", exit_code),
        paste("Log:", log_path)
      ))
      warning("Plot failed; continuing to next task: ", basename(target), call. = FALSE)
    }

    results[[i]] <- data.frame(
      TaskID = task_id,
      Group = as.character(row$Group),
      Target = target,
      Status = status,
      ExitCode = as.integer(exit_code),
      ElapsedSeconds = elapsed,
      Log = log_path,
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  summary_path <- file.path(
    status_dir,
    paste0("render_summary_", paste(group, collapse = "_"), ".csv")
  )
  utils::write.csv(result_df, summary_path, row.names = FALSE)
  message("Render summary: ", summary_path)
  message("Success/current: ", sum(result_df$Status %in% c("success", "skipped_current")),
    "; failed: ", sum(result_df$Status == "failed"),
    "; timeout: ", sum(result_df$Status == "timeout"))
  invisible(result_df)
}

finalize_plot_queue <- function() {
  meta <- .list_plot_meta(pipeline_cfg$plot_queue_dir)
  if (nrow(meta) == 0L) stop("No plot tasks found.")
  current <- vapply(seq_len(nrow(meta)), function(i) {
    task <- as.character(meta$TaskRDS[i])
    target <- as.character(meta$Target[i])
    if (!file.exists(task) || !file.exists(target)) return(FALSE)
    dep_mtime <- .render_dependency_mtime(task)
    file.info(target)$mtime >= dep_mtime && file.info(target)$size >= 100L
  }, logical(1))
  status <- data.frame(meta, RenderedCurrent = current, stringsAsFactors = FALSE)
  out <- file.path(pipeline_cfg$plot_queue_dir, "plot_render_completeness.csv")
  utils::write.csv(status, out, row.names = FALSE)
  marker <- file.path(pipeline_cfg$out_base, "RUN_PLOTS_COMPLETED.ok")
  if (all(current)) {
    writeLines(c(
      paste0("All queued plots rendered: ", Sys.time()),
      paste0("NPlots: ", nrow(status))
    ), marker)
    message("All queued plots are current. Marker: ", marker)
  } else {
    unlink(marker)
    warning(sum(!current), " plot task(s) remain missing/stale. See: ", out, call. = FALSE)
  }
  invisible(status)
}
