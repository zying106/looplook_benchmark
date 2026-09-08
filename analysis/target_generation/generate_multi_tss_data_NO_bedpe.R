#!/usr/bin/env Rscript

# ==============================================================================
# looplook BRD4 multi-TSS-window data-generation workflow (R39-compatible)
#
# Purpose
# -------
# Generate the eight benchmark objects for each promoter/TSS window:
#
#   res, res2,                        (annotation - uses expression for conflict resolution)
#   refined_res, refined_res2,        (expression refinement)
#   cr, cr2,                          (expression → chromatin refinement)
#   cr_only, cr2_only                 (annotation → chromatin refinement, no expression reclassification)
#
# Improvements over the four duplicated scripts
# ----------------------------------------------
# 1. The consensus BEDPE is generated once and reused by every TSS window.
# 2. Only `tss_region` changes among 1/2/5/10 kb analyses.
# 3. The same `tss_region` is passed to chromatin refinement.
# 4. Redundant chromatin validation is removed from expression refinement.
# 5. R39 defaults that affect interpretation are explicitly fixed.
# 6. Input MD5s, package versions and parameters are saved with every RData.
# 7. Every plot stored in the eight benchmark objects is exported automatically
#    as a layered PNG/PDF gallery with TSV and HTML indexes.
#
# Default output files remain compatible with the downstream benchmark:
#
#   v0.14/tss1000/tss_1000_res_chromatin.RData
#   v0.14/tss2000/tss_2000_res_chromatin.RData
#   v0.14/tss5000/tss_5000_res_chromatin.RData
#   v0.14/tss10000/tss_10000_res_chromatin.RData
#
# Run:
#   Rscript clean_brd4_r21_generate_multi_tss_R39_visual_export.R
#
# Optional subset:
#   Rscript clean_brd4_r21_generate_multi_tss_R39_visual_export.R --tss=1000,5000
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Library path and dependencies
# ------------------------------------------------------------------------------

.required_env_path <- function(name) {
  value <- trimws(Sys.getenv(name, ""))
  if (!nzchar(value)) stop("Set environment variable ", name, " to an absolute path.", call. = FALSE)
  path.expand(value)
}
.r_lib <- trimws(Sys.getenv("LOOPLOOK_R_LIB", ""))
if (nzchar(.r_lib)) .libPaths(c(path.expand(.r_lib), .libPaths()))

required_packages <- c(
  "looplook",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(looplook)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})


# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

cfg <- list(
  version_tag = "v0.15",
  project_name = "141_brd4_arv825",
  
  # Input/output roots
  data_dir = .required_env_path("LOOPLOOK_TARGET_DATA_DIR"),
  out_dir = .required_env_path("LOOPLOOK_TARGET_OUT_DIR"),
  
  # Replicate loop files
  loop_files = c(
    "GSM6836729_LPS141_1.filt.intra.loop_counts_fixed.bedpe",
    "GSM6836730_LPS141_2.filt.intra.loop_counts_fixed.bedpe"
  ),
  
  # Main inputs
  expr_file = "TPM_clean.txt",
  target_file = "lps141_folsl2_peaks.narrowPeak",
  
  # Control samples used both for initial annotation conflict resolution
  # and formal expression refinement.
  sample_columns = c(
    "lps141_dmso_1",
    "lps141_dmso_2"
  ),
  
  # TSS-window sensitivity analysis
  tss_windows = c(
    1000L,
    2000L,
    5000L,
    10000L
  ),
  
  # Consensus-loop parameters
  loop_gap = 500L,
  min_consensus = 2L,
  score_col = NULL,
  min_raw_score = 2,
  min_score = 3,
  
  # "warn" preserves the behavior of the existing scripts.
  # For an initial strict QC run, change to "error".
  chaining_policy = "warn",
  
  blacklist_species = "hg38",
  roi_file = "141_H3K27ac_peaks.narrowPeak",
  roi_mode = "both",
  
  # Annotation parameters
  species = "hg38",
  anchor_merge_gap = 0L,
  annotation_anchor_gap = -1L,
  annotation_anchor_min_overlap = 1L,
  annotation_anchor_min_frac = 0,
  min_expr = 0,
  conflict_strategy = "biotype_first",
  co_dominance_ratio = 0.1,
  hub_percentile = 0.95,
  hub_metric = "unique_contacts",
  
  # Expression-refinement parameters
  expression_threshold = 0.5,
  expression_threshold_mode = "absolute",
  expression_unit = "TPM",
  reclassify_by_expression = TRUE,
  
  # Chromatin-refinement parameters
  chromatin_anchor_gap = 200L,
  chromatin_anchor_min_overlap = 100L,
  bw_ratio_threshold = 3,
  recompute_targets = TRUE,
  
  # Chromatin BED files (named list: key = signal type, value = filename relative to data_dir)
  # Used by refine_loop_anchors_by_chromatin() to compute peak-level proximity signals.
  chromatin_beds = list(
    H3K4me1  = "lps141_H3K4me1_peaks.narrowPeak",
    H3K27ac  = "141_H3K27ac_peaks.narrowPeak",
    H3K4me3  = "lps141_H3K4me3_peaks.narrowPeak",
    H3K27me3 = "1ps141_h3k27me3_broad_peaks.broadPeak",
    ATAC     = "AT_70.ATAC_peaks.narrowPeak"
  ),
  
  # Chromatin bigWig files (named list: key = signal type, value = filename relative to data_dir)
  # Used by refine_loop_anchors_by_chromatin() for continuous signal ratio calculation.
  chromatin_bw = list(
    H3K4me1 = "lps141_H3K4me1_subtract_nonegative.bw",
    H3K4me3 = "lps141_H3K4me3_subtract_nonegative.bw"
  ),
  
  # Automatic visualization export
  export_visualizations = TRUE,
  visualization_dir_name = "visualization_export",
  visualization_formats = c("pdf"),
  visualization_width = 8,
  visualization_height = 7,
  visualization_dpi = 300L,
  visualization_max_depth = 14L,
  overwrite_visualizations = TRUE,
  
  # Execution
  overwrite_consensus = FALSE,
  overwrite_tss_results = TRUE,
  write_output = TRUE,
  quiet = FALSE
)


# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------

log_message <- function(...) {
  message(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    ...
  )
}

normalize_existing_path <- function(path, label) {
  path <- path.expand(path)

  if (!file.exists(path)) {
    stop(
      label,
      " not found: ",
      path,
      call. = FALSE
    )
  }

  normalizePath(
    path,
    mustWork = TRUE
  )
}

safe_md5 <- function(path) {
  if (length(path) != 1L ||
      is.na(path) ||
      !file.exists(path)) {
    return(NA_character_)
  }

  unname(
    as.character(
      tools::md5sum(path)
    )
  )
}

get_script_path <- function() {
  args <- commandArgs(
    trailingOnly = FALSE
  )

  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )

  if (length(file_arg) > 0L) {
    return(
      normalizePath(
        sub("^--file=", "", file_arg[1L]),
        mustWork = FALSE
      )
    )
  }

  NA_character_
}

parse_tss_argument <- function(default_windows) {
  args <- commandArgs(
    trailingOnly = TRUE
  )

  tss_arg <- grep(
    "^--tss=",
    args,
    value = TRUE
  )

  if (length(tss_arg) == 0L) {
    return(default_windows)
  }

  if (length(tss_arg) > 1L) {
    stop(
      "Provide at most one --tss= argument.",
      call. = FALSE
    )
  }

  raw_values <- strsplit(
    sub("^--tss=", "", tss_arg[1L]),
    ",",
    fixed = TRUE
  )[[1L]]

  windows <- suppressWarnings(
    as.integer(
      trimws(raw_values)
    )
  )

  if (length(windows) == 0L ||
      anyNA(windows) ||
      any(windows <= 0L)) {
    stop(
      "`--tss` must contain positive integers, e.g. --tss=1000,5000.",
      call. = FALSE
    )
  }

  unique(windows)
}

validate_expression_samples <- function(expr_path, sample_columns) {
  expr_header <- tryCatch(
    names(
      utils::read.delim(
        expr_path,
        header = TRUE,
        nrows = 1L,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    ),
    error = function(e) {
      stop(
        "Unable to read expression-matrix header: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  missing_samples <- setdiff(
    sample_columns,
    expr_header
  )

  if (length(missing_samples) > 0L) {
    stop(
      "Expression matrix is missing sample column(s): ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

validate_tss_window <- function(x) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0 ||
      x != floor(x)) {
    stop(
      "Each TSS window must be one finite positive integer.",
      call. = FALSE
    )
  }

  as.integer(x)
}

save_named_objects <- function(objects, file) {
  if (!is.list(objects) ||
      is.null(names(objects)) ||
      any(names(objects) == "")) {
    stop(
      "`objects` must be a fully named list.",
      call. = FALSE
    )
  }

  save_env <- list2env(
    objects,
    parent = emptyenv()
  )

  save(
    list = names(objects),
    file = file,
    envir = save_env,
    version = 3
  )
}


# ------------------------------------------------------------------------------
# 2b. Automatic visualization export helpers
# ------------------------------------------------------------------------------

sanitize_file_component <- function(x, fallback = "plot") {
  x <- as.character(x)[1L]

  if (is.na(x) || !nzchar(x)) {
    x <- fallback
  }

  x_ascii <- suppressWarnings(
    iconv(
      x,
      from = "",
      to = "ASCII//TRANSLIT"
    )
  )

  if (is.na(x_ascii) || !nzchar(x_ascii)) {
    x_ascii <- fallback
  }

  x_ascii <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x_ascii
  )

  x_ascii <- gsub(
    "_+",
    "_",
    x_ascii
  )

  x_ascii <- gsub(
    "^[_\\.]+|[_\\.]+$",
    "",
    x_ascii
  )

  if (!nzchar(x_ascii)) {
    x_ascii <- fallback
  }

  substr(
    x_ascii,
    1L,
    100L
  )
}

is_exportable_plot <- function(x) {
  inherits(
    x,
    c(
      "ggplot",
      "patchwork",
      "grob",
      "gTree",
      "gtable",
      "recordedplot",
      "trellis",
      "htmlwidget"
    )
  )
}

plot_object_type <- function(x) {
  if (inherits(x, "htmlwidget")) {
    return("htmlwidget")
  }

  if (inherits(x, "recordedplot")) {
    return("recordedplot")
  }

  if (inherits(x, "trellis")) {
    return("trellis")
  }

  if (inherits(x, c("grob", "gTree", "gtable"))) {
    return("grid")
  }

  if (inherits(x, c("ggplot", "patchwork"))) {
    return("ggplot")
  }

  "unknown"
}

collect_plot_objects <- function(
    x,
    path = "root",
    depth = 0L,
    max_depth = 14L
) {
  if (is_exportable_plot(x)) {
    return(
      list(
        list(
          path = path,
          plot = x,
          type = plot_object_type(x)
        )
      )
    )
  }

  if (depth >= max_depth ||
      !is.list(x) ||
      is.data.frame(x)) {
    return(list())
  }

  # Do not recurse through common genomic containers even if a class happens
  # to expose list-like behavior.
  if (inherits(
      x,
      c(
        "GRanges",
        "GRangesList",
        "GInteractions",
        "SummarizedExperiment",
        "DataFrame"
      )
  )) {
    return(list())
  }

  child_names <- names(x)

  if (is.null(child_names)) {
    child_names <- rep(
      "",
      length(x)
    )
  }

  output <- list()

  for (i in seq_along(x)) {
    child <- x[[i]]

    # Atomic vectors, matrices and S4 genomic objects cannot contain the plot
    # objects produced by looplook and are skipped to avoid scanning large data.
    if (!is_exportable_plot(child) &&
        (!is.list(child) ||
         is.data.frame(child))) {
      next
    }

    child_name <- child_names[i]

    if (is.na(child_name) ||
        !nzchar(child_name)) {
      child_name <- paste0(
        "item_",
        i
      )
    }

    child_path <- paste(
      path,
      child_name,
      sep = "/"
    )

    output <- c(
      output,
      collect_plot_objects(
        child,
        path = child_path,
        depth = depth + 1L,
        max_depth = max_depth
      )
    )
  }

  output
}

draw_exportable_plot <- function(plot_object) {
  if (inherits(
      plot_object,
      c(
        "ggplot",
        "patchwork",
        "trellis"
      )
  )) {
    print(
      plot_object
    )
    return(
      invisible(TRUE)
    )
  }

  if (inherits(
      plot_object,
      c(
        "grob",
        "gTree",
        "gtable"
      )
  )) {
    grid::grid.newpage()
    grid::grid.draw(
      plot_object
    )
    return(
      invisible(TRUE)
    )
  }

  if (inherits(
      plot_object,
      "recordedplot"
  )) {
    grDevices::replayPlot(
      plot_object
    )
    return(
      invisible(TRUE)
    )
  }

  stop(
    "Unsupported plot class: ",
    paste(
      class(plot_object),
      collapse = ", "
    ),
    call. = FALSE
  )
}

save_plot_file <- function(
    plot_object,
    filename,
    format,
    width,
    height,
    dpi
) {
  format <- tolower(
    format
  )

  if (inherits(
      plot_object,
      "htmlwidget"
  )) {
    if (format != "html") {
      return(
        list(
          saved = FALSE,
          message = "HTML widget is exported only as HTML."
        )
      )
    }

    if (!requireNamespace(
        "htmlwidgets",
        quietly = TRUE
    )) {
      return(
        list(
          saved = FALSE,
          message = "Package 'htmlwidgets' is unavailable."
        )
      )
    }

    htmlwidgets::saveWidget(
      widget = plot_object,
      file = filename,
      selfcontained = FALSE
    )

    return(
      list(
        saved = TRUE,
        message = ""
      )
    )
  }

  if (format == "png") {
    grDevices::png(
      filename = filename,
      width = round(width * dpi),
      height = round(height * dpi),
      units = "px",
      res = dpi,
      bg = "white"
    )

    device_open <- TRUE

    tryCatch(
      {
        draw_exportable_plot(
          plot_object
        )
      },
      finally = {
        if (device_open) {
          grDevices::dev.off()
        }
      }
    )

    return(
      list(
        saved = TRUE,
        message = ""
      )
    )
  }

  if (format == "pdf") {
    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      onefile = FALSE,
      useDingbats = FALSE
    )

    device_open <- TRUE

    tryCatch(
      {
        draw_exportable_plot(
          plot_object
        )
      },
      finally = {
        if (device_open) {
          grDevices::dev.off()
        }
      }
    )

    return(
      list(
        saved = TRUE,
        message = ""
      )
    )
  }

  list(
    saved = FALSE,
    message = paste0(
      "Unsupported visualization format: ",
      format
    )
  )
}

relative_to_root <- function(path, root) {
  if (length(path) != 1L ||
      is.na(path) ||
      !nzchar(path)) {
    return(NA_character_)
  }

  root_norm <- normalizePath(
    root,
    winslash = "/",
    mustWork = FALSE
  )

  path_norm <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )

  root_prefix <- paste0(
    root_norm,
    "/"
  )

  if (startsWith(
      path_norm,
      root_prefix
  )) {
    return(
      substring(
        path_norm,
        nchar(root_prefix) + 1L
      )
    )
  }

  path_norm
}

html_escape <- function(x) {
  x <- as.character(x)

  x <- gsub(
    "&",
    "&amp;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "<",
    "&lt;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    ">",
    "&gt;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "\"",
    "&quot;",
    x,
    fixed = TRUE
  )

  x
}

write_visualization_gallery <- function(
    manifest,
    root_dir,
    title
) {
  index_path <- file.path(
    root_dir,
    "index.html"
  )

  successful <- manifest[
    manifest$Status %in%
      c(
        "exported",
        "partially_exported"
      ),
    ,
    drop = FALSE
  ]

  cards <- character()

  if (nrow(successful) == 0L) {
    cards <- paste0(
      "<p class='empty'>No exportable plot objects were found.</p>"
    )
  } else {
    for (i in seq_len(
        nrow(successful)
    )) {
      row <- successful[i, , drop = FALSE]

      png_rel <- relative_to_root(
        row$PNG,
        root_dir
      )

      pdf_rel <- relative_to_root(
        row$PDF,
        root_dir
      )

      html_rel <- relative_to_root(
        row$HTML,
        root_dir
      )

      preview <- if (!is.na(png_rel) &&
                     nzchar(png_rel) &&
                     file.exists(row$PNG)) {
        paste0(
          "<a href='",
          html_escape(png_rel),
          "'><img loading='lazy' src='",
          html_escape(png_rel),
          "' alt='",
          html_escape(row$Plot_Path),
          "'></a>"
        )
      } else {
        "<div class='no-preview'>No PNG preview</div>"
      }

      links <- character()

      if (!is.na(png_rel) &&
          nzchar(png_rel) &&
          file.exists(row$PNG)) {
        links <- c(
          links,
          paste0(
            "<a href='",
            html_escape(png_rel),
            "'>PNG</a>"
          )
        )
      }

      if (!is.na(pdf_rel) &&
          nzchar(pdf_rel) &&
          file.exists(row$PDF)) {
        links <- c(
          links,
          paste0(
            "<a href='",
            html_escape(pdf_rel),
            "'>PDF</a>"
          )
        )
      }

      if (!is.na(html_rel) &&
          nzchar(html_rel) &&
          file.exists(row$HTML)) {
        links <- c(
          links,
          paste0(
            "<a href='",
            html_escape(html_rel),
            "'>Interactive HTML</a>"
          )
        )
      }

      cards <- c(
        cards,
        paste0(
          "<article class='plot-card'>",
          preview,
          "<div class='meta'>",
          "<div class='stage'>",
          html_escape(row$Stage),
          " · ",
          html_escape(row$Hop),
          "</div>",
          "<h2>",
          html_escape(row$Plot_Path),
          "</h2>",
          "<div class='sub'>",
          html_escape(row$Object),
          " · ",
          html_escape(row$Plot_Type),
          "</div>",
          "<div class='links'>",
          paste(
            links,
            collapse = " · "
          ),
          "</div>",
          "</div>",
          "</article>"
        )
      )
    }
  }

  html <- c(
    "<!doctype html>",
    "<html lang='en'>",
    "<head>",
    "<meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'>",
    paste0(
      "<title>",
      html_escape(title),
      "</title>"
    ),
    "<style>",
    "body{font-family:Arial,sans-serif;margin:0;background:#f5f6f8;color:#1d2433}",
    "header{padding:24px 28px;background:#fff;border-bottom:1px solid #dfe3e8}",
    "h1{font-size:24px;margin:0 0 6px}",
    "header p{margin:0;color:#5f6877}",
    ".grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:18px;padding:22px}",
    ".plot-card{background:#fff;border:1px solid #dfe3e8;border-radius:12px;overflow:hidden}",
    ".plot-card img{width:100%;height:250px;object-fit:contain;background:#fff;display:block}",
    ".no-preview{height:250px;display:flex;align-items:center;justify-content:center;color:#777;background:#fafafa}",
    ".meta{padding:14px 16px 18px}",
    ".stage{font-size:12px;color:#56657a;text-transform:uppercase;letter-spacing:.04em}",
    "h2{font-size:16px;line-height:1.35;margin:7px 0;overflow-wrap:anywhere}",
    ".sub{font-size:13px;color:#687386}",
    ".links{margin-top:12px;font-size:14px}",
    "a{color:#2459a9;text-decoration:none}",
    "a:hover{text-decoration:underline}",
    ".empty{padding:30px}",
    "</style>",
    "</head>",
    "<body>",
    "<header>",
    paste0(
      "<h1>",
      html_escape(title),
      "</h1>"
    ),
    paste0(
      "<p>",
      nrow(successful),
      " plot object(s) exported. Open the PNG for quick review or PDF for vector output.</p>"
    ),
    "</header>",
    "<main class='grid'>",
    cards,
    "</main>",
    "</body>",
    "</html>"
  )

  writeLines(
    html,
    index_path,
    useBytes = TRUE
  )

  index_path
}

benchmark_visualization_layout <- function() {
  data.frame(
    Object = c(
      "res",
      "res2",
      "refined_res",
      "refined_res2",
      "cr",
      "cr2",
      "cr_only",
      "cr2_only"
    ),
    Stage = c(
      "01_annotation",
      "01_annotation",
      "02_expression",
      "02_expression",
       "03_expression_chromatin",
       "03_expression_chromatin",
       "04_chromatin_refinement_no_expr_reclass",
       "04_chromatin_refinement_no_expr_reclass"
    ),
    Hop = c(
      "hop0",
      "hop1",
      "hop0",
      "hop1",
      "hop0",
      "hop1",
      "hop0",
      "hop1"
    ),
    stringsAsFactors = FALSE
  )
}

export_benchmark_visualizations <- function(
    objects,
    root_dir,
    tss_window,
    cfg
) {
  empty_manifest <- data.frame(
    TSS_Window = integer(),
    Object = character(),
    Stage = character(),
    Hop = character(),
    Plot_Path = character(),
    Plot_Type = character(),
    PNG = character(),
    PDF = character(),
    HTML = character(),
    Status = character(),
    Message = character(),
    stringsAsFactors = FALSE
  )

  if (!isTRUE(
      cfg$export_visualizations
  )) {
    return(
      list(
        root = NA_character_,
        index = NA_character_,
        manifest_path = NA_character_,
        manifest = empty_manifest
      )
    )
  }

  if (isTRUE(
      cfg$overwrite_visualizations
  ) &&
      dir.exists(root_dir)) {
    unlink(
      root_dir,
      recursive = TRUE,
      force = TRUE
    )
  }

  dir.create(
    root_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  layout <- benchmark_visualization_layout()
  rows <- list()
  row_i <- 0L

  for (layout_i in seq_len(
      nrow(layout)
  )) {
    object_name <- layout$Object[layout_i]
    stage_name <- layout$Stage[layout_i]
    hop_name <- layout$Hop[layout_i]

    if (!object_name %in% names(objects) ||
        is.null(objects[[object_name]])) {
      next
    }

    object_dir <- file.path(
      root_dir,
      stage_name,
      hop_name
    )

    dir.create(
      object_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    plot_nodes <- collect_plot_objects(
      objects[[object_name]],
      path = object_name,
      max_depth = cfg$visualization_max_depth
    )

    if (length(plot_nodes) == 0L) {
      next
    }

    for (plot_i in seq_along(
        plot_nodes
    )) {
      node <- plot_nodes[[plot_i]]
      plot_object <- node$plot
      plot_type <- node$type

      plot_slug <- sanitize_file_component(
        node$path,
        fallback = paste0(
          "plot_",
          plot_i
        )
      )

      file_stem <- sprintf(
        "%03d_%s",
        plot_i,
        plot_slug
      )

      requested_formats <- unique(
        tolower(
          cfg$visualization_formats
        )
      )

      if (plot_type == "htmlwidget") {
        requested_formats <- "html"
      } else {
        requested_formats <- intersect(
          requested_formats,
          c(
            "png",
            "pdf"
          )
        )
      }

      output_paths <- c(
        PNG = NA_character_,
        PDF = NA_character_,
        HTML = NA_character_
      )

      messages <- character()
      saved_any <- FALSE
      failed_any <- FALSE

      for (format in requested_formats) {
        filename <- file.path(
          object_dir,
          paste0(
            file_stem,
            ".",
            format
          )
        )

        result <- tryCatch(
          save_plot_file(
            plot_object = plot_object,
            filename = filename,
            format = format,
            width = cfg$visualization_width,
            height = cfg$visualization_height,
            dpi = cfg$visualization_dpi
          ),
          error = function(e) {
            list(
              saved = FALSE,
              message = conditionMessage(e)
            )
          }
        )

        if (isTRUE(
            result$saved
        )) {
          output_paths[
            toupper(format)
          ] <- filename

          saved_any <- TRUE
        } else {
          failed_any <- TRUE
          messages <- c(
            messages,
            paste0(
              format,
              ": ",
              result$message
            )
          )
        }
      }

      status <- if (saved_any &&
                    failed_any) {
        "partially_exported"
      } else if (saved_any) {
        "exported"
      } else {
        "failed"
      }

      row_i <- row_i + 1L

      rows[[row_i]] <- data.frame(
        TSS_Window = tss_window,
        Object = object_name,
        Stage = stage_name,
        Hop = hop_name,
        Plot_Path = node$path,
        Plot_Type = plot_type,
        PNG = unname(
          output_paths["PNG"]
        ),
        PDF = unname(
          output_paths["PDF"]
        ),
        HTML = unname(
          output_paths["HTML"]
        ),
        Status = status,
        Message = paste(
          unique(messages),
          collapse = " | "
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  manifest <- if (length(rows) > 0L) {
    do.call(
      rbind,
      rows
    )
  } else {
    empty_manifest
  }

  manifest_path <- file.path(
    root_dir,
    "Visualization_Export_Index.tsv"
  )

  utils::write.table(
    manifest,
    file = manifest_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  index_path <- write_visualization_gallery(
    manifest = manifest,
    root_dir = root_dir,
    title = paste0(
      "looplook benchmark visualization gallery — TSS ±",
      tss_window,
      " bp"
    )
  )

  list(
    root = root_dir,
    index = index_path,
    manifest_path = manifest_path,
    manifest = manifest
  )
}

export_visualizations_from_rdata <- function(
    rdata_path,
    tss_out_dir,
    tss_window,
    cfg
) {
  load_env <- new.env(
    parent = emptyenv()
  )

  loaded_names <- load(
    rdata_path,
    envir = load_env
  )

  benchmark_names <- benchmark_visualization_layout()$Object

  available_names <- intersect(
    benchmark_names,
    loaded_names
  )

  objects <- mget(
    available_names,
    envir = load_env,
    inherits = FALSE
  )

  export_benchmark_visualizations(
    objects = objects,
    root_dir = file.path(
      tss_out_dir,
      cfg$visualization_dir_name
    ),
    tss_window = tss_window,
    cfg = cfg
  )
}

write_global_visualization_overview <- function(
    run_summary,
    version_dir
) {
  overview_path <- file.path(
    version_dir,
    "Visualization_Overview.html"
  )

  links <- character()

  for (i in seq_len(
      nrow(run_summary)
  )) {
    index_path <- run_summary$Visualization_Index[i]

    if (is.na(index_path) ||
        !nzchar(index_path) ||
        !file.exists(index_path)) {
      next
    }

    rel_index <- relative_to_root(
      index_path,
      version_dir
    )

    links <- c(
      links,
      paste0(
        "<li><a href='",
        html_escape(rel_index),
        "'>TSS ±",
        run_summary$TSS_Window[i],
        " bp visualization gallery</a></li>"
      )
    )
  }

  if (length(links) == 0L) {
    links <- "<li>No visualization galleries were generated.</li>"
  }

  html <- c(
    "<!doctype html>",
    "<html lang='en'>",
    "<head>",
    "<meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'>",
    "<title>looplook TSS-window visualization overview</title>",
    "<style>",
    "body{font-family:Arial,sans-serif;max-width:850px;margin:36px auto;padding:0 22px;color:#1d2433}",
    "h1{font-size:26px}",
    "li{margin:14px 0;font-size:17px}",
    "a{color:#2459a9;text-decoration:none}",
    "a:hover{text-decoration:underline}",
    "</style>",
    "</head>",
    "<body>",
    "<h1>looplook TSS-window visualization overview</h1>",
    "<p>Open a TSS-specific gallery to review all plots exported from annotation and refinement objects.</p>",
    "<ul>",
    links,
    "</ul>",
    "</body>",
    "</html>"
  )

  writeLines(
    html,
    overview_path,
    useBytes = TRUE
  )

  overview_path
}


make_run_manifest <- function(
    tss_window,
    tss_region,
    consensus_bedpe,
    consensus_n,
    paths,
    cfg,
    script_path
) {
  package_names <- c(
    "looplook",
    "GenomicRanges",
    "InteractionSet",
    "ChIPseeker",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db"
  )

  package_versions <- vapply(
    package_names,
    function(pkg) {
      as.character(
        tryCatch(
          utils::packageVersion(pkg),
          error = function(e) NA_character_
        )
      )
    },
    character(1)
  )

  list(
    run_time = Sys.time(),
    script_path = script_path,
    script_md5 = safe_md5(script_path),

    tss_window = tss_window,
    tss_region = tss_region,

    consensus_bedpe = consensus_bedpe,
    consensus_bedpe_md5 = safe_md5(consensus_bedpe),
    consensus_loop_n = consensus_n,

    input_paths = paths,
    input_md5 = vapply(
      paths,
      safe_md5,
      character(1)
    ),

    package_versions = package_versions,

    consensus_parameters = list(
      gap = cfg$loop_gap,
      mode = "consensus",
      min_consensus = cfg$min_consensus,
      score_col = cfg$score_col,
      min_raw_score = cfg$min_raw_score,
      min_score = cfg$min_score,
      chaining_policy = cfg$chaining_policy,
      blacklist_species = cfg$blacklist_species,
      roi_mode = cfg$roi_mode
    ),

    annotation_parameters = list(
      tss_region = tss_region,
      anchor_merge_gap = cfg$anchor_merge_gap,
      min_expr = cfg$min_expr,
      conflict_strategy = cfg$conflict_strategy,
      co_dominance_ratio = cfg$co_dominance_ratio,
      anchor_gap = cfg$annotation_anchor_gap,
      anchor_min_overlap = cfg$annotation_anchor_min_overlap,
      anchor_min_frac = cfg$annotation_anchor_min_frac,
      hub_percentile = cfg$hub_percentile,
      hub_metric = cfg$hub_metric
    ),

    expression_parameters = list(
      threshold = cfg$expression_threshold,
      threshold_mode = cfg$expression_threshold_mode,
      unit_type = cfg$expression_unit,
      reclassify_by_expression = cfg$reclassify_by_expression,
      chromatin_beds = "not supplied; formal chromatin refinement performed separately"
    ),

    visualization_parameters = list(
      enabled = cfg$export_visualizations,
      directory_name = cfg$visualization_dir_name,
      formats = cfg$visualization_formats,
      width = cfg$visualization_width,
      height = cfg$visualization_height,
      dpi = cfg$visualization_dpi,
      max_depth = cfg$visualization_max_depth,
      overwrite = cfg$overwrite_visualizations
    ),

    chromatin_parameters = list(
      tss_region = tss_region,
      anchor_gap = cfg$chromatin_anchor_gap,
      anchor_min_overlap = cfg$chromatin_anchor_min_overlap,
      bw_ratio_threshold = cfg$bw_ratio_threshold,
      recompute_targets = cfg$recompute_targets,
      hub_metric = cfg$hub_metric
    )
  )
}


# ------------------------------------------------------------------------------
# 3. Resolve and validate paths
# ------------------------------------------------------------------------------

dir.create(
  cfg$out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

version_dir <- file.path(
  cfg$out_dir,
  cfg$version_tag
)

common_dir <- file.path(
  version_dir,
  "common"
)

dir.create(
  common_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

loop_files <- vapply(
  file.path(
    cfg$data_dir,
    cfg$loop_files
  ),
  normalize_existing_path,
  character(1),
  label = "Loop file"
)

expr_path <- normalize_existing_path(
  file.path(
    cfg$data_dir,
    cfg$expr_file
  ),
  "Expression matrix"
)

target_file <- normalize_existing_path(
  file.path(
    cfg$data_dir,
    cfg$target_file
  ),
  "Target peak BED"
)

roi_file <- normalize_existing_path(
  file.path(
    cfg$data_dir,
    cfg$roi_file
  ),
  "Region-of-interest BED"
)

chromatin_beds <- lapply(cfg$chromatin_beds, function(fname) {
  normalize_existing_path(file.path(cfg$data_dir, fname), fname)
})

chromatin_bw <- lapply(cfg$chromatin_bw, function(fname) {
  normalize_existing_path(file.path(cfg$data_dir, fname), fname)
})

validate_expression_samples(
  expr_path,
  cfg$sample_columns
)

tss_windows <- parse_tss_argument(
  cfg$tss_windows
)

tss_windows <- vapply(
  tss_windows,
  validate_tss_window,
  integer(1)
)

script_path <- get_script_path()

log_message(
  "looplook version: ",
  as.character(
    packageVersion("looplook")
  )
)

log_message(
  "TSS windows: ",
  paste(tss_windows, collapse = ", "),
  " bp"
)


# ------------------------------------------------------------------------------
# 4. Generate one common consensus BEDPE
# ------------------------------------------------------------------------------

consensus_bedpe <- file.path(
  common_dir,
  "consensus_loops_gap500.bedpe"
)

if (!cfg$overwrite_consensus &&
    file.exists(consensus_bedpe)) {
  log_message(
    "Reusing existing consensus BEDPE: ",
    consensus_bedpe
  )

  consensus_result <- NULL
  consensus_lines <- tryCatch(
    readLines(consensus_bedpe, warn = FALSE),
    error = function(e) character(0))
  consensus_n <- sum(nzchar(consensus_lines) & !grepl("^#", consensus_lines))
  log_message("Reused consensus loops: ", consensus_n)
} else {
  log_message(
    "Generating common consensus BEDPE..."
  )

  consensus_result <- consolidate_chromatin_loops(
    files = loop_files,
    gap = cfg$loop_gap,
    mode = "consensus",
    min_consensus = cfg$min_consensus,
    score_col = cfg$score_col,
    min_raw_score = cfg$min_raw_score,
    min_score = cfg$min_score,
    chaining_policy = cfg$chaining_policy,
    blacklist_species = cfg$blacklist_species,
    region_of_interest = roi_file,
    roi_mode = cfg$roi_mode,
    out_file = consensus_bedpe,
    write_output = TRUE,
    quiet = cfg$quiet
  )

  consensus_n <- length(
    consensus_result
  )

  if (!file.exists(consensus_bedpe)) {
    stop(
      "Consensus BEDPE was not created: ",
      consensus_bedpe,
      call. = FALSE
    )
  }

  if (consensus_n == 0L) {
    stop(
      "Consensus loop set is empty. Check score-column detection, ",
      "score thresholds, ROI filtering and input BEDPE files.",
      call. = FALSE
    )
  }

  log_message(
    "Consensus loops: ",
    consensus_n,
    " | MD5=",
    safe_md5(consensus_bedpe)
  )
}


# ------------------------------------------------------------------------------
# 5. Per-window workflow
# ------------------------------------------------------------------------------

run_one_tss_window <- function(tss_window) {
  tss_window <- validate_tss_window(
    tss_window
  )

  tss_region_use <- c(
    -tss_window,
    tss_window
  )

  tss_tag <- paste0(
    "tss",
    tss_window
  )

  tss_out_dir <- file.path(
    version_dir,
    tss_tag
  )

  rdata_path <- file.path(
    tss_out_dir,
    sprintf(
      "tss_%d_res_chromatin.RData",
      tss_window
    )
  )

  if (!cfg$overwrite_tss_results &&
      file.exists(rdata_path)) {
    log_message(
      "[",
      tss_tag,
      "] Reusing existing RData: ",
      rdata_path
    )

    visualization_output <- if (isTRUE(
        cfg$export_visualizations
    )) {
      log_message(
        "[",
        tss_tag,
        "] Exporting visualizations from reused RData"
      )

      export_visualizations_from_rdata(
        rdata_path = rdata_path,
        tss_out_dir = tss_out_dir,
        tss_window = tss_window,
        cfg = cfg
      )
    } else {
      list(
        root = NA_character_,
        index = NA_character_,
        manifest_path = NA_character_,
        manifest = NULL
      )
    }

    return(
      data.frame(
        TSS_Window = tss_window,
        TSS_Region = paste(
          tss_region_use,
          collapse = ","
        ),
        Consensus_BEDPE = consensus_bedpe,
        Consensus_MD5 = safe_md5(consensus_bedpe),
        RData = rdata_path,
        RData_MD5 = safe_md5(rdata_path),
        Manifest = file.path(
          tss_out_dir,
          sprintf(
            "tss_%d_run_manifest.rds",
            tss_window
          )
        ),
        Visualization_Dir = visualization_output$root,
        Visualization_Index = visualization_output$index,
        Visualization_Manifest = visualization_output$manifest_path,
        Visualization_Plots = if (is.null(
            visualization_output$manifest
        )) {
          NA_integer_
        } else {
          sum(
            visualization_output$manifest$Status %in%
              c(
                "exported",
                "partially_exported"
              )
          )
        },
        Status = "reused",
        stringsAsFactors = FALSE
      )
    )
  }

  dir.create(
    tss_out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  annotation_hop0_dir <- file.path(
    tss_out_dir,
    "anno_hop0"
  )

  annotation_hop1_dir <- file.path(
    tss_out_dir,
    "anno_hop1"
  )

  refined_hop0_dir <- file.path(
    tss_out_dir,
    "refined"
  )

  refined_hop1_dir <- file.path(
    tss_out_dir,
    "refined_hop1"
  )

  chromatin_hop0_dir <- file.path(
    tss_out_dir,
    "chromatin_hop0"
  )

  chromatin_hop1_dir <- file.path(
    tss_out_dir,
    "chromatin_hop1"
  )

  chromatin_only_hop0_dir <- file.path(
    tss_out_dir,
    "chromatin_only_hop0"
  )

  chromatin_only_hop1_dir <- file.path(
    tss_out_dir,
    "chromatin_only_hop1"
  )

  dir.create(chromatin_hop0_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(chromatin_hop1_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(chromatin_only_hop0_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(chromatin_only_hop1_dir, recursive = TRUE, showWarnings = FALSE)

  log_message(
    "[",
    tss_tag,
    "] Annotation hop=0"
  )

  res <- annotate_peaks_and_loops(
    bedpe_file = consensus_bedpe,
    target_bed = target_file,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    org_db = "org.Hs.eg.db",
    species = cfg$species,
    tss_region = tss_region_use,
    anchor_merge_gap = cfg$anchor_merge_gap,
    out_dir = annotation_hop0_dir,
    expr_matrix_file = expr_path,
    sample_columns = cfg$sample_columns,
    project_name = cfg$project_name,
    neighbor_hop = 0L,
    hub_percentile = cfg$hub_percentile,
    hub_metric = cfg$hub_metric,
    min_expr = cfg$min_expr,
    conflict_strategy = cfg$conflict_strategy,
    co_dominance_ratio = cfg$co_dominance_ratio,
    anchor_gap = cfg$annotation_anchor_gap,
    anchor_min_overlap = cfg$annotation_anchor_min_overlap,
    anchor_min_frac = cfg$annotation_anchor_min_frac,
    write_output = cfg$write_output,
    quiet = cfg$quiet
  )

  log_message(
    "[",
    tss_tag,
    "] Annotation hop=1"
  )

  res2 <- annotate_peaks_and_loops(
    bedpe_file = consensus_bedpe,
    target_bed = target_file,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    org_db = "org.Hs.eg.db",
    species = cfg$species,
    tss_region = tss_region_use,
    anchor_merge_gap = cfg$anchor_merge_gap,
    out_dir = annotation_hop1_dir,
    expr_matrix_file = expr_path,
    sample_columns = cfg$sample_columns,
    project_name = cfg$project_name,
    neighbor_hop = 1L,
    hub_percentile = cfg$hub_percentile,
    hub_metric = cfg$hub_metric,
    min_expr = cfg$min_expr,
    conflict_strategy = cfg$conflict_strategy,
    co_dominance_ratio = cfg$co_dominance_ratio,
    anchor_gap = cfg$annotation_anchor_gap,
    anchor_min_overlap = cfg$annotation_anchor_min_overlap,
    anchor_min_frac = cfg$annotation_anchor_min_frac,
    write_output = cfg$write_output,
    quiet = cfg$quiet
  )

  log_message(
    "[",
    tss_tag,
    "] Expression refinement hop=0"
  )

  refined_res <- refine_loop_anchors_by_expression(
    annotation_res = res,
    expr_matrix_file = expr_path,
    sample_columns = cfg$sample_columns,
    threshold = cfg$expression_threshold,
    threshold_mode = cfg$expression_threshold_mode,
    unit_type = cfg$expression_unit,
    species = cfg$species,
    out_dir = refined_hop0_dir,
    project_name = cfg$project_name,
    reclassify_by_expression = cfg$reclassify_by_expression,
    hub_percentile = cfg$hub_percentile,
    hub_metric = cfg$hub_metric,

    # Deliberately omitted:
    # chromatin_beds
    #
    # Expression refinement uses chromatin_beds only for an optional
    # validation table. Formal chromatin reclassification is performed below.
    chromatin_beds = list(),

    write_output = cfg$write_output,
    quiet = cfg$quiet,
    allow_rerefine = FALSE
  )

  log_message(
    "[",
    tss_tag,
    "] Expression refinement hop=1"
  )

  refined_res2 <- refine_loop_anchors_by_expression(
    annotation_res = res2,
    expr_matrix_file = expr_path,
    sample_columns = cfg$sample_columns,
    threshold = cfg$expression_threshold,
    threshold_mode = cfg$expression_threshold_mode,
    unit_type = cfg$expression_unit,
    species = cfg$species,
    out_dir = refined_hop1_dir,
    project_name = cfg$project_name,
    reclassify_by_expression = cfg$reclassify_by_expression,
    hub_percentile = cfg$hub_percentile,
    hub_metric = cfg$hub_metric,
    chromatin_beds = list(),
    write_output = cfg$write_output,
    quiet = cfg$quiet,
    allow_rerefine = FALSE
  )

  log_message(
    "[",
    tss_tag,
    "] Expression -> chromatin refinement hop=0"
  )

  cr <- refine_loop_anchors_by_chromatin(
    annotation_res = refined_res,
    chromatin_beds = chromatin_beds,
    anchor_gap = cfg$chromatin_anchor_gap,
    anchor_min_overlap = cfg$chromatin_anchor_min_overlap,
    species = cfg$species,

    # Critical: use the same TSS window as the upstream annotation.
    tss_region = tss_region_use,

    out_dir = chromatin_hop0_dir,
    project_name = cfg$project_name,
    candidate_types = NULL,
    recompute_targets = cfg$recompute_targets,
    write_output = cfg$write_output,
    quiet = cfg$quiet,
    chromatin_bw = chromatin_bw,
    bw_ratio_threshold = cfg$bw_ratio_threshold,
    enhancer_bed = NULL,
    hub_metric = cfg$hub_metric
  )

  log_message(
    "[",
    tss_tag,
    "] Expression -> chromatin refinement hop=1"
  )

  cr2 <- refine_loop_anchors_by_chromatin(
    annotation_res = refined_res2,
    chromatin_beds = chromatin_beds,
    anchor_gap = cfg$chromatin_anchor_gap,
    anchor_min_overlap = cfg$chromatin_anchor_min_overlap,
    species = cfg$species,
    tss_region = tss_region_use,
    out_dir = chromatin_hop1_dir,
    project_name = cfg$project_name,
    candidate_types = NULL,
    recompute_targets = cfg$recompute_targets,
    write_output = cfg$write_output,
    quiet = cfg$quiet,
    chromatin_bw = chromatin_bw,
    bw_ratio_threshold = cfg$bw_ratio_threshold,
    enhancer_bed = NULL,
    hub_metric = cfg$hub_metric
  )

  log_message(
    "[",
    tss_tag,
    "] Chromatin-only refinement hop=0"
  )

  cr_only <- refine_loop_anchors_by_chromatin(
    annotation_res = res,
    chromatin_beds = chromatin_beds,
    anchor_gap = cfg$chromatin_anchor_gap,
    anchor_min_overlap = cfg$chromatin_anchor_min_overlap,
    species = cfg$species,
    tss_region = tss_region_use,
    out_dir = chromatin_only_hop0_dir,
    project_name = cfg$project_name,
    candidate_types = NULL,
    recompute_targets = cfg$recompute_targets,
    write_output = cfg$write_output,
    quiet = cfg$quiet,
    chromatin_bw = chromatin_bw,
    bw_ratio_threshold = cfg$bw_ratio_threshold,
    enhancer_bed = NULL,
    hub_metric = cfg$hub_metric
  )

  log_message(
    "[",
    tss_tag,
    "] Chromatin refinement (no expression reclassification) hop=1"
  )

  cr2_only <- refine_loop_anchors_by_chromatin(
    annotation_res = res2,
    chromatin_beds = chromatin_beds,
    anchor_gap = cfg$chromatin_anchor_gap,
    anchor_min_overlap = cfg$chromatin_anchor_min_overlap,
    species = cfg$species,
    tss_region = tss_region_use,
    out_dir = chromatin_only_hop1_dir,
    project_name = cfg$project_name,
    candidate_types = NULL,
    recompute_targets = cfg$recompute_targets,
    write_output = cfg$write_output,
    quiet = cfg$quiet,
    chromatin_bw = chromatin_bw,
    bw_ratio_threshold = cfg$bw_ratio_threshold,
    enhancer_bed = NULL,
    hub_metric = cfg$hub_metric
  )

  input_paths <- c(
    loop_rep1 = loop_files[1L],
    loop_rep2 = loop_files[2L],
    expression = expr_path,
    target_peaks = target_file,
    roi_h3k27ac = roi_file,
    chromatin_beds = unlist(
      chromatin_beds,
      use.names = TRUE
    ),
    chromatin_bigwig = unlist(
      chromatin_bw,
      use.names = TRUE
    )
  )

  run_manifest <- make_run_manifest(
    tss_window = tss_window,
    tss_region = tss_region_use,
    consensus_bedpe = consensus_bedpe,
    consensus_n = consensus_n,
    paths = input_paths,
    cfg = cfg,
    script_path = script_path
  )

  benchmark_objects <- list(
    res = res,
    refined_res = refined_res,
    cr = cr,
    cr_only = cr_only,
    res2 = res2,
    refined_res2 = refined_res2,
    cr2 = cr2,
    cr2_only = cr2_only
  )

  visualization_output <- if (isTRUE(
      cfg$export_visualizations
  )) {
    log_message(
      "[",
      tss_tag,
      "] Exporting all stored visualizations"
    )

    export_benchmark_visualizations(
      objects = benchmark_objects,
      root_dir = file.path(
        tss_out_dir,
        cfg$visualization_dir_name
      ),
      tss_window = tss_window,
      cfg = cfg
    )
  } else {
    list(
      root = NA_character_,
      index = NA_character_,
      manifest_path = NA_character_,
      manifest = NULL
    )
  }

  objects_to_save <- c(
    benchmark_objects,
    list(
      run_manifest = run_manifest,
      visualization_manifest =
        visualization_output$manifest
    )
  )

  save_named_objects(
    objects_to_save,
    rdata_path
  )

  manifest_path <- file.path(
    tss_out_dir,
    sprintf(
      "tss_%d_run_manifest.rds",
      tss_window
    )
  )

  saveRDS(
    run_manifest,
    manifest_path,
    version = 3
  )

  session_info_path <- file.path(
    tss_out_dir,
    "sessionInfo.txt"
  )

  writeLines(
    capture.output(
      sessionInfo()
    ),
    session_info_path
  )

  log_message(
    "[",
    tss_tag,
    "] Saved: ",
    rdata_path,
    " | MD5=",
    safe_md5(rdata_path)
  )

  rm(
    res,
    res2,
    refined_res,
    refined_res2,
    cr,
    cr2,
    cr_only,
    cr2_only
  )

  invisible(
    gc()
  )

  data.frame(
    TSS_Window = tss_window,
    TSS_Region = paste(
      tss_region_use,
      collapse = ","
    ),
    Consensus_BEDPE = consensus_bedpe,
    Consensus_MD5 = safe_md5(consensus_bedpe),
    RData = rdata_path,
    RData_MD5 = safe_md5(rdata_path),
    Manifest = manifest_path,
    Visualization_Dir = visualization_output$root,
    Visualization_Index = visualization_output$index,
    Visualization_Manifest = visualization_output$manifest_path,
    Visualization_Plots = if (is.null(
        visualization_output$manifest
    )) {
      NA_integer_
    } else {
      sum(
        visualization_output$manifest$Status %in%
          c(
            "exported",
            "partially_exported"
          )
      )
    },
    Status = "completed",
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 6. Run requested windows and export summary
# ------------------------------------------------------------------------------

run_summary_list <- vector(
  "list",
  length(tss_windows)
)

for (i in seq_along(tss_windows)) {
  tss_window <- tss_windows[i]

  log_message(
    "============================================================"
  )

  log_message(
    "Starting TSS window +/-",
    tss_window,
    " bp (",
    i,
    "/",
    length(tss_windows),
    ")"
  )

  run_summary_list[[i]] <- run_one_tss_window(
    tss_window
  )
}

run_summary <- do.call(
  rbind,
  run_summary_list
)

summary_path <- file.path(
  version_dir,
  "TSS_Window_Run_Summary.tsv"
)

utils::write.table(
  run_summary,
  file = summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

visualization_overview_path <- if (isTRUE(
    cfg$export_visualizations
)) {
  write_global_visualization_overview(
    run_summary = run_summary,
    version_dir = version_dir
  )
} else {
  NA_character_
}

log_message(
  "All requested TSS windows completed."
)

log_message(
  "Summary: ",
  summary_path
)

if (!is.na(
    visualization_overview_path
)) {
  log_message(
    "Visualization overview: ",
    visualization_overview_path
  )
}

print(
  run_summary,
  row.names = FALSE
)
