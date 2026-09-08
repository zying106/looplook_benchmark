#!/usr/bin/env Rscript

# Replot the primary hop0 convergence panel from existing GSEA iteration CSVs.
# This script does not run GSEA and does not modify the benchmark master/cache.
#
# Usage:
#   Rscript plot_convergence_cma_hop0.R /absolute/path/to/res_v13
#
# Optional arguments:
#   2: sample size (default: 200)
#   3: expected iterations per mode (default: 300)
#   4: output PDF path
#
# Alternatively, set LOOPLOOK_OUT_BASE instead of providing argument 1.

args <- commandArgs(trailingOnly = TRUE)

out_base <- if (length(args) >= 1L && nzchar(args[1L])) {
  path.expand(args[1L])
} else {
  path.expand(Sys.getenv("LOOPLOOK_OUT_BASE", ""))
}

if (!nzchar(out_base)) {
  stop(
    "Provide the res_v13 directory as argument 1 or set LOOPLOOK_OUT_BASE.",
    call. = FALSE
  )
}

sample_size <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
expected_iterations <- if (length(args) >= 3L) as.integer(args[3L]) else 300L

if (is.na(sample_size) || sample_size <= 0L) {
  stop("Sample size must be a positive integer.", call. = FALSE)
}
if (is.na(expected_iterations) || expected_iterations <= 0L) {
  stop("Expected iterations must be a positive integer.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.", call. = FALSE)
}

sim_dir <- file.path(
  out_base,
  paste0("size", sample_size),
  "tmp_sim"
)

if (!dir.exists(sim_dir)) {
  stop("GSEA iteration directory not found: ", sim_dir, call. = FALSE)
}

output_dir <- file.path(out_base, "convergence_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- if (length(args) >= 4L && nzchar(args[4L])) {
  path.expand(args[4L])
} else {
  file.path(
    output_dir,
    sprintf("Convergence_CMA_hop0_size%d.pdf", sample_size)
  )
}
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

source_csv <- file.path(
  dirname(output_pdf),
  sprintf("Convergence_CMA_hop0_size%d_source_data.csv", sample_size)
)
summary_csv <- file.path(
  dirname(output_pdf),
  sprintf("Convergence_CMA_hop0_size%d_summary.csv", sample_size)
)

sim_files <- list.files(
  sim_dir,
  pattern = "^sim_.*\\.csv$",
  full.names = TRUE
)

if (length(sim_files) == 0L) {
  stop("No sim_*.csv files found in: ", sim_dir, call. = FALSE)
}

mode_from_file <- function(path) {
  sub("\\.csv$", "", sub("^sim_", "", basename(path)))
}

is_primary_hop0_mode <- function(mode) {
  grepl(
    "^(anno|refined|chrom|chrom_only)_(all|promoter)_(F|T)$",
    mode
  )
}

sim_modes <- vapply(sim_files, mode_from_file, character(1))
keep <- vapply(sim_modes, is_primary_hop0_mode, logical(1))
sim_files <- sim_files[keep]
sim_modes <- sim_modes[keep]

if (length(sim_files) == 0L) {
  stop(
    "No primary hop0 mode files were found. Expected names such as ",
    "sim_anno_all_F.csv or sim_chrom_promoter_T.csv.",
    call. = FALSE
  )
}

if (length(sim_files) != 16L) {
  warning(
    sprintf(
      "Expected 16 primary hop0 modes but found %d. The available modes will be plotted.",
      length(sim_files)
    ),
    call. = FALSE
  )
}

make_cma <- function(x) {
  valid <- is.finite(x)
  cumulative_n <- cumsum(valid)
  cumulative_sum <- cumsum(ifelse(valid, x, 0))
  out <- cumulative_sum / pmax(cumulative_n, 1L)
  out[cumulative_n == 0L] <- NA_real_
  out
}

read_one_mode <- function(path, mode) {
  dat <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(e) {
      warning("Failed to read ", path, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(dat) || nrow(dat) == 0L) return(NULL)
  if (!"NES" %in% names(dat)) {
    warning("Skipping file without an NES column: ", path, call. = FALSE)
    return(NULL)
  }

  if ("ID" %in% names(dat)) {
    dat <- dat[as.character(dat$ID) == "looplook", , drop = FALSE]
  }
  if (nrow(dat) == 0L) {
    warning("Skipping file without looplook rows: ", path, call. = FALSE)
    return(NULL)
  }

  if ("Iteration" %in% names(dat)) {
    iteration_value <- suppressWarnings(as.integer(dat$Iteration))
    ord <- order(iteration_value, na.last = TRUE)
    dat <- dat[ord, , drop = FALSE]
  }

  dat <- utils::head(dat, expected_iterations)
  nes <- suppressWarnings(as.numeric(dat$NES))

  data.frame(
    Iteration = seq_along(nes),
    NES = nes,
    CMA = make_cma(nes),
    Mode = mode,
    SourceFile = normalizePath(path, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

rows <- Map(read_one_mode, sim_files, sim_modes)
rows <- rows[!vapply(rows, is.null, logical(1))]
if (length(rows) == 0L) {
  stop("None of the hop0 CSV files contained usable NES values.", call. = FALSE)
}

conv_df <- do.call(rbind, rows)
rownames(conv_df) <- NULL

conv_df$Pipeline <- ifelse(
  grepl("^anno_", conv_df$Mode), "Basic",
  ifelse(
    grepl("^refined_", conv_df$Mode), "E-refined",
    ifelse(
      grepl("^chrom_only_", conv_df$Mode), "C-refined",
      ifelse(grepl("^chrom_", conv_df$Mode), "I-refined", "Other")
    )
  )
)
conv_df$Pipeline <- factor(
  conv_df$Pipeline,
  levels = c("Basic", "E-refined", "C-refined", "I-refined")
)

mode_split <- split(conv_df, conv_df$Mode)
summary_rows <- lapply(mode_split, function(x) {
  last_valid <- tail(x$CMA[is.finite(x$CMA)], 1L)
  data.frame(
    Mode = x$Mode[1L],
    Pipeline = as.character(x$Pipeline[1L]),
    ObservedIterations = nrow(x),
    FiniteNES = sum(is.finite(x$NES)),
    ExpectedIterations = expected_iterations,
    Complete = nrow(x) >= expected_iterations,
    FinalCMA = if (length(last_valid) == 1L) last_valid else NA_real_,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
summary_df <- summary_df[order(summary_df$Pipeline, summary_df$Mode), , drop = FALSE]

utils::write.csv(conv_df, source_csv, row.names = FALSE)
utils::write.csv(summary_df, summary_csv, row.names = FALSE)

pipeline_colors <- c(
  "Basic" = "#1F77B4",
  "E-refined" = "#9467BD",
  "C-refined" = "#FFA500",
  "I-refined" = "#E64B35"
)

p <- ggplot2::ggplot(
  conv_df,
  ggplot2::aes(
    x = Iteration,
    y = CMA,
    group = Mode,
    color = Pipeline
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted",
    color = "grey60",
    linewidth = 0.3
  ) +
  ggplot2::geom_line(linewidth = 0.3, alpha = 0.55) +
  ggplot2::scale_color_manual(values = pipeline_colors, drop = FALSE) +
  ggplot2::labs(
    title = sprintf(
      "CMA convergence of primary hop0 modes (sample size = %d)",
      sample_size
    ),
    subtitle = sprintf(
      "%d modes; up to %d GSEA iterations per mode",
      length(unique(conv_df$Mode)),
      expected_iterations
    ),
    x = "Iteration",
    y = "Cumulative moving average of NES",
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 9, color = "grey40", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3)
  )

ggplot2::ggsave(
  filename = output_pdf,
  plot = p,
  width = 8,
  height = 4,
  units = "in"
)

message("Primary hop0 modes plotted: ", length(unique(conv_df$Mode)))
message("PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Source data: ", normalizePath(source_csv, mustWork = FALSE))
message("Summary: ", normalizePath(summary_csv, mustWork = FALSE))

