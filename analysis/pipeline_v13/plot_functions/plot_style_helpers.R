# Helpers for v12 queued-plot customization.

`%||%` <- function(x, y) if (is.null(x)) y else x

merge_lists <- function(x, y) {
  if (is.null(y)) return(x)
  for (nm in names(y)) x[[nm]] <- y[[nm]]
  x
}

infer_plot_family <- function(filename) {
  x <- tolower(basename(filename))
  if (grepl("violin|boxplot|effect.?size|distribution", x)) return("distribution")
  if (grepl("ridge", x)) return("ridge")
  if (grepl("distance", x)) return("distance")
  if (grepl("go_|enrich|gsea", x)) return("enrichment")
  if (grepl("venn|upset|overlap", x)) return("overlap")
  if (grepl("summary|ranking", x)) return("summary")
  "defaults"
}

resolve_plot_settings <- function(target, family = NULL, style = plot_style) {
  fam <- family %||% infer_plot_family(target)
  out <- style$defaults %||% list()
  if (!identical(fam, "defaults") && !is.null(style[[fam]])) {
    out <- merge_lists(out, style[[fam]])
  }
  for (rule in style$overrides %||% list()) {
    if (!is.null(rule$pattern) && grepl(rule$pattern, basename(target), ignore.case = TRUE)) {
      out <- merge_lists(out, rule$settings %||% list())
    }
  }
  out
}

.geom_class <- function(layer) class(layer$geom)[1L]
.stat_class <- function(layer) class(layer$stat)[1L]

is_line_layer <- function(layer) {
  .geom_class(layer) %in% c("GeomLine", "GeomPath", "GeomSegment")
}

is_summary_layer <- function(layer) {
  grepl("StatSummary", .stat_class(layer), fixed = TRUE) ||
    grepl("summary", paste(names(layer$stat_params), collapse = " "), ignore.case = TRUE)
}

is_mean_layer <- function(layer) {
  txt <- paste(capture.output(str(layer$stat_params, max.level = 1)), collapse = " ")
  is_summary_layer(layer) && grepl("mean", txt, ignore.case = TRUE)
}

is_median_layer <- function(layer) {
  txt <- paste(capture.output(str(layer$stat_params, max.level = 1)), collapse = " ")
  is_summary_layer(layer) && grepl("median", txt, ignore.case = TRUE)
}

adjust_violin_layers <- function(p, adjust = NULL) {
  if (is.null(adjust) || !is.finite(adjust) || adjust <= 0) return(p)
  for (i in seq_along(p$layers)) {
    if (.geom_class(p$layers[[i]]) == "GeomViolin") {
      p$layers[[i]]$stat_params$adjust <- adjust
    }
  }
  p
}

strip_optional_layers <- function(p, settings) {
  keep <- rep(TRUE, length(p$layers))
  if (identical(settings$paired_line, FALSE)) {
    keep <- keep & !vapply(p$layers, is_line_layer, logical(1))
  }
  if (identical(settings$show_mean, FALSE)) {
    keep <- keep & !vapply(p$layers, is_mean_layer, logical(1))
  }
  if (identical(settings$show_median, FALSE)) {
    keep <- keep & !vapply(p$layers, is_median_layer, logical(1))
  }
  p$layers <- p$layers[keep]
  p
}

apply_publication_style <- function(p, settings) {
  p <- strip_optional_layers(p, settings)
  p <- adjust_violin_layers(p, settings$violin_adjust)

  base_size <- settings$base_size %||% 10
  base_family <- settings$base_family %||% ""
  p <- p + ggplot2::theme_classic(base_size = base_size, base_family = base_family)

  th <- list()
  if (!is.null(settings$legend_position)) th$legend.position <- settings$legend_position
  if (identical(settings$remove_grid, TRUE)) {
    th$panel.grid.major <- ggplot2::element_blank()
    th$panel.grid.minor <- ggplot2::element_blank()
  }
  if (!is.null(settings$x_text_angle)) {
    th$axis.text.x <- ggplot2::element_text(angle = settings$x_text_angle, hjust = 1)
  }
  if (length(th)) p <- p + do.call(ggplot2::theme, th)

  if (!is.null(settings$y_limits) && length(settings$y_limits) == 2L && all(is.finite(settings$y_limits))) {
    # coord_cartesian avoids deleting violin density tails / statistical rows.
    p <- p + ggplot2::coord_cartesian(ylim = settings$y_limits)
  }
  p
}

customize_queued_plot <- function(task, style = plot_style) {
  settings <- resolve_plot_settings(task$Target, task$StyleFamily %||% NULL, style)
  p <- apply_publication_style(task$Plot, settings)

  save_args <- task$GgsaveArgs %||% list()
  if (!is.null(settings$width)) save_args$width <- settings$width
  if (!is.null(settings$height)) save_args$height <- settings$height

  list(plot = p, ggsave_args = save_args, settings = settings)
}

extract_plot_ready_data <- function(p) {
  plot_data <- p$data
  layer_data <- lapply(p$layers, function(layer) {
    d <- layer$data
    if (inherits(d, "waiver")) NULL else d
  })
  names(layer_data) <- paste0("layer_", seq_along(layer_data))
  list(plot_data = plot_data, layer_data = layer_data)
}
