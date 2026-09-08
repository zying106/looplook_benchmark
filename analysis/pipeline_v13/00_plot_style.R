# looplook v12 centralized figure style
# Edit this file for publication styling. Changes here do NOT invalidate GSEA caches.

plot_style <- list(
  defaults = list(
    width = NULL,          # NULL keeps the size requested by the analysis code
    height = NULL,
    base_size = 10,
    base_family = "sans",
    legend_position = "top",
    remove_grid = TRUE,
    paired_line = TRUE,
    show_mean = TRUE,
    show_median = TRUE,
    violin_adjust = NULL,
    y_limits = NULL,       # e.g. c(-2.5, 1); applied with coord_cartesian, not scale limits
    x_text_angle = NULL
  ),

  # User-preferred clean distribution style. These are only applied when a
  # matching plot contains the relevant layers; otherwise they are harmless.
  distribution = list(
    width = 3.5,
    height = 6.6,
    paired_line = FALSE,
    show_mean = FALSE,
    show_median = FALSE,
    violin_adjust = 1.2
  ),

  ridge = list(
    width = 7.0,
    height = 6.0
  ),

  distance = list(
    width = 8.0,
    height = 5.5,
    paired_line = FALSE
  ),

  enrichment = list(
    width = 8.0,
    height = 6.0
  ),

  overlap = list(
    width = 7.0,
    height = 6.0
  ),

  summary = list(
    width = 8.0,
    height = 6.0
  ),

  # Optional filename-specific overrides. First matching rule wins after
  # family/default settings are merged.
  overrides = list(
    list(pattern = "LFC.*Violin|Violin.*LFC", settings = list(
      width = 3.5, height = 6.6,
      paired_line = FALSE, show_mean = FALSE, show_median = FALSE,
      violin_adjust = 1.2
    )),
    list(pattern = "NES.*Boxplot|Effect.*Size", settings = list(
      paired_line = FALSE,
      show_mean = FALSE,
      show_median = FALSE
    ))
  )
)

invisible(plot_style)
