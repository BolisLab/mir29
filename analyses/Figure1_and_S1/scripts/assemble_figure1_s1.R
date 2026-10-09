# Assemble the current nine-panel Supplementary Figure S1.
# Run after all S1 calculation chunks in Figure1_and_S1.Rmd.
# This file changes only layout, numbering, and exports, not source estimates.

s1_required_objects <- paste0("s1_panel_", letters[1:9])
s1_evaluation_environment <- environment()
stopifnot(all(vapply(s1_required_objects, function(name) {
  exists(name, envir = s1_evaluation_environment, inherits = TRUE)
}, logical(1))))

# The relocated panels already carry standalone tags. Use their untagged cores
# when a composite supplies aligned outer labels.
s1_display_panels <- list(
  panel_a = s1_panel_a,
  panel_b = s1_panel_b,
  panel_c = s1_panel_c_core,
  panel_d = s1_panel_d_core,
  panel_e = s1_panel_e_core,
  panel_f = s1_panel_f,
  panel_g = s1_panel_g,
  panel_h = s1_panel_h,
  panel_i = s1_panel_i
)

s1_panel_dimensions <- data.frame(
  name = names(s1_display_panels),
  width = c(8.8, 6.0, 9.5, 7.2, 10.0, 8.4, 8.4, 8.8, 9.4),
  height = c(5.8, 5.8, 4.8, 6.2, 5.3, 4.7, 4.7, 5.3, 5.3),
  stringsAsFactors = FALSE
)

for (index in seq_len(nrow(s1_panel_dimensions))) {
  specification <- s1_panel_dimensions[index, ]
  plot_object <- s1_display_panels[[specification$name]]
  for (format in c("svg", "png")) {
    ggsave(
      filename = file.path(s1_output_dir, paste0(specification$name, ".", format)),
      plot = plot_object,
      device = if (format == "svg") svglite::svglite else ragg::agg_png,
      width = specification$width, height = specification$height,
      units = "in", dpi = 300, bg = "white", limitsize = FALSE
    )
  }
}

# The A4 overview places panels drawn at different native device sizes next
# to each other. Use their high-resolution PNG previews here to preserve all
# fixed-size text/heatmap geometry. Individual SVG exports remain editable.
s1_scaled_previews <- lapply(names(s1_display_panels), function(name) {
  image_path <- file.path(s1_output_dir, paste0(name, ".png"))
  if (name == "panel_d") {
    # Remove only the empty margins around the fixed-size clock heatmap.
    image <- png::readPNG(image_path)
    occupied <- image[, , 1] < 0.99 | image[, , 2] < 0.99 | image[, , 3] < 0.99
    rows <- which(rowSums(occupied) > 0)
    columns <- which(colSums(occupied) > 0)
    row_range <- seq(max(1, min(rows) - 12), min(nrow(image), max(rows) + 12))
    column_range <- seq(max(1, min(columns) - 12), min(ncol(image), max(columns) + 12))
    image <- image[row_range, column_range, , drop = FALSE]
  } else {
    image <- png::readPNG(image_path, native = TRUE)
  }
  ggdraw() + draw_grob(grid::rasterGrob(image, interpolate = TRUE))
})
names(s1_scaled_previews) <- names(s1_display_panels)

s1_labeled_row <- function(panel_ids, relative_widths) {
  plot_grid(
    plotlist = s1_scaled_previews[paste0("panel_", panel_ids)],
    nrow = 1, rel_widths = relative_widths,
    labels = panel_ids, label_size = 14, label_fontface = "plain",
    label_x = 0.005, label_y = 0.995, hjust = 0, vjust = 1
  )
}

s1_row_1 <- s1_labeled_row(c("a", "b"), c(1.30, 0.95))
s1_row_2 <- s1_labeled_row(c("c", "d", "e"), c(1.40, 0.82, 1.33))
s1_row_3 <- s1_labeled_row(c("f", "g"), c(1.02, 1.12))
s1_row_4 <- s1_labeled_row(c("h", "i"), c(1.15, 1.20))
s1_composite_tight <- plot_grid(
  s1_row_1, s1_row_2, s1_row_3, s1_row_4,
  ncol = 1, rel_heights = c(1.06, 0.86, 0.90, 1.13)
)
s1_composite_a4_body <- plot_grid(
  s1_composite_tight, NULL, ncol = 1, rel_heights = c(0.90, 0.10)
)
s1_composite_a4 <- ggdraw() +
  draw_plot(s1_composite_a4_body, x = 0.035, y = 0.005,
            width = 0.93, height = 0.985)


ggsave(
  file.path(s1_output_dir, "supplementary_figure_s1_analytical.pdf"),
  s1_composite_a4, width = 8.2677, height = 11.6929, units = "in",
  device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
)
ggsave(
  file.path(s1_output_dir, "supplementary_figure_s1_analytical.png"),
  s1_composite_a4, width = 8.2677, height = 11.6929, units = "in",
  device = ragg::agg_png, dpi = 300, bg = "white", limitsize = FALSE
)

# Store rendered grobs rather than ggplot evaluation environments, avoiding
# accidental serialization of the large RNA analysis workspace.
s1_panel_grobs <- lapply(s1_display_panels, ggplotGrob)
saveRDS(
  s1_panel_grobs,
  file.path(s1_output_dir, "supplementary_figure_s1_panel_grobs.rds"),
  compress = "gzip"
)
write.table(
  s1_panel_dimensions,
  file.path(s1_output_dir, "current_panel_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("All nine S1 panels exported with the 2 October 2026 numbering.")
