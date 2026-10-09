# Portable reconstruction of final Figure 3l from the frozen target taxonomy.
# The original taxonomy-building source and target-score table are retained in
# provenance/historical_code/source_legacy/ and data/main/alluvial/.

figure3_build_alluvial <- function(results) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Install ggplot2 and ggalluvial to rebuild Figure 3l.")
  }

  taxo_plot <- results$taxonomy_plotted
  required <- c(
    "Gene_lbl", "Fine_lbl", "Macro_lbl", "Freq", "flow_group"
  )
  stopifnot(all(required %in% names(taxo_plot)))

  p0 <- ggplot2::ggplot(
    taxo_plot,
    ggplot2::aes(
      y = Freq, axis1 = Gene_lbl, axis2 = Fine_lbl, axis3 = Macro_lbl
    )
  ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = flow_group, alpha = flow_group),
      width = 0.10, knot.pos = 0.40, decreasing = FALSE,
      lode.guidance = "forward", aes.bind = "flows", discern = TRUE,
      color = NA
    ) +
    ggalluvial::geom_stratum(
      width = 0.10, fill = "grey97", color = "grey35",
      linewidth = 0.25, discern = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        Background = "#C8C8C8", `DNA methylation` = "magenta",
        `Chromatin regulation` = "#7B61FF"
      ),
      breaks = c("DNA methylation", "Chromatin regulation"), name = NULL
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        Background = 0.28, `DNA methylation` = 0.95,
        `Chromatin regulation` = 0.95
      ),
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(
      limits = c("Gene", "Fine", "Macro"), expand = c(0.02, 0.02)
    ) +
    ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.width = grid::unit(0.45, "cm"),
      plot.title = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 8)),
      plot.margin = ggplot2::margin(8, 75, 8, 24)
    ) +
    ggplot2::labs(title = "miR-29 targeted genes")

  strata <- ggplot2::ggplot_build(p0)$data[[2]]
  epi <- strata[
    strata$x == 3 & as.character(strata$stratum) == "M: Epigenetics",
    , drop = FALSE
  ]
  stopifnot(nrow(epi) == 1L)
  n_target <- nrow(taxo_plot)

  p0 +
    ggplot2::annotate(
      "rect", xmin = 2.95, xmax = 3.05,
      ymin = epi$ymin, ymax = epi$ymax, fill = "black"
    ) +
    ggplot2::annotate(
      "text", x = 3.07, y = epi$y, label = "Epigenetics",
      hjust = 0, fontface = "bold", size = 3.2
    ) +
    ggplot2::annotate(
      "text", x = 3.08,
      y = seq(epi$y - 12, epi$y - 57, length.out = 6),
      label = c("TET3", "TET1", "TET2", "TDG", "DNMT3A", "DNMT3B"),
      hjust = 0, size = 3.1
    ) +
    ggplot2::annotate(
      "segment", x = 0.82, xend = 0.82, y = 1, yend = n_target + 2,
      linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed")
    ) +
    ggplot2::annotate(
      "text", x = 0.72, y = n_target / 2,
      label = "Higher predicted repression", angle = 90, size = 3.1
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.9, 3.05), ylim = c(0, n_target + 2), clip = "off"
    )
}
