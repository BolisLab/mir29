export_all_panels <- function() {
  panels <- list(
    Figure3_b = list(arrangeGrob(grobs = growth_panels, ncol = 2), 9, 3.7),
    Figure3_c = list(de_plot, 8, 3.5), Figure3_d = list(pca_plot, 6, 4.5),
    Figure3_e = list(arrangeGrob(grobs = main_clock_panels, ncol = 2), 9, 4.5),
    Figure3_f_pan = list(pan_plot, 8, 1.8),
    Figure3_f_HF_HGF = list(arrangeGrob(hf_plot, hgf_plot, ncol = 2), 6.5, 5),
    Figure3_f_HGF_vertical = list(hgf_vertical_plot, 2.5, 4.5),
    Figure3_g = list(main_gsea_plot, 4, 5.5),
    Figure3_h_TF = list(tf_heatmap, 7, 2.6), Figure3_h_PROGENy = list(progeny_heatmap, 7, 5.4),
    Figure3_i = list(gene_plot, 7.5, 6.5),
    Figure3_j = list(arrangeGrob(grobs = network_panels, ncol = 6), 12, 5),
    Figure3_k = list(arrangeGrob(umap_plot, projection_plot, ncol = 2), 10, 4.5),
    Figure3_l = list(alluvial_plot, 7, 5.5), Figure3_m = list(regulator_plot, 7, 3.5),
    FigureS3_a = list(jaccard_plot, 8, 4.2),
    FigureS3_b = list(arrangeGrob(grobs = supplementary_clock_panels, ncol = 3), 12, 4.6),
    FigureS3_c = list(coverage_plot, 8.5, 4.5), FigureS3_d = list(contribution_plot, 13, 3.5),
    FigureS3_e = list(arrangeGrob(grobs = enrichment_plots, ncol = 6), 13, 3.7),
    FigureS3_f = list(overlap_plot, 6.2, 4.8),
    FigureS3_g = list(arrangeGrob(grobs = loto_plots, ncol = 7), 14, 3.7),
    FigureS3_h = list(activity_plot, 14, 3.5),
    FigureS3_i = list(arrangeGrob(marker_plot, umap_plot, ncol = 2, widths = c(2.2, 1)), 12, 4.2),
    FigureS3_j = list(target_plot, 10, 4))
  for (name in names(panels)) {
    x <- panels[[name]]
    png_path <- file.path("figures/reconstructed", paste0(name, ".png"))
    pdf_path <- file.path("figures/reconstructed", paste0(name, ".pdf"))
    ggsave(png_path, x[[1]], width = x[[2]], height = x[[3]], dpi = 240, bg = "white")
    if (name %in% c("Figure3_k", "FigureS3_i")) {
      # Dense UMAP layers are rasterized, as in the accepted Scanpy export.
      grDevices::pdf(pdf_path, width = x[[2]], height = x[[3]])
      grid.newpage(); grid.raster(png::readPNG(png_path)); dev.off()
    } else {
      ggsave(pdf_path, x[[1]], width = x[[2]], height = x[[3]], device = cairo_pdf, bg = "white")
    }
  }
  invisible(names(panels))
}
