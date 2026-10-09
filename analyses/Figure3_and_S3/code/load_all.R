# Run from the portable project root; only bundled relative paths are used.
modules <- c("common", "upstream_methods", "signature_analysis", "enrichment_methods", "signature_plots",
             "clock_panels", "activity_panels", "other_panels", "singlecell_panels", "figure3_alluvial_panel")
for (module in modules) source(file.path("code", paste0(module, ".R")), local = environment())
alluvial_plot <- figure3_build_alluvial(alluvial)
capture.output(sessionInfo(), file = "validation/sessionInfo.txt")
