# Run from the root of this bundle. No network or original workstation paths are required.
for (directory in c("tables", "validation", "figures/reconstructed")) dir.create(directory, recursive = TRUE,
  showWarnings = FALSE)
for (module in c("00_setup.R", "01_main_analysis.R", "02_main_plots.R", "03_supplementary_analysis.R",
  "04_supplementary_plots.R")) {
  message("Running ", module)
  source(file.path("code", module), local = TRUE)
}
jsonlite::write_json(CHECKS, "validation/numerical_checks.json", pretty = TRUE, auto_unbox = TRUE)
writeLines(vapply(CHECKS, function(x) paste(if (x$pass) "PASS" else "FAIL", x$check, x$detail, sep = " | "),
  character(1)), "validation/numerical_checks.txt")
writeLines(capture.output(sessionInfo()), "validation/sessionInfo.txt")
saveRDS(list(tcga_counts = tcga_counts, tcga_ranks = tcga_ranks, loo_tcga = loo_tcga, ccle_counts = ccle_counts,
  ccle_ranks = ccle_ranks, loo_ccle = loo_ccle, rna_boxes = rna_boxes, pair_counts = pair_counts, loo_paired = loo_paired,
  reference_qc = reference_qc, reference_distributions = reference_distributions, fidelity_age = fidelity_age,
  fidelity_age_meta = fidelity_age_meta), "data/FigureS5_panel_data.rds")
message("Completed ", length(CHECKS), " numerical checks.")
