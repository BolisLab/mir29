# One entry point. Execute from the bundle root; source files are English only.
source("code/00_setup.R")
writeLines("check\tstatus\tmaximum_absolute_error", "validation/numeric_checks.tsv")
modules <- c("01_primary.R", "02_clocks.R", "02b_clock_decomposition.R", "03_fidelity.R",
             "04a_wnt_locus.R", "04b_wnt_responses.R", "04c_wnt_statistics.R", "04d_wnt_circos.R",
             "05_htert_hpf.R", "05b_htert_clocks.R", "06a_statistics.R",
             "06b_statistics_sensitivity.R", "06c_clock_statistics.R", "07_validate_assets.R",
             "08_current_annotations.R")
for (module in modules) {
  message("Running ", module)
  source(file.path("code", module), local = FALSE)
}
capture.output(sessionInfo(), file = "validation/sessionInfo.txt")
writeLines(names(panels), "validation/reconstructed_panels.txt")
