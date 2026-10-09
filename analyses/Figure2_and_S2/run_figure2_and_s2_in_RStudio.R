# Open Figure2_main.Rproj, then copy/paste these lines into RStudio or source
# this script. All analysis inputs are in this portable folder.
stopifnot(file.exists("Figure2_and_S2.Rmd"))
options(figure2.main.root = normalizePath(getwd(), mustWork = TRUE))
source("code/figure2_mirna_phylogeny.R")
source("code/figure2_methylation.R")
figure2_plot_objects <- list(panel_a = panel_a_tree, panel_f = panel_f,
                            panel_g = panel_g, panel_h = panel_h,
                            panel_i = panel_i, panel_j = panel_j, panel_k = panel_k)
figure2_s2_root <- "."
s2_export <- TRUE
source("code/figure2_s2_build.R")
stopifnot(all(figure2_mirna_checks$Status == "PASS"),
          all(figure2_methylation_checks$Status == "PASS"))
message("Done. Main exports: figures/generated; S2 exports: figures/supplementary_s2.")
# To reload later without running the analyses:
# load("figures/supplementary_s2/Figure2_and_S2_workspace.RData")
# p_S2h
# figure2_plot_objects$panel_k
