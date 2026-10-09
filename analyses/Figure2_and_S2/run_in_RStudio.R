# Figure 2 main — one-click analytical rebuild for RStudio
#
# Open Figure2_main.Rproj first, then run:
#   source("run_in_RStudio.R")

options(stringsAsFactors = FALSE)
options(figure2.main.root = normalizePath(getwd(), mustWork = TRUE))

source(file.path("code", "figure2_mirna_phylogeny.R"))
source(file.path("code", "figure2_methylation.R"))

figure2_plot_objects <- list(
  panel_a = panel_a_tree,
  panel_f = panel_f,
  panel_g = panel_g,
  panel_h = panel_h,
  panel_i = panel_i,
  panel_j = panel_j,
  panel_k = panel_k
)

dir.create(file.path("figures", "generated"), recursive = TRUE, showWarnings = FALSE)

saveRDS(
  figure2_plot_objects,
  file.path("figures", "generated", "figure2_plot_objects.rds"),
  compress = "gzip"
)

save(
  panel_a_tree, panel_f, panel_g, panel_h, panel_i, panel_j, panel_k,
  figure2_mirna_checks, figure2_methylation_checks,
  file = file.path("figures", "generated", "figure2_workspace.RData"),
  compress = "gzip"
)

stopifnot(
  all(figure2_mirna_checks$Status == "PASS"),
  all(figure2_methylation_checks$Status == "PASS")
)

message(
  "Figure 2 analytical rebuild complete.\n",
  "Plots and workspaces: ",
  normalizePath(file.path("figures", "generated"), mustWork = TRUE)
)
