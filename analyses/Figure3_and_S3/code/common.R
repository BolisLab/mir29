# Shared display conventions and local input access.
required <- c("ggplot2", "gridExtra", "matrixStats", "minpack.lm", "circlize", "ggforce", "ggalluvial")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages before rendering: ", paste(missing, collapse = ", "))
options(device = function(...) grDevices::pdf(file = NULL, ...))
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})
read_input <- function(name) readRDS(file.path("data/verified", paste0(name, ".rds")))
theme_figure <- function(size = 10) {
  theme_classic(base_size = size, base_family = "Helvetica") +
    theme(plot.title = element_text(size = size, face = "italic"),
          legend.position = "bottom", plot.margin = margin(7, 9, 7, 7))
}
draw_panel <- function(x) {
  if (inherits(x, "ggplot")) print(x) else { grid.newpage(); grid.draw(x) }
}
time_hours <- c(72, 96, 168, 192, 264, 288)
time_positions <- c(1, 2, 4, 5, 7, 8)
time_colours <- c("#E87466", "#D72822", "#4DBAB2", "#286E6A", "#8C73B8", "#4F3492")
arm_colours <- c(INH = "#328A37", inhCTRL = "#A3CDA7", MIM = "#9B42B7", mimCTRL = "#DBAFE8", CTRL = "#666666")
arm_lines <- c(INH = "solid", inhCTRL = "dashed", MIM = "solid", mimCTRL = "dashed", CTRL = "solid")
map_time <- function(hours) approx(time_hours, time_positions, xout = hours, rule = 2)$y
signature_key <- data.frame(
  letter = c("B", "C", "D", "E", "F", "G"),
  label = c("PanTissue-UP", "PanTissue-DN", "HF-UP", "HF-DN", "HGF-UP", "HGF-DN"),
  id = c("pan_aging_up", "pan_aging_dn", "new_panfib_up_R1", "new_panfib_dn_R1", "new_hgf_up_R1", "new_hgf_dn_R1"),
  stringsAsFactors = FALSE
)
final_letter <- function(id) signature_key$letter[match(id, signature_key$id)]
final_label <- function(id) signature_key$label[match(id, signature_key$id)]
support <- readRDS("data/supplementary/signature_support_data.rds")
multilevel <- readRDS("data/supplementary/LOTO_multilevel_results.rds")
core <- readRDS("data/supplementary/core_candidates.rds")
dir.create("validation", showWarnings = FALSE)
dir.create("figures/reconstructed", recursive = TRUE, showWarnings = FALSE)
