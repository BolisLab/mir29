# Figure S2 panels d-g: portable, frozen-data reconstruction
#
# Usage from the portable Figure2_main project root:
#   figure2_s2_root <- "."
#   s2_export <- FALSE
#   source("code/figure2_s2_defg.R")
#
# Objects exposed by this module:
#   p_S2d, p_S2e, p_S2f, p_S2g
#   s2_defg_data, s2_defg_checks, s2_defg_dimensions
#
# Analytical values are read from frozen tables. This module does not rerun
# correlations, methylation preprocessing, sequence classification, or any
# bootstrap. Exports are optional and are disabled by default.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gghalves)
  library(cowplot)
  library(grid)
})

if (!exists("figure2_s2_root", inherits = TRUE)) {
  figure2_s2_root <- "."
}
if (!exists("s2_export", inherits = TRUE)) {
  s2_export <- FALSE
}
s2_export <- isTRUE(s2_export)

figure2_s2_root <- normalizePath(
  path.expand(as.character(figure2_s2_root)[1]),
  mustWork = TRUE
)
s2_defg_data_dir <- file.path(
  figure2_s2_root, "data", "figure2_s2", "defg"
)
s2_defg_export_dir <- file.path(
  figure2_s2_root, "figures", "supplementary_s2"
)

s2_defg_paths <- c(
  S2d_coverage = file.path(
    s2_defg_data_dir, "S2d_methylation_coverage_source.tsv"
  ),
  S2e_target = file.path(
    s2_defg_data_dir, "S2e_tissue_target_summary.tsv"
  ),
  S2e_points = file.path(
    s2_defg_data_dir, "S2e_tissue_display_points.tsv"
  ),
  S2e_distribution = file.path(
    s2_defg_data_dir, "S2e_tissue_full_distributions.tsv.gz"
  ),
  S2f_beta_species = file.path(
    s2_defg_data_dir, "S2f_beta_species_summary.tsv"
  ),
  S2f_beta_tissue_provenance = file.path(
    s2_defg_data_dir, "S2f_beta_tissue_summary_provenance.tsv"
  ),
  S2g_rho_species = file.path(
    s2_defg_data_dir, "S2g_rho_species_summary.tsv"
  )
)

missing_s2_defg <- s2_defg_paths[!file.exists(s2_defg_paths)]
if (length(missing_s2_defg)) {
  stop(
    "Missing frozen Figure S2 d-g input(s):\n  ",
    paste(missing_s2_defg, collapse = "\n  ")
  )
}

read_frozen_tsv <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "")
  )
}

read_frozen_tsv_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "")
  )
}

abbreviate_species <- function(x) {
  vapply(strsplit(x, " +"), function(parts) {
    if (length(parts) < 2L) return(parts[[1]])
    paste0(substr(parts[[1]], 1, 1), ". ", paste(parts[-1], collapse = " "))
  }, character(1))
}

pretty_tissue_compact <- function(x) {
  out <- x
  out[out == "BoneMarrow"] <- "B. marrow"
  out[out == "CerebralCortex"] <- "Cortex"
  out[out == "LymphNode"] <- "L. node"
  out
}

clade_colours <- c(
  Primates = "#477DB3",
  Chiroptera = "#7856A6",
  Carnivora = "#C95C78",
  Artiodactyla = "#E89B4C",
  Cetacea = "#42A7A1"
)

sign_colours <- c(
  Negative = "#72BF98",
  Positive = "#6D3E91"
)

sequence_status_colours <- c(
  `Intact CpG` = "#72BF98",
  `Disrupted/non-conserved` = "#D95F5F"
)

target_fill <- "#FFD92F"
target_probe <- "cg10501210"
probe_universe_n <- 37554L

theme_s2_defg <- theme_classic(base_family = "Helvetica", base_size = 7.1) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.36),
    plot.margin = margin(3, 3, 3, 3)
  )

# -----------------------------------------------------------------------------
# S2d. Methylation species/sample atlas (36 species; 5 macroclades).
#
# `n / tissues` reports profiles and tissue strata meeting the upstream n >= 10
# rule. It intentionally describes coverage, including the two Human strata
# with no age variance; those two strata are omitted from S2e because no rho can
# be estimated there.
# -----------------------------------------------------------------------------

s2d_coverage <- read_frozen_tsv(s2_defg_paths[["S2d_coverage"]])
s2d_coverage$clade <- factor(
  s2d_coverage$clade,
  levels = names(clade_colours)
)
s2d_coverage <- s2d_coverage[order(s2d_coverage$species_order), , drop = FALSE]

stopifnot(
  nrow(s2d_coverage) == 36L,
  identical(as.integer(table(s2d_coverage$clade)), c(5L, 10L, 6L, 9L, 6L)),
  sum(s2d_coverage$n_filtered) == 5774L,
  sum(s2d_coverage$n_eligible) == 5660L,
  sum(s2d_coverage$n_eligible_strata) == 67L,
  sum(s2d_coverage$n_target) == 5617L,
  sum(s2d_coverage$n_target_strata) == 65L
)

s2d_coverage$y <- nrow(s2d_coverage) - seq_len(nrow(s2d_coverage)) + 1
s2d_coverage$stats_label <- paste0(
  format(s2d_coverage$n_eligible, big.mark = ",", trim = TRUE),
  " / ",
  s2d_coverage$n_eligible_strata
)

s2d_log_range <- range(log10(s2d_coverage$n_eligible))
s2d_coverage$bar_start <- 0.545
s2d_coverage$bar_end <- s2d_coverage$bar_start +
  0.185 * (log10(s2d_coverage$n_eligible) - s2d_log_range[[1]]) /
  diff(s2d_log_range)

s2d_clade_blocks <- do.call(
  rbind,
  lapply(split(s2d_coverage, s2d_coverage$clade), function(z) {
    data.frame(
      clade = as.character(z$clade[[1]]),
      ymin = min(z$y) - 0.39,
      ymax = max(z$y) + 0.39,
      ymid = mean(range(z$y)),
      stringsAsFactors = FALSE
    )
  })
)
s2d_clade_blocks$clade <- factor(
  s2d_clade_blocks$clade,
  levels = names(clade_colours)
)

p_S2d <- ggplot(s2d_coverage, aes(y = y)) +
  geom_segment(
    data = s2d_clade_blocks,
    aes(x = 0.300, xend = 0.300, y = ymin, yend = ymax, colour = clade),
    inherit.aes = FALSE,
    linewidth = 2.05,
    lineend = "butt",
    show.legend = FALSE
  ) +
  geom_text(
    data = s2d_clade_blocks,
    aes(x = 0.315, y = ymid, label = as.character(clade), colour = clade),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.12,
    show.legend = FALSE
  ) +
  geom_segment(
    aes(x = bar_start, xend = bar_end, yend = y, colour = clade),
    linewidth = 2.45,
    lineend = "round",
    alpha = 0.58,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = bar_end, fill = clade),
    shape = 21,
    size = 1.85,
    stroke = 0.31,
    colour = "#222222",
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = 0.006, label = species_label),
    hjust = 0,
    family = "Helvetica",
    fontface = "italic",
    size = 2.16
  ) +
  geom_text(
    aes(x = 0.985, label = stats_label),
    hjust = 1,
    family = "Helvetica",
    size = 2.08
  ) +
  annotate(
    "text", x = 0.006, y = 37.30, label = "Species",
    hjust = 0, family = "Helvetica", fontface = "bold", size = 2.40
  ) +
  annotate(
    "text", x = 0.300, y = 37.30, label = "Macroclade",
    hjust = 0, family = "Helvetica", fontface = "bold", size = 1.90
  ) +
  annotate(
    "text", x = 0.505, y = 37.30, label = "Sample coverage (log)",
    hjust = 0, family = "Helvetica", fontface = "bold", size = 1.65
  ) +
  annotate(
    "text", x = 0.985, y = 37.30, label = "n / tissues",
    hjust = 1, family = "Helvetica", fontface = "bold", size = 1.85
  ) +
  scale_colour_manual(values = clade_colours, drop = FALSE) +
  scale_fill_manual(values = clade_colours, drop = FALSE) +
  coord_cartesian(
    xlim = c(0, 1), ylim = c(0.45, 37.65), clip = "off"
  ) +
  theme_void(base_family = "Helvetica", base_size = 7.1) +
  theme(plot.margin = margin(2, 4, 2, 2))

# -----------------------------------------------------------------------------
# S2e. Tissue-resolved distribution of age-correlation coefficients.
#
# Each half-violin/boxplot is the complete frozen distribution of 37,554
# tissue-specific CpG Spearman rho values. The yellow point is cg10501210.
# Human Keratinocyte and Placenta are deliberately omitted: although each met
# n >= 10, age was invariant and neither a tissue-wide rho distribution nor a
# target rho could be estimated. Exactly 41 informative strata remain.
# -----------------------------------------------------------------------------

s2e_target_all <- read_frozen_tsv(s2_defg_paths[["S2e_target"]])
s2e_points_all <- read_frozen_tsv(s2_defg_paths[["S2e_points"]])
s2e_distribution_all <- read_frozen_tsv_gz(
  s2_defg_paths[["S2e_distribution"]]
)

s2e_target <- s2e_target_all[
  s2e_target_all$informative & is.finite(s2e_target_all$target_rho),
  , drop = FALSE
]
s2e_valid_keys <- s2e_target$tissue_key
s2e_points <- s2e_points_all[
  s2e_points_all$tissue_key %in% s2e_valid_keys &
    is.finite(s2e_points_all$rho),
  , drop = FALSE
]
s2e_distribution <- s2e_distribution_all[
  s2e_distribution_all$tissue_key %in% s2e_valid_keys &
    is.finite(s2e_distribution_all$rho),
  , drop = FALSE
]
rm(s2e_distribution_all)

s2e_species_rows <- list(
  top = c(
    "Phoca vitulina", "Zalophus californianus", "Bos taurus", "Ovis aries",
    "Delphinapterus leucas", "Lagenorhynchus obliquidens", "Orcinus orca",
    "Tursiops truncatus"
  ),
  middle = c(
    "Chlorocebus sabaeus", "Macaca mulatta", "Papio hamadryas"
  ),
  bottom = "Homo sapiens"
)

s2e_expected_species <- unname(unlist(s2e_species_rows, use.names = FALSE))
stopifnot(
  nrow(s2e_target_all) == 43L,
  nrow(s2e_target) == 41L,
  setequal(unique(s2e_target$species), s2e_expected_species),
  sum(s2e_target$species %in% s2e_species_rows$top) == 17L,
  sum(s2e_target$species %in% s2e_species_rows$middle) == 11L,
  sum(s2e_target$species %in% s2e_species_rows$bottom) == 13L,
  nrow(s2e_distribution) == probe_universe_n * 41L,
  length(unique(s2e_distribution$probe)) == probe_universe_n,
  target_probe == "cg10501210"
)

make_s2e_row <- function(species_in_row) {
  meta_parts <- lapply(species_in_row, function(sp) {
    z <- s2e_target[s2e_target$species == sp, , drop = FALSE]
    z[order(tolower(z$tissue)), , drop = FALSE]
  })
  meta <- do.call(rbind, meta_parts)
  rownames(meta) <- NULL
  meta$x <- seq_len(nrow(meta))
  x_lookup <- setNames(meta$x, meta$tissue_key)

  dat <- s2e_distribution[
    s2e_distribution$tissue_key %in% meta$tissue_key,
    , drop = FALSE
  ]
  dat$x <- unname(x_lookup[dat$tissue_key])

  pts <- s2e_points[
    s2e_points$tissue_key %in% meta$tissue_key,
    , drop = FALSE
  ]
  pts$x <- unname(x_lookup[pts$tissue_key])
  pts$display_index <- ave(
    seq_len(nrow(pts)), pts$tissue_key, FUN = seq_along
  )
  pts$x_jitter <- pts$x + 0.150 +
    0.038 * sin(pts$display_index * 2.399963)
  pts$sign <- ifelse(pts$rho < 0, "Negative", "Positive")

  meta$x_target <- meta$x + 0.225

  species_positions <- do.call(rbind, lapply(species_in_row, function(sp) {
    xx <- meta$x[meta$species == sp]
    data.frame(
      species = abbreviate_species(sp),
      x = mean(range(xx)),
      stringsAsFactors = FALSE
    )
  }))

  species_n <- as.integer(table(factor(meta$species, levels = species_in_row)))
  boundaries <- head(cumsum(species_n), -1L) + 0.5
  x_labels <- setNames(
    paste0(
      pretty_tissue_compact(meta$tissue),
      "\nn=", meta$n_samples
    ),
    meta$x
  )

  ggplot(dat, aes(x = x, y = rho, group = tissue_key)) +
    geom_hline(yintercept = 0, colour = "grey58", linewidth = 0.25) +
    geom_segment(
      data = data.frame(x = boundaries),
      aes(x = x, xend = x, y = -1, yend = 1.02),
      inherit.aes = FALSE,
      colour = "grey70",
      linewidth = 0.30
    ) +
    gghalves::geom_half_violin(
      side = "r",
      fill = NA,
      colour = "black",
      linetype = "dotted",
      linewidth = 0.25,
      trim = FALSE,
      scale = "width",
      width = 0.61
    ) +
    gghalves::geom_half_boxplot(
      side = "l",
      center = FALSE,
      errorbar.draw = TRUE,
      width = 0.27,
      outlier.shape = NA,
      fill = "grey15",
      colour = "black",
      linewidth = 0.29
    ) +
    geom_point(
      data = pts,
      aes(x = x_jitter, y = rho, colour = sign),
      inherit.aes = FALSE,
      size = 0.16,
      alpha = 0.72
    ) +
    geom_point(
      data = meta,
      aes(x = x_target, y = target_rho),
      inherit.aes = FALSE,
      shape = 21,
      fill = target_fill,
      colour = "black",
      stroke = 0.34,
      size = 1.42
    ) +
    geom_text(
      data = species_positions,
      aes(x = x, y = 1.105, label = species),
      inherit.aes = FALSE,
      family = "Helvetica",
      fontface = "italic",
      size = 2.05
    ) +
    scale_colour_manual(values = sign_colours, guide = "none") +
    scale_x_continuous(
      breaks = meta$x,
      labels = x_labels,
      limits = c(0.45, nrow(meta) + 0.55),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1", "-0.5", "0", "+0.5", "+1")
    ) +
    coord_cartesian(ylim = c(-1, 1), clip = "off") +
    labs(
      x = NULL,
      y = expression("CpG methylation ~ age (" * rho * ")")
    ) +
    theme_s2_defg +
    theme(
      axis.text.x = element_text(
        angle = 54, hjust = 1, vjust = 1,
        size = 4.80, lineheight = 0.88
      ),
      axis.text.y = element_text(size = 4.95),
      axis.title.y = element_text(size = 5.70, margin = margin(r = 2)),
      plot.margin = margin(8, 2, 1, 4)
    )
}

s2e_row_plots <- lapply(s2e_species_rows, make_s2e_row)
p_S2e <- cowplot::plot_grid(
  plotlist = s2e_row_plots,
  ncol = 1,
  rel_heights = c(0.96, 0.96, 1.08),
  align = "v",
  axis = "lr"
)

# -----------------------------------------------------------------------------
# S2f. Mean normalized beta signal by mapped-sequence status.
#
# Frozen upstream aggregation represented by each point:
#   1. retain confidence-filtered samples in tissue strata with n >= 10;
#   2. calculate the arithmetic mean cg10501210 beta across finite samples
#      separately within every eligible species-tissue stratum;
#   3. calculate an unweighted arithmetic mean of those tissue means within
#      species. Thus each species contributes one point and its tissues receive
#      equal weight, irrespective of tissue sample count.
# -----------------------------------------------------------------------------

s2f_beta <- read_frozen_tsv(s2_defg_paths[["S2f_beta_species"]])
s2f_beta$conservation_status <- factor(
  s2f_beta$conservation_status,
  levels = names(sequence_status_colours)
)
s2f_status_n <- table(s2f_beta$conservation_status)

stopifnot(
  nrow(s2f_beta) == 66L,
  identical(as.integer(s2f_status_n), c(41L, 25L)),
  all(is.finite(s2f_beta$beta_tissue_mean_to_species_mean)),
  all(s2f_beta$beta_tissue_mean_to_species_mean >= 0 &
        s2f_beta$beta_tissue_mean_to_species_mean <= 1)
)

s2_status_labels <- c(
  `Intact CpG` = "Intact CpG\nSp=41",
  `Disrupted/non-conserved` = "Disrupted\nSp=25"
)

p_S2f <- ggplot(
  s2f_beta,
  aes(x = conservation_status, y = beta_tissue_mean_to_species_mean)
) +
  geom_violin(
    aes(fill = conservation_status),
    alpha = 0.18,
    colour = "black",
    linewidth = 0.31,
    trim = TRUE,
    width = 0.82
  ) +
  geom_boxplot(
    width = 0.21,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.34
  ) +
  geom_point(
    aes(fill = conservation_status),
    shape = 21,
    colour = "black",
    stroke = 0.20,
    size = 0.92,
    alpha = 0.84,
    position = position_jitter(width = 0.115, height = 0, seed = 29024)
  ) +
  scale_fill_manual(values = sequence_status_colours, guide = "none") +
  scale_x_discrete(labels = s2_status_labels) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    x = NULL,
    y = "Normalized probe beta signal"
  ) +
  theme_s2_defg +
  theme(
    axis.text.x = element_text(size = 5.05, lineheight = 0.90),
    axis.text.y = element_text(size = 5.15),
    axis.title.y = element_text(size = 5.70, margin = margin(r = 2)),
    plot.margin = margin(3, 2, 2, 3)
  )

# -----------------------------------------------------------------------------
# S2g. Species-level age association by mapped-sequence status.
#
# Frozen upstream aggregation represented by each point:
#   1. calculate Spearman rho between cg10501210 beta and chronological age
#      separately within every eligible species-tissue stratum (n >= 10);
#   2. take the unweighted arithmetic mean of the finite tissue-level rho values
#      within species. No Fisher-z transformation is used at this tissue-to-
#      species step. Point size is the sum of samples across contributing tissue
#      strata (`target_effective_n`).
# -----------------------------------------------------------------------------

s2g_rho_all <- read_frozen_tsv(s2_defg_paths[["S2g_rho_species"]])
s2g_rho <- s2g_rho_all[
  s2g_rho_all$included & is.finite(s2g_rho_all$rho),
  , drop = FALSE
]
s2g_rho$conservation_status <- factor(
  s2g_rho$conservation_status,
  levels = names(sequence_status_colours)
)
s2g_status_n <- table(s2g_rho$conservation_status)

stopifnot(
  nrow(s2g_rho) == 66L,
  identical(as.integer(s2g_status_n), c(41L, 25L)),
  all(s2g_rho$rho >= -1 & s2g_rho$rho <= 1),
  all(s2g_rho$target_effective_n > 0)
)

p_S2g <- ggplot(s2g_rho, aes(x = conservation_status, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.28) +
  geom_violin(
    aes(fill = conservation_status),
    alpha = 0.18,
    colour = "black",
    linewidth = 0.31,
    trim = FALSE,
    width = 0.82
  ) +
  geom_boxplot(
    width = 0.21,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.34
  ) +
  geom_point(
    aes(
      fill = conservation_status,
      size = pmax(target_effective_n, 1)
    ),
    shape = 21,
    colour = "black",
    stroke = 0.19,
    alpha = 0.84,
    position = position_jitter(width = 0.115, height = 0, seed = 29012)
  ) +
  scale_fill_manual(values = sequence_status_colours, guide = "none") +
  scale_size_continuous(
    trans = "sqrt",
    range = c(0.62, 2.10),
    breaks = c(100, 1000, 3000),
    labels = c("100", "1,000", "3,000"),
    name = "samples"
  ) +
  scale_x_discrete(labels = s2_status_labels) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    labels = c("-1", "-0.5", "0", "+0.5", "+1")
  ) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    x = NULL,
    y = expression("Species aggregate Spearman " * rho)
  ) +
  guides(
    size = guide_legend(
      title.position = "top",
      override.aes = list(fill = "white")
    )
  ) +
  theme_s2_defg +
  theme(
    axis.text.x = element_text(size = 5.05, lineheight = 0.90),
    axis.text.y = element_text(size = 5.15),
    axis.title.y = element_text(size = 5.70, margin = margin(r = 2)),
    legend.position = "right",
    legend.title = element_text(size = 5.2),
    legend.text = element_text(size = 4.8),
    legend.key.height = unit(0.17, "in"),
    legend.key.width = unit(0.13, "in"),
    plot.margin = margin(3, 2, 2, 2)
  )

# Recommended native export geometry for the approved A4 S2 layout. Panel
# letters are deliberately absent and remain controlled by the assembly layer.
s2_defg_dimensions <- data.frame(
  panel = c("S2d", "S2e", "S2f", "S2g"),
  width_in = c(2.75, 5.15, 1.45, 2.25),
  height_in = c(5.25, 4.65, 1.70, 1.70),
  stringsAsFactors = FALSE
)

s2_defg_aggregation <- c(
  S2d = paste0(
    "Coverage n is the number of confidence/age-filtered profiles in ",
    "species-tissue strata with n >= 10; the tissue count uses the same rule."
  ),
  S2e = paste0(
    "Each tissue distribution contains all 37,554 frozen CpG-age Spearman ",
    "rho values; yellow is cg10501210. Human age-invariant strata are omitted."
  ),
  S2f = paste0(
    "For each species, cg10501210 beta was averaged across finite samples ",
    "within each eligible tissue and those tissue means were then averaged ",
    "arithmetically with equal tissue weight."
  ),
  S2g = paste0(
    "For each species, tissue-level cg10501210 Spearman rho values were ",
    "averaged arithmetically with equal tissue weight (no Fisher-z at the ",
    "tissue-to-species step)."
  )
)

s2_defg_data <- list(
  S2d_coverage = s2d_coverage,
  S2e_target_all = s2e_target_all,
  S2e_target_plotted = s2e_target,
  S2e_display_points = s2e_points,
  S2e_full_distributions = s2e_distribution,
  S2e_species_rows = s2e_species_rows,
  S2f_beta_species = s2f_beta,
  S2g_rho_species_all = s2g_rho_all,
  S2g_rho_species_plotted = s2g_rho,
  source_paths = s2_defg_paths,
  aggregation = s2_defg_aggregation
)

s2_defg_omitted_strata <- sort(paste(
  s2e_target_all$species[!s2e_target_all$informative],
  s2e_target_all$tissue[!s2e_target_all$informative],
  sep = "||"
))

s2_defg_checks <- data.frame(
  check = c(
    "all_frozen_inputs_present",
    "S2d_species_n",
    "S2d_macroclade_n",
    "S2d_eligible_profile_n",
    "S2d_eligible_tissue_strata_n",
    "S2e_species_n",
    "S2e_informative_tissue_strata_n",
    "S2e_omitted_strata_are_human_age_invariant",
    "S2e_full_rho_rows",
    "S2e_cpgs_per_tissue",
    "S2f_species_n",
    "S2f_intact_disrupted_n",
    "S2g_species_n",
    "S2g_intact_disrupted_n"
  ),
  status = c(
    all(file.exists(s2_defg_paths)),
    nrow(s2d_coverage) == 36L,
    nlevels(droplevels(s2d_coverage$clade)) == 5L,
    sum(s2d_coverage$n_eligible) == 5660L,
    sum(s2d_coverage$n_eligible_strata) == 67L,
    length(unique(s2e_target$species)) == 12L,
    nrow(s2e_target) == 41L,
    identical(
      s2_defg_omitted_strata,
      c("Homo sapiens||Keratinocyte", "Homo sapiens||Placenta")
    ),
    nrow(s2e_distribution) == probe_universe_n * 41L,
    length(unique(s2e_distribution$probe)) == probe_universe_n,
    nrow(s2f_beta) == 66L,
    identical(as.integer(s2f_status_n), c(41L, 25L)),
    nrow(s2g_rho) == 66L,
    identical(as.integer(s2g_status_n), c(41L, 25L))
  ),
  observed = c(
    paste(sum(file.exists(s2_defg_paths)), "of", length(s2_defg_paths)),
    as.character(nrow(s2d_coverage)),
    as.character(nlevels(droplevels(s2d_coverage$clade))),
    as.character(sum(s2d_coverage$n_eligible)),
    as.character(sum(s2d_coverage$n_eligible_strata)),
    as.character(length(unique(s2e_target$species))),
    as.character(nrow(s2e_target)),
    paste(s2_defg_omitted_strata, collapse = "; "),
    as.character(nrow(s2e_distribution)),
    as.character(length(unique(s2e_distribution$probe))),
    as.character(nrow(s2f_beta)),
    paste(as.integer(s2f_status_n), collapse = "/"),
    as.character(nrow(s2g_rho)),
    paste(as.integer(s2g_status_n), collapse = "/")
  ),
  expected = c(
    paste(length(s2_defg_paths), "of", length(s2_defg_paths)),
    "36", "5", "5660", "67", "12", "41",
    "Homo sapiens||Keratinocyte; Homo sapiens||Placenta",
    as.character(probe_universe_n * 41L),
    as.character(probe_universe_n),
    "66", "41/25", "66", "41/25"
  ),
  stringsAsFactors = FALSE
)

stopifnot(all(s2_defg_checks$status))

if (s2_export) {
  dir.create(s2_defg_export_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Package 'svglite' is required when s2_export is TRUE.")
  }

  s2_plots_for_export <- list(
    S2d = p_S2d,
    S2e = p_S2e,
    S2f = p_S2f,
    S2g = p_S2g
  )

  for (i in seq_len(nrow(s2_defg_dimensions))) {
    stem <- s2_defg_dimensions$panel[[i]]
    plot_i <- s2_plots_for_export[[stem]]
    width_i <- s2_defg_dimensions$width_in[[i]]
    height_i <- s2_defg_dimensions$height_in[[i]]

    ggsave(
      file.path(s2_defg_export_dir, paste0(stem, ".pdf")),
      plot_i,
      width = width_i,
      height = height_i,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
    ggsave(
      file.path(s2_defg_export_dir, paste0(stem, ".svg")),
      plot_i,
      width = width_i,
      height = height_i,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
    ggsave(
      file.path(s2_defg_export_dir, paste0(stem, ".png")),
      plot_i,
      width = width_i,
      height = height_i,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
}
