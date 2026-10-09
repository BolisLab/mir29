###############################################################################
# MIR29 Supplementary Figure S1
#
# Reconstructs the six original data-driven panels, now S1a, b, f, g, h, i,
# from unchanged frozen inputs. The relocated S1c, d, e are built in the
# corresponding sections of Figure1_and_S1.Rmd. The complete notebook then
# calls assemble_figure1_s1.R to export the current nine-panel supplementary.
# No raw GTEx matrix or donor-bootstrap refitting is needed for this helper.
###############################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(scales)
  library(grid)
})

s1_project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(s1_project_dir, "data", "figure1_s1"))) {
  stop(
    "Set the RStudio working directory to the folder containing Figure1_and_S1.Rmd ",
    "before running this section."
  )
}

s1_data_dir <- file.path(s1_project_dir, "data", "figure1_s1")
s1_output_dir <- file.path(s1_project_dir, "outputs", "supplementary_s1")
s1_asset_dir <- file.path(s1_data_dir, "assets")
dir.create(s1_output_dir, recursive = TRUE, showWarnings = FALSE)

s1_read_tsv <- function(name) {
  path <- file.path(s1_data_dir, name)
  if (!file.exists(path)) stop("Missing Supplementary Figure S1 input: ", path)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"")
}

s1_atlas <- s1_read_tsv("gtex_30tissue_rna_sample_donor_atlas.tsv")
s1_funnel <- s1_read_tsv("gtex_7tissue_multiomic_funnel.tsv")
s1_age30 <- s1_read_tsv(
  "gtex_LINC_PINT_30tissue_donor_aware_age_sensitivity.tsv"
)
s1_age30_meta <- s1_read_tsv(
  "gtex_LINC_PINT_30tissue_REML_cluster_bootstrap_summary.tsv"
)
s1_paired <- s1_read_tsv(
  "gtex_LINC_PINT_paired_partial_spearman_tissue.tsv"
)
s1_paired_meta <- s1_read_tsv(
  "gtex_LINC_PINT_paired_partial_spearman_REML_meta.tsv"
)
s1_upset_membership <- s1_read_tsv("rz_upset_gene_membership.tsv")
s1_upset_counts <- s1_read_tsv("rz_upset_intersection_counts.tsv")
s1_arz_scores <- s1_read_tsv("ARZ_multiset_score_heatmap_source_data.tsv")
s1_mir29a <- s1_read_tsv("gtex_mir29a_tissue_associations_BH30.tsv")
s1_mir29_meta <- s1_read_tsv(
  "gtex_host_mature_mir29_pan_tissue_summaries.tsv"
)
s1_method4_targets <- s1_read_tsv(
  "mirna_method4_target_tissue_positions.tsv"
)
s1_method4 <- readRDS(file.path(s1_data_dir, "mirna_method4_matrix.rds"))

# Fail loudly if a different analysis snapshot is substituted by accident.
s1_expected_patterns <- c(
  "A only" = 174L, "R only" = 522L, "Z only" = 168L,
  "A∩R" = 3L, "A∩Z" = 1L, "R∩Z" = 5L, "A∩R∩Z" = 1L
)
s1_observed_patterns <- setNames(
  as.integer(s1_upset_counts$n), s1_upset_counts$pattern
)[names(s1_expected_patterns)]

s1_paired_tissues_expected <- c(
  "Breast", "Colon", "Lung", "Muscle", "Ovary", "Prostate", "Testis"
)
s1_paired_predictors_expected <- c("AGE", "DNAmAge", "AgeAccel")
s1_paired_grid_observed <- with(
  s1_paired,
  sort(paste(tissue, predictor, sep = "::"))
)
s1_paired_grid_expected <- sort(as.vector(outer(
  s1_paired_tissues_expected,
  s1_paired_predictors_expected,
  paste,
  sep = "::"
)))
s1_paired_unique_n <- unique(s1_paired[c("tissue", "n")])

stopifnot(
  nrow(s1_atlas) == 30L,
  sum(s1_atlas$n_rna_sample_rows) == 19136L,
  length(unique(s1_atlas$tissue)) == 30L,
  sum(s1_atlas$selected_for_7tissue_dnam_analysis) == 7L,
  identical(as.integer(s1_funnel$sample_records), c(3750L, 883L, 551L)),
  identical(as.integer(s1_funnel$unique_donors_across_tissues), c(927L, 423L, 348L)),
  nrow(s1_age30) == 30L,
  nrow(s1_age30_meta) == 2L,
  nrow(s1_paired) == 21L,
  nrow(s1_paired_meta) == 3L,
  identical(s1_paired_grid_observed, s1_paired_grid_expected),
  identical(sort(unique(s1_paired_meta$predictor)), sort(s1_paired_predictors_expected)),
  sum(s1_paired_unique_n$n) == 551L,
  all(s1_paired$bootstrap_replicates == 5000L),
  all(grepl("REML", s1_paired_meta$meta_method, fixed = TRUE)),
  all(s1_paired_meta$ci_low < s1_paired_meta$meta_rho),
  all(s1_paired_meta$meta_rho < s1_paired_meta$ci_high),
  all(s1_paired_meta$ci_low > 0),
  s1_paired_meta$meta_rho[s1_paired_meta$predictor == "DNAmAge"] >
    s1_paired_meta$meta_rho[s1_paired_meta$predictor == "AGE"],
  identical(unname(s1_observed_patterns), unname(s1_expected_patterns)),
  length(unique(s1_arz_scores$gene_id)) == 10L,
  nrow(s1_arz_scores) == 30L,
  sum(s1_upset_membership$pattern == "A∩R∩Z") == 1L,
  s1_upset_membership$symbol[s1_upset_membership$pattern == "A∩R∩Z"] ==
    "LINC-PINT",
  nrow(subset(s1_mir29a, host == "LINC-PINT")) == 30L,
  sum(subset(s1_mir29a, host == "LINC-PINT")$mir29a_tissue_BH_FDR < 0.05) == 23L,
  identical(dim(s1_method4), c(884L, 30L))
)

s1_colors <- c(
  age = "#E72F37",
  rejuvenation = "#A44297",
  zygotic = "#3559A7",
  sample = "#A8BBC6",
  donor = "#3E586B",
  dnam = "#B5152B",
  pair = "#F2D300",
  mir29a = "#F1C51D",
  mir29b = "#8DCCD7",
  mir29c = "#EF6476",
  positive_bg = "#FCF2F3",
  negative_bg = "#EDF6F8"
)

s1_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "#202020"),
      axis.title = element_text(color = "#202020"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(4, 4, 4, 4)
    )
}

s1_add_optional_asset <- function(plot, filename, x, y, width, height) {
  asset_path <- file.path(s1_asset_dir, filename)
  if (!file.exists(asset_path)) return(plot)
  image_data <- png::readPNG(asset_path, native = TRUE)
  image_grob <- rasterGrob(image_data, interpolate = TRUE)
  ggdraw(plot) + draw_grob(image_grob, x = x, y = y, width = width, height = height)
}

s1_draw_png <- function(filename, x, y, width, height) {
  asset_path <- file.path(s1_asset_dir, filename)
  image_data <- png::readPNG(asset_path, native = TRUE)
  draw_grob(
    rasterGrob(image_data, interpolate = TRUE),
    x = x, y = y, width = width, height = height
  )
}

###############################################################################
# Panel a: complete 30-tissue RNA sample/donor atlas and modality tracks
###############################################################################

s1_atlas <- s1_atlas %>% arrange(tissue)
s1_atlas$tissue_factor <- factor(s1_atlas$tissue, levels = s1_atlas$tissue)

s1_atlas_bars <- s1_atlas %>%
  select(tissue_factor, n_rna_sample_rows, n_unique_donors) %>%
  pivot_longer(
    cols = c(n_rna_sample_rows, n_unique_donors),
    names_to = "measure", values_to = "count"
  ) %>%
  mutate(
    measure = factor(
      measure,
      levels = c("n_rna_sample_rows", "n_unique_donors"),
      labels = c("RNA-seq samples", "Unique donors")
    )
  )

s1_panel_a_bars <- ggplot(
  s1_atlas_bars,
  aes(x = tissue_factor, y = count, fill = measure)
) +
  geom_hline(
    yintercept = c(500, 1500, 2500),
    linetype = "dotted", linewidth = 0.35, color = "grey25"
  ) +
  geom_col(position = position_dodge(width = 0.78), width = 0.74) +
  scale_fill_manual(
    values = c(
      "RNA-seq samples" = s1_colors[["sample"]],
      "Unique donors" = s1_colors[["donor"]]
    ),
    guide = guide_legend(
      title = NULL, direction = "vertical", byrow = TRUE,
      keyheight = unit(3.0, "mm"), keywidth = unit(2.5, "mm")
    )
  ) +
  scale_y_continuous(
    breaks = c(0, 500, 1500, 2500),
    labels = label_number(big.mark = ","),
    limits = c(0, 3350), expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "number of samples and unique donors") +
  s1_theme(8.6) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 8.2),
    legend.position = "inside", legend.position.inside = c(0.015, 0.94),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 8.2),
    legend.key.spacing.y = unit(0, "pt"),
    plot.margin = margin(3, 5, 0, 13)
  )

s1_availability <- bind_rows(
  data.frame(
    tissue_factor = s1_atlas$tissue_factor,
    track = "RNA-Seq", status = "RNA available", stringsAsFactors = FALSE
  ),
  data.frame(
    tissue_factor = s1_atlas$tissue_factor,
    track = "DNA-Methylation",
    status = ifelse(
      s1_atlas$selected_for_7tissue_dnam_analysis,
      "Selected DNAm tissue",
      ifelse(s1_atlas$has_dnam_clock_profiles, "Other DNAm tissue", NA_character_)
    ),
    stringsAsFactors = FALSE
  )
) %>%
  filter(!is.na(status)) %>%
  mutate(
    track = factor(track, levels = c("RNA-Seq", "DNA-Methylation")),
    status = factor(
      status,
      levels = c("RNA available", "Other DNAm tissue", "Selected DNAm tissue")
    )
  )

s1_panel_a_tracks <- ggplot(
  s1_availability,
  aes(x = tissue_factor, y = track, fill = status)
) +
  geom_hline(yintercept = c(1, 2), linewidth = 0.25, color = "grey75") +
  geom_tile(width = 0.34, height = 0.45, color = NA) +
  scale_fill_manual(
    values = c(
      "RNA available" = "#202020",
      "Other DNAm tissue" = "#EDC5D2",
      "Selected DNAm tissue" = s1_colors[["dnam"]]
    ), guide = "none"
  ) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  s1_theme(8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_text(size = 5.4, hjust = 0),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
    panel.grid = element_blank(),
    plot.margin = margin(0, 5, 2, 13)
  )

s1_panel_a_core <- plot_grid(
  s1_panel_a_bars, s1_panel_a_tracks,
  ncol = 1, rel_heights = c(3.8, 1.4), align = "v", axis = "lr"
)
s1_panel_a <- s1_add_optional_asset(
  s1_panel_a_core, "render_all_tissues_wheel.png",
  x = 0.78, y = 0.71, width = 0.19, height = 0.26
)

###############################################################################
# Panel b: seven-tissue modality availability and strict sample pairing
###############################################################################

s1_funnel <- s1_funnel %>%
  arrange(stage_order) %>%
  mutate(
    label = c("RNA-seq samples", "DNAm profiles", "RNA-DNAm pairs"),
    label = factor(label, levels = rev(label)),
    value_label = c(
      "3750 samples / 927 donors",
      "883 samples / 423 donors",
      "551 samples / 348 donors"
    )
  )

s1_panel_b_top <- ggplot(
  s1_funnel,
  aes(x = sample_records, y = label, fill = label)
) +
  geom_col(width = 0.50) +
  geom_text(
    data = filter(s1_funnel, stage_order == 1),
    aes(label = value_label), hjust = 1.03, size = 2.75
  ) +
  geom_text(
    data = filter(s1_funnel, stage_order != 1),
    aes(label = value_label), hjust = -0.03, size = 2.75
  ) +
  scale_fill_manual(
    values = c(
      "RNA-seq samples" = s1_colors[["sample"]],
      "DNAm profiles" = s1_colors[["dnam"]],
      "RNA-DNAm pairs" = s1_colors[["pair"]]
    ), guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 4850), breaks = c(0, 1000, 2000, 3000),
    labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL, title = "tissues profiled for DNAm, RNA (N = 7)") +
  s1_theme(8.8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8.8),
    plot.title = element_text(size = 7.2, face = "plain", hjust = 0),
    plot.margin = margin(2, 3, 0, 3)
  )

s1_selected_atlas <- s1_atlas %>%
  filter(selected_for_7tissue_dnam_analysis) %>%
  arrange(match(tissue, c("Breast", "Colon", "Lung", "Muscle", "Ovary", "Prostate", "Testis")))

# The first lower-panel count is not the complete RNA atlas. It is the set of
# RNA biospecimen cores that match a DNAm core in any tissue before enforcing
# same-tissue pairing; the explicit label prevents it being mistaken for 3,750.
s1_selected_long <- s1_selected_atlas %>%
  select(
    tissue,
    `RNA core matches` = n_rna_core_matches_any_tissue,
    `DNAm samples` = n_dnam_clock_profiles,
    `RNA-DNAm pairs` = n_strict_same_biospecimen_paired
  ) %>%
  pivot_longer(-tissue, names_to = "measure", values_to = "count") %>%
  mutate(
    tissue = factor(
      tissue,
      levels = rev(c("Breast", "Colon", "Lung", "Muscle", "Ovary", "Prostate", "Testis"))
    ),
    measure = factor(
      measure,
      levels = c("RNA core matches", "DNAm samples", "RNA-DNAm pairs")
    )
  )

s1_panel_b_bottom <- ggplot(
  s1_selected_long,
  aes(x = count, y = tissue, fill = measure)
) +
  geom_col(position = position_dodge(width = 0.70), width = 0.23) +
  scale_fill_manual(
    values = c(
      "RNA core matches" = s1_colors[["sample"]],
      "DNAm samples" = s1_colors[["dnam"]],
      "RNA-DNAm pairs" = s1_colors[["pair"]]
    ),
    guide = guide_legend(title = NULL, nrow = 1, byrow = TRUE)
  ) +
  scale_x_continuous(
    limits = c(0, 250), breaks = seq(0, 200, 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.35, 1.15))) +
  labs(x = "number of samples", y = NULL) +
  s1_theme(8.6) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 8.5),
    legend.position = "bottom", legend.text = element_text(size = 4.6),
    legend.key.width = unit(2.6, "mm"), legend.key.height = unit(2.2, "mm"),
    plot.margin = margin(8, 3, 0, 3)
  )

if (all(file.exists(file.path(
  s1_asset_dir,
  c(
    "render_breast.png", "render_colon.png", "render_lung.png",
    "render_skeletal_muscle.png", "render_skeletal_muscle_bone.png",
    "render_ovary_uterus.png", "render_prostate.png", "render_testis.png"
  )
)))) {
  s1_panel_b_bottom <- ggdraw(s1_panel_b_bottom) +
    s1_draw_png("render_breast.png", x = 0.145, y = 0.855, width = 0.052, height = 0.085) +
    s1_draw_png("render_colon.png", x = 0.262, y = 0.850, width = 0.052, height = 0.090) +
    s1_draw_png("render_lung.png", x = 0.384, y = 0.855, width = 0.056, height = 0.085) +
    s1_draw_png("render_skeletal_muscle_bone.png", x = 0.500, y = 0.855, width = 0.058, height = 0.085) +
    s1_draw_png("render_skeletal_muscle.png", x = 0.510, y = 0.860, width = 0.044, height = 0.073) +
    s1_draw_png("render_ovary_uterus.png", x = 0.622, y = 0.850, width = 0.055, height = 0.090) +
    s1_draw_png("render_prostate.png", x = 0.741, y = 0.850, width = 0.055, height = 0.090) +
    s1_draw_png("render_testis.png", x = 0.858, y = 0.850, width = 0.055, height = 0.090)
}
s1_panel_b <- plot_grid(
  s1_panel_b_top, s1_panel_b_bottom,
  ncol = 1, rel_heights = c(0.72, 1.35), align = "v", axis = "lr"
)

###############################################################################
# Panel f: 30-tissue LINC-PINT chronological-age associations
###############################################################################

s1_age_tissues <- sort(s1_age30$tissue)
s1_age30_long <- bind_rows(
  s1_age30 %>% transmute(
    label = tissue, estimator = "all biospecimens",
    rho = sample_row_partial_rho,
    ci_low = sample_row_partial_cluster_boot_ci_low,
    ci_high = sample_row_partial_cluster_boot_ci_high
  ),
  s1_age30 %>% transmute(
    label = tissue, estimator = "donor aggregated",
    rho = donor_aggregate_partial_rho,
    ci_low = donor_aggregate_partial_boot_ci_low,
    ci_high = donor_aggregate_partial_boot_ci_high
  )
) %>%
  mutate(
    estimator = factor(
      estimator, levels = c("all biospecimens", "donor aggregated")
    ),
    x_base = match(label, s1_age_tissues),
    x = x_base + ifelse(estimator == "all biospecimens", -0.10, 0.10)
  )

s1_age30_meta_plot <- s1_age30_meta %>%
  transmute(
    estimator = factor(
      ifelse(estimator == "Sample rows", "all biospecimens", "donor aggregated"),
      levels = c("all biospecimens", "donor aggregated")
    ),
    rho = meta_rho, ci_low, ci_high,
    x = length(s1_age_tissues) + ifelse(
      estimator == "all biospecimens", 0.86, 1.14
    )
  )

s1_age30_diamonds <- bind_rows(lapply(seq_len(nrow(s1_age30_meta_plot)), function(i) {
  row <- s1_age30_meta_plot[i, ]
  data.frame(
    estimator = row$estimator,
    group = i,
    x = c(row$x, row$x - 0.18, row$x, row$x + 0.18),
    y = c(row$ci_low, row$rho, row$ci_high, row$rho)
  )
}))

s1_panel_f <- ggplot() +
  annotate(
    "rect", xmin = 0.35, xmax = length(s1_age_tissues) + 1.55,
    ymin = -0.50, ymax = 0, fill = s1_colors[["negative_bg"]]
  ) +
  annotate(
    "rect", xmin = 0.35, xmax = length(s1_age_tissues) + 1.55,
    ymin = 0, ymax = 0.70, fill = s1_colors[["positive_bg"]]
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.35) +
  geom_linerange(
    data = s1_age30_long,
    aes(x = x, ymin = ci_low, ymax = ci_high, color = estimator),
    linewidth = 0.40
  ) +
  geom_point(
    data = s1_age30_long,
    aes(x = x, y = rho, color = estimator, shape = estimator),
    size = 1.25
  ) +
  geom_polygon(
    data = s1_age30_diamonds,
    aes(x = x, y = y, group = group, fill = estimator),
    color = NA
  ) +
  scale_color_manual(
    values = c(
      "all biospecimens" = s1_colors[["sample"]],
      "donor aggregated" = s1_colors[["donor"]]
    ), name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "all biospecimens" = s1_colors[["sample"]],
      "donor aggregated" = s1_colors[["donor"]]
    ), guide = "none"
  ) +
  scale_shape_manual(values = c("all biospecimens" = 16, "donor aggregated" = 15)) +
  scale_x_continuous(
    breaks = seq_len(length(s1_age_tissues) + 1L),
    labels = c(s1_age_tissues, "meta"),
    limits = c(0.35, length(s1_age_tissues) + 1.55),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    position = "right", breaks = c(-0.5, 0, 0.5),
    limits = c(-0.50, 0.70), expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL, y = "partial Spearman correlation (rho)",
    title = "LINC-PINT ~ AGE (all RNA-seq samples)"
  ) +
  s1_theme(8.1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.8),
    axis.title.y.right = element_text(size = 8.5),
    plot.title = element_text(size = 9.5, face = "plain"),
    legend.position = "inside", legend.position.inside = c(0.91, 0.08),
    legend.justification = c(1, 0),
    legend.text = element_text(size = 7.0),
    legend.key.height = unit(2.5, "mm"),
    plot.margin = margin(3, 8, 3, 4)
  ) +
  guides(shape = "none", color = guide_legend(override.aes = list(shape = c(16, 15))))

###############################################################################
# Panel g: LINC-PINT with age, DNAmAge, and age acceleration in matched data
###############################################################################

s1_paired_tissues <- c("Breast", "Colon", "Lung", "Muscle", "Ovary", "Prostate", "Testis")
s1_predictors <- c("AGE", "DNAmAge", "AgeAccel")
s1_paired_plot <- s1_paired %>%
  filter(predictor %in% s1_predictors) %>%
  mutate(
    predictor = factor(predictor, levels = s1_predictors),
    r = rho,
    y = match(tissue, rev(s1_paired_tissues)) + 1
  )
s1_paired_meta_plot <- s1_paired_meta %>%
  filter(predictor %in% s1_predictors) %>%
  mutate(
    predictor = factor(predictor, levels = s1_predictors),
    r = meta_rho,
    y = 1
  )

s1_predictor_colors <- c(
  AGE = "#6E6E6E", DNAmAge = "#C4312E", AgeAccel = "#4E69C3"
)
s1_paired_diamonds <- bind_rows(lapply(seq_len(nrow(s1_paired_meta_plot)), function(i) {
  row <- s1_paired_meta_plot[i, ]
  data.frame(
    predictor = row$predictor, group = i,
    x = c(row$ci_low, row$r, row$ci_high, row$r),
    y = c(1, 1.21, 1, 0.79)
  )
}))

s1_paired_ranges <- bind_rows(
  s1_paired_plot %>% select(predictor, ci_low, ci_high),
  s1_paired_meta_plot %>% select(predictor, ci_low, ci_high)
) %>%
  group_by(predictor) %>%
  summarise(
    min_x = min(-0.25, min(ci_low, na.rm = TRUE) - 0.04),
    max_x = max(0.70, max(ci_high, na.rm = TRUE) + 0.04),
    .groups = "drop"
  )
s1_paired_background <- bind_rows(
  s1_paired_ranges %>% transmute(predictor, xmin = min_x, xmax = 0, fill = "negative"),
  s1_paired_ranges %>% transmute(predictor, xmin = 0, xmax = max_x, fill = "positive")
)

s1_panel_g <- ggplot() +
  geom_rect(
    data = s1_paired_background,
    aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = 8.55, fill = fill),
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6E9AAD", linewidth = 0.35) +
  geom_hline(yintercept = 1.5, linewidth = 0.55, color = "black") +
  geom_errorbarh(
    data = s1_paired_plot,
    aes(y = y, xmin = ci_low, xmax = ci_high),
    height = 0, linewidth = 0.45, color = "#202020"
  ) +
  geom_point(
    data = s1_paired_plot,
    aes(x = r, y = y, color = predictor), size = 1.75
  ) +
  geom_polygon(
    data = s1_paired_diamonds,
    aes(x = x, y = y, group = group, fill = predictor), color = NA
  ) +
  facet_wrap(
    ~ predictor, nrow = 1, scales = "free_x",
    labeller = as_labeller(c(
      AGE = "LINC-PINT ~ AGE",
      DNAmAge = "LINC-PINT ~ DNAmAge",
      AgeAccel = "LINC-PINT ~ AgeAcc"
    ))
  ) +
  scale_fill_manual(
    values = c(
      negative = s1_colors[["negative_bg"]],
      positive = s1_colors[["positive_bg"]],
      AGE = s1_predictor_colors[["AGE"]],
      DNAmAge = s1_predictor_colors[["DNAmAge"]],
      AgeAccel = s1_predictor_colors[["AgeAccel"]]
    ), guide = "none"
  ) +
  scale_color_manual(values = s1_predictor_colors, guide = "none") +
  scale_y_continuous(
    breaks = 1:8,
    labels = c("meta", rev(s1_paired_tissues)),
    limits = c(0.5, 8.55), expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "partial Spearman's correlation (rho)", y = NULL) +
  s1_theme(8.5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 7.3, face = "plain"),
    axis.text.y = element_text(size = 8.3),
    axis.title.x = element_text(size = 8.5),
    panel.spacing.x = unit(2.0, "mm"),
    plot.margin = margin(4, 4, 2, 3)
  )

###############################################################################
# Panel h: exact fixed-set intersections and raw A/R/Z statistics
###############################################################################

s1_pattern_order <- c(
  "A only", "R only", "Z only", "A∩R", "A∩Z", "R∩Z", "A∩R∩Z"
)
s1_pattern_label_code <- c(
  "A only" = "A~only",
  "R only" = "R~only",
  "Z only" = "Z~only",
  "A∩R" = "A %intersection% R",
  "A∩Z" = "A %intersection% Z",
  "R∩Z" = "R %intersection% Z",
  "A∩R∩Z" = "A %intersection% R %intersection% Z"
)
s1_signature_colors <- c(
  A = s1_colors[["age"]],
  R = s1_colors[["rejuvenation"]],
  Z = s1_colors[["zygotic"]]
)
s1_pattern_colors <- c(
  "A only" = s1_signature_colors[["A"]],
  "R only" = s1_signature_colors[["R"]],
  "Z only" = s1_signature_colors[["Z"]],
  "A∩R" = "#E2E2E2", "A∩Z" = "#E2E2E2", "R∩Z" = "#E2E2E2",
  "A∩R∩Z" = "#263545"
)

s1_upset_counts_plot <- s1_upset_counts %>%
  mutate(pattern = factor(pattern, levels = s1_pattern_order)) %>%
  arrange(pattern)

s1_small_labels <- data.frame(
  pattern = factor(c("A∩R", "A∩Z", "R∩Z", "A∩R∩Z"), levels = s1_pattern_order),
  label = c(
    "CLEC2B, STC1, RASSF8-AS1",
    "USP50",
    "INPP4B, BICC1, DEPDC7, SLC41A2, SCML1",
    "LINC-PINT"
  ),
  stringsAsFactors = FALSE
)

# RASSF8-AS1 is a display-only alias present in the final Illustrator artwork.
# The frozen manuscript annotation table has no symbol for ENSG00000246695; the
# stable Ensembl identifier is retained in source data and this override is
# therefore explicit rather than silently changing the analytical identifier.
s1_annotation_override <- data.frame(
  gene_id = "ENSG00000246695",
  frozen_symbol = NA_character_,
  final_display_label = "RASSF8-AS1",
  stringsAsFactors = FALSE
)

s1_panel_h_bar <- ggplot(
  s1_upset_counts_plot,
  aes(x = pattern, y = n, fill = pattern)
) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", linewidth = 0.35) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.25, size = 2.7) +
  geom_text(
    data = s1_small_labels,
    aes(x = pattern, y = 8, label = label),
    inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.3,
    size = 2.45, color = "#202020"
  ) +
  scale_fill_manual(values = s1_pattern_colors, guide = "none") +
  scale_y_sqrt(
    breaks = c(0, 5, 25, 100, 500),
    limits = c(0, 650), expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = NULL, y = "Genes (square-root scale)") +
  s1_theme(8.4) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8.6),
    plot.margin = margin(4, 3, 0, 4)
  )

s1_upset_matrix <- expand.grid(
  pattern = s1_pattern_order,
  signature = c("Z", "R", "A"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    present = mapply(function(pattern, signature) {
      if (pattern == "A only") return(signature == "A")
      if (pattern == "R only") return(signature == "R")
      if (pattern == "Z only") return(signature == "Z")
      signature %in% strsplit(pattern, "∩", fixed = TRUE)[[1]]
    }, pattern, signature),
    pattern = factor(pattern, levels = s1_pattern_order),
    signature = factor(signature, levels = c("Z", "R", "A"))
  )
s1_upset_lines <- s1_upset_matrix %>% filter(present) %>% group_by(pattern) %>%
  filter(n() > 1L) %>% ungroup()

s1_panel_h_matrix <- ggplot(s1_upset_matrix, aes(pattern, signature)) +
  geom_point(color = "grey82", size = 2.0) +
  geom_line(
    data = s1_upset_lines,
    aes(group = pattern), color = "grey35", linewidth = 0.45
  ) +
  geom_point(
    data = filter(s1_upset_matrix, present),
    aes(color = signature), size = 2.3
  ) +
  scale_color_manual(values = s1_signature_colors, guide = "none") +
  scale_x_discrete(labels = function(x) parse(text = unname(s1_pattern_label_code[x]))) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8.2, base_family = "sans") +
  theme(
    axis.text.x = element_text(
      angle = -48, hjust = 0, vjust = 1, size = 6.2, color = "#202020"
    ),
    axis.text.y = element_text(face = "bold", size = 7.5, color = "#202020"),
    plot.margin = margin(0, 3, 3, 4)
  )

s1_panel_h_upset <- plot_grid(
  s1_panel_h_bar, s1_panel_h_matrix,
  ncol = 1, rel_heights = c(3.3, 1.1), align = "v", axis = "lr"
)

s1_gene_order <- c(
  "CLEC2B", "STC1", "RASSF8-AS1", "USP50", "INPP4B",
  "BICC1", "DEPDC7", "SLC41A2", "SCML1", "LINC-PINT"
)
s1_arz_plot <- s1_arz_scores %>%
  mutate(
    display_gene = ifelse(gene_id == "ENSG00000246695", "RASSF8-AS1", display_gene),
    display_gene = factor(display_gene, levels = rev(s1_gene_order)),
    signature = factor(signature, levels = c("A", "R", "Z"))
  ) %>%
  group_by(signature) %>%
  mutate(
    symmetric_limit = max(abs(raw_score), na.rm = TRUE),
    display_score = raw_score / symmetric_limit,
    membership_alpha = ifelse(fixed_set_member, 1, 0.16)
  ) %>%
  ungroup()

s1_break_positions <- c(1.5, 6.5, 7.5)
s1_panel_h_heatmap <- ggplot(
  s1_arz_plot,
  aes(x = signature, y = display_gene)
) +
  geom_tile(
    aes(fill = display_score, alpha = membership_alpha),
    width = 0.92, height = 0.92, color = "grey88", linewidth = 0.25
  ) +
  geom_hline(
    yintercept = s1_break_positions,
    linetype = "dotted", color = "grey30", linewidth = 0.35
  ) +
  geom_tile(
    data = filter(s1_arz_plot, linc_pint),
    fill = NA, color = "black", linewidth = 0.7,
    width = 0.92, height = 0.92
  ) +
  scale_fill_gradient2(
    low = "#3559A7", mid = "white", high = "#E54853",
    midpoint = 0, limits = c(-1, 1), oob = squish,
    name = NULL
  ) +
  scale_alpha_identity() +
  scale_x_discrete(
    position = "top", labels = c(A = "A", R = "R", Z = "Z")
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8.3, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 9.2, color = "white"),
    axis.text.y = element_text(size = 8.0, color = "#202020"),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(18, "mm"), legend.key.height = unit(2.2, "mm"),
    legend.text = element_text(size = 6.5),
    plot.margin = margin(12, 1, 0, 2)
  ) +
  guides(fill = guide_colorbar(
    barwidth = unit(15, "mm"), barheight = unit(1.8, "mm"),
    ticks = FALSE, label = FALSE
  ))

s1_membership_gene <- s1_arz_plot %>%
  distinct(display_gene, gene_id, linc_pint) %>%
  left_join(
    s1_upset_membership %>% select(gene_id, A, R, Z),
    by = "gene_id"
  )
s1_membership_long <- s1_membership_gene %>%
  pivot_longer(c(A, R, Z), names_to = "signature", values_to = "present") %>%
  mutate(signature = factor(signature, levels = c("A", "R", "Z")))
s1_membership_lines <- s1_membership_long %>%
  filter(present) %>% group_by(display_gene) %>% filter(n() > 1) %>% ungroup()

s1_panel_h_membership <- ggplot(
  s1_membership_long,
  aes(x = signature, y = display_gene)
) +
  geom_point(color = "grey80", size = 1.5) +
  geom_line(
    data = s1_membership_lines,
    aes(group = display_gene), color = "grey35", linewidth = 0.35
  ) +
  geom_point(
    data = filter(s1_membership_long, present),
    aes(color = signature), size = 1.9
  ) +
  geom_hline(
    yintercept = s1_break_positions,
    linetype = "dotted", color = "grey30", linewidth = 0.35
  ) +
  scale_color_manual(values = s1_signature_colors, guide = "none") +
  scale_x_discrete(position = "top", labels = c(A = "", R = "", Z = "")) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8.2) +
  theme(plot.margin = margin(12, 2, 17, 0))

s1_panel_h_right <- plot_grid(
  s1_panel_h_heatmap, s1_panel_h_membership,
  nrow = 1, rel_widths = c(2.0, 0.80), align = "h", axis = "tb"
)
s1_panel_h <- plot_grid(
  s1_panel_h_upset, s1_panel_h_right,
  nrow = 1, rel_widths = c(1.15, 1.00), align = "h"
)

###############################################################################
# Panel i: LINC-PINT/miR-29a coupling and mature miR-29 age distributions
###############################################################################

s1_mir29a_linc <- s1_mir29a %>%
  filter(host == "LINC-PINT", mirna == "miR-29a-3p") %>%
  arrange(tissue) %>%
  mutate(
    tissue = factor(tissue, levels = rev(sort(unique(tissue)))),
    significance = factor(
      ifelse(mir29a_tissue_BH_FDR < 0.05, "BH FDR < 0.05", "BH FDR >= 0.05"),
      levels = c("BH FDR < 0.05", "BH FDR >= 0.05")
    )
  )
s1_host_limit <- max(abs(s1_mir29a_linc$primary_partial_rho), na.rm = TRUE)
s1_host_limit <- max(0.60, ceiling(s1_host_limit * 10) / 10)

s1_panel_i_host <- ggplot(
  s1_mir29a_linc,
  aes(x = 1, y = tissue)
) +
  geom_point(
    aes(
      fill = primary_partial_rho,
      size = abs(primary_partial_rho),
      alpha = significance
    ),
    shape = 21, color = "grey25", stroke = 0.28
  ) +
  scale_fill_gradient2(
    low = "#3C67A7", mid = "white", high = "#D55D67",
    midpoint = 0, limits = c(-s1_host_limit, s1_host_limit),
    oob = squish, name = "partial Spearman"
  ) +
  scale_size_continuous(range = c(0.8, 4.0), limits = c(0, s1_host_limit), guide = "none") +
  scale_alpha_manual(
    values = c("BH FDR < 0.05" = 1, "BH FDR >= 0.05" = 0.20),
    guide = "none"
  ) +
  scale_x_continuous(breaks = 1, labels = "") +
  labs(x = "miR-29a ~ LINC-PINT", y = NULL) +
  s1_theme(7.0) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_text(size = 4.7),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 5.2),
    legend.position = "top",
    legend.title = element_text(size = 5.4, face = "plain"),
    legend.text = element_text(size = 4.8),
    legend.key.height = unit(2.3, "mm"),
    legend.key.width = unit(8, "mm"),
    plot.margin = margin(0, 0, 2, 1)
  )

s1_method4_long <- as.data.frame(s1_method4, check.names = FALSE) %>%
  mutate(feature_id = rownames(s1_method4)) %>%
  pivot_longer(-feature_id, names_to = "tissue", values_to = "coefficient") %>%
  filter(is.finite(coefficient))

s1_mir29_labels <- c("miR-29a", "miR-29b", "miR-29c")
s1_mir29_colors <- c(
  "miR-29a" = s1_colors[["mir29a"]],
  "miR-29b" = s1_colors[["mir29b"]],
  "miR-29c" = s1_colors[["mir29c"]]
)
s1_method4_background <- crossing(
  s1_method4_long,
  mirna = factor(s1_mir29_labels, levels = s1_mir29_labels)
) %>%
  mutate(tissue = factor(tissue, levels = rev(colnames(s1_method4))))
s1_method4_targets <- s1_method4_targets %>%
  mutate(
    mirna = factor(mirna, levels = s1_mir29_labels),
    tissue = factor(tissue, levels = rev(colnames(s1_method4)))
  )

s1_make_method4_panel <- function(target_name) {
  ggplot(
    filter(s1_method4_background, mirna == target_name),
    aes(x = coefficient, y = tissue, group = tissue)
  ) +
    annotate(
      "rect", xmin = -0.75, xmax = 0, ymin = -Inf, ymax = Inf,
      fill = s1_colors[["negative_bg"]]
    ) +
    annotate(
      "rect", xmin = 0, xmax = 0.80, ymin = -Inf, ymax = Inf,
      fill = s1_colors[["positive_bg"]]
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.28, color = "grey55") +
    geom_boxplot(
      width = 0.53, fill = "#7A8386", color = "#4D5658",
      linewidth = 0.24, outlier.shape = 16, outlier.size = 0.22,
      outlier.alpha = 0.35
    ) +
    geom_point(
      data = filter(s1_method4_targets, mirna == target_name),
      aes(x = coefficient, y = tissue), inherit.aes = FALSE,
      shape = 21, size = 1.35, stroke = 0.28,
      fill = s1_mir29_colors[[target_name]], color = "black"
    ) +
    scale_x_continuous(
      breaks = c(-0.5, 0, 0.5), expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = c(-0.75, 0.80)) +
    labs(x = NULL, y = NULL, title = paste0(target_name, " ~ AGE")) +
    s1_theme(6.9) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(size = 5.2),
      plot.title = element_text(size = 6.2, face = "plain", hjust = 0.5),
      plot.margin = margin(10, 1, 2, 1)
    )
}

s1_panel_i_mir29a <- s1_make_method4_panel("miR-29a")
s1_panel_i_mir29b <- s1_make_method4_panel("miR-29b")
s1_panel_i_mir29c <- s1_make_method4_panel("miR-29c")
s1_panel_i_body <- plot_grid(
  s1_panel_i_host,
  s1_panel_i_mir29a, s1_panel_i_mir29b, s1_panel_i_mir29c,
  nrow = 1, rel_widths = c(0.78, 1.13, 1.13, 1.13), align = "h", axis = "tb"
)
s1_panel_i <- plot_grid(
  s1_panel_i_body,
  ggdraw() + draw_label(
    "mature miR-29 ~ age across healthy tissues (Spearman rho)",
    x = 0.63, y = 0.64, size = 5.5
  ),
  ncol = 1, rel_heights = c(1, 0.075)
)

s1_checks <- data.frame(
  Quantity = c(
    "GTEx tissues", "RNA-seq sample rows", "Selected paired tissues",
    "RNA / DNAm / strict-pair records", "Thirty-tissue estimators",
    "Paired partial-Spearman outcomes", "A/R/Z exclusive counts",
    "Genes shared by at least two signatures", "Triple-intersection genes",
    "LINC-PINT-miR-29a tissues", "LINC-PINT-miR-29a BH discoveries",
    "Method-4 matrix"
  ),
  Observed = c(
    nrow(s1_atlas), sum(s1_atlas$n_rna_sample_rows),
    sum(s1_atlas$selected_for_7tissue_dnam_analysis),
    paste(s1_funnel$sample_records, collapse = " / "),
    paste(nrow(s1_age30), "x", nrow(s1_age30_meta)),
    nrow(s1_paired_plot), paste(s1_observed_patterns, collapse = ", "),
    length(unique(s1_arz_scores$gene_id)),
    sum(s1_upset_membership$pattern == "A∩R∩Z"),
    nrow(s1_mir29a_linc), sum(s1_mir29a_linc$mir29a_tissue_BH_FDR < 0.05),
    paste(dim(s1_method4), collapse = " x ")
  ),
  Status = "PASS",
  stringsAsFactors = FALSE
)

message(
  "Panels S1a, b, f, g, h, i rebuilt from unchanged frozen inputs."
)
