## Supplementary Figure S2a-c - portable plotting module
##
## The module consumes only frozen, analysis-ready tables stored inside the
## Figure2_main bundle. It does not recalculate correlations or infer locus
## complements from expression data. Source this file with the working
## directory set to Figure2_main, or set `figure2_s2_root` first.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gghalves)
  library(readr)
  library(tibble)
})

options(stringsAsFactors = FALSE)

if (!exists("figure2_s2_root", inherits = FALSE)) {
  figure2_s2_root <- getOption("figure2_s2_root", ".")
}
if (!exists("s2_export", inherits = FALSE)) {
  s2_export <- isTRUE(getOption("s2_export", FALSE))
}

figure2_s2_root <- normalizePath(
  path.expand(as.character(figure2_s2_root[[1]])),
  mustWork = FALSE
)
s2_export <- isTRUE(s2_export)

s2_abc_data_dir <- file.path(
  figure2_s2_root, "data", "figure2_s2", "abc"
)
s2_abc_figure_dir <- file.path(
  figure2_s2_root, "figures", "supplementary_s2"
)

if (!dir.exists(s2_abc_data_dir)) {
  stop(
    "Cannot find the portable S2a-c data directory: ", s2_abc_data_dir,
    ". Set figure2_s2_root to the Figure2_main bundle root."
  )
}

s2_abc_dimensions <- list(
  S2a = c(width = 9.00, height = 2.90),
  S2b = c(width = 3.80, height = 2.90),
  S2c = c(width = 14.80, height = 4.35)
)

read_s2_tsv <- function(filename) {
  path <- file.path(s2_abc_data_dir, filename)
  if (!file.exists(path)) stop("Missing portable S2 input: ", path)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

format_signed_tick <- function(x) {
  ifelse(
    x > 0,
    paste0("+", formatC(x, format = "f", digits = 1)),
    ifelse(x < 0, formatC(x, format = "f", digits = 1), "0")
  )
}

## -------------------------------------------------------------------------
## S2a - exact miR-29 locus complement across the Figure 2 phylogeny
## -------------------------------------------------------------------------

s2a_source <- read_s2_tsv("S2a_mir29_complement_source.tsv") %>%
  arrange(species_order, match(row_id, c(
    "I_1", "I_2", "II_1", "II_2", "III_1", "III_2", "S_1"
  )))

s2a_identity_levels <- c(
  "miR-29a", "miR-29b", "miR-29c",
  "miR-29β1", "miR-29β2", "miR-29OR"
)
s2a_identity_colours <- c(
  "miR-29a" = "#F9C515",
  "miR-29b" = "#8FD3DC",
  "miR-29c" = "#F46872",
  "miR-29β1" = "#4169C1",
  "miR-29β2" = "#2E8B73",
  "miR-29OR" = "#9A9A9A"
)
s2a_identity_labels <- expression(
  "miR-29a", "miR-29b", "miR-29c",
  "miR-29" * beta * "1", "miR-29" * beta * "2", "miR-29 variant"
)

s2a_row_positions <- tibble::tribble(
  ~row_id, ~y,
  "I_1",    7.75,
  "I_2",    7.05,
  "II_1",   5.55,
  "II_2",   4.85,
  "III_1",  3.35,
  "III_2",  2.65,
  "S_1",    1.05
)

s2a_species <- s2a_source %>%
  distinct(species_order, species_full, species_short, Genome) %>%
  arrange(species_order)

s2a_plot_data <- s2a_source %>%
  left_join(s2a_row_positions, by = "row_id") %>%
  mutate(
    identity = factor(identity, levels = s2a_identity_levels),
    duplicated = as.logical(duplicated),
    present = as.logical(present)
  )

s2a_block_labels <- tibble::tribble(
  ~block_label,   ~roman, ~y,   ~fontface, ~circle_fill, ~circle_text,
  "LINC-PINT",   "I",    7.40, "italic",  "#4B4B4B",   "white",
  "MIR29B2CHG",  "II",   5.20, "italic",  "#777777",   "white",
  "duplication", "III",  3.00, "plain",   "#B5B5B5",   "#303030",
  "singleton",   "",     1.05, "plain",   NA_character_, NA_character_
)

s2a_validation <- list(
  species_59 = nrow(s2a_species) == 59L,
  seven_rows_per_species = nrow(s2a_plot_data) == 59L * 7L &&
    all(table(s2a_plot_data$Genome) == 7L),
  exact_endpoints = identical(s2a_species$Genome[[1]], "Callorhinchus_milii") &&
    identical(tail(s2a_species$Genome, 1), "Homo_sapiens"),
  no_source_tree_outgroup = !any(s2a_species$Genome == "Saccoglossus_kowalevskii"),
  singleton_or_six = sum(
    s2a_plot_data$row_id == "S_1" &
      as.character(s2a_plot_data$identity) == "miR-29OR",
    na.rm = TRUE
  ) == 6L,
  singleton_beta2_three = sum(
    s2a_plot_data$row_id == "S_1" &
      as.character(s2a_plot_data$identity) == "miR-29β2",
    na.rm = TRUE
  ) == 3L,
  cluster_iii_six_members = sum(
    s2a_plot_data$block == "Cluster III" & s2a_plot_data$present,
    na.rm = TRUE
  ) == 6L,
  fifteen_duplicate_markers = sum(s2a_plot_data$duplicated, na.rm = TRUE) == 15L
)
stopifnot(all(unlist(s2a_validation, use.names = FALSE)))

p_S2a <- ggplot(s2a_plot_data, aes(x = species_order, y = y)) +
  geom_hline(
    yintercept = s2a_row_positions$y,
    linewidth = 0.20, colour = "#E9E9E9"
  ) +
  geom_vline(
    xintercept = s2a_species$species_order,
    linewidth = 0.17, colour = "#F0F0F0"
  ) +
  geom_point(
    shape = 21, size = 1.38, stroke = 0.19,
    fill = "white", colour = "#D2D2D2"
  ) +
  geom_point(
    data = filter(s2a_plot_data, present),
    aes(fill = identity),
    shape = 21, size = 1.70, stroke = 0.22, colour = "#555555"
  ) +
  geom_point(
    data = filter(s2a_plot_data, duplicated),
    aes(shape = "intra-cluster duplication"),
    size = 0.52, colour = "black"
  ) +
  geom_text(
    data = filter(s2a_block_labels, fontface == "italic"),
    aes(x = -2.70, y = y, label = block_label),
    inherit.aes = FALSE, hjust = 1,
    family = "Helvetica", fontface = "italic", size = 2.55,
    colour = "#303030"
  ) +
  geom_text(
    data = filter(s2a_block_labels, fontface == "plain"),
    aes(x = -2.70, y = y, label = block_label),
    inherit.aes = FALSE, hjust = 1,
    family = "Helvetica", size = 2.55, colour = "#303030"
  ) +
  geom_point(
    data = filter(s2a_block_labels, roman != ""),
    aes(x = -2.10, y = y, fill = circle_fill),
    inherit.aes = FALSE, shape = 21, size = 3.10,
    stroke = 0, show.legend = FALSE
  ) +
  geom_text(
    data = filter(s2a_block_labels, roman != ""),
    aes(x = -2.10, y = y, label = roman, colour = circle_text),
    inherit.aes = FALSE, family = "Helvetica", fontface = "bold",
    size = 1.92, show.legend = FALSE
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(s2a_identity_colours, stats::setNames(
      s2a_block_labels$circle_fill[!is.na(s2a_block_labels$circle_fill)],
      s2a_block_labels$circle_fill[!is.na(s2a_block_labels$circle_fill)]
    )),
    breaks = s2a_identity_levels,
    labels = s2a_identity_labels,
    drop = FALSE,
    guide = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, size = 2.30, colour = "#555555")
    )
  ) +
  scale_colour_identity() +
  scale_shape_manual(
    name = NULL,
    values = c("intra-cluster duplication" = 16),
    guide = guide_legend(
      order = 2,
      override.aes = list(shape = 16, size = 1.45, colour = "black")
    )
  ) +
  scale_x_continuous(
    breaks = s2a_species$species_order,
    labels = s2a_species$species_short,
    limits = c(-7.3, 59.6),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(limits = c(0.55, 8.10), expand = expansion(mult = c(0, 0))) +
  labs(title = "miR-29 locus complement", x = NULL, y = NULL) +
  theme_minimal(base_family = "Helvetica", base_size = 7.0) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, face = "italic",
      colour = "#202020", size = 4.05, margin = margin(t = 1.4)
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 5.65, colour = "#333333"),
    legend.spacing.x = grid::unit(1.5, "pt"),
    legend.key.width = grid::unit(7.2, "pt"),
    legend.key.height = grid::unit(5.4, "pt"),
    plot.title = element_text(
      face = "bold", size = 7.4, colour = "#171717",
      margin = margin(b = 0.4)
    ),
    plot.margin = margin(2.5, 3.5, 2.5, 3.5)
  )

## -------------------------------------------------------------------------
## S2b - compact miRNA-age sample/stratum coverage table
## -------------------------------------------------------------------------

s2b_coverage <- read_s2_tsv("S2b_mirna_coverage_source.tsv") %>%
  arrange(species_order) %>%
  mutate(
    n = as.integer(n),
    n_tissues = as.integer(n_tissues),
    y = if_else(species_key == "human", 0.62, 13 - as.numeric(species_order)),
    bar_start = 0.315,
    bar_end = bar_start + 0.400 *
      (log10(n) - min(log10(n))) /
      (max(log10(n)) - min(log10(n))),
    sample_tissue_label = paste0(
      format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
      " / ", n_tissues
    )
  )

s2b_provenance <- read_s2_tsv("S2b_dataset_provenance.tsv") %>%
  arrange(match(species_key, s2b_coverage$species_key))

s2b_validation <- list(
  twelve_species = nrow(s2b_coverage) == 12L,
  total_samples_17635 = sum(s2b_coverage$n) == 17635L,
  total_strata_68 = sum(s2b_coverage$n_tissues) == 68L,
  nonhuman_samples_1403 = sum(s2b_coverage$n[s2b_coverage$species_key != "human"]) == 1403L,
  human_gtex_16232_30 = with(
    filter(s2b_coverage, species_key == "human"),
    length(n) == 1L && n == 16232L && n_tissues == 30L
  ),
  macaca_brain_only_14_1 = with(
    filter(s2b_coverage, species_key == "macaca"),
    length(n) == 1L && n == 14L && n_tissues == 1L && tissue_names == "Brain"
  ),
  complete_provenance = setequal(s2b_coverage$species_key, s2b_provenance$species_key)
)
stopifnot(all(unlist(s2b_validation, use.names = FALSE)))

p_S2b <- ggplot(s2b_coverage, aes(y = y)) +
  geom_hline(yintercept = 1.40, linetype = "dotted", linewidth = 0.38, colour = "#777777") +
  geom_hline(yintercept = 0.08, linewidth = 0.32, colour = "#8A8A8A") +
  geom_vline(xintercept = c(0.285, 0.755), linewidth = 0.26, colour = "#B8B8B8") +
  geom_segment(
    aes(x = bar_start, xend = bar_end, yend = y),
    linewidth = 3.15, lineend = "round", colour = "#D1D1D1"
  ) +
  geom_point(
    aes(x = bar_end), shape = 21, size = 2.05,
    stroke = 0.30, fill = "#656565", colour = "#353535"
  ) +
  geom_text(
    aes(x = 0.015, label = species_label),
    hjust = 0, fontface = "italic", family = "Helvetica", size = 2.90
  ) +
  geom_text(
    aes(x = 0.985, label = sample_tissue_label),
    hjust = 1, fontface = "italic", family = "Helvetica", size = 2.78
  ) +
  annotate(
    "text", x = 0.015, y = 13.05, label = "species",
    hjust = 0, fontface = "bold", family = "Helvetica", size = 2.60
  ) +
  annotate(
    "text", x = 0.315, y = 13.05, label = "samples (log)",
    hjust = 0, fontface = "bold", family = "Helvetica", size = 2.60
  ) +
  annotate(
    "text", x = 0.985, y = 13.05, label = "samples/tissues",
    hjust = 1, fontface = "bold", family = "Helvetica", size = 2.60
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 13.35), clip = "off") +
  theme_void(base_family = "Helvetica", base_size = 7.0) +
  theme(plot.margin = margin(2.5, 3.5, 2.5, 3.5))

## -------------------------------------------------------------------------
## S2c - tissue-resolved mature-miRNA age-association distributions
## -------------------------------------------------------------------------

s2c_distribution <- read_s2_tsv("S2c_mirna_tissue_distribution_source.tsv") %>%
  mutate(
    x = as.numeric(x),
    rho = as.numeric(rho),
    x_points = x + 0.120
  )

s2c_highlights <- read_s2_tsv("S2c_mirna_tissue_highlights.tsv") %>%
  mutate(
    x = as.numeric(x),
    x_highlight = as.numeric(x_highlight),
    rho = as.numeric(rho),
    highlight_class = factor(
      highlight_class,
      levels = c("miR-29b", "miR-29c", "miR-29beta1", "miR-29a")
    )
  )

s2c_tissue_index <- read_s2_tsv("S2c_mirna_tissue_index.tsv") %>%
  arrange(species_order, tissue_order) %>%
  mutate(x = as.numeric(x), n_samples = as.integer(n_samples))

s2c_species_spans <- s2c_tissue_index %>%
  group_by(species_key, species_label, species_order) %>%
  summarise(
    xmin = min(x), xmax = max(x), xmid = mean(c(xmin, xmax)),
    .groups = "drop"
  ) %>%
  arrange(species_order)

s2c_separators <- head(s2c_species_spans$xmax + 0.5, -1L)
s2c_x_labels <- paste0(
  s2c_tissue_index$tissue_label,
  "<br><span style='color:#A8A8A8'>n=",
  s2c_tissue_index$n_samples,
  "</span>"
)
names(s2c_x_labels) <- as.character(s2c_tissue_index$x)

s2c_highlight_colours <- c(
  "miR-29b" = "#A3D9DE",
  "miR-29c" = "#F36275",
  "miR-29beta1" = "#3155C6",
  "miR-29a" = "#F8C717"
)

s2c_validation <- list(
  five_multitissue_species = identical(
    s2c_species_spans$species_key,
    c("zfish", "furzeri", "goat", "mouse", "rat")
  ),
  exact_strata = nrow(s2c_tissue_index) == 32L,
  exact_strata_per_species = identical(
    as.integer(table(factor(
      s2c_tissue_index$species_key,
      levels = c("zfish", "furzeri", "goat", "mouse", "rat")
    ))),
    c(2L, 4L, 2L, 13L, 11L)
  ),
  exact_distribution_rows = nrow(s2c_distribution) == 23772L,
  exact_highlight_rows = nrow(s2c_highlights) == 102L,
  no_human_profiles = !any(s2c_distribution$species_key == "human"),
  x_index_complete = setequal(unique(s2c_distribution$x), s2c_tissue_index$x),
  finite_rho = all(is.finite(s2c_distribution$rho)),
  source_highlight_colours_unchanged = all(
    unique(s2c_highlights$highlight_color) %in% unname(s2c_highlight_colours)
  )
)
stopifnot(all(unlist(s2c_validation, use.names = FALSE)))

set.seed(290629)
p_S2c <- ggplot(s2c_distribution, aes(x = x, y = rho, group = x)) +
  geom_segment(
    data = s2c_tissue_index,
    aes(x = x, xend = x, y = -1, yend = 1),
    inherit.aes = FALSE, linetype = "dotted",
    linewidth = 0.20, colour = "#3E3E3E"
  ) +
  geom_vline(
    xintercept = s2c_separators, linetype = "dashed",
    linewidth = 0.38, colour = "#888888"
  ) +
  gghalves::geom_half_boxplot(
    side = "l", width = 0.40, outlier.shape = NA,
    fill = "#5A5A5A", colour = "#202020", linewidth = 0.34,
    show.legend = FALSE
  ) +
  gghalves::geom_half_violin(
    side = "r", scale = "width", trim = TRUE, width = 0.58,
    fill = NA, colour = "#444444", linetype = "dotted",
    linewidth = 0.25, adjust = 1.05, show.legend = FALSE
  ) +
  geom_point(
    data = filter(s2c_distribution, sign == "negative"),
    aes(x = x_points, y = rho), inherit.aes = FALSE,
    colour = "#33417F", shape = 16, size = 0.20, alpha = 0.76,
    position = position_jitter(width = 0.019, height = 0, seed = 290629),
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(s2c_distribution, sign == "positive"),
    aes(x = x_points, y = rho), inherit.aes = FALSE,
    colour = "#EE3137", shape = 16, size = 0.20, alpha = 0.76,
    position = position_jitter(width = 0.019, height = 0, seed = 290630),
    show.legend = FALSE
  ) +
  geom_point(
    data = s2c_highlights,
    aes(x = x_highlight, y = rho, fill = highlight_class),
    inherit.aes = FALSE, shape = 21, colour = "#555555",
    stroke = 0.30, size = 1.62
  ) +
  geom_segment(
    data = s2c_species_spans,
    aes(x = xmin + 0.08, xend = xmax - 0.08, y = 1.105, yend = 1.105),
    inherit.aes = FALSE, linewidth = 0.32, colour = "#737373"
  ) +
  geom_text(
    data = s2c_species_spans,
    aes(x = xmid, y = 1.205, label = species_label),
    inherit.aes = FALSE, family = "Helvetica", fontface = "italic",
    size = 2.82
  ) +
  scale_fill_manual(
    values = s2c_highlight_colours,
    breaks = names(s2c_highlight_colours),
    labels = expression(
      "miR-29b", "miR-29c", "miR-29" * beta * "1", "miR-29a"
    ),
    name = NULL,
    guide = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, size = 2.55, colour = "#555555")
    )
  ) +
  scale_x_continuous(
    breaks = s2c_tissue_index$x,
    labels = s2c_x_labels[as.character(s2c_tissue_index$x)],
    limits = c(0.45, max(s2c_tissue_index$x) + 0.65),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = format_signed_tick,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(-1, 1.29), clip = "off") +
  labs(
    x = NULL,
    y = expression("miRNA abundance ~ Age (" * rho * ")")
  ) +
  theme_bw(base_family = "Helvetica", base_size = 7.0) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#454545", fill = NA, linewidth = 0.45),
    axis.text.x = ggtext::element_markdown(
      angle = 62, hjust = 1, vjust = 1,
      size = 4.75, colour = "#202020", lineheight = 0.90
    ),
    axis.text.y = element_text(size = 6.15, colour = "#202020"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.30),
    axis.title.y = element_text(size = 7.35, margin = margin(r = 4)),
    legend.position = "bottom",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(7.0, "pt"),
    legend.key.height = grid::unit(6.0, "pt"),
    legend.spacing.x = grid::unit(1.4, "pt"),
    legend.text = element_text(size = 5.5),
    plot.margin = margin(11, 5, 2, 5)
  )

## -------------------------------------------------------------------------
## Public objects, checks, provenance, and optional exports
## -------------------------------------------------------------------------

s2_abc_data <- list(
  S2a = list(
    locus_matrix = s2a_plot_data,
    species_order = s2a_species,
    identity_colours = s2a_identity_colours
  ),
  S2b = list(
    coverage = s2b_coverage,
    dataset_provenance = s2b_provenance
  ),
  S2c = list(
    distributions = s2c_distribution,
    highlights = s2c_highlights,
    tissue_index = s2c_tissue_index,
    species_spans = s2c_species_spans
  )
)

s2_abc_checks <- dplyr::bind_rows(
  tibble::tibble(
    panel = "S2a",
    check = names(s2a_validation),
    status = as.logical(unlist(s2a_validation, use.names = FALSE)),
    detail = c(
      "59 displayed vertebrates in the frozen anticlockwise order.",
      "413 matrix cells: seven locus rows for every species.",
      "Order runs from Callorhinchus milii to Homo sapiens.",
      "The Saccoglossus source-tree outgroup is not displayed.",
      "Six ancestral singleton miR-29-variant calls.",
      "Three singleton miR-29beta2 calls.",
      "Three teleosts carry both members of Cluster III.",
      "Fifteen Cluster II cells carry an intra-cluster duplication marker."
    )
  ),
  tibble::tibble(
    panel = "S2b",
    check = names(s2b_validation),
    status = as.logical(unlist(s2b_validation, use.names = FALSE)),
    detail = c(
      "Twelve rows including the contextual human GTEx reference.",
      "17,635 displayed profiles including GTEx.",
      "68 displayed analysis tissues/strata including GTEx.",
      "1,403 profiles across the 11 non-human Figure 2F datasets.",
      paste0(
        "Human contextual GTEx RNA-atlas reference: 16,232 profiles across ",
        "30 tissues; this is not a count of directly measured, paired ",
        "mature-miRNA libraries."
      ),
      "Accepted macaque analysis: brain only, n=14.",
      "Every coverage row has a matching repository/provenance record."
    )
  ),
  tibble::tibble(
    panel = "S2c",
    check = names(s2c_validation),
    status = as.logical(unlist(s2c_validation, use.names = FALSE)),
    detail = c(
      "D. rerio, N. furzeri, C. hircus, M. musculus and R. norvegicus.",
      "Thirty-two exact species-by-tissue strata.",
      "Strata per species: 2, 4, 2, 13 and 11.",
      "23,772 finite mature-miRNA tissue-level Spearman coefficients.",
      "102 frozen miR-29 highlight points.",
      "Human profiles are not part of this non-human tissue strip.",
      "Every plotted coefficient maps to the frozen tissue index.",
      "No non-finite Spearman coefficient is plotted.",
      "The source highlight palette is reproduced without alteration."
    )
  )
)

s2_abc_notes <- list(
  dimensions_inches = s2_abc_dimensions,
  provenance = list(
    S2a = paste0(
      "Frozen 59-vertebrate anticlockwise locus matrix; complement derives ",
      "from the curated phylogeny annotation, not from expression datasets."
    ),
    S2b = paste0(
      "Non-human entries summarize the Figure 2F species datasets. The ",
      "H. sapiens entry (16,232 profiles; 30 tissues) is the broader GTEx ",
      "RNA-atlas reference from Figure 1/S1 and is displayed for context. ",
      "It is not a count of directly measured, paired mature-miRNA libraries ",
      "and is not included in the non-human Figure 2F/G aggregation."
    ),
    S2c = paste0(
      "Precomputed tissue-specific Spearman coefficients for the five ",
      "multi-tissue non-human species; distributions show all finite mature-",
      "miRNA coefficients and coloured circles mark the frozen miR-29 calls."
    )
  ),
  illustrator_manual_elements = c(
    S2a = "Panel letter is intentionally omitted.",
    S2b = "Panel letter is intentionally omitted.",
    S2c = paste0(
      "Animal silhouettes are intentionally not recreated in R. Add the ",
      "five established Illustrator silhouettes above the species labels; ",
      "the plot leaves top clearance and contains the corresponding anchors."
    )
  )
)

save_s2_abc_panel <- function(plot, stem, dimensions, dpi = 450) {
  dir.create(s2_abc_figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(s2_abc_figure_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(s2_abc_figure_dir, paste0(stem, ".svg"))
  png_path <- file.path(s2_abc_figure_dir, paste0(stem, ".png"))

  ggplot2::ggsave(
    pdf_path, plot = plot,
    width = dimensions[["width"]], height = dimensions[["height"]],
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Package 'svglite' is required when s2_export = TRUE.")
  }
  ggplot2::ggsave(
    svg_path, plot = plot,
    width = dimensions[["width"]], height = dimensions[["height"]],
    units = "in", device = svglite::svglite, bg = "white"
  )
  ggplot2::ggsave(
    png_path, plot = plot,
    width = dimensions[["width"]], height = dimensions[["height"]],
    units = "in", dpi = dpi, bg = "white"
  )

  invisible(c(pdf = pdf_path, svg = svg_path, png = png_path))
}

s2_abc_export_paths <- NULL
if (s2_export) {
  s2_abc_export_paths <- list(
    S2a = save_s2_abc_panel(p_S2a, "S2a", s2_abc_dimensions$S2a),
    S2b = save_s2_abc_panel(p_S2b, "S2b", s2_abc_dimensions$S2b),
    S2c = save_s2_abc_panel(p_S2c, "S2c", s2_abc_dimensions$S2c)
  )
}

invisible(list(
  p_S2a = p_S2a,
  p_S2b = p_S2b,
  p_S2c = p_S2c,
  s2_abc_data = s2_abc_data,
  s2_abc_checks = s2_abc_checks,
  exports = s2_abc_export_paths
))
