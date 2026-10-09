# Figure 2 main — panels h–k (pan-mammalian DNA methylation)
#
# RStudio use:
#   1. Open the portable Figure2_main directory as the working directory.
#   2. Source this file, or copy/paste it into the console.
#
# The script reads only data/methylation/figure2_methylation_portable_input.rds.
# It never reads Marco's Desktop/Dropbox and never downloads data.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gghalves)
})

required_input <- file.path("data", "methylation", "figure2_methylation_portable_input.rds")

find_figure2_root <- function() {
  starts <- unique(c(
    normalizePath(getwd(), mustWork = FALSE),
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      editor_path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
      if (nzchar(editor_path)) dirname(normalizePath(editor_path, mustWork = FALSE)) else character()
    } else {
      character()
    }
  ))
  for (start in starts) {
    candidate <- start
    for (i in 0:6) {
      if (file.exists(file.path(candidate, required_input))) return(candidate)
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop(
    "Cannot locate Figure2_main. Set the RStudio working directory to the ",
    "portable Figure2_main folder and rerun this script."
  )
}

figure2_root <- find_figure2_root()
input_path <- file.path(figure2_root, required_input)
figure_dir <- file.path(figure2_root, "figures", "generated")
data_dir <- file.path(figure2_root, "data", "methylation")
validation_dir <- file.path(figure2_root, "validation", "methylation")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

methylation_input <- readRDS(input_path)
target_probe <- methylation_input$metadata$target_probe
tissue_n_min <- methylation_input$metadata$tissue_n_min
clade_species <- methylation_input$clade_species
human <- "Homo sapiens"
probe_ids <- rownames(methylation_input$rho_by_species_tissue[[human]])
fisher_epsilon <- 1e-12

stopifnot(
  target_probe == "cg10501210",
  tissue_n_min == 10L,
  length(probe_ids) == 37554L,
  target_probe %in% probe_ids,
  length(unlist(clade_species, use.names = FALSE)) == 36L
)

arithmetic_row_mean <- function(x) {
  x <- as.matrix(x)
  n_valid <- rowSums(is.finite(x))
  out <- rowMeans(x, na.rm = TRUE)
  out[n_valid == 0L] <- NA_real_
  names(out) <- rownames(x)
  out
}

fisher_row_mean <- function(x) {
  x <- as.matrix(x)
  valid <- is.finite(x)
  n_valid <- rowSums(valid)
  x[valid & x >= 1] <- 1 - fisher_epsilon
  x[valid & x <= -1] <- -1 + fisher_epsilon
  out <- tanh(rowMeans(atanh(x), na.rm = TRUE))
  out[n_valid == 0L] <- NA_real_
  names(out) <- rownames(x)
  out
}

# ---------------------------------------------------------------------------
# Reconstruct the locked hierarchy
# ---------------------------------------------------------------------------
# Tissue-level entries are cached Spearman rho values.  Within each species,
# eligible tissues are assigned equal weight and averaged directly in rho
# space.  Every aggregation above species uses equal-unit Fisher z.

all_species <- unname(unlist(clade_species, use.names = FALSE))
species_rho <- matrix(
  NA_real_,
  nrow = length(probe_ids),
  ncol = length(all_species),
  dimnames = list(probe_ids, all_species)
)
for (sp in all_species) {
  species_rho[, sp] <- arithmetic_row_mean(
    methylation_input$rho_by_species_tissue[[sp]][probe_ids, , drop = FALSE]
  )
}

myotis <- methylation_input$special_collapses$Myotis
tursiops <- methylation_input$special_collapses$Tursiops
other_chiroptera <- setdiff(clade_species$Chiroptera, myotis)
other_cetacea <- setdiff(clade_species$Cetacea, tursiops)
primates_nonhuman <- setdiff(clade_species$Primates, human)

myotis_rho <- fisher_row_mean(species_rho[, myotis, drop = FALSE])
tursiops_rho <- fisher_row_mean(species_rho[, tursiops, drop = FALSE])

chiroptera_units <- cbind(
  species_rho[, other_chiroptera, drop = FALSE],
  Myotis = myotis_rho
)
cetacea_units <- cbind(
  species_rho[, other_cetacea, drop = FALSE],
  Tursiops = tursiops_rho
)

panel_j_rho <- cbind(
  Primates = fisher_row_mean(species_rho[, primates_nonhuman, drop = FALSE]),
  `mu-Chiroptera` = fisher_row_mean(chiroptera_units),
  Carnivora = fisher_row_mean(species_rho[, clade_species$Carnivora, drop = FALSE]),
  Artiodactyla = fisher_row_mean(species_rho[, clade_species$Artiodactyla, drop = FALSE]),
  Cetacea = fisher_row_mean(cetacea_units)
)

# Panel k differs from j only in the primate input: Human is included here.
panel_k_clades <- cbind(
  Primates = fisher_row_mean(species_rho[, clade_species$Primates, drop = FALSE]),
  Chiroptera = panel_j_rho[, "mu-Chiroptera"],
  Carnivora = panel_j_rho[, "Carnivora"],
  Artiodactyla = panel_j_rho[, "Artiodactyla"],
  Cetacea = panel_j_rho[, "Cetacea"]
)
panel_k_rho <- fisher_row_mean(panel_k_clades)

sample_annotation <- methylation_input$sample_annotation
raw_n <- function(species) sum(sample_annotation$SpeciesLatinName %in% species)

panel_j_counts <- data.frame(
  group = colnames(panel_j_rho),
  display_group = c("Primates", "µ-Chiroptera", "Carnivora", "Artiodactyla", "Cetacea"),
  n = c(
    raw_n(primates_nonhuman),
    raw_n(clade_species$Chiroptera),
    raw_n(clade_species$Carnivora),
    raw_n(clade_species$Artiodactyla),
    raw_n(clade_species$Cetacea)
  ),
  species_n = c(4L, 10L, 6L, 9L, 6L),
  stringsAsFactors = FALSE
)

target_rank <- sum(panel_k_rho[is.finite(panel_k_rho)] < panel_k_rho[[target_probe]]) + 1L
target_k <- unname(panel_k_rho[[target_probe]])

stopifnot(
  raw_n(human) == 1126L,
  identical(panel_j_counts$n, c(944L, 347L, 1051L, 1344L, 962L)),
  sum(panel_j_counts$n) + raw_n(human) == 5774L,
  sum(panel_j_counts$species_n) + 1L == 36L,
  abs(species_rho[target_probe, human] - (-0.509146370103404)) < 1e-12,
  abs(target_k - (-0.617636318578119)) < 1e-10,
  target_rank == 1L
)

# ---------------------------------------------------------------------------
# Reusable plotting helpers
# ---------------------------------------------------------------------------

negative_colour <- "#64C99B"
positive_colour <- "#6537A6"
target_colour <- "#FFEB00"

distribution_long <- function(x, display_names = colnames(x)) {
  stopifnot(length(display_names) == ncol(x))
  out <- data.frame(
    probe = rep(rownames(x), times = ncol(x)),
    group = rep(display_names, each = nrow(x)),
    rho = as.vector(x),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$rho), , drop = FALSE]
  out$group <- factor(out$group, levels = display_names)
  out$x <- as.numeric(out$group)
  set.seed(2902)
  out$x_point <- out$x + 0.19 + stats::runif(nrow(out), -0.022, 0.022)
  out$direction <- ifelse(out$rho >= 0, "positive", "negative")
  out$is_target <- out$probe == target_probe
  out
}

half_distribution_plot <- function(x, display_names, y_limits, x_labels = display_names) {
  plot_df <- distribution_long(x, display_names)
  target_df <- plot_df[plot_df$is_target, , drop = FALSE]

  p <- ggplot(plot_df, aes(x = x, y = rho, group = group)) +
    geom_vline(xintercept = seq_along(display_names), colour = "grey92", linewidth = 0.45) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    gghalves::geom_half_violin(
      side = "r", fill = NA, colour = "black", linewidth = 0.45,
      linetype = "dotted", trim = FALSE, scale = "width", width = 0.82
    ) +
    gghalves::geom_half_boxplot(
      side = "l", center = FALSE, errorbar.draw = TRUE,
      width = 0.34, outlier.shape = NA, fill = "grey25",
      colour = "black", linewidth = 0.45
    )

  point_layer <- if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(
      data = plot_df,
      mapping = aes(x = x_point, y = rho, fill = direction),
      inherit.aes = FALSE, shape = 21, colour = "white", stroke = 0.05,
      size = 0.38, alpha = 0.78,
      raster.dpi = 600
    )
  } else {
    geom_point(
      data = plot_df,
      mapping = aes(x = x_point, y = rho, fill = direction),
      inherit.aes = FALSE, shape = 21, colour = "white", stroke = 0.05,
      size = 0.38, alpha = 0.45
    )
  }

  p +
    point_layer +
    geom_point(
      data = target_df,
      aes(x = x + 0.27, y = rho), inherit.aes = FALSE,
      shape = 21, fill = target_colour, colour = "black",
      size = 3.1, stroke = 0.8
    ) +
    scale_fill_manual(values = c(negative = negative_colour, positive = positive_colour)) +
    scale_x_continuous(
      breaks = seq_along(display_names), labels = x_labels,
      expand = expansion(mult = c(0.08, 0.11))
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    labs(x = NULL, y = expression("DNA methylation ~ Age (" * rho * ")")) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.line.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.x = element_text(colour = "black", size = 10),
      axis.text.y = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black", margin = margin(r = 7)),
      plot.margin = margin(7, 10, 7, 7)
    )
}

# ---------------------------------------------------------------------------
# Panel h — Human methylation/age distribution
# ---------------------------------------------------------------------------

panel_h_matrix <- matrix(
  species_rho[, human], ncol = 1,
  dimnames = list(probe_ids, "H. sapiens")
)
p_F2H <- half_distribution_plot(
  panel_h_matrix,
  display_names = "H. sapiens",
  y_limits = c(-0.7, 0.75),
  x_labels = "H. sapiens\nn = 1126"
)

# ---------------------------------------------------------------------------
# Panel i — sequence inventory supporting the manual alignment artwork
# ---------------------------------------------------------------------------
# The Illustrator panel also contains silhouettes, the DNA cartoon, an
# extended orthologue alignment and manual target-position annotation.  Those
# elements are preserved in reference/panel_i_reference.png.  The plot below
# is the reproducible analytical component available from the supplied source
# workbook: the 50-mer sequence inventory, grouped by mapped CpG status.

alignment <- methylation_input$alignment
alignment <- alignment[!is.na(alignment$sequence) & nzchar(alignment$sequence), , drop = FALSE]
alignment$sequence_clean <- gsub("[^ACGT-]", "-", alignment$sequence)
max_sequence_length <- max(nchar(alignment$sequence_clean))
alignment$status <- ifelse(
  is.na(alignment$conserved), "Status unavailable",
  ifelse(alignment$conserved, "Intact CpG", "Disrupted/non-conserved")
)
alignment$status <- factor(
  alignment$status,
  levels = c("Disrupted/non-conserved", "Intact CpG", "Status unavailable")
)
alignment <- alignment[order(alignment$status, alignment$species), , drop = FALSE]
alignment$row <- rev(seq_len(nrow(alignment)))

sequence_tiles <- do.call(rbind, lapply(seq_len(nrow(alignment)), function(i) {
  chars <- strsplit(alignment$sequence_clean[[i]], "", fixed = TRUE)[[1]]
  data.frame(
    species = alignment$species[[i]], status = alignment$status[[i]],
    row = alignment$row[[i]], position = seq_along(chars), nucleotide = chars,
    stringsAsFactors = FALSE
  )
}))
sequence_tiles$nucleotide <- factor(
  sequence_tiles$nucleotide,
  levels = c("A", "C", "G", "T", "-")
)

p_F2I <- ggplot(sequence_tiles, aes(x = position, y = row, fill = nucleotide)) +
  geom_tile(width = 1.02, height = 1.02) +
  facet_grid(status ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(A = "#F2BFC1", C = "#BDD2FA", G = "#BFE8C4", T = "#FFE2A5", `-` = "white"),
    drop = FALSE
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 10, 20, 30, 40, 50)) +
  labs(
    x = "Position in supplied cg10501210 probe-sequence inventory",
    y = NULL, fill = "Base"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), strip.text.y = element_text(angle = 0),
    legend.position = "top", plot.margin = margin(6, 7, 6, 7)
  )

# ---------------------------------------------------------------------------
# Panel j — five non-Human macrogroup distributions
# ---------------------------------------------------------------------------

j_x_labels <- paste0(
  panel_j_counts$display_group,
  "\nn = ", panel_j_counts$n,
  "\nSp = ", panel_j_counts$species_n
)
p_F2J <- half_distribution_plot(
  panel_j_rho,
  display_names = panel_j_counts$group,
  y_limits = c(-1, 1),
  x_labels = j_x_labels
) +
  annotate("segment", x = 0.75, xend = 2.25, y = 0.96, yend = 0.96, linewidth = 0.55) +
  annotate("text", x = 1.5, y = 0.99, label = "Euarchontoglires", fontface = "italic", vjust = 0) +
  annotate("segment", x = 2.75, xend = 5.25, y = 0.96, yend = 0.96, linewidth = 0.55) +
  annotate("text", x = 4, y = 0.99, label = "Laurasiatheria", fontface = "italic", vjust = 0) +
  coord_cartesian(ylim = c(-1, 1), clip = "off")

# ---------------------------------------------------------------------------
# Panel k — equal-clade Boreoeutheria rank curve
# ---------------------------------------------------------------------------

k_order <- order(panel_k_rho, na.last = NA)
panel_k_df <- data.frame(
  probe = names(panel_k_rho)[k_order],
  rho = unname(panel_k_rho[k_order]),
  rank = seq_along(k_order),
  stringsAsFactors = FALSE
)
panel_k_df$rank_fraction <- (panel_k_df$rank - 1) / (nrow(panel_k_df) - 1)
panel_k_df$direction <- ifelse(panel_k_df$rho >= 0, "positive", "negative")
target_k_df <- panel_k_df[panel_k_df$probe == target_probe, , drop = FALSE]

p_F2K <- ggplot(panel_k_df, aes(x = rank_fraction, y = rho, colour = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_path(linewidth = 1.15, lineend = "round") +
  geom_point(size = 0.25, alpha = 0.65) +
  geom_point(
    data = target_k_df, aes(x = rank_fraction, y = rho),
    inherit.aes = FALSE, shape = 21, fill = target_colour,
    colour = "black", size = 3.2, stroke = 0.8
  ) +
  annotate(
    "text", x = 0.16, y = target_k + 0.09,
    label = sprintf("MIR29B2CHG\nρ = %.3f", target_k),
    hjust = 0, fontface = "italic", size = 3.7
  ) +
  annotate(
    "segment", x = 0.15, xend = target_k_df$rank_fraction,
    y = target_k + 0.07, yend = target_k + 0.015,
    arrow = arrow(length = grid::unit(0.10, "cm")), linewidth = 0.4
  ) +
  scale_colour_manual(values = c(negative = negative_colour, positive = positive_colour)) +
  scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(limits = c(-0.72, 0.82), breaks = c(-0.4, 0, 0.4, 0.8)) +
  labs(
    x = expression(bar(rho)), y = NULL,
    title = "Boreoeutheria",
    subtitle = "n = 5774; Sp = 36; Human included in Primates"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none", axis.text.y = element_text(colour = "black"),
    plot.title = element_text(face = "italic", hjust = 0.5),
    plot.subtitle = element_text(colour = "grey40", size = 8.5, hjust = 0.5),
    plot.margin = margin(7, 9, 7, 7)
  )

# Stable names consumed by Figure2_main.Rmd when this module is sourced.
panel_h <- p_F2H
panel_i <- p_F2I
panel_j <- p_F2J
panel_k <- p_F2K

figure2_methylation_summary <- list(
  target_probe = target_probe,
  human = list(n = raw_n(human), rho = unname(species_rho[target_probe, human])),
  panel_j = panel_j_counts,
  boreoeutheria = list(
    n = raw_n(all_species), species_n = length(all_species),
    rho = target_k, rank = target_rank, rank_denominator = sum(is.finite(panel_k_rho)),
    human_included_in_primates = TRUE
  )
)
figure2_methylation_check_pass <- c(
  probe_universe_37554 = length(probe_ids) == 37554L,
  human_n_1126 = raw_n(human) == 1126L,
  total_n_5774 = raw_n(all_species) == 5774L,
  species_n_36 = length(all_species) == 36L,
  panel_k_human_included = isTRUE(methylation_input$metadata$human_in_panel_k),
  panel_k_target_rho_matches = abs(target_k - (-0.617636318578119)) < 1e-10,
  panel_k_target_rank_1 = target_rank == 1L
)
figure2_methylation_checks <- data.frame(
  Quantity = c(
    "Probe universe", "Human filtered samples", "All-panel filtered samples",
    "Species in K", "Human included in K Primates", "K target rho",
    "K target ascending rank"
  ),
  Observed = c(
    length(probe_ids), raw_n(human), raw_n(all_species), length(all_species),
    isTRUE(methylation_input$metadata$human_in_panel_k),
    sprintf("%.12f", target_k), target_rank
  ),
  Expected = c(37554L, 1126L, 5774L, 36L, TRUE, "-0.617636318578", 1L),
  Status = ifelse(figure2_methylation_check_pass, "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
stopifnot(all(figure2_methylation_checks$Status == "PASS"))

# ---------------------------------------------------------------------------
# Export data, plots and a reusable R object cache
# ---------------------------------------------------------------------------

target_summary <- rbind(
  data.frame(panel = "H", group = "Homo sapiens", rho = species_rho[target_probe, human],
             n = 1126L, species_n = 1L, human_included = TRUE),
  data.frame(panel = "J", group = panel_j_counts$display_group,
             rho = as.numeric(panel_j_rho[target_probe, ]), n = panel_j_counts$n,
             species_n = panel_j_counts$species_n, human_included = FALSE),
  data.frame(panel = "K", group = "Boreoeutheria",
             rho = target_k, n = 5774L, species_n = 36L, human_included = TRUE)
)
write.table(
  target_summary,
  file.path(data_dir, "figure2_methylation_target_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  figure2_methylation_checks,
  file.path(validation_dir, "figure2_methylation_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(
  list(
    metadata = methylation_input$metadata,
    species_rho = species_rho,
    panel_h_rho = panel_h_matrix,
    panel_j_rho = panel_j_rho,
    panel_j_counts = panel_j_counts,
    panel_k_clades = panel_k_clades,
    panel_k_rho = panel_k_rho,
    panel_k_target_rank = target_rank,
    target_summary = target_summary,
    checks = figure2_methylation_checks,
    summary = figure2_methylation_summary,
    sequence_tiles = sequence_tiles,
    plots = list(panel_h = panel_h, panel_i = panel_i, panel_j = panel_j, panel_k = panel_k)
  ),
  file.path(data_dir, "figure2_methylation_derived_objects.rds"),
  compress = "xz"
)

ggsave(file.path(figure_dir, "Figure2H_human_methylation_age.pdf"), p_F2H,
       width = 2.4, height = 4.5, units = "in", device = cairo_pdf)
ggsave(file.path(figure_dir, "Figure2H_human_methylation_age.png"), p_F2H,
       width = 2.4, height = 4.5, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "Figure2I_cg10501210_sequence_inventory.pdf"), p_F2I,
       width = 6.1, height = 4.2, units = "in", device = cairo_pdf)
ggsave(file.path(figure_dir, "Figure2I_cg10501210_sequence_inventory.png"), p_F2I,
       width = 6.1, height = 4.2, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "Figure2J_methylation_macroclades.pdf"), p_F2J,
       width = 9.1, height = 4.8, units = "in", device = cairo_pdf)
ggsave(file.path(figure_dir, "Figure2J_methylation_macroclades.png"), p_F2J,
       width = 9.1, height = 4.8, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "Figure2K_boreoeutheria_rank.pdf"), p_F2K,
       width = 3.2, height = 4.5, units = "in", device = cairo_pdf)
ggsave(file.path(figure_dir, "Figure2K_boreoeutheria_rank.png"), p_F2K,
       width = 3.2, height = 4.5, units = "in", dpi = 300, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(figure_dir, "Figure2H_human_methylation_age.svg"), p_F2H,
         width = 2.4, height = 4.5, units = "in", device = svglite::svglite)
  ggsave(file.path(figure_dir, "Figure2I_cg10501210_sequence_inventory.svg"), p_F2I,
         width = 6.1, height = 4.2, units = "in", device = svglite::svglite)
  ggsave(file.path(figure_dir, "Figure2J_methylation_macroclades.svg"), p_F2J,
         width = 9.1, height = 4.8, units = "in", device = svglite::svglite)
  ggsave(file.path(figure_dir, "Figure2K_boreoeutheria_rank.svg"), p_F2K,
         width = 3.2, height = 4.5, units = "in", device = svglite::svglite)
}

if (requireNamespace("cowplot", quietly = TRUE)) {
  p_F2HJK <- cowplot::plot_grid(
    p_F2H, p_F2J, p_F2K, nrow = 1,
    rel_widths = c(1.0, 3.8, 1.35), align = "h", axis = "tb"
  )
  ggsave(file.path(figure_dir, "Figure2_HJK_methylation_combined.pdf"), p_F2HJK,
         width = 14, height = 4.9, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, "Figure2_HJK_methylation_combined.png"), p_F2HJK,
         width = 14, height = 4.9, units = "in", dpi = 300, bg = "white")
}

validation_lines <- c(
  "Figure 2 methylation module validation",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Portable input: ", normalizePath(input_path)),
  paste0("Probe universe: ", length(probe_ids)),
  paste0("Tissue eligibility threshold: n >= ", tissue_n_min),
  "Aggregation: tissue Spearman rho -> arithmetic tissue mean within species -> Fisher-z above species",
  sprintf("Panel H: Human n=%d, target rho=%.12f", raw_n(human), species_rho[target_probe, human]),
  paste0(
    "Panel J: ",
    paste(sprintf("%s n=%d Sp=%d rho=%.6f", panel_j_counts$display_group,
                  panel_j_counts$n, panel_j_counts$species_n,
                  panel_j_rho[target_probe, ]), collapse = "; ")
  ),
  sprintf("Panel K: n=%d, Sp=%d, Human included=TRUE, target rho=%.12f, ascending rank=%d/%d",
          raw_n(all_species), length(all_species), target_k, target_rank,
          sum(is.finite(panel_k_rho))),
  paste0("Panel I sequence rows available from source workbook: ", nrow(alignment)),
  paste0("Panel I caveat: the exact extended alignment/artwork is imported from reference/panel_i_reference.png; ",
         "the generated plot is the reproducible supplied 50-mer sequence inventory."),
  "VALIDATION STATUS: PASS"
)
writeLines(validation_lines, file.path(validation_dir, "figure2_methylation_validation.txt"))
writeLines(capture.output(sessionInfo()), file.path(validation_dir, "figure2_methylation_sessionInfo.txt"))

message("Figure 2 methylation panels regenerated in: ", figure_dir)
message(sprintf("Panel K audit: cg10501210 rho = %.12f; rank = %d/%d",
                target_k, target_rank, sum(is.finite(panel_k_rho))))
