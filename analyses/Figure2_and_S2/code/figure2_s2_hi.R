# Supplementary Figure 2: revised h and compact leave-one-clade-out inset.
# Uses audited coefficients and percentile bootstrap intervals unchanged.
# Source in RStudio or run with Rscript. No command-line arguments required.
suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

if (!exists("figure2_s2_root")) figure2_s2_root <- "."
if (!exists("s2_export")) s2_export <- FALSE
output_dir <- file.path(figure2_s2_root, "figures/supplementary_s2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
human_path <- file.path(figure2_s2_root, "data/figure2_s2/hi/human_bootstrap.tsv")
loo_path <- file.path(figure2_s2_root, "data/figure2_s2/hi/leave_one_out.tsv")
human <- read.delim(human_path, check.names = FALSE, stringsAsFactors = FALSE)
loo <- read.delim(loo_path, check.names = FALSE, stringsAsFactors = FALSE)

# Match the user's current human panel: the 13 estimable tissues only.
# Keratinocyte and Placenta remain in the full source table, with NA rho.
human_tissues <- human[human$level == "tissue" & is.finite(human$rho_cached), ]
aggregate <- human[human$level == "aggregate", ]
pretty_names <- c(BoneMarrow = "Bone marrow", CerebralCortex = "Cerebral cortex", LymphNode = "Lymph node")
human_tissues$tissue_label <- human_tissues$tissue
named <- human_tissues$tissue %in% names(pretty_names)
human_tissues$tissue_label[named] <- unname(pretty_names[human_tissues$tissue[named]])
human_tissues <- human_tissues[order(tolower(human_tissues$tissue_label), method = "radix"), ]
stopifnot(
  nrow(human_tissues) == 13L,
  nrow(aggregate) == 1L,
  human_tissues$tissue_label[[1]] == "Adipose",
  tail(human_tissues$tissue_label, 1) == "Spleen",
  all(is.finite(human_tissues$ci_low)), all(is.finite(human_tissues$ci_high)),
  abs(aggregate$rho_cached - mean(human_tissues$rho_cached)) < 1e-12
)
human_tissues$y <- seq(nrow(human_tissues) + 1L, 2L)
human_tissues$label <- sprintf("%s (n=%d)", human_tissues$tissue_label, human_tissues$n_complete)
aggregate$y <- 0.65
aggregate$label <- "Aggregate"
forest <- rbind(
  human_tissues[, c("level", "tissue", "y", "label", "n_complete", "rho_cached", "ci_low", "ci_high")],
  aggregate[, c("level", "tissue", "y", "label", "n_complete", "rho_cached", "ci_low", "ci_high")]
)

p_h <- ggplot(forest, aes(x = rho_cached, y = y)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.20, colour = "#333333") +
  geom_hline(yintercept = 1.30, linewidth = 0.16, colour = "#CCD6D9") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = y), colour = "#6B6B6B", linewidth = 0.22) +
  geom_point(data = subset(forest, level == "tissue"), aes(size = n_complete),
             shape = 21, fill = "#72BF98", colour = "#333333", stroke = 0.18) +
  geom_point(data = subset(forest, level == "aggregate"),
             shape = 23, fill = "#FFD92F", colour = "#222222", stroke = 0.22, size = 1.5) +
  scale_size_continuous(range = c(0.80, 1.60), guide = "none") +
  scale_x_continuous(breaks = c(-1, -0.5, 0), labels = c("-1", "-0.5", "0"),
                     limits = c(-1, 0.10), expand = expansion(mult = c(0.015, 0))) +
  scale_y_continuous(breaks = forest$y, labels = forest$label,
                     limits = c(0.05, 14.65), expand = expansion(mult = c(0, 0))) +
  labs(x = expression("Spearman's " * rho), y = NULL) +
  theme_classic(base_family = "Helvetica", base_size = 5) +
  theme(
    panel.background = element_rect(fill = "#F5FAFC", colour = NA),
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.22, colour = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.18),
    axis.ticks.length = grid::unit(0.6, "mm"),
    axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = 1.4)),
    axis.text.x = element_text(size = 4.8, colour = "black", margin = margin(t = 1)),
    axis.title.x = element_text(size = 5.2, margin = margin(t = 2)),
    plot.margin = margin(0.5, 1, 0.5, 0.5)
  )

# Compact clade summary: the available space is ~10 mm high. Explicit values
# and ranks retain the result without compressing a second correlation axis.
clade_order <- c("Primates", "Chiroptera", "Carnivora", "Artiodactyla", "Cetacea")
clade_colours <- c(Primates = "#4E86BD", Chiroptera = "#8463AD", Carnivora = "#D55B88",
                   Artiodactyla = "#E69B52", Cetacea = "#42A8A0")
clades <- loo[loo$omission_type == "clade", ]
clades <- clades[match(clade_order, clades$omitted), ]
stopifnot(
  nrow(clades) == 5L,
  identical(clades$omitted, clade_order),
  identical(as.integer(clades$global_rank_ascending_min), c(1L, 1L, 1L, 2L, 2L)),
  all(clades$global_rank_denominator == 37554L)
)
# Align the five summaries to the centres of the clade groups in the user's
# existing species plot (PDF page coordinates minus the inset x = 326 pt).
clades$x <- c(43.15, 88.45, 133.90, 177.55, 220.60)
clades$rho_label <- sprintf("\u03c1 = %.3f", clades$global_rho)
clades$rank_label <- paste("rank", clades$global_rank_ascending_min)
p_i_clade <- ggplot(clades, aes(x = x)) +
  annotate("text", x = 0.53, y = 3.27, label = "Leave one clade out", hjust = 0,
           size = 5.4 / ggplot2::.pt, fontface = "bold", family = "Helvetica", colour = "#333333") +
  geom_segment(aes(x = x - 17, xend = x + 17, y = 2.64, yend = 2.64, colour = omitted),
               linewidth = 0.35) +
  geom_text(aes(y = 2.05, label = omitted, colour = omitted), size = 5.0 / ggplot2::.pt, family = "Helvetica") +
  geom_text(aes(y = 1.21, label = rho_label), size = 5.0 / ggplot2::.pt, family = "Helvetica", colour = "#222222") +
  geom_text(aes(y = 0.40, label = rank_label), size = 4.7 / ggplot2::.pt, family = "Helvetica", colour = "#555555") +
  scale_colour_manual(values = clade_colours, guide = "none") +
  scale_x_continuous(limits = c(0, 238), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3.7), expand = c(0, 0)) +
  theme_void(base_family = "Helvetica") + theme(plot.margin = margin(0, 1, 0, 1))

save_panel <- function(plot, stem, width_mm, height_mm) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".svg")), plot,
         width = width_mm, height = height_mm, units = "mm", device = svglite::svglite, bg = "white")
}
# Species order is inherited from the frozen table, which follows main K.
species_loo <- loo[loo$omission_type == "species", ]
species_loo$omitted_clade <- factor(species_loo$omitted_clade, levels = clade_order)
species_loo$short_name <- sub("^([A-Z])[a-z]+ ", "\\1. ", species_loo$omitted)
species_loo$short_name <- sub("^C[.] lupus familiaris$", "C. l. familiaris", species_loo$short_name)
species_loo$short_name <- sub("^O[.] rosmarus divergens$", "O. r. divergens", species_loo$short_name)
species_loo$short_name <- sub("^C[.] taurinus albojubatus$", "C. t. albojubatus", species_loo$short_name)
species_loo$short_name <- factor(species_loo$short_name, levels = species_loo$short_name)
baseline <- unique(loo$baseline_global_rho)
stopifnot(length(baseline) == 1L, abs(baseline + 0.617636318578119) < 1e-12)
p_S2i_species <- ggplot(species_loo, aes(short_name, global_rho, colour = omitted_clade)) +
  geom_hline(yintercept = baseline, linetype = "dashed", colour = "grey45", linewidth = .2) +
  geom_segment(aes(xend = short_name, y = baseline, yend = global_rho), linewidth = .3) +
  geom_point(aes(fill = omitted_clade), shape = 21, colour = "#333333", stroke = .2, size = 1.1) +
  facet_grid(~omitted_clade, scales = "free_x", space = "fixed") +
  scale_colour_manual(values = clade_colours, guide = "none") +
  scale_fill_manual(values = clade_colours, guide = "none") +
  scale_y_continuous(breaks = c(-.64, -.62, -.60),
                     labels = scales::label_number(accuracy = .01)) +
  coord_cartesian(ylim = c(-.648, -.587)) +
  labs(x = NULL, y = expression(atop("Spearman's " * rho, "(leave one species out)"))) +
  theme_classic(base_size = 5, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 4),
        axis.text = element_text(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = .2),
        axis.ticks = element_line(linewidth = .15),
        axis.ticks.length = grid::unit(.6, "mm"),
        strip.background = element_blank(), strip.text = element_text(size = 4.8),
        panel.spacing = grid::unit(1, "mm"),
        plot.title = element_text(size = 5.4, face = "bold"),
        plot.margin = margin(1, 1, 1, 1))
p_S2h <- p_h
p_S2i_clades <- p_i_clade
p_S2i <- cowplot::plot_grid(p_S2i_clades, p_S2i_species, ncol = 1,
                           rel_heights = c(10, 39))

# Verify the entire leave-one-out analysis directly from the portable
# 37,554-probe tissue matrices, including the balanced Myotis/Tursiops units.
s2_check_loo <- function() {
  x <- readRDS(file.path(figure2_s2_root, "data/methylation/figure2_methylation_portable_input.rds"))
  fm <- function(m) {
    m[!is.finite(m)] <- NA_real_
    m[] <- pmax(-1 + 1e-12, pmin(1 - 1e-12, m))
    tanh(rowMeans(atanh(m), na.rm = TRUE))
  }
  sm <- do.call(cbind, lapply(x$rho_by_species_tissue, function(m) rowMeans(m, na.rm = TRUE)))
  build <- function(m) {
    z <- lapply(names(x$clade_species), function(cl) {
      members <- intersect(x$clade_species[[cl]], colnames(m))
      groups <- x$special_collapses
      grouped <- intersect(unlist(groups), members)
      units <- m[, setdiff(members, grouped), drop = FALSE]
      for (g in groups) {
        keep <- intersect(g, members)
        if (length(keep)) units <- cbind(units, fm(m[, keep, drop = FALSE]))
      }
      fm(units)
    })
    out <- do.call(cbind, z)
    colnames(out) <- names(x$clade_species)
    rownames(out) <- rownames(m)
    out
  }
  all_clades <- build(sm)
  result <- lapply(seq_len(nrow(loo)), function(i) {
    z <- if (loo$omission_type[i] == "species")
      build(sm[, colnames(sm) != loo$omitted[i], drop = FALSE]) else
      all_clades[, colnames(all_clades) != loo$omitted[i], drop = FALSE]
    k <- fm(z)
    at <- match("cg10501210", rownames(z))
    data.frame(omitted = loo$omitted[i], rho = k[at],
               global_rank = rank(k, ties.method = "min", na.last = "keep")[at],
               clade_rank = if (loo$omission_type[i] == "species")
                 rank(z[, loo$omitted_clade[i]], ties.method = "min", na.last = "keep")[at] else NA_real_)
  })
  result <- do.call(rbind, result)
  stopifnot(max(abs(result$rho - loo$global_rho)) < 1e-12,
            all(result$global_rank == loo$global_rank_ascending_min),
            all(result$clade_rank[1:36] == loo$clade_rank_ascending_min[1:36]))
  result
}
loo_verified <- s2_check_loo()
# Reproduce the 95% percentile intervals from the compact one-probe sample
# extract, without loading the original multi-gigabyte beta matrix.
s2_recompute_human_bootstrap <- function(reps = 2000L, seed = 290105L) {
  a <- readRDS(file.path(figure2_s2_root, "data/figure2_s2/hi/cg10501210_sample_annotation_beta.rds"))
  a <- a[a$SpeciesLatinName == "Homo sapiens", , drop = FALSE]
  counts <- table(a$Tissue)
  eligible <- names(counts[counts >= 10L])
  draws <- matrix(NA_real_, reps, length(eligible), dimnames = list(NULL, eligible))
  point <- setNames(rep(NA_real_, length(eligible)), eligible)
  set.seed(seed)
  for (t in eligible) {
    d <- a[a$Tissue == t & is.finite(a$Age) & is.finite(a$target_beta), ]
    if (nrow(d) >= 3L && length(unique(d$Age)) >= 2L && length(unique(d$target_beta)) >= 2L)
      point[t] <- cor(d$Age, d$target_beta, method = "spearman")
    if (nrow(d) >= 20L && length(unique(d$Age)) >= 5L && length(unique(d$target_beta)) >= 5L) {
      for (b in seq_len(reps)) {
        idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
        draws[b, t] <- suppressWarnings(cor(d$Age[idx], d$target_beta[idx], method = "spearman"))
      }
    }
  }
  keep <- names(point)[is.finite(point)]
  qs <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < .9 * reps) return(c(NA_real_, NA_real_))
    unname(quantile(x, c(.025, .975), type = 6))
  }
  cis <- t(apply(draws, 2, qs))
  tab <- data.frame(tissue = eligible, rho = unname(point), ci_low = cis[, 1], ci_high = cis[, 2])
  agg_draws <- rowMeans(draws[, keep, drop = FALSE])
  ag_ci <- qs(agg_draws)
  tab <- rbind(tab, data.frame(tissue = "Human unweighted tissue mean", rho = mean(point, na.rm = TRUE),
                               ci_low = ag_ci[1], ci_high = ag_ci[2]))
  original <- human[match(tab$tissue, human$tissue), ]
  stopifnot(isTRUE(all.equal(tab$rho, original$rho_cached, tolerance = 1e-12)),
            isTRUE(all.equal(tab$ci_low, original$ci_low, tolerance = 1e-12)),
            isTRUE(all.equal(tab$ci_high, original$ci_high, tolerance = 1e-12)))
  list(summary = tab, draws = draws, aggregate_draws = agg_draws, seed = seed, reps = reps)
}
human_bootstrap_recomputed <- s2_recompute_human_bootstrap()
s2_hi_data <- list(human_source = human, human_display = forest,
                   human_bootstrap_recomputed = human_bootstrap_recomputed,
                   leave_one_out = loo, leave_one_out_recomputed = loo_verified,
                   clade_display = clades, species_display = species_loo)
s2_hi_checks <- data.frame(
  check = c("Human: 13 estimable tissues", "Human: alphabetical tissue order",
            "Human: aggregate is the arithmetic tissue mean", "K: 41 leave-one-out recomputations match",
            "LOSO: global rank 1 in all 36 omissions", "LOCO: global ranks 1, 1, 1, 2, 2",
            "Human: all 2000-replicate bootstrap CIs reproduced from sample data"),
  status = rep("PASS", 7))
if (isTRUE(s2_export)) {
  save_panel(p_S2h, "S2h_human_tissues", 40, 36)
  save_panel(p_S2i_clades, "S2i_clade_strip", 84, 10)
  save_panel(p_S2i_species, "S2i_species_omissions", 84, 39)
  save_panel(p_S2i, "S2i_leave_one_out", 84, 49)
  write.table(forest, file.path(output_dir, "S2h_plot_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(loo, file.path(output_dir, "S2i_plot_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
