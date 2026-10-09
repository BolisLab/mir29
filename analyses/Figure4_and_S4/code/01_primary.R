# Figure 4 primary methylome panels. Frozen checkpoint replay.
# ---- 4a
pca_data <- cache("pca", function() {
  meta <- metadata()
  x <- as.matrix(read_input("pan_mval_all.june2026.RDS"))[, meta$sample, drop = FALSE]
  design <- model.matrix(~ TIME * GROUP, meta)
  x <- limma::removeBatchEffect(x, batch = meta$EXPERIMENT, design = design)
  x <- limma::removeBatchEffect(x, batch = meta$BATCH, design = design)
  variance <- matrixStats::rowVars(x, na.rm = TRUE)
  x <- x[order(variance, decreasing = TRUE)[seq_len(20000)], , drop = FALSE]
  x <- x[rowSums(!is.finite(x)) == 0, , drop = FALSE]
  fit <- prcomp(t(x), center = TRUE, scale. = TRUE)
  percent <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  # A PCA axis sign is arbitrary. Orient INH to the right, as in the artwork.
  scores <- fit$x[, 1:2]
  if (mean(scores[meta$GROUP == "INH", 1]) < mean(scores[meta$GROUP == "MIM", 1])) scores[, 1] <- -scores[, 1]
  list(samples = cbind(meta, scores), percent = percent,
       probes = rownames(x), rotation = fit$rotation[, 1:2])
})
print(pca_data$percent[1:4])
save_table(pca_data$samples, "4a_PCA_scores")
z <- pca_data$samples
low <- c(CTRL = "#ededed", INH = "#B9DAB5", MIM = "#E4C5E0")
high <- c(CTRL = "#515151", INH = "#00736D", MIM = "#841452")
z$fill <- vapply(seq_len(nrow(z)), function(i) colorRampPalette(c(low[z$GROUP[i]], high[z$GROUP[i]]))(6)[z$TP2[i]], character(1))
z$shape <- ifelse(z$GROUP == "CTRL", ifelse(z$EXPERIMENT == "inh", 22, 23), 21)
p <- ggplot(z, aes(PC1, PC2)) +
  stat_ellipse(aes(group = GROUP, color = GROUP), type = "norm", level = .99,
               linewidth = .35, show.legend = FALSE) +
  geom_point(aes(fill = fill, shape = shape), size = 2.2, stroke = .25) +
  scale_fill_identity() + scale_shape_identity() + scale_color_manual(values = COLORS) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_data$percent[1]),
       y = sprintf("PC2 (%.1f%%)", pca_data$percent[2])) + theme_figure() +
  theme(aspect.ratio = .6)
save_panel(p, "4a_PCA", 5.2, 3.6)


# ---- 4b
# unrestricted panT1 count table in the same file.
# 1. Select temporal probes within each condition using Time BH < .05.
# 2. For each timepoint-vs-CTRL-T1 contrast, recompute BH only within that set.
# 3. Count positive/negative effects at this restricted BH < .05.
# This is a historical, data-selected test family, NOT genome-wide FDR.
counts <- cache("panel_b_historical_restricted_counts", function() {
  r <- read_input("methy_res.notcollapsed_2026.RDS")
  sets <- lapply(r$Time, function(z) rownames(z)[which(z$adj.P.Val < .05)])
  out <- bind_rows(lapply(names(r$panT1C), function(nm) {
    group <- toupper(sub(".*\\.", "", nm)); if (group == "CTR") group <- "CTRL"
    ids <- sets[[paste0("TIME_", tolower(group))]]
    z <- r$panT1C[[nm]][ids, , drop = FALSE]
    q <- p.adjust(z$P.Value, "BH")
    data.frame(GROUP = group, TP2 = as.integer(sub(".*T([1-6])vs.*", "\\1", nm)),
               n_temporal_set = length(ids), hyper = sum(q < .05 & z$logFC > 0, na.rm = TRUE),
               hypo = sum(q < .05 & z$logFC < 0, na.rm = TRUE))
  }))
  bind_rows(out, data.frame(GROUP = "CTRL", TP2 = 1, n_temporal_set = length(sets$TIME_ctrl), hyper = 0, hypo = 0))
})
save_table(counts, "panel_b_historical_temporal_subset_counts")
print(counts)
z <- counts %>% pivot_longer(c(hyper, hypo), names_to = "direction", values_to = "n") %>%
  mutate(GROUP = factor(GROUP, c("CTRL", "INH", "MIM")),
         y = 7 - TP2 + c(CTRL = .20, INH = 0, MIM = -.20)[GROUP])
parts <- lapply(c("hypo", "hyper"), function(direction) {
  p <- ggplot(z[z$direction == direction, ], aes(n, y, fill = GROUP)) +
    geom_col(orientation = "y", width = .17, color = "#444444", linewidth = .2) +
    scale_fill_manual(values = COLORS) + scale_y_continuous(breaks = 1:6, labels = paste0("T", 6:1)) +
    labs(x = if (direction == "hypo") "Hypomethylated" else "Hypermethylated", y = NULL) +
    theme_figure() + theme(legend.position = "none")
  if (direction == "hypo") p <- p + scale_x_reverse(labels = scales::label_number(scale = 1e-3, suffix = "k"))
  else p <- p + scale_x_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
})
save_panel(cowplot::plot_grid(plotlist = parts, nrow = 1, align = "h"), "4b_temporal_DMP_counts", 4.4, 3.5)


# ---- 4c
r <- list(counts = input("primary/limma_contrast_counts.rds"))
d <- subset(r$counts, section == "Time")
d$GROUP <- c("CTRL", "INH", "MIM")[match(d$contrast, c("TIME_ctrl", "TIME_inh", "TIME_mim"))]
d$total <- d$hyper + d$hypo; d$ratio <- ifelse(d$hypo > 0, d$hyper / d$hypo, NA)
save_table(d, "panel_c_temporal_DMP_totals")
bars <- d %>% pivot_longer(c(hyper, hypo), names_to = "direction", values_to = "count")
p1 <- ggplot(bars, aes(count, factor(GROUP, c("MIM", "INH", "CTRL")), fill = direction)) +
  geom_col(width = .65, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c(hyper = "#DC3237", hypo = "#4878BA")) +
  geom_text(data = d, aes(total, GROUP, label = format(total, big.mark = ",")),
            inherit.aes = FALSE, hjust = -.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, .3))) +
  labs(x = "DMPs", y = NULL) + theme_figure()
p2 <- ggplot(subset(d, is.finite(ratio)), aes(ratio, factor(GROUP, c("MIM", "INH", "CTRL")), color = GROUP)) +
  geom_segment(aes(x = 1, xend = ratio, yend = GROUP), linewidth = .45) +
  geom_point(size = 2) + geom_text(aes(label = sprintf("%.1fx", ratio)), hjust = -.3, size = 3) +
  scale_color_manual(values = COLORS) + scale_x_log10(limits = c(.9, 20)) +
  labs(x = "Hyper / hypo", y = NULL) + theme_figure() + theme(legend.position = "none")
save_panel(cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(3, 1)), "4c_temporal_DMP_summary", 5.5, 2.2)



# ---- 4d
# historical actual-time mixed model and its 90% confidence intervals.
source("code/upstream/lib_global_methylation.R")
baseline_time_hours <- 72
time_values <- c(72, 96, 168, 192, 264, 288)
confidence_level <- .90
plot_group_levels <- c("CTRL", "Inhibitor", "Mimic")
sample_global_metrics <- cache("global_raw_collapsed_metrics", function() {
  b <- as.matrix(read_input("pan_beta_all.june2026.RDS"))
  ids <- base_probe(rownames(b))
  # Mean of available replicate probes, on beta scale (not inverse-logit of M).
  sums <- rowsum(replace(b, is.na(b), 0), ids, reorder = FALSE)
  counts <- rowsum(1L * !is.na(b), ids, reorder = FALSE)
  b <- sums / counts
  rm(sums, counts); gc()
  complete <- rowSums(is.finite(b)) == ncol(b)
  x <- b[complete, , drop = FALSE]
  saveRDS(list(n_collapsed = nrow(b), n_complete = sum(complete),
               retained_probes = rownames(x)), "data/primary/global_mean_probe_universe.rds")
  meta <- metadata()
  out <- data.frame(sample = colnames(x), mean_beta = colMeans(x))
  out <- left_join(out, meta, by = "sample")
  out$mGROUP <- factor(out$mGROUP, c("inh_CTRL", "inh_INH", "mim_CTRL", "mim_MIM"))
  out$TIME_HOURS <- out$TP3
  out$TIME_factor <- factor(out$TP3, levels = time_values)
  out$TIME_FROM_T1_DAYS <- (out$TP3 - 72) / 24
  out$BLOCK <- interaction(out$EXPERIMENT, out$BATCH, drop = TRUE)
  out
})
fit <- cache("global_mean_mixed_model", function() fit_own_baseline_metric(sample_global_metrics, "mean_beta"))
save_table(sample_global_metrics, "panel_d_sample_metrics")
save_table(fit$trajectory, "panel_d_model_adjusted_delta_beta")
save_table(fit$slope_contrasts, "panel_d_slope_contrasts")
z <- fit$trajectory
z$x <- (z$TIME_HOURS - 72) / 24
colors <- c(CTRL = COLORS["CTRL"], Inhibitor = COLORS["INH"], Mimic = COLORS["MIM"])
names(colors) <- plot_group_levels
p <- ggplot(z, aes(x, estimate * 1000, color = GROUP, fill = GROUP, group = GROUP)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey65", linewidth = .3) +
  geom_ribbon(aes(ymin = lower.CL * 1000, ymax = upper.CL * 1000), alpha = .12, color = NA) +
  geom_line(linetype = "dashed", linewidth = .25) +
  geom_errorbar(aes(ymin = lower.CL * 1000, ymax = upper.CL * 1000), width = .12, linewidth = .3) +
  geom_smooth(method = "lm", formula = y ~ 0 + x + I(x^2), se = FALSE, linewidth = .8) +
  geom_point(size = 1.6) + scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = (time_values - 72) / 24, labels = paste0("T", 1:6)) +
  labs(x = NULL, y = expression(Delta~"mean methylation ("*10^3*")")) + theme_figure()
save_panel(p, "4d_global_methylation", 5.5, 2.5)


# ---- 4e_f_h
d <- readRDS("data/primary/hgf_density_adjusted_no_reflection.rds") %>% group_by(GROUP, TP2, beta) %>%
  summarise(mean = mean(density), sd = sd(density), .groups = "drop")
p <- ggplot(d, aes(beta, mean, color = GROUP, fill = GROUP)) +
  geom_ribbon(aes(ymin = pmax(0, mean - sd), ymax = mean + sd), alpha = .1, color = NA) +
  geom_line(linewidth = .55) + facet_wrap(~paste0("T", TP2), nrow = 1) +
  scale_color_manual(values = COLORS) + scale_fill_manual(values = COLORS) +
  labs(x = expression(beta), y = "Density") + theme_figure()
save_panel(p, "4e_beta_distributions", 5.0, 2.9)
metrics <- readRDS("data/primary/hgf_compartment_metrics_adjusted_verified_annotation.rds")
d <- subset(metrics, comp == "All_CpGs" & TP2 == 6) %>% pivot_longer(c(entropy, FAE))
d$value[d$name == "FAE"] <- 100 * d$value[d$name == "FAE"]
p <- ggplot(d, aes(GROUP, value, fill = GROUP)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .55, linewidth = .3) +
  geom_point(position = position_jitter(width = .09, height = 0, seed = 42), size = 1.3, stroke = 0) +
  facet_wrap(~name, scales = "free_y", nrow = 1) + scale_fill_manual(values = COLORS) +
  labs(x = NULL, y = NULL, title = "T6") + theme_figure() + theme(legend.position = "none")
save_panel(p, "4f_entropy_FAE", 4.6, 2.6)

compartment_order <- c("All_CpGs", "OpenSea", "Shelf", "Shore", "Island", "Promoter_CGI")
heatmap_data <- metrics %>% group_by(comp, GROUP, TP2) %>%
  summarise(entropy = mean(entropy), mean_beta = mean(mean_beta), .groups = "drop") %>%
  pivot_longer(c(entropy, mean_beta), names_to = "metric") %>% group_by(metric, comp) %>%
  mutate(z = as.numeric(scale(value))) %>% ungroup() %>%
  mutate(x = TP2 + 6.6 * (as.integer(GROUP) - 1),
         y = match(comp, rev(compartment_order)))
save_table(heatmap_data, "panel_j_heatmap_values")
pal <- scico::scico(256, palette = "vik")
p <- ggplot(heatmap_data, aes(x, y, fill = z)) +
  geom_tile(width = .98, height = .96, color = "black", linewidth = .3) +
  facet_wrap(~factor(metric, c("entropy", "mean_beta"), c("Entropy", "Mean beta")), ncol = 1) +
  scale_fill_gradient2(low = pal[1], mid = "white", high = pal[256], limits = c(-2, 2), oob = scales::squish) +
  scale_x_continuous(breaks = c(1:6, 1:6 + 6.6, 1:6 + 13.2), labels = rep(1:6, 3), expand = expansion(add = .5)) +
  scale_y_continuous(breaks = 1:6, labels = gsub("_", " ", rev(compartment_order)), expand = expansion(add = .5)) +
  labs(x = "CTRL                       INH                       MIM", y = NULL, fill = "Row z") +
  theme_figure() + theme(axis.line = element_blank(), axis.ticks = element_blank())
save_panel(p, "4h_compartment_heatmaps", 5.5, 4.1)


# ---- 4g
# but chronological AGE for the density panel. Both scatter candidates are
# exported and labeled explicitly; this ambiguity must not be hidden.
source("code/upstream/lib_gtex_density.R")
suppressPackageStartupMessages(library(mgcv))
annotations <- NULL # Only needed for optional upstream density regeneration.
metrics <- cache("gtex_sample_metrics", function() {
  b <- read_input("gtex.methylation.RDS")
  out <- bind_rows(lapply(names(annotations), function(tissue) {
    ann <- annotations[[tissue]]
    ids <- intersect(rownames(ann), colnames(b))
    bind_rows(lapply(ids, function(id) {
      v <- as.numeric(b[, id]); v[v < 0] <- 0; v[v > 1] <- 1
      p <- hist(v, breaks = seq(0, 1, .05), plot = FALSE)$counts
      p <- p / sum(p) + 1e-8
      data.frame(sample = id, tissue, AGE = ann[id, "AGE"], DNAmAge = ann[id, "DNAmAge"],
                 entropy = -sum(p * log(p)), FAE = mean(v < .2 | v > .8, na.rm = TRUE),
                 mean_beta = mean(v, na.rm = TRUE))
    }))
  }))
  rm(b); gc()
  out
})
save_table(metrics, "panel_g_GTEx_sample_metrics")
make_windows <- function(d, age_column) {
  d$age <- d[[age_column]]; d <- d[is.finite(d$age), ]
  design <- model.matrix(~age, d)
  d$entropy_adj <- as.numeric(limma::removeBatchEffect(matrix(d$entropy, 1), batch = d$tissue, design = design))
  z <- qlogis(pmax(1e-6, pmin(1-1e-6, d$FAE)))
  d$FAE_adj <- plogis(as.numeric(limma::removeBatchEffect(matrix(z, 1), batch = d$tissue, design = design)))
  breaks <- seq(min(d$age), max(d$age), length.out = 21)
  windows <- bind_rows(lapply(split(d, d$tissue), function(tis) {
    tis <- tis[order(tis$age), ]
    if (nrow(tis) < 15) return(NULL)
    bind_rows(lapply(seq(1, nrow(tis)-14, by = 5), function(i) {
      x <- tis[i:(i+14), ]
      ac <- mean(x$age)
      data.frame(tissue = tis$tissue[1], age = ac, age_bin = max(1, min(20, findInterval(ac, breaks, rightmost.closed = TRUE))),
                 entropy = mean(x$entropy_adj), FAE = mean(x$FAE_adj))
    }))
  }))
  keep <- names(which(table(windows$age_bin) >= 6))
  windows[windows$age_bin %in% keep, ]
}
plot_scatter <- function(windows, age_column) {
  z <- pivot_longer(windows, c(entropy, FAE), names_to = "metric", values_to = "value")
  ggplot(z, aes(age, value)) + geom_point(color = "grey45", alpha = .5, size = .65, stroke = 0) +
    geom_smooth(aes(color = metric, fill = metric), method = "gam", formula = y ~ s(x, k = 3),
                se = TRUE, linewidth = .7) +
    facet_wrap(~factor(metric, c("entropy", "FAE")), nrow = 1, scales = "free_y") +
    scale_color_manual(values = c(entropy = "#D8343B", FAE = "#386FA8")) +
    scale_fill_manual(values = c(entropy = "#D8343B", FAE = "#386FA8")) +
    labs(x = if (age_column == "AGE") "Chronological age (years)" else "DNAm age (years)", y = NULL) +
    theme_figure() + theme(legend.position = "none")
}
for (age_column in "DNAmAge") {
  windows <- make_windows(metrics, age_column)
  save_table(windows, paste0("panel_g_windows_", age_column))
  # The original requested k=1, which mgcv raises to its minimum k=3.
  p <- plot_scatter(windows, age_column)
  save_panel(p, paste0("panel_g_top_", age_column), 5.4, 2.6)
}
density_data <- cache("gtex_meta_density", function() {
  b <- as.matrix(read_input("gtex.methylation.RDS"))
  set.seed(20260917)
  meta_curves_fast_absolute(b, annotations, win_size = 15, step_size = 5,
                           target_n = 20000, nbins_beta = 256, n_age_bins = 20,
                           smooth_sigma = .010, min_windows_per_bin = 2)
})
save_table(density_data, "panel_g_beta_density_by_chronological_age")
p <- ggplot(density_data, aes(beta, density, color = age_mid, group = age_bin)) +
  geom_line(linewidth = .55) + scale_color_viridis_c(option = "plasma", name = "Age") +
  coord_cartesian(xlim = c(0, 1)) + labs(x = expression("Methylation "*beta), y = "Density") + theme_figure()
save_panel(p, "panel_g_bottom_age_density", 5.4, 2.5)
top <- plot_scatter(make_windows(metrics, "DNAmAge"), "DNAmAge")
save_panel(cowplot::plot_grid(top, p, ncol = 1, rel_heights = c(1, .85)),
           "4g_GTEx_context", 5.4, 4.8)


# ---- 4i_left
out <- input("primary/panel_k_compartment_burden.rds")$results
sets <- input("primary/panel_k_compartment_burden.rds")$feature_sets
out$y <- match(out$comp, rev(names(sets))) + c(CTRL = .17, INH = 0, MIM = -.17)[out$GROUP]
p <- ggplot(out, aes(net, y, color = GROUP, shape = GROUP)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#E9EFFA", alpha = .55) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F8DADD", alpha = .55) +
  geom_vline(xintercept = 0, linewidth = .3, color = "grey50") +
  geom_segment(aes(x = 0, xend = net, yend = y), linewidth = .3) +
  geom_point(aes(size = burden)) + scale_color_manual(values = COLORS) +
  scale_shape_manual(values = c(CTRL = 16, INH = 15, MIM = 17)) +
  scale_size_continuous(range = c(1, 3.5), guide = "none") +
  scale_y_continuous(breaks = 1:5, labels = rev(paste0(names(sets), "\nn = ",
                           format(lengths(sets), big.mark = ",", trim = TRUE)))) +
  labs(x = "Net temporal DMP burden (%)", y = NULL) + theme_figure()
save_panel(p, "4i_left_compartment_burden", 4.5, 3.4)


# ---- 4j_left
# Historical overlap semantics include probe strand and exact collapsed
# repClass/repFamily strings; multi-label entries are not silently expanded.
suppressPackageStartupMessages(library(GenomicRanges))
membership <- NULL # Stored membership was already used in the per-sample checkpoint.
spec <- data.frame(label = c("LTR", "ERVL", "ERV1", "ERVK", "LINE", "L1", "L2", "DNA/hAT", "Charlie", "Tip100", "SINE/Alu"),
                   field = c("repClass", rep("repFamily", 3), "repClass", rep("repFamily", 2), "repClass", rep("repFamily", 2), "repClass"),
                   block = c(rep("LTR", 4), rep("LINE", 3), rep("DNA/hAT", 3), "SINE/Alu"))
values <- list("LTR", c("ERVL", "ERVL-MaLR"), "ERV1", "ERVK", "LINE", "L1", "L2", "DNA", "hAT-Charlie", "hAT-Tip100", "SINE")
repeat_data <- cache("repeat_beta_metrics_historical_specification", function() {
  b <- as.matrix(read_input("pan_beta_all.adj.june2026.RDS"))
  bind_rows(lapply(seq_len(nrow(spec)), function(i) {
    probes <- intersect(membership$probe[membership[[spec$field[i]]] %in% values[[i]]], rownames(b))
    if (length(probes) < 50) return(NULL)
    data.frame(sample = colnames(b), label = spec$label[i], block = spec$block[i],
               n_probes = length(probes), beta = colMeans(b[probes, , drop = FALSE], na.rm = TRUE))
  })) %>% left_join(metadata(), by = "sample")
})
z <- repeat_data %>% group_by(label, block, GROUP, TP2) %>%
  summarise(value = mean(beta), n_probes = dplyr::first(n_probes), .groups = "drop") %>%
  group_by(label) %>% mutate(z = as.numeric(scale(value))) %>% ungroup() %>%
  mutate(x = TP2 + 6.6 * (as.integer(GROUP) - 1),
         y = length(spec$label) + 1 - match(label, spec$label))
save_table(z, "4j_left_repeat_heatmap")
save_table(unique(z[, c("label", "n_probes")]), "panel_i_left_repeat_probe_counts")
pal <- scico::scico(256, palette = "vik")
p <- ggplot(z, aes(x, y, fill = z)) +
  geom_tile(width = .98, height = .95, color = "black", linewidth = .3) +
  scale_fill_gradient2(low = pal[1], mid = "white", high = pal[256], limits = c(-2, 2), oob = scales::squish) +
  scale_x_continuous(breaks = c(1:6, 1:6 + 6.6, 1:6 + 13.2), labels = rep(1:6, 3), expand = expansion(add = .5)) +
  scale_y_continuous(breaks = 1:length(spec$label), labels = rev(spec$label), expand = expansion(add = .5)) +
  labs(x = "CTRL                       INH                       MIM", y = NULL, fill = "Row z") +
  theme_figure() + theme(axis.line = element_blank(), axis.ticks = element_blank())
save_panel(p, "4j_left_repeat_heatmap", 5.0, 4.1)


# ---- 4l
suppressPackageStartupMessages(library(mgcv))
objects <- cache("ctrl_mover_trajectories", function() {
  meta <- metadata()
  meta$bio_rep <- factor(stringr::str_extract(meta$SHORT_ID, "P\\d+|exp\\d+"))
  if (anyNA(meta$bio_rep)) stop("Historical biological replicate parser did not match all samples.")
  r <- readRDS("data/primary/results_notcollapsed.rds")$Time$TIME_ctrl
  ids <- rownames(r)[which(r$adj.P.Val < .05 & abs(r$logFC) > 0)]
  direction <- setNames(ifelse(r[ids, "logFC"] > 0, "Hyper", "Hypo"), ids)
  b <- as.matrix(read_input("pan_beta_all.adj.june2026.RDS"))[ids, meta$sample, drop = FALSE]
  # Subsetting before conversion to long format is exactly equivalent to the
  # original post-pivot subset and avoids materializing tens of millions of rows.
  for (condition in levels(meta$GROUP)) {
    idx <- which(meta$GROUP == condition)
    baseline <- rowMeans(b[, meta$GROUP == condition & meta$TP2 == 1, drop = FALSE], na.rm = TRUE)
    b[, idx] <- b[, idx, drop = FALSE] - baseline
  }
  dat <- bind_rows(lapply(c("Hypo", "Hyper"), function(d) {
    data.frame(sample = colnames(b), direction = d,
               value = colMeans(b[names(direction)[direction == d], , drop = FALSE], na.rm = TRUE))
  })) %>% left_join(meta, by = "sample") %>%
    group_by(bio_rep, direction, GROUP, TP2) %>% summarise(rep_mean = mean(value), .groups = "drop")
  fits <- dat %>% group_by(direction, GROUP) %>% group_modify(~ {
    fit <- gam(rep_mean ~ s(TP2, k = 4) + s(bio_rep, bs = "re"), data = .x, method = "REML", select = TRUE)
    grid <- data.frame(TP2 = seq(1, 6, length.out = 200), bio_rep = .x$bio_rep[1])
    pred <- predict(fit, grid, se.fit = TRUE, exclude = "s(bio_rep)")
    data.frame(TP2 = grid$TP2, fit = pred$fit, lower = pred$fit - qnorm(.95) * pred$se.fit,
               upper = pred$fit + qnorm(.95) * pred$se.fit)
  }) %>% ungroup()
  list(probes = direction, replicate_means = dat, fitted = fits)
})
save_table(objects$replicate_means, "panel_l_replicate_means")
save_table(objects$fitted, "panel_l_GAM_fitted_curves")
p <- ggplot(objects$replicate_means, aes(TP2, rep_mean, color = GROUP)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = .3, color = "grey60") +
  geom_line(aes(group = interaction(GROUP, bio_rep)), alpha = .12, linewidth = .35) +
  geom_point(aes(shape = bio_rep), size = 1.1) +
  geom_ribbon(data = objects$fitted, aes(x = TP2, ymin = lower, ymax = upper, fill = GROUP),
              inherit.aes = FALSE, alpha = .15) +
  geom_line(data = objects$fitted, aes(y = fit), linewidth = .7) +
  facet_wrap(~factor(direction, c("Hypo", "Hyper")), scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  scale_color_manual(values = COLORS) + scale_fill_manual(values = COLORS) +
  labs(x = NULL, y = expression(Delta*beta~"from own T1")) + theme_figure()
save_panel(p, "4l_CTRL_temporal_movers", 5.8, 3.8)


# ---- 4m_S4f
# P < .05 supports direction; interaction significance is not a class gate.
r <- NULL # Needed only for the optional upstream cache calculation.
standardize <- function(x, type = "DNAm") {
  x <- as.data.frame(x)
  if (type == "DNAm") data.frame(id = rownames(x), slope = x$logFC, p = x$P.Value, q = x$adj.P.Val)
  else data.frame(id = rownames(x), slope = x$log2FoldChange, p = x$pvalue, q = x$padj)
}
pair <- function(ctrl, trt, type, treatment) {
  x <- merge(standardize(ctrl, type), standardize(trt, type), by = "id", suffixes = c("_ctrl", "_trt"))
  sig <- function(x) is.finite(x) & x < .05
  x <- x[(sig(x$q_ctrl) | sig(x$q_trt)) & is.finite(x$slope_ctrl) & is.finite(x$slope_trt), ]
  x$layer <- type; x$treatment <- treatment
  support_c <- sig(x$p_ctrl) & abs(x$slope_ctrl) > 1e-12
  support_t <- sig(x$p_trt) & abs(x$slope_trt) > 1e-12
  same <- sign(x$slope_ctrl) == sign(x$slope_trt)
  larger <- abs(x$slope_trt) > abs(x$slope_ctrl) + 1e-12
  x$change <- (abs(x$slope_trt) - abs(x$slope_ctrl)) / pmax(abs(x$slope_trt), abs(x$slope_ctrl), 1e-12)
  x$detailed_class <- case_when(
    support_c & !support_t ~ "Attenuation_lost",
    !support_c & support_t & same & larger ~ "Amplification_emergent_concordant",
    !support_c & support_t ~ "Recruitment",
    support_c & support_t & !same ~ "Rerouting",
    support_c & support_t & same & x$change < -.2 ~ "Attenuation_weaker",
    support_c & support_t & same & x$change > .2 ~ "Amplification_shared",
    support_c & support_t & same & abs(x$change) <= .2 ~ "Preserved_shared",
    TRUE ~ "Unclassified")
  x$class <- case_when(grepl("^Attenuation", x$detailed_class) ~ "Attenuated",
                       grepl("^Amplification", x$detailed_class) ~ "Amplified",
                       grepl("^Preserved", x$detailed_class) ~ "Preserved",
                       x$detailed_class %in% c("Recruitment", "Rerouting") ~ "Rewired", TRUE ~ "Unclassified")
  x
}
geometries <- cache("temporal_slope_geometries", function() {
  out <- list(pair(r$Time$TIME_ctrl, r$Time$TIME_inh, "DNAm", "INH"),
              pair(r$Time$TIME_ctrl, r$Time$TIME_mim, "DNAm", "MIM"))
  for (arm in c("inh", "mim")) {
    z <- readRDS(file.path(Sys.getenv("FIG4_RNA_DIR", "external_inputs"), paste0("res_", arm, ".RDS")))
    saveRDS(z[c("res_time_ctrl", "res_time_trt", "res_interaction")], paste0("data/primary/rna_", arm, "_temporal_results.rds"))
    out[[length(out) + 1]] <- pair(z$res_time_ctrl, z$res_time_trt, "RNA", toupper(arm))
  }
  bind_rows(out)
})
stopifnot(!any(geometries$class == "Unclassified"))
summary <- geometries %>% count(layer, treatment, class) %>% group_by(layer, treatment) %>%
  mutate(total = sum(n), percent = 100 * n / total) %>% ungroup()
save_table(summary, "S4f_trajectory_outcomes")
print(summary)
scatter_stats <- geometries %>% filter(layer == "DNAm") %>% group_by(treatment) %>%
  summarise(n = n(), slope_ratio = coef(lm(slope_trt ~ 0 + slope_ctrl))[1], .groups = "drop")
save_table(scatter_stats, "panel_m_slope_statistics")
print(scatter_stats)
plots <- lapply(c("INH", "MIM"), function(arm) {
  x <- subset(geometries, layer == "DNAm" & treatment == arm)
  stat <- subset(scatter_stats, treatment == arm)
  main <- ggplot(x, aes(slope_ctrl, slope_trt)) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey65", linewidth = .3) +
    geom_vline(xintercept = 0, linetype = 3, color = "grey65", linewidth = .3) +
    geom_abline(slope = 1, linetype = 2, linewidth = .35) +
    geom_point(color = COLORS[arm], alpha = .12, size = .25, stroke = 0) +
    geom_abline(slope = stat$slope_ratio, color = COLORS[arm], linewidth = .6) +
    annotate("text", x = .45, y = -.43, label = sprintf("SR = %.2f", stat$slope_ratio), hjust = 1, color = COLORS[arm], size = 3) +
    coord_equal(xlim = c(-.5, .5), ylim = c(-.5, .5)) + labs(x = "Slope CTRL", y = paste("Slope", arm)) + theme_figure()
  top <- ggplot(x %>% pivot_longer(c(slope_ctrl, slope_trt)), aes(value, factor(name, c("slope_trt", "slope_ctrl")))) +
    geom_boxplot(width = .55, outlier.shape = NA, fill = COLORS[arm], alpha = .35, linewidth = .3) +
    coord_cartesian(xlim = c(-.5, .5)) + labs(title = paste0("n = ", format(nrow(x), big.mark = ",")), x = NULL, y = NULL) +
    theme_figure(8) + theme(axis.text = element_blank(), axis.ticks = element_blank())
  cowplot::plot_grid(top, main, ncol = 1, rel_heights = c(.5, 2.6), align = "v")
})
save_panel(cowplot::plot_grid(plotlist = plots, nrow = 1), "4m_temporal_slope_scatter", 6.2, 3.9)
summary$class <- factor(summary$class, c("Attenuated", "Preserved", "Rewired", "Amplified"))
p <- ggplot(summary, aes(percent, factor(treatment, c("MIM", "INH")), fill = class)) +
  geom_col(width = .6, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(percent >= 5, round(percent), "")),
            position = position_stack(vjust = .5, reverse = TRUE), size = 3, color = "white") +
  geom_text(data = unique(summary[, c("layer", "treatment", "total")]),
            aes(x = 102, y = factor(treatment, c("MIM", "INH")),
                label = format(total, big.mark = ",")), inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_x_continuous(limits = c(0, 124), breaks = c(0, 50, 100)) +
  facet_wrap(~layer, ncol = 1) + scale_fill_manual(values = c(Attenuated = "#6F809F", Preserved = "#E7E7E7", Rewired = "#989C9D", Amplified = "#A73558")) +
  labs(x = "Temporal trajectory outcomes (%)", y = NULL) + theme_figure()
save_panel(p, "S4f_trajectory_outcomes", 5.2, 3.0)


# ---- 4n
# The original script left MVALS externally assigned. This reconstruction
# explicitly records the candidate matrix; the audit compares its geometry.
source("code/upstream/lib_methylation_walk.R")
TIME_LEVELS <- 1:6
BASELINE_TIME <- 1L
CONDITIONS <- c("inh_CTRL", "inh_INH", "mim_CTRL", "mim_MIM")
PLOT_SCALE <- 1000
SMOOTH_CURVE_N <- 300L
variant <- Sys.getenv("FIG4_WALK_MATRIX", "pan_mval_all.adj_26june2026.RDS")
key <- paste0("walk_", sub("\\.RDS$", "", variant))
obj <- cache(key, function() {
  x <- as.matrix(read_input(variant))
  meta <- metadata(); meta <- meta[match(colnames(x), meta$sample), ]
  meta$condition_walk <- meta$mGROUP; meta$time_walk <- meta$TP2
  probes <- select_top_variable(x, colnames(x), 50000)
  out <- make_walk(x[probes, , drop = FALSE], meta, "top50k_all", "Top 50,000 variable CpGs", "RMS")
  out$matrix_file <- variant; out$probes <- probes
  out
})
save_table(obj$metrics, paste0("panel_o_metrics_", key))
save_table(obj$path_df, paste0("panel_o_coordinates_", key))
colors <- c("Matched CTRL" = COLORS[["CTRL"]], "miR-29 inhibitor" = COLORS[["INH"]], "miR-29 mimic" = COLORS[["MIM"]])
# The accepted final panel uses the original unrotated coordinates.
orient <- function(d) d
d <- orient(obj$path_df); curves <- orient(obj$curve_df)
max_r <- max(sqrt(d$x_plot^2 + d$y_plot^2)) * 1.12
radii <- pretty(c(0, max_r), n = 4); radii <- radii[radii > 0 & radii <= max_r]
rings <- bind_rows(lapply(radii, function(r) data.frame(r, x = r * cos(seq(0, 2*pi, length.out = 361)),
                                                       y = r * sin(seq(0, 2*pi, length.out = 361)))))
p <- ggplot() + geom_path(data = rings, aes(x, y, group = r), color = "grey75", linewidth = .25, linetype = 3) +
  geom_path(data = d, aes(x_plot, y_plot, color = arm, group = condition), linewidth = .3, alpha = .5) +
  geom_path(data = curves, aes(x_plot, y_plot, color = arm, group = condition), linewidth = .7) +
  geom_point(data = d, aes(x_plot, y_plot, color = arm), size = 2) +
  geom_text(data = subset(d, time %in% c(1, 6)), aes(x_plot, y_plot, label = time_label, color = arm), nudge_y = max_r*.065, size = 3) +
  facet_wrap(~experiment, nrow = 1) + coord_equal(xlim = c(-max_r,max_r), ylim = c(-max_r,max_r)) +
  scale_x_continuous(breaks = c(-1000, 0, 1000)) +
  scale_color_manual(values = colors) + labs(x = expression("Dimension 1 (RMS "*10^3*")"), y = "Dimension 2") +
  theme_figure() + theme(panel.spacing = grid::unit(8, "mm"), legend.text = element_text(size = 8))
save_panel(p, "4n_methylation_walk", 6.0, 3.5)
