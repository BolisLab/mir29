# S4j, m and n: independent reconstruction and numerical checks.
hp_input <- function(name) input(paste0("htert_hpf/", name, ".rds"))
hp_check <- function(name, observed, reference, tolerance = 1e-10) assert_close(observed, reference, name, tolerance)
retention <- function(delta) {
  stopifnot(all(is.finite(delta)), "T6" %in% colnames(delta))
  ref <- delta[, "T6"]
  as.numeric(crossprod(ref, delta) / sum(ref^2))
}
scatter_stats <- function(d) data.frame(n = nrow(d),
  retention = sum(d$effect_T6 * d$effect_T8) / sum(d$effect_T6^2),
  pearson_r = cor(d$effect_T6, d$effect_T8),
  spearman_r = cor(d$effect_T6, d$effect_T8, method = "spearman"))

temporal <- hp_input("htert_temporal_tables")
historical_camera <- hp_input("htert_historical_enrichment")
counts <- do.call(rbind, lapply(names(temporal), function(condition) {
  tab <- temporal[[condition]]
  data.frame(condition, hyper = sum(tab$q < .05 & tab$slope > 0),
             hypo = sum(tab$q < .05 & tab$slope < 0), n_tested = nrow(tab))
}))
counts$total <- counts$hyper + counts$hypo
counts$hyper_percent <- ifelse(counts$total > 0, 100 * counts$hyper / counts$total, NA)
hp_check("hTERT inhibitor hyper/hypo", unlist(counts[counts$condition == "Inhibitor", c("hyper", "hypo")]), c(3207, 425))
hp_check("hTERT CTRL and Mimic significant DMP counts", counts$total[counts$condition != "Inhibitor"], c(0, 0))
save_table(counts, "panel_p_DMP_counts")

common <- Reduce(intersect, lapply(temporal, `[[`, "probe_id"))
tables_common <- lapply(temporal, function(x) x[x$probe_id %in% common, ])
inh <- tables_common$Inhibitor
sets <- list(Hypo = inh$probe_id[inh$q < .05 & inh$slope < 0],
             Hyper = inh$probe_id[inh$q < .05 & inh$slope > 0])
stopifnot(setequal(sets$Hyper, historical_camera$signatures$Hyper),
          setequal(sets$Hypo, historical_camera$signatures$Hypo))
camera <- do.call(rbind, lapply(c("CTRL", "Mimic"), function(target) {
  tab <- tables_common[[target]]
  stats <- setNames(tab$statistic, tab$probe_id)
  result <- limma::cameraPR(stats, index = sets, use.ranks = FALSE,
                           inter.gene.cor = .01, sort = FALSE)
  data.frame(target, signature = rownames(result), n = result$NGenes,
             direction = result$Direction, p = result$PValue, q = p.adjust(result$PValue, "BH"))
}))
old <- historical_camera$results$projection_results
idx <- match(paste(camera$target, camera$signature), paste(old$target, old$signature))
hp_check("CAMERA p matches historical", camera$p, old$camera_p_value[idx])
hp_check("CAMERA within-target BH matches historical", camera$q, old$FDR[idx])
camera$minus_log10_q <- -log10(camera$q)
save_table(camera, "panel_p_CAMERA_statistics")
rank <- temporal$Inhibitor[order(temporal$Inhibitor$statistic, decreasing = TRUE), ]
rank$rank <- seq_len(nrow(rank))
rank$direction <- ifelse(rank$q < .05, ifelse(rank$slope > 0, "Hyper", "Hypo"), "Other")
hp_check("hTERT rank universe", nrow(rank), 937687)
save_table(rank[rank$direction != "Other", ], "panel_p_significant_CpGs")
htert_panel <- list(counts = counts, rank = rank, camera = camera, signatures = sets,
             settings = list(q_cutoff = .05, camera_inter_probe_correlation = .01,
                             camera_use_ranks = FALSE, BH_family = "Hyper/Hypo within each target"))

rna <- hp_input("hpf_rna_inputs")
sel <- rna$selection
ids <- sel$gene_id[is.finite(sel$pvalue) & sel$pvalue < .05]
stopifnot(setequal(ids, rna$exact$core_ids))
rna_curves <- do.call(rbind, lapply(names(rna$exact$universe_ids), function(universe) {
  ids <- rna$exact$universe_ids[[universe]]
  data.frame(universe, timepoint = 1:8, n_features = length(ids),
             retention = retention(rna$delta[ids, , drop = FALSE]))
}))
idx <- match(paste(rna_curves$universe, rna_curves$timepoint),
             paste(rna$exact$nested_curves$universe, rna$exact$nested_curves$timepoint))
hp_check("RNA all seven nested-universe curves", rna_curves$retention, rna$exact$nested_curves$retention[idx])
rna_scatter <- data.frame(feature_id = rna$exact$core_ids,
                          effect_T6 = rna$delta[rna$exact$core_ids, "T6"],
                          effect_T8 = rna$delta[rna$exact$core_ids, "T8"])
idx <- match(rna_scatter$feature_id, rna$exact$scatter_data$feature_id)
hp_check("RNA observed T6 coordinates", rna_scatter$effect_T6, rna$exact$scatter_data$effect_T6[idx])
hp_check("RNA observed T8 coordinates", rna_scatter$effect_T8, rna$exact$scatter_data$effect_T8[idx])

dnam <- hp_input("hpf_dnam_inputs")
dnam_sel <- hp_input("hpf_dnam_selection_audit")
ids <- dnam_sel$probe[is.finite(dnam_sel$T6_offset_FDR) & dnam_sel$T6_offset_FDR < .05]
stopifnot(setequal(ids, dnam$exact$exact_core_ids))
dnam_curves <- dnam$exact$nested_curves[, c("universe", "timepoint", "n_features", "retention")]
dec <- dnam$decision$curves$nested_responsive
dec <- dec[dec$method == "RAW_T6_UNWEIGHTED", ]
idx <- match(paste(dnam_curves$universe, dnam_curves$timepoint), paste(dec$universe, dec$timepoint))
hp_check("DNAm curves match the final collapsed-CpG checkpoint", dnam_curves$retention, dec$retention[idx])
dnam_scatter <- dnam$exact$scatter_data
# These are two historical exports, not bit-identical numerical objects.
# Preserve both: their T8 retention differs by 4.505e-6 (0.00045 percentage points).
hp_check("DNAm T8 curve/scatter checkpoint agreement within 1e-5", scatter_stats(dnam_scatter)$retention,
              dnam$exact$exact_core_curve$retention[dnam$exact$exact_core_curve$timepoint == 8], tolerance = 1e-5)

make_envelope <- function(curves) do.call(rbind, lapply(1:8, function(timepoint) {
  v <- curves$retention[curves$timepoint == timepoint]
  data.frame(timepoint, minimum = min(v), maximum = max(v), median = median(v), n_universes = length(v))
}))
rna_envelope <- make_envelope(rna_curves)
dnam_envelope <- make_envelope(dnam_curves)
hp_check("RNA envelope minima/maxima", c(rna_envelope$minimum, rna_envelope$maximum),
              c(rna$exact$envelope_curve$minimum_retention, rna$exact$envelope_curve$maximum_retention))
hp_check("DNAm envelope minima/maxima", c(dnam_envelope$minimum, dnam_envelope$maximum),
              c(dnam$exact$envelope_curve$minimum_retention, dnam$exact$envelope_curve$maximum_retention))
core_curves <- rbind(transform(rna_curves[rna_curves$universe == "EXACT_CORE", ], omic = "RNA"),
                     transform(dnam_curves[dnam_curves$universe == "PANEL_P_CORE_EXACT", ], omic = "DNAm"))
stats <- rbind(cbind(omic = "RNA", scatter_stats(rna_scatter)), cbind(omic = "DNAm", scatter_stats(dnam_scatter)))
stats$retention_percent <- 100 * stats$retention
save_table(core_curves, "panel_r_core_curves")
save_table(rbind(transform(rna_envelope, omic = "RNA"), transform(dnam_envelope, omic = "DNAm")), "panel_r_envelopes")
save_table(rbind(transform(rna_curves, omic = "RNA"), transform(dnam_curves, omic = "DNAm")), "panel_r_all_universe_curves")
save_table(stats, "panel_s_statistics")
save_table(rna_scatter, "panel_s_RNA_observed")
save_table(dnam_scatter, "panel_s_DNAm_observed")
hpf_panel <- list(core = core_curves, RNA_envelope = rna_envelope, DNAm_envelope = dnam_envelope,
             RNA_scatter = rna_scatter, DNAm_scatter = dnam_scatter, statistics = stats,
             RNA_universes = rna$exact$universe_ids, DNAm_universes = dnam$decision$universe_ids)

colors <- c(Hyper = "#D9483B", Hypo = "#3F6FB6", RNA = "#248636", DNAm = "#3F7099")
theme_panel <- function() theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = .4), axis.ticks = element_line(linewidth = .35),
        axis.text = element_text(color = "black"), plot.title = element_text(size = 11, hjust = .5),
        plot.margin = margin(7, 10, 7, 7), legend.title = element_blank())
p <- htert_panel
c <- p$counts
c$y <- match(c$condition, c("Mimic", "Inhibitor", "CTRL"))
bars <- rbind(data.frame(y = c$y, xmin = 0, xmax = c$hyper, direction = "Hyper"),
              data.frame(y = c$y, xmin = c$hyper, xmax = c$total, direction = "Hypo"))
bars <- bars[bars$xmax > bars$xmin, ]
counts <- ggplot() + geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = y - .27,
                                               ymax = y + .27, fill = direction), color = NA) +
  scale_fill_manual(values = colors[c("Hyper", "Hypo")], guide = "none") +
  geom_text(data = c[c$total == 0, ], aes(x = 100, y = y, label = "0"), hjust = 0, size = 3.6) +
  annotate("text", x = 3207 / 2, y = 2, label = "88%", color = "white", size = 4.1) +
  annotate("text", x = 3207 + 425 / 2, y = 2.62, label = "12%", size = 3.8) +
  annotate("segment", x = 3207 + 425 / 2, xend = 3207 + 425 / 2, y = 2.4, yend = 2.09, linewidth = .35) +
  annotate("text", x = 3750, y = 2, label = "3,632", hjust = 0, size = 3.7) +
  scale_y_continuous(breaks = 1:3, labels = c("Mimic", "Inhibitor", "CTRL"), limits = c(.5, 3.5)) +
  scale_x_continuous(breaks = seq(0, 4000, 1000), labels = c("0", "1k", "2k", "3k", "4k"),
                     limits = c(0, 4650), expand = expansion(mult = 0)) +
  labs(x = "Temporal DMPs (FDR < 0.05)", y = NULL, title = "hTERT HGF") + theme_panel() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

d <- p$rank
background <- d[d$direction == "Other", ]
# Uniform rank thinning affects only the dense grey visual background.
# The complete ranked input and every significant CpG are retained in data/.
background <- background[unique(round(seq(1, nrow(background), length.out = 5000))), ]
selected <- d[d$direction != "Other", ]
thresholds <- c(min(selected$statistic[selected$direction == "Hyper"]),
                max(selected$statistic[selected$direction == "Hypo"]))
ranking <- ggplot() +
  geom_hline(yintercept = 0, color = "grey65", linewidth = .35) +
  geom_point(data = background, aes(rank, statistic), shape = 16, stroke = 0, size = .45, color = "grey70") +
  geom_point(data = selected, aes(rank, statistic, color = direction), shape = 16, stroke = 0, size = .85) +
  geom_hline(yintercept = thresholds[1], color = colors["Hyper"], linetype = "dashed", linewidth = .35) +
  geom_hline(yintercept = thresholds[2], color = colors["Hypo"], linetype = "dashed", linewidth = .35) +
  scale_color_manual(values = colors[c("Hyper", "Hypo")], guide = "none") +
  annotate("text", x = nrow(d) * .10, y = max(d$statistic) * .85, label = "Hyper", hjust = 0, size = 3.5) +
  annotate("text", x = nrow(d) * .9, y = min(d$statistic) * .85, label = "Hypo", hjust = 1, size = 3.5) +
  scale_x_continuous(breaks = c(1, 500000, nrow(d)), labels = c("1", "500k", "938k"),
                     expand = expansion(mult = c(.03, .03))) +
  labs(x = "CpG rank", y = "Inhibitor ~ time (moderated t)", title = "Temporal CpG signatures") + theme_panel()

e <- p$camera
e$y <- c(CTRL_Hyper = 4, CTRL_Hypo = 3, Mimic_Hyper = 1.6, Mimic_Hypo = .6)[paste(e$target, e$signature, sep = "_")]
e$color <- ifelse(e$q < .05, ifelse(e$direction == "Up", colors["Hyper"], colors["Hypo"]), "#D0D0D0")
e$label <- ifelse(e$q < .05, sprintf("q = %.3f", e$q), "n.s.")
enrichment <- ggplot(e, aes(x = minus_log10_q, y = y)) +
  geom_col(aes(fill = color), orientation = "y", width = .62, color = NA) +
  scale_fill_identity() +
  geom_vline(xintercept = -log10(.05), linetype = "dashed", linewidth = .4, color = "grey55") +
  geom_hline(yintercept = 2.3, linewidth = .4) +
  geom_text(aes(x = minus_log10_q + .12, label = label), hjust = 0, size = 3.6) +
  scale_y_continuous(breaks = c(.6, 1.6, 3, 4), labels = c("hypo", "hyper", "hypo", "hyper"),
                     limits = c(0, 4.9), expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(0, 3.3), expand = expansion(mult = 0)) +
  annotate("text", x = .02, y = 4.65, label = "CTRL", hjust = 0, size = 3.3, color = "grey35") +
  annotate("text", x = .02, y = 2.02, label = "Mimic", hjust = 0, size = 3.3, color = "#923D8E") +
  labs(x = expression(-log[10]("BH FDR")), y = NULL, title = "Probe-set enrichment (CAMERA)") + theme_panel() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
panel_p <- plot_grid(counts, ranking, enrichment, nrow = 1, rel_widths = c(1.0, 1.15, 1.3),
                     labels = c("", "", ""), label_size = 15, align = "h")
save_panel(counts, "S4j_hTERT_DMP_counts", 3.5, 2.15)
save_panel(ranking, "S4j_hTERT_CpG_ranking", 3.0, 4.0)
save_panel(enrichment, "S4j_hTERT_CAMERA", 4.4, 3.3)
save_panel(panel_p, "S4j_hTERT", 11.0, 3.75)

hp <- hpf_panel
retention_plot <- function(omic) {
  data <- hp$core[hp$core$omic == omic, ]
  data <- data[order(data$timepoint), ]
  env <- hp[[paste0(omic, "_envelope")]]
  clr <- colors[omic]
  labs_end <- data[data$timepoint > 6, ]
  labs_end$label <- sprintf("%d%%", round(100 * labs_end$retention))
  ggplot(data, aes(timepoint, retention)) +
    annotate("rect", xmin = .55, xmax = 6, ymin = -Inf, ymax = Inf, fill = "#E4F3E1") +
    annotate("rect", xmin = 6, xmax = 8.55, ymin = -Inf, ymax = Inf, fill = "#EFEFEF") +
    geom_ribbon(data = env, aes(x = timepoint, ymin = minimum, ymax = maximum),
                fill = clr, alpha = .19, inherit.aes = FALSE) +
    geom_hline(yintercept = c(0, 1), color = "grey60", linetype = "dotted", linewidth = .4) +
    geom_vline(xintercept = 6, linetype = "dashed", color = "grey40", linewidth = .4) +
    geom_line(data = data[data$timepoint <= 6, ], color = clr, linewidth = .7) +
    geom_line(data = data[data$timepoint >= 6, ], color = "grey30", linetype = "dashed", linewidth = .6) +
    geom_point(data = data[data$timepoint <= 6, ], fill = clr, color = "white", shape = 21, size = 2.4, stroke = .35) +
    geom_point(data = data[data$timepoint == 6, ], fill = clr, color = "grey25", shape = 21, size = 2.5, stroke = .4) +
    geom_point(data = data[data$timepoint > 6, ], fill = "white", color = "grey30", shape = 21, size = 2.9, stroke = .6) +
    geom_text(data = labs_end, aes(label = label), nudge_y = if (omic == "RNA") .09 else .13, size = 3.5) +
    scale_x_continuous(breaks = 1:8, labels = paste0("T", 1:8), limits = c(.55, 8.55), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(-.25, 1.18), expand = c(0, 0)) +
    labs(x = "treated                          withdrawal", y = paste("Retention of the", omic, "state"),
          title = paste("T6 responsive core (", omic, ")", sep = "")) +
    theme_panel() + theme(axis.title.x = element_text(size = 9), panel.border = element_rect(fill = NA, color = "grey55", linewidth = .35))
}
panel_r <- plot_grid(retention_plot("RNA"), retention_plot("DNAm"), nrow = 1,
                     labels = c("", ""), label_size = 15, align = "hv", axis = "tblr")
save_panel(panel_r, "S4m_HPF_retention", 7.0, 3.05)

scatter_plot <- function(omic) {
  data <- hp[[paste0(omic, "_scatter")]]
  stats <- hp$statistics[hp$statistics$omic == omic, ]
  clr <- colors[omic]
  lim <- 1.08 * max(abs(c(data$effect_T6, data$effect_T8)))
  ggplot(data, aes(effect_T6, effect_T8)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", linewidth = .35) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = .35) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey35", linewidth = .5) +
    geom_point(shape = 16, stroke = 0, size = if (omic == "RNA") .75 else .55,
               color = clr, alpha = if (omic == "RNA") .68 else .50) +
    geom_abline(intercept = 0, slope = stats$retention, color = clr, linewidth = .7) +
    annotate("text", x = -.91 * lim, y = .81 * lim,
             label = sprintf("r = %.2f", stats$pearson_r), color = clr, hjust = 0, fontface = "italic", size = 3.8) +
    annotate("text", x = .9 * lim, y = -.83 * lim,
             label = sprintf("res = %d%%", round(stats$retention_percent)), color = clr, hjust = 1, size = 3.5) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(title = paste0("Residual effect (", omic, ")"), subtitle = paste0("n = ", format(stats$n, big.mark = ",")),
          x = if (omic == "RNA") expression("Observed effect at T6 ("*Delta*"VST)") else expression("Observed effect at T6 ("*Delta*"M)"),
          y = if (omic == "RNA") expression("Observed effect at T8 ("*Delta*"VST)") else expression("Observed effect at T8 ("*Delta*"M)")) +
    theme_panel() + theme(plot.title = element_text(hjust = 0), plot.subtitle = element_text(face = "italic", size = 10))
}
panel_s <- plot_grid(scatter_plot("RNA"), scatter_plot("DNAm"), nrow = 1,
                     labels = c("", ""), label_size = 15, align = "hv", axis = "tblr")
save_panel(panel_s, "S4n_HPF_residual_effect", 7.0, 3.35)
