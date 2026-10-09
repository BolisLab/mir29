# Final Figure S5: retained panels only. Illustrator-only decorations are not redrawn.
stopifnot(requireNamespace("ggrastr", quietly = TRUE))
tissue_labels <- c(breast = "Breast", central_nervous_system = "CNS", haematopoietic_and_lymphoid_tissue = "Blood / lymphoid",
  large_intestine = "Intestine", liver = "Liver", lung = "Lung", oesophagus = "Esophagus", ovary = "Ovary",
  pancreas = "Pancreas", skin = "Skin", stomach = "Stomach")
supp_theme <- function(size = 10) panel_theme(size) + theme(plot.title = element_text(size = size, face = "plain"))

# ---- S5a
d <- rbind(data.frame(Cancer = tcga_counts$Cancer, Type = "Primary", N = tcga_counts$Primary), data.frame(Cancer = tcga_counts$Cancer,
  Type = "Metastatic", N = tcga_counts$Metastatic))
d$Cancer <- factor(d$Cancer, levels = rev(a_cancers))
d$Type <- factor(d$Type, levels = c("Primary", "Metastatic"))
ann <- tcga_counts
ann$Cancer <- factor(ann$Cancer, levels = rev(a_cancers))
ann$Label <- ifelse(ann$Metastatic > 0, paste0(ann$Primary, " / ", ann$Metastatic), ann$Primary)
p <- ggplot(d, aes(N, Cancer, fill = Type)) + geom_col(width = 0.7) + geom_text(data = ann, inherit.aes = FALSE,
  mapping = aes(x = N, y = Cancer, label = Label), hjust = -0.12, size = 2.6, colour = "#477CB4") +
  scale_fill_manual(values = c(Primary = "#4C7FB8", Metastatic = "#DE983C")) + scale_x_continuous(expand = expansion(mult = c(0,
  0.28))) + labs(title = "TCGA: miRNA and DNAmAge", x = "Patients", y = NULL) + supp_theme(9)
save_panel(p, "S5a_TCGA_cohort", 5, 7)

# ---- S5b
tcga_rank_plot <- function(d, ylabel) {
  h <- d[!is.na(d$Feature), ]
  background <- d[is.na(d$Feature), ]
  h$Label <- c(`miR-29a` = "29a", `miR-29b` = "29b", `miR-29c` = "29c", `pan-miR-29` = "Σ(miR-29)")[h$Feature]
  zero <- max(d$Rank[d$rho <= 0])
  ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#EAEDF7") + annotate("rect",
    xmin = zero, xmax = Inf, ymin = 0, ymax = Inf, fill = "#FAF0F0") + geom_hline(yintercept = 0,
    linewidth = 0.4) + geom_vline(xintercept = zero, linetype = "dashed", linewidth = 0.3) + geom_point(data = background,
    aes(Rank, rho, colour = rho > 0, size = 0.35 + 4 * abs(rho)), alpha = 0.88) + scale_colour_manual(values = c(`FALSE` = "#446F9E",
    `TRUE` = "#D65242"), guide = "none") + scale_size_identity() + geom_point(data = h, aes(Rank,
    rho, fill = Feature, shape = Feature), size = 3, stroke = 0.3) + scale_fill_manual(values = COLORS,
    guide = "none") + scale_shape_manual(values = c(`miR-29a` = 21, `miR-29b` = 21, `miR-29c` = 21,
    `pan-miR-29` = 24), guide = "none") + ggrepel::geom_text_repel(data = h, aes(Rank, rho, label = Label),
    nudge_x = -320, direction = "y", hjust = 1, size = 3.2, seed = 29, min.segment.length = 0, box.padding = 0.4,
    max.overlaps = Inf) + annotate("text", x = 2040, y = -0.214, label = "N = 8,428", hjust = 1,
    size = 3.3, fontface = "italic") + scale_x_continuous(limits = c(-300, 2100), breaks = c(0, 500,
    1000, 1500), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.225, 0.265), breaks = seq(-0.2,
    0.2, 0.1), expand = c(0, 0)) + labs(title = "Cross-cancer meta-rank", x = "miRNA meta-rank",
    y = ylabel) + supp_theme(10) + theme(aspect.ratio = 2.53, panel.border = element_rect(fill = NA,
    colour = "grey25", linewidth = 0.35), axis.line = element_blank())
}
p <- cowplot::plot_grid(tcga_rank_plot(tcga_ranks$AGE, "miRNA expression ~ chronological age [ρ]"),
  tcga_rank_plot(tcga_ranks$Acceleration, "miRNA expression ~ age acceleration (P) [ρ]"), nrow = 1,
  align = "hv", axis = "tblr")
save_panel(p, "S5b_TCGA_age_and_acceleration", 7.2, 8.5)

# ---- S5c/e/i: horizontal segments show movement from the full estimate, not uncertainty
loo_plot_final <- function(d, title, labels = NULL, colour = "#396896", xlabel = "Equal-stratum Fisher mean of Spearman ρ",
  full_label = "All cancers") {
  if (is.null(labels))
    labels <- setNames(d$Excluded[-1], d$Excluded[-1])
  yl <- c(None = full_label, setNames(paste("Omit", labels[d$Excluded[-1]]), d$Excluded[-1]))
  d$Row <- factor(d$Excluded, levels = rev(d$Excluded))
  ggplot(d, aes(rho, Row)) + geom_vline(xintercept = d$Full_rho[1], linetype = "dashed", colour = "grey50",
    linewidth = 0.4) + geom_segment(aes(x = Full_rho, xend = rho, yend = Row), colour = colour, linewidth = 0.45,
    alpha = 0.7) + geom_point(aes(shape = Excluded == "None"), colour = colour, fill = "white", size = 2.5) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 23), guide = "none") + scale_y_discrete(labels = yl) +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) + labs(title = title, x = xlabel,
    y = NULL) + supp_theme(9)
}
p <- loo_plot_final(loo_tcga, "Leave-one-cancer-out: Σ(miR-29) ~ DNAmAge (P)")
save_panel(p, "S5c_TCGA_leave_one_out", 6, 7.5)

# ---- S5d
d <- ccle_counts
d$Tissue <- factor(d$Tissue, levels = b_tissues)
p <- ggplot(d, aes(Tissue, N)) + geom_col(fill = "#58A397", width = 0.7) + geom_text(aes(label = N),
  vjust = -0.4, size = 3, colour = "#438F88") + scale_x_discrete(labels = tissue_labels) + scale_y_continuous(expand = expansion(mult = c(0,
  0.12))) + labs(x = NULL, y = "Eligible cell lines") + supp_theme(10) + theme(axis.text.x = element_text(angle = 90,
  hjust = 1, vjust = 0.5))
save_panel(p, "S5d_CCLE_cohort", 6, 3.4)

# ---- S5e
p <- loo_plot_final(loo_ccle, "Leave-one-lineage-out: Σ(miR-29) ~ DNAmAge (S)", tissue_labels, "#45998C",
  full_label = "All lineages")
save_panel(p, "S5e_CCLE_leave_one_out", 6, 4.5)

# ---- S5f
ccle_rank_plot <- function(d, title) {
  h <- d[!is.na(d$Feature), ]
  h$Label <- c(`miR-29a` = "29a", `miR-29b` = "29b", `miR-29c` = "29c", `pan-miR-29` = "Σ(miR-29)")[h$Feature]
  ggplot(d, aes(Rank, rho)) + geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) + ggrastr::rasterise(geom_point(colour = "#B6C5CA",
    size = 0.6), dpi = 600, dev = "cairo") + geom_point(data = h, aes(fill = Feature, shape = Feature),
    size = 2.6, colour = "grey35", stroke = 0.25) + scale_fill_manual(values = COLORS, guide = "none") +
    scale_shape_manual(values = c(`miR-29a` = 21, `miR-29b` = 21, `miR-29c` = 21, `pan-miR-29` = 23),
      guide = "none") + ggrepel::geom_text_repel(data = h, aes(label = Label), seed = 29, nudge_x = -215,
    direction = "y", size = 2.7, min.segment.length = 0, max.overlaps = Inf) + scale_x_continuous(breaks = c(0,
    250, 500, 750)) + labs(title = title, x = "miRNA meta-rank", y = "Spearman ρ") + supp_theme(10)
}
titles <- c(Age = "miRNA ~ chronological age", Skin = "miRNA ~ DNAmAge (S)", Horvath = "miRNA ~ DNAmAge (P)",
  HorvathAcceleration = "miRNA ~ age acceleration (P)")
p <- cowplot::plot_grid(plotlist = lapply(names(titles), function(m) ccle_rank_plot(ccle_ranks[[m]],
  titles[[m]])), ncol = 2, align = "hv", axis = "tblr")
save_panel(p, "S5f_CCLE_four_ranks", 8, 7.4)

# ---- S5g
s <- rna_boxes$summary
o <- rna_boxes$outliers
h <- rna_boxes$highlights
for (obj in c("s", "o", "h")) {
  z <- get(obj)
  z$Row <- 31 - match(z$Unit, a_cancers)
  assign(obj, z)
}
h$Row <- h$Row + seq(-0.16, 0.16, length.out = 7)[match(h$Gene, c(MRC_GENES, "MIR29B2CHG", "LINC-PINT"))]
p <- ggplot(s, aes(x = Row, group = Unit)) + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
  linewidth = 0.35) + geom_boxplot(aes(ymin = Low, lower = Q25, middle = Median, upper = Q75, ymax = High),
  stat = "identity", orientation = "x", width = 0.64, fill = "#626868", colour = "#444848", linewidth = 0.25) +
  ggrastr::rasterise(geom_point(data = o, aes(y = rho), size = 0.16, alpha = 0.45, colour = "#6C7479"),
    dpi = 600, dev = "cairo") + geom_point(data = h, aes(x = Row, y = rho, fill = Feature), shape = 21,
  size = 1.9, colour = "grey35", stroke = 0.22) + scale_fill_manual(values = c(`MRC genes` = "#E8D833",
  MIR29B2CHG = "#E2A1A1", `LINC-PINT` = "#AE8EB7")) + scale_x_continuous(breaks = 1:30, labels = rev(a_cancers),
  limits = c(0.4, 30.6), expand = c(0, 0)) + scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1,
  1), expand = expansion(mult = 0.015)) + coord_flip() + labs(x = NULL, y = "Gene expression ~ DNAmAge (P) [ρ]",
  fill = NULL) + supp_theme(9) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
  legend.position = "top", legend.text = element_text(size = 8))
save_panel(p, "S5g_TCGA_RNA_distributions", 6.7, 7.7)

# ---- S5h
d <- pair_counts
d$Cancer <- factor(d$Cancer, levels = rev(CANCERS16))
p <- ggplot(d, aes(N_pairs, Cancer)) + geom_col(fill = "#9E6CAF", width = 0.7) + geom_text(aes(label = N_pairs),
  hjust = -0.15, size = 3, colour = "#86589A") + scale_x_continuous(expand = expansion(mult = c(0,
  0.2))) + labs(title = "Matched tumor-normal pairs", x = "Paired patients", y = NULL) + supp_theme(10)
save_panel(p, "S5h_paired_cohort", 4.5, 5.2)

# ---- S5i
p <- loo_plot_final(loo_paired, "Leave-one-cancer-out: ΔΣ(miR-29) ~ ΔMRC", colour = "#8139A4", xlabel = "Pearson r (paired cancer means)")
save_panel(p, "S5i_paired_leave_one_out", 5.8, 5.8)

# ---- S5j
d <- reference_qc
d$Cancer <- factor(d$Cancer, levels = rev(CANCERS15))
d$Nlabel <- paste0("N = ", d$N_normal_reference, ifelse(d$Cancer == "KICH", "*", ""))
p0 <- ggplot(d, aes(N_stable_CpGs/1000, Cancer)) + geom_col(fill = "#5A91A3", width = 0.7) + geom_text(aes(label = sprintf("%.1fk",
  N_stable_CpGs/1000)), hjust = -0.15, size = 2.8) + geom_text(aes(x = 4, label = Nlabel), hjust = 0,
  colour = "white", size = 2.8) + scale_x_continuous(limits = c(0, 365), breaks = c(0, 100, 200, 300),
  expand = c(0, 0)) + labs(x = "Retained CpGs (thousands)", y = NULL, title = "Normal-tissue references") +
  supp_theme(10)
ref_summary_plot <- function(metric, cutoff, xlabel) {
  z <- reference_distributions[reference_distributions$Metric == metric, ]
  z$Cancer <- factor(z$Cancer, levels = rev(CANCERS15))
  ggplot(z, aes(Median, Cancer)) + geom_vline(xintercept = cutoff, colour = "#C84242", linetype = "dashed",
    linewidth = 0.4) + geom_segment(aes(x = P05, xend = P95, yend = Cancer), colour = "#AAC3CC",
    linewidth = 0.4) + geom_segment(aes(x = Q25, xend = Q75, yend = Cancer), colour = "#5AA38C",
    linewidth = 1.3) + geom_point(size = 1.3, colour = "#526C6B") + scale_x_continuous(limits = c(0,
    cutoff * 1.08)) + labs(x = xlabel, y = NULL, title = paste("Cutoff:", cutoff)) + supp_theme(10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
}
p <- cowplot::plot_grid(p0, ref_summary_plot("normal_sd_beta", 0.05, "SD(β)"), ref_summary_plot("age_delta_beta_decade",
  0.01, "|Δβ| / decade"), nrow = 1, rel_widths = c(1.85, 1, 1), align = "h", axis = "tb")
save_panel(p, "S5j_reference_QC", 10, 6)

# ---- S5k
d <- subset(fidelity_age, Contrast == "Age acceleration ~ fidelity")
m <- subset(fidelity_age_meta, Contrast == "Age acceleration ~ fidelity")
d$Row <- factor(d$Cancer, levels = rev(c("Pooled", CANCERS15)))
m$Row <- factor("Pooled", levels = levels(d$Row))
p <- ggplot(d, aes(rho, Row)) + annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#E7EDF5") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FAE9E9") + geom_vline(xintercept = 0,
  linetype = "dashed", colour = "grey55", linewidth = 0.35) + geom_segment(aes(x = low, xend = high,
  yend = Row), linewidth = 0.45, colour = "grey50") + geom_point(aes(colour = rho > 0, alpha = q <
  0.05), size = 2) + geom_segment(data = m, aes(x = low, xend = high, yend = Row), linewidth = 0.7) +
  geom_point(data = m, shape = 18, size = 3) + scale_colour_manual(values = c(`FALSE` = "#426FAD",
  `TRUE` = "#D95C64"), guide = "none") + scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1),
  guide = "none") + scale_y_discrete(labels = c(Pooled = "Pooled", setNames(CANCERS15, CANCERS15))) +
  labs(title = sprintf("Age acceleration (P) ~ fidelity loss\nPooled q = %.3f", m$q), x = "Age-adjusted Spearman ρ",
    y = NULL) + supp_theme(10)
save_panel(p, "S5k_acceleration_fidelity", 5.3, 6.1)
