# Figure 5 drawing layer. All values are read from the validated panel tables.  Graphical assembly
# is regenerated, not inferred from pixels in Illustrator.
if (!exists("PANELS")) PANELS <- readRDS(file.path(BUNDLE, "data", "Figure5_panel_data.rds"))
COLORS <- c(`miR-29a` = "#E8D833", `miR-29b` = "#91C7D3", `miR-29c` = "#E74337", `pan-miR-29` = "#111111")
panel_theme <- function(size = 10) theme_classic(base_size = size, base_family = "sans") + theme(plot.title = element_text(size = size +
  1, face = "bold", margin = margin(b = 8)), axis.line = element_line(linewidth = 0.35, colour = "#4A4A4A"),
  axis.ticks = element_line(linewidth = 0.3, colour = "#4A4A4A"), plot.margin = margin(8, 10, 8, 10),
  legend.position = "bottom", legend.title = element_blank())
save_panel <- function(plot, name, width, height) {
  ggsave(file.path(PDF_DIR, paste0(name, ".pdf")), plot, width = width, height = height, device = cairo_pdf,
    bg = "white", symbolfamily = grDevices::cairoSymbolFont("Arial", usePUA = FALSE))
  ggsave(file.path(PDF_DIR, paste0(name, ".png")), plot, width = width, height = height, dpi = 180,
    bg = "white")
}
combine_panels <- function(plots, widths, title) {
  # Align panel regions exactly even when one graph has a bottom legend and another has top-axis
  # labels. Shared cancer rows must have identical height.
  grobs <- lapply(plots, ggplotGrob)
  above <- lapply(grobs, function(g) sum(g$heights[seq_len(g$layout$t[g$layout$name == "panel"] - 1)]))
  below <- lapply(grobs, function(g) sum(g$heights[seq.int(g$layout$b[g$layout$name == "panel"] + 1,
    length(g$heights))]))
  max_above <- do.call(grid::unit.pmax, above)
  max_below <- do.call(grid::unit.pmax, below)
  grobs <- lapply(seq_along(grobs), function(i) {
    g <- gtable::gtable_add_rows(grobs[[i]], max_above - above[[i]], pos = 0)
    gtable::gtable_add_rows(g, max_below - below[[i]], pos = -1)
  })
  body <- cowplot::plot_grid(plotlist = grobs, nrow = 1, rel_widths = widths)
  heading <- cowplot::ggdraw() + cowplot::draw_label(title, x = 0.015, hjust = 0, fontface = "bold",
    size = 14)
  cowplot::plot_grid(heading, body, ncol = 1, rel_heights = c(0.075, 1))
}
rank_plot <- function(d, ylab, title, ylim = NULL) {
  highlighted <- d[!is.na(d$Feature), ]
  background <- d[is.na(d$Feature), ]
  zerorank <- max(c(1, d$Rank[d$rho <= 0]))
  p <- ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#FAE9E9") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#E7EDF5") + geom_hline(yintercept = 0,
    linewidth = 0.35) + geom_vline(xintercept = zerorank, linetype = "dashed", linewidth = 0.3, color = "grey60") +
    geom_point(data = background, aes(Rank, rho, colour = rho > 0), size = 0.5) + scale_colour_manual(values = c(`FALSE` = "#426FAD",
    `TRUE` = "#D95C64"), guide = "none") + geom_point(data = highlighted, aes(Rank, rho, fill = Feature,
    shape = Feature), size = 2.5, colour = "#555555", stroke = 0.3) + scale_fill_manual(values = COLORS,
    guide = "none") + scale_shape_manual(values = c(`miR-29a` = 21, `miR-29b` = 21, `miR-29c` = 21,
    `pan-miR-29` = 24), guide = "none") + ggrepel::geom_text_repel(data = highlighted, aes(Rank,
    rho, label = Feature), size = 2.7, direction = "y", nudge_x = -0.22 * max(d$Rank), hjust = 1,
    seed = 29, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.18, segment.size = 0.25) +
    labs(x = "miRNA meta-rank", y = ylab, title = title) + panel_theme(9)
  if (!is.null(ylim))
    p <- p + coord_cartesian(ylim = ylim)
  p
}
bubble_plot <- function(d, summary, row_names, metric_order, metric_labels, labels = row_names, title = NULL,
  summary_label = "Pan-cancer") {
  n <- length(row_names)
  d$row <- n + 1 - match(d$Row, row_names)
  d$x <- match(d$Metric, metric_order)
  summary$x <- match(summary$Metric, metric_order)
  summary$row <- -1.15
  display <- rbind(d[, c("row", "x", "rho")], summary[, c("row", "x", "rho")])
  ggplot(display, aes(x, row)) + annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -1.6, ymax = -0.65,
    fill = "#F0F0F0", colour = "grey40", linewidth = 0.3) + geom_hline(yintercept = 0.05, linewidth = 0.25,
    linetype = "dashed", colour = "grey65") + geom_vline(xintercept = 1:3, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(size = abs(rho), fill = rho), shape = 21, colour = "grey45", stroke = 0.25) +
    scale_size_area(max_size = 5, limits = c(0, max(abs(display$rho)) + 0.01), guide = "none") +
    scale_fill_gradient2(low = "#426FAD", mid = "white", high = "#D65359", midpoint = 0, limits = c(-0.55,
      0.55), oob = scales::squish, guide = "none") + geom_text(data = summary, aes(label = sprintf("%.2f",
    rho)), nudge_y = -0.9, size = 2.8) + scale_x_continuous(breaks = 1:3, labels = metric_labels,
    position = "top", limits = c(0.5, 3.5), expand = c(0, 0)) + scale_y_continuous(breaks = c(-1.15,
    seq_len(n)), labels = c(summary_label, rev(labels)), limits = c(-2.7, n + 0.7), expand = c(0,
    0)) + labs(x = NULL, y = NULL, title = title) + panel_theme(9) + theme(axis.line = element_blank(),
    axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0), axis.text.y = element_text(size = 8.5),
    plot.margin = margin(8, 4, 8, 5))
}

# a: identical cancer ordering across bubble and distribution plots.
a <- PANELS$a
ab <- a$bubbles
ab$Row <- ab$Cancer
pa_bubble <- bubble_plot(ab, a$summary, a$cancers, c("AGE", "Horvath", "Acceleration"), c("Age", "DNAmAge (P)",
  "Age acc. (P)"), title = "Σ(miR-29)")
ad <- a$distribution
ad$row <- length(a$cancers) + 1 - match(ad$Cancer, a$cancers)
pa_box <- ggplot(ad, aes(rho, row, group = Cancer)) + geom_vline(xintercept = 0, linetype = "dashed",
  colour = "grey60", linewidth = 0.3) + geom_boxplot(orientation = "y", width = 0.65, fill = "#626868",
  colour = "#444848", linewidth = 0.25, outlier.size = 0.18, outlier.alpha = 0.5) + geom_point(data = ad[!is.na(ad$Feature),
  ], aes(fill = Feature), shape = 21, size = 1.7, stroke = 0.25, colour = "#666666") + scale_fill_manual(values = COLORS,
  breaks = names(COLORS)[1:3]) + scale_y_continuous(limits = c(-2.7, length(a$cancers) + 0.7), expand = c(0,
  0)) + scale_x_continuous(breaks = c(-0.5, 0, 0.5)) + labs(x = "miRNA expression ~ DNAmAge [ρ]",
  y = NULL, title = "Within-cancer correlations") + panel_theme(9) + theme(axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.text = element_text(size = 8),
  legend.key.width = grid::unit(0.65, "cm"))
pa_rank <- rank_plot(a$ranking, "miRNA ~ DNAmAge [ρ]", "Cross-cancer meta-rank")
pa <- combine_panels(list(pa_bubble, pa_box, pa_rank), c(1.15, 1.6, 1.05), "a | TCGA miR-29 and epigenetic age")
save_panel(pa, "Figure5a_TCGA_miRNA_age", 9.3, 7.1)

# b: chronological age, Horvath age, Skin&Blood age; ranking uses Skin&Blood acceleration.
b <- PANELS$b
bb <- b$bubbles
bb$Row <- bb$Tissue
tissue_labels <- c("breast", "CNS", "blood", "intestine", "liver", "lung", "esophagus", "ovary", "pancreas",
  "skin", "stomach")
pb_bubble <- bubble_plot(bb, b$summary, b$tissues, c("Age", "Horvath", "Skin"), c("Age", "DNAmAge (P)",
  "DNAmAge (S)"), labels = tissue_labels, title = "Σ(miR-29)", summary_label = "Pan-lineage")
pb_rank <- rank_plot(b$ranking, "miRNA ~ age acceleration (S) [ρ]", "Skin & Blood acceleration")
pb <- combine_panels(list(pb_bubble, pb_rank), c(1.1, 1), "b | Cancer cell lines: 11 lineages")
save_panel(pb, "Figure5b_CCLE_miRNA_age", 6.8, 4.8)

# c: actual gene meta-ranks and the negative-tail selection.
c <- PANELS$c
cr <- c$ranking
show_genes <- c("MIR29B2CHG", "LINC-PINT", MRC_GENES)
ch <- cr[cr$SYMBOL %in% show_genes & !is.na(cr$SYMBOL), ]
ch <- ch[order(ch$SYMBOL, -abs(ch$rho)), ]
ch <- ch[!duplicated(ch$SYMBOL), ]
ch$Colour <- ifelse(ch$SYMBOL %in% MRC_GENES, "#D36E93", ifelse(ch$SYMBOL == "LINC-PINT", "#AE8EB7",
  "#E2A1A1"))
pc <- ggplot(cr, aes(MetaRank, rho)) + annotate("rect", xmin = 0, xmax = max(c$selected$MetaRank), ymin = -Inf,
  ymax = c$threshold, fill = "#F5E696", alpha = 0.7) + geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_line(linewidth = 0.4, colour = "#333333") + geom_point(data = ch, aes(fill = Colour), shape = 21,
  size = 2.6, stroke = 0.35, colour = "#555555") + scale_fill_identity() + ggrepel::geom_text_repel(data = ch,
  aes(label = SYMBOL), seed = 29, size = 3.2, fontface = "italic", nudge_x = ifelse(ch$rho < 0, 8500,
    -7500), direction = "y", min.segment.length = 0, segment.size = 0.25, box.padding = 0.35, max.overlaps = Inf) +
  annotate("text", x = 16000, y = -0.25, label = "185 genes\nbelow mean - 3 SD", size = 3.4, hjust = 0) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000, 40000), labels = c("0", "10k", "20k", "30k",
    "40k")) + labs(x = "Gene meta-rank", y = "Gene expression ~ DNAmAge [meta-ρ]", title = "c | Pan-cancer transcriptome ranking") +
  panel_theme(10)
save_panel(pc, "Figure5c_RNA_age_ranking", 5.1, 4.6)

# Generic compact forest. Both the uncertainty intervals and pooling are unchanged.
forest_plot <- function(d, meta, title, xlab, limits, show_labels = TRUE, expression_panel = FALSE) {
  rows <- if (expression_panel)
    CANCERS16 else CANCERS15
  d$row <- length(rows) + 1 - match(d$Cancer, rows)
  d$alpha <- if (expression_panel)
    ifelse(d$FDR < 0.05, 1, 0.35) else ifelse(d$CI_low * d$CI_high > 0, 1, 0.4)
  d$colour <- ifelse(d$Effect >= 0, "#BD3D4B", "#355EA6")
  meta$row <- c(-0.5, -1.6)
  diamonds <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    data.frame(x = c(meta$CI_low[i], meta$Effect[i], meta$CI_high[i], meta$Effect[i]), y = meta$row[i] +
      c(0, 0.25, 0, -0.25), group = i, colour = ifelse(meta$Effect[i] >= 0, "#BD3D4B", "#355EA6"))
  }))
  p <- ggplot() + annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F2F8FB") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FCF5F6") + geom_vline(xintercept = 0,
    linetype = "dashed", colour = "grey65", linewidth = 0.35) + geom_segment(data = d, aes(x = CI_low,
    xend = CI_high, y = row, yend = row), colour = "#777777", linewidth = 0.4) + geom_point(data = d,
    aes(Effect, row, colour = colour, alpha = alpha), size = 1.8) + geom_hline(yintercept = 0.35,
    linewidth = 0.35) + geom_polygon(data = diamonds, aes(x, y, group = group, fill = colour), linewidth = 0.25) +
    scale_colour_identity() + scale_fill_identity() + scale_alpha_identity() + scale_y_continuous(breaks = c(-1.6,
    -0.5, seq_along(rows)), labels = c("Organ", "Histotype", rev(rows)), limits = c(-2.3, length(rows) +
    0.7), expand = c(0, 0)) + coord_cartesian(xlim = limits) + labs(x = xlab, y = NULL, title = title) +
    panel_theme(10) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 11, face = "plain"))
  if (!show_labels)
    p <- p + theme(axis.text.y = element_blank())
  p
}
ed <- PANELS$e$effects
em <- PANELS$e$meta
ed$CI_low <- ed$Lower
ed$CI_high <- ed$Upper
pe1 <- forest_plot(subset(ed, Readout == "pan-miR-29"), subset(em, Readout == "pan-miR-29"), "Σ(miR-29)",
  "Paired difference (SD units)", c(-3.2, 2.1), TRUE, TRUE)
pe2 <- forest_plot(subset(ed, Readout == "ERC"), subset(em, Readout == "ERC"), "MRC", "Paired difference (SD units)",
  c(-2, 4.1), FALSE, TRUE)
pe <- combine_panels(list(pe1, pe2), c(1, 1), "e | Matched tumor-normal expression")
save_panel(pe, "Figure5e_paired_expression", 6.7, 5.4)

# f: original weighted KDE; Pearson coefficient and p-value are kept together.
f <- PANELS$f
pf <- ggplot() + geom_contour_filled(data = f$kde$grid, aes(x, y, z = density_scaled), breaks = c(0.05,
  0.15, 0.3, 0.5, 0.75, Inf), alpha = 0.85) + scale_fill_manual(values = c("#F4F8F4", "#E1EEDB", "#C9DFCF",
  "#ACCED0", "#86B8C2"), guide = "none") + geom_smooth(data = f$means, aes(Effect_miRNA, Effect_MRC),
  method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, colour = "#555555") + geom_point(data = f$means,
  aes(Effect_miRNA, Effect_MRC), colour = "#D84147", size = 1.9) + ggrepel::geom_text_repel(data = f$means,
  aes(Effect_miRNA, Effect_MRC, label = Cancer), seed = 29, size = 3, segment.size = 0.25, min.segment.length = 0,
  box.padding = 0.23, max.overlaps = Inf) + annotate("text", x = 1.8, y = 4.2, hjust = 1, label = sprintf("Pearson r = %.2f\np = %.4f",
  unname(f$pearson$estimate), f$pearson$p.value), size = 3.4) + coord_cartesian(xlim = c(-3.5, 2.1),
  ylim = c(-2.3, 4.6)) + labs(x = "Σ(miR-29) Δ(Tumor - Normal)", y = "MRC difference (Tumor - Normal)",
  title = "f | Matched tumor-normal pairs") + panel_theme(10)
save_panel(pf, "Figure5f_paired_scatter_Pearson", 5.3, 5)

# h: each forest retains its original (different) assay-availability cohort.
hl <- PANELS$h$left
hl$Effect <- hl$rho
hr <- PANELS$h$right
hr$Effect <- hr$rho
hm1 <- PANELS$h$left_meta
hm1$Effect <- hm1$rho
hm2 <- PANELS$h$right_meta
hm2$Effect <- hm2$rho
ph1 <- forest_plot(hl, hm1, "Σ(miR-29) ~ ΔM", "Age-adjusted Spearman's ρ", c(-0.55, 0.55))
ph2 <- forest_plot(hr, hm2, "MRC ~ ΔM", "Age-adjusted Spearman's ρ", c(-0.55, 0.85), FALSE)
ph <- combine_panels(list(ph1, ph2), c(1, 1), "h | Epigenetic fidelity loss")
save_panel(ph, "Figure5h_fidelity_associations", 6.7, 5.4)

saveRDS(list(a = pa, b = pb, c = pc, e = pe, f = pf, h = ph), file.path(BUNDLE, "data", "Figure5_ggplots.rds"))
