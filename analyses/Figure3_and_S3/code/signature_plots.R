# The y axis retains correlations or the historical normalized enrichment score.
curve_plot <- function(d, title, reverse_rank = FALSE, vertical = FALSE) {
  d <- d[is.finite(d$r), ]
  d <- d[order(if (reverse_rank) d$r else -d$r, d$gene), ]
  d$x <- if (vertical) 0 else seq_len(nrow(d)) / nrow(d)
  d$tail <- ifelse(d$UP, "UP", ifelse(d$DN, "DN", "Other"))
  ggplot(d, aes(x, r)) +
    geom_point(aes(colour = tail), size = .7, shape = 16) +
    scale_colour_manual(values = c(UP = "#ED1B24", DN = "#294DA3", Other = "black")) +
    scale_y_continuous(breaks = c(-1, 0, 1), expand = expansion(mult = .08)) +
    labs(x = if (vertical) NULL else "Gene rank (fraction)", y = "mRNA ~ age correlation", title = title) +
    theme_figure() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
hf_plot <- curve_plot(hf_curve, "HF: D / E")
hgf_plot <- curve_plot(hgf_curve, "HGF: F / G", reverse_rank = TRUE)
hgf_vertical_plot <- curve_plot(hgf_curve, "HGF: F / G", vertical = TRUE)
pan <- authentic$Pan
pan$tail <- ifelse(pan$SYMBOL %in% signature_sets$pan_aging_up, "UP",
                   ifelse(pan$SYMBOL %in% signature_sets$pan_aging_dn, "DN", "Other"))
pan_plot <- ggplot(pan, aes(r, 0, colour = tail)) + geom_point(size = .7) +
  scale_colour_manual(values = c(UP = "#ED1B24", DN = "#294DA3", Other = "black")) +
  labs(x = "Historical pan-tissue age correlation", y = NULL, title = "B: PanTissue-UP   C: PanTissue-DN") +
  theme_figure() + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Compress only the central NES interval, matching the approved figure.
compress_nes <- scales::trans_new("central_NES_compression",
  transform = function(y) ifelse(abs(y) <= .75, y / 10, sign(y) * (.075 + abs(y) - .75)),
  inverse = function(y) ifelse(abs(y) <= .075, y * 10, sign(y) * (.75 + abs(y) - .075)))
gsea_plot <- function(statistics, title = "miR-29 inhibition over time") {
  background <- authentic$background
  highlights <- authentic$highlights
  ids <- rownames(highlights)
  stopifnot(all(ids %in% statistics$signature))
  highlights$NES <- statistics$NES[match(ids, statistics$signature)]
  highlights$letter <- final_letter(ids)
  ggplot(background, aes(x, NES)) +
    geom_point(aes(colour = color, alpha = alpha), size = .28) +
    geom_point(data = highlights, aes(x, NES, fill = fill), colour = "#202020", shape = 21, size = 5.1, stroke = .5) +
    geom_text(data = highlights, aes(x, NES, label = letter), colour = "white", size = 3.1, fontface = "bold") +
    scale_colour_identity() + scale_fill_identity() + scale_alpha_identity() +
    scale_y_continuous(trans = compress_nes, breaks = c(-2.5, -2, -1.5, -1, 0, 1, 1.5, 2, 2.5)) +
    labs(x = NULL, y = "Normalized enrichment score (NES)", title = title) + theme_figure() +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
main_gsea_plot <- gsea_plot(main_statistics)
loto_plots <- lapply(definitions$tissues, function(tissue) {
  d <- loto_statistics[loto_statistics$variant == paste0("omit_", tissue), ]
  d$NES <- d$NES_simple
  label <- unname(support$tissue_labels[tissue])
  gsea_plot(d, paste("Omit", label)) + theme_figure(8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
})

# Sampling coverage and descriptive tissue contributions.
md <- support$source_metadata
donor_labels <- setNames(paste0(support$source_donors$donor, "\n", support$source_donors$SEX), support$source_donors$donor)
coverage_plot <- ggplot(md, aes(donor, tissue, fill = AGE)) +
  geom_tile(width = .9, height = .88, colour = "grey50", linewidth = .25) +
  geom_text(aes(label = AGE, colour = text_colour), size = 3) + scale_colour_identity() +
  scale_fill_gradient(low = "white", high = "#C81D25", limits = c(25, 89), name = "Age (years)") +
  scale_x_discrete(labels = donor_labels, drop = FALSE) + scale_y_discrete(drop = FALSE) +
  labs(x = "Individual / sex (ordered by age)", y = NULL) + coord_fixed() + theme_figure(9)
contribution <- support$tissue_partial_correlations
contribution <- contribution[is.finite(contribution$r), ]
contribution$signature <- factor(contribution$signature, levels = c("HF-UP", "HF-DN"))
contribution$tissue <- factor(contribution$tissue, levels = unname(support$tissue_labels[definitions$tissues]))
contribution_means <- aggregate(r ~ tissue + signature, contribution, mean)
write.csv(contribution_means, "validation/tissue_mean_partial_Pearson.csv", row.names = FALSE)
contribution_plot <- ggplot(contribution, aes(signature, r, fill = signature)) +
  geom_hline(yintercept = 0, colour = "grey85", linetype = 2, linewidth = .25) +
  geom_violin(trim = TRUE, colour = NA, alpha = .35, scale = "width") +
  geom_boxplot(width = .12, fill = "white", linewidth = .25, outlier.shape = NA) +
  geom_point(aes(colour = signature), position = position_jitter(width = .10, seed = 29), size = .4) +
  geom_text(data = contribution_means, aes(y = 1.16, label = sprintf("%.2f", r), colour = signature), size = 2.8, show.legend = FALSE) +
  scale_fill_manual(values = c("HF-UP" = "#D86D70", "HF-DN" = "#829CC8")) +
  scale_colour_manual(values = c("HF-UP" = "#A92C23", "HF-DN" = "#315995")) +
  facet_wrap(~tissue, nrow = 1) + scale_y_continuous(breaks = seq(-1, 1, .5)) +
  labs(x = NULL, y = "mRNA ~ age partial Pearson r") + theme_figure(9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank())

# Reconstruct weighted running enrichment without altering duplicate-symbol
# first-match semantics in the frozen 25,987-entry Wald vector.
running_enrichment <- function(statistic, genes) {
  hits <- sort(match(intersect(genes, names(statistic)), names(statistic)))
  step <- rep(-1 / (length(statistic) - length(hits)), length(statistic))
  step[hits] <- abs(statistic[hits]) / sum(abs(statistic[hits]))
  es <- c(0, cumsum(step))
  display_indices <- sort(unique(c(0, length(statistic), hits, pmax(hits - 1, 0), seq(0, length(statistic), by = 25))))
  list(curve = data.frame(rank = display_indices, ES = es[display_indices + 1L]), hits = hits,
       ES = if (max(es) > -min(es)) max(es) else min(es))
}
plot_enrichment <- function(id) {
  rank <- support$main_rank
  z <- running_enrichment(rank, signature_sets[[id]])
  s <- main_statistics[main_statistics$signature == id, ]
  stopifnot(abs(z$ES - s$ES) < 1e-8)
  colour <- if (z$ES > 0) "#D83C60" else "#315995"
  p <- ggplot(z$curve, aes(rank, ES)) +
    geom_area(fill = colour, alpha = .15) + geom_line(colour = colour, linewidth = .55) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = .2) +
    scale_x_continuous(limits = c(0, length(rank)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-.65, .65), breaks = c(-.5, 0, .5)) +
    labs(title = paste(final_label(id), final_letter(id)), x = NULL, y = "ES",
         subtitle = sprintf("NES = %+.2f\nFDR = %.2g", s$NES, s$q_multilevel_global_4721)) +
    theme_figure(9) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  hits <- ggplot(data.frame(rank = z$hits), aes(rank)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1), linewidth = .15) +
    scale_x_continuous(limits = c(0, length(rank)), expand = c(0, 0)) + theme_void()
  indices <- sort(unique(c(1, length(rank), which(diff(sign(rank)) != 0), which(diff(sign(rank)) != 0) + 1, seq(1, length(rank), by = 25))))
  metric <- data.frame(rank = indices, statistic = as.numeric(rank[indices]))
  strip <- ggplot(metric, aes(rank, statistic)) +
    geom_area(aes(y = pmax(statistic, 0)), fill = "#CE7D81") +
    geom_area(aes(y = pmin(statistic, 0)), fill = "#859BC1") +
    scale_x_continuous(limits = c(0, length(rank)), breaks = c(0, 10000, 20000, length(rank)), labels = c("0", "10k", "20k", "26k"), expand = c(0, 0)) +
    labs(x = paste0("Rank (", s$n_defined, " defined genes)"), y = "Wald") + theme_figure(8)
  grobs <- lapply(list(p, hits, strip), ggplotGrob)
  widths <- do.call(unit.pmax, lapply(grobs, `[[`, "widths"))
  grobs <- lapply(grobs, function(g) { g$widths <- widths; g })
  arrangeGrob(grobs = grobs, ncol = 1, heights = c(3, .3, 1.25))
}
ordered_signatures <- signature_key$id[c(3, 4, 5, 6, 1, 2)]
enrichment_plots <- lapply(ordered_signatures, plot_enrichment)

# The final S3f contains only full-definition overlap, not leading-edge overlap.
overlap <- expand.grid(a = ordered_signatures, b = ordered_signatures, stringsAsFactors = FALSE)
overlap$n <- mapply(function(a, b) length(intersect(signature_sets[[a]], signature_sets[[b]])), overlap$a, overlap$b)
overlap$jaccard <- mapply(function(a, b) length(intersect(signature_sets[[a]], signature_sets[[b]])) / length(union(signature_sets[[a]], signature_sets[[b]])), overlap$a, overlap$b)
overlap$jaccard[overlap$a == overlap$b] <- NA_real_
overlap$a <- factor(overlap$a, levels = ordered_signatures)
overlap$b <- factor(overlap$b, levels = rev(ordered_signatures))
overlap_plot <- ggplot(overlap, aes(a, b, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = .3) +
  geom_text(aes(label = n), size = 3) +
  scale_fill_gradient(low = "white", high = "#23616D", na.value = "#D2D2D2", name = "Jaccard") +
  scale_x_discrete(labels = final_label) + scale_y_discrete(labels = final_label) + coord_fixed() +
  labs(x = NULL, y = NULL) + theme_figure(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
