# Reconstruct matched-control matrices and row dendrograms.
tf <- read_input("tf")
progeny <- read_input("progeny")
heatmap_checkpoint <- read_input("heatmaps")
matched_control <- function(scores, metadata) {
  scores <- scores[, rownames(metadata), drop = FALSE]
  cell <- paste(metadata$GROUP, metadata$TIME, sep = "_")
  means <- vapply(unique(cell), function(g) rowMeans(scores[, cell == g, drop = FALSE]), numeric(nrow(scores)))
  rownames(means) <- rownames(scores)
  result <- cbind(means[, paste0("INH_", 1:6)] - means[, paste0("inhCTRL_", 1:6)],
                  means[, paste0("MIM_", 1:6)] - means[, paste0("mimCTRL_", 1:6)])
  colnames(result) <- c(paste0("INH_", 1:6), paste0("MIM_", 1:6))
  result
}
tf_matrix <- matched_control(tf$primary$raw_scores, tf$metadata)
progeny_matrix <- matched_control(progeny$scores, progeny$metadata)
stopifnot(max(abs(tf_matrix[rownames(tf$main), ] - tf$main)) < 1e-10,
          max(abs(progeny_matrix - progeny$main)) < 1e-10)
palette <- heatmap_checkpoint$palette
colour_activity <- circlize::colorRamp2(palette$values, palette$colours, space = "LAB")

dendrogram_segments <- function(tree) {
  labels <- labels(tree)
  segments <- list()
  xcoord <- function(h) -.18 - 2.45 * h / attr(tree, "height")
  recurse <- function(node) {
    if (is.leaf(node)) return(list(x = xcoord(0), y = length(labels) + 1 - match(attr(node, "label"), labels)))
    children <- lapply(node, recurse)
    x <- xcoord(attr(node, "height"))
    y <- vapply(children, `[[`, numeric(1), "y")
    for (child in children) segments[[length(segments) + 1L]] <<- data.frame(x = x, xend = child$x, y = child$y, yend = child$y)
    segments[[length(segments) + 1L]] <<- data.frame(x = x, xend = x, y = min(y), yend = max(y))
    list(x = x, y = mean(y))
  }
  recurse(tree)
  do.call(rbind, segments)
}
plot_activity_heatmap <- function(name) {
  m <- if (name == "TF") tf_matrix else progeny_matrix
  frozen <- heatmap_checkpoint$trees[[name]]
  m <- m[rownames(frozen$matrix), colnames(frozen$matrix), drop = FALSE]
  hc <- hclust(dist(m), method = "complete")
  stopifnot(isTRUE(all.equal(hc$height, frozen$hclust$height, tolerance = 1e-10)),
            identical(hc$merge, frozen$hclust$merge))
  tree <- frozen$dendrogram
  rows <- labels(tree)
  d <- expand.grid(row = rows, condition = colnames(m), stringsAsFactors = FALSE)
  d$value <- as.vector(m[rows, , drop = FALSE])
  d$hex <- colour_activity(d$value)
  d$x <- match(d$condition, colnames(m))
  d$x <- d$x + ifelse(d$x > 6, .4, 0)
  d$y <- length(rows) + 1 - match(d$row, rows)
  branch <- dendrogram_segments(tree)
  row_labels <- sub("TP53", "p53", rows)
  row_labels <- sub("JAK.STAT", "JAK/STAT", row_labels, fixed = TRUE)
  ggplot(d, aes(x, y)) +
    geom_segment(data = branch, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE, linewidth = .25) +
    geom_tile(aes(fill = hex), colour = "#303030", linewidth = .22) + scale_fill_identity() +
    annotate("text", x = 13.2, y = length(rows):1, label = row_labels, hjust = 0, size = 3.2) +
    annotate("text", x = c(3.5, 9.9), y = length(rows) + 1, label = c("inh", "mimic"), fontface = "italic") +
    scale_x_continuous(breaks = c(1:6, 7:12 + .4), labels = rep(1:6, 2)) +
    coord_fixed(xlim = c(-2.9, 17.3), ylim = c(0, length(rows) + 1.5), clip = "off") +
    theme_void() + theme(axis.text.x = element_text(size = 8), plot.margin = margin(5, 5, 5, 5))
}
tf_heatmap <- plot_activity_heatmap("TF")
progeny_heatmap <- plot_activity_heatmap("PROGENy")

activity <- read_input("activity")
activity_order <- c("p53", "EGFR", "MAPK", "ELK1", "MYC", "ATF4")
d <- activity$means
f <- activity$fits
d$feature <- factor(d$feature, levels = activity_order)
f$feature <- factor(f$feature, levels = activity_order)
d$x <- map_time(d$hours)
f$x <- map_time(f$hours)
activity_plot <- ggplot(d, aes(x, mean, group = group, colour = group)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = group), colour = NA, alpha = .12) +
  geom_line(aes(linetype = group), linewidth = .25) +
  geom_line(data = f, aes(x, fitted, linetype = group), linewidth = .85) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .12, linewidth = .25) +
  geom_point(size = 1.4) + facet_wrap(~feature, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = arm_colours) + scale_fill_manual(values = arm_colours, guide = "none") +
  scale_linetype_manual(values = arm_lines) +
  scale_x_continuous(breaks = time_positions, labels = paste0("T", 1:6)) +
  labs(x = NULL, y = "Inferred activity (score)", colour = NULL, linetype = NULL) + theme_figure(9)
write.csv(activity$tests, "validation/selected_activity_statistics.csv", row.names = FALSE)
