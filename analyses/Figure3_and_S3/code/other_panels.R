# Growth, PCA and selected-gene displays from analysis-ready measurements.
growth <- readRDS("data/main/panel_b_audit.rds")
source("code/growth_data.R")
growth_plot <- function(d, colour, title, raw) {
  d$group <- as.character(d$group)
  long <- normalize_growth(raw)
  check <- mean_sem(long[is.finite(long$y), ], "y", c("group", "time"))
  keys <- function(x) paste(x$group, x$time)
  stopifnot(max(abs(check$mean[match(keys(d), keys(check))] - d$mean)) < 1e-10)
  fit <- fit_growth(long)
  cols <- setNames(c("grey50", colour), c("ctrl", setdiff(unique(d$group), "ctrl")))
  ggplot(d, aes(time, mean, colour = group, group = group)) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 4, linewidth = .3) +
    geom_point(size = 1.8) + geom_line(linewidth = .25, linetype = 2) +
    geom_line(data = fit, aes(y = yhat), linewidth = .65) +
    scale_colour_manual(values = cols) + labs(x = "Time (h)", y = "Relative cell growth (%)", title = title, colour = NULL) + theme_figure()
}
growth_panels <- list(growth_plot(growth$summary_inh, "#328A37", "Inhibitor", growth_raw$INH),
                     growth_plot(growth$summary_mim, "#923D8E", "Mimic", growth_raw$MIM))

# Recount differential expression directly from the twelve audited contrasts.
de <- read_input("differential_expression")
de_counts <- t(vapply(de$results, function(z) c(up = sum(z$padj < .05 & z$log2FoldChange > 0, na.rm = TRUE),
  down = sum(z$padj < .05 & z$log2FoldChange < 0, na.rm = TRUE)), numeric(2)))
stopifnot(all(unname(de_counts) == unname(de$stored_counts)))
de_display <- data.frame(arm = rep(c("INH", "MIM"), each = 6), time = rep(1:6, 2), de_counts)
de_long <- rbind(transform(de_display, count = up, direction = "UP"), transform(de_display, count = -down, direction = "DN"))
de_plot <- ggplot(de_long, aes(time, count, fill = direction)) + geom_col(width = .65) +
  geom_hline(yintercept = 0, linewidth = .3) + facet_wrap(~arm, nrow = 1) +
  scale_fill_manual(values = c(UP = "#E5262F", DN = "#284FA0")) +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  scale_y_continuous(labels = abs) + labs(x = NULL, y = "Differentially expressed genes", fill = NULL) + theme_figure()
write.csv(de_display, "validation/differential_expression_counts.csv", row.names = FALSE)

# The frozen PCA was fitted to the 20,000 most variable corrected VST genes.
pca_checkpoint <- readRDS("data/main/panel_d_audit.rds")
expression_input <- read_input("expression")
top_genes <- names(sort(matrixStats::rowVars(expression_input$vst_corrected, useNames = TRUE), decreasing = TRUE))[1:20000]
stopifnot(identical(top_genes, pca_checkpoint$top_gene_ids))
pca_recomputed <- prcomp(t(expression_input$vst_corrected[top_genes, ]), center = TRUE, scale. = TRUE)
stopifnot(max(abs(pca_recomputed$sdev - pca_checkpoint$pca$sdev)) < 1e-8)
pca_data <- pca_checkpoint$coordinates
pca_data$fill <- "grey50"
for (group in unique(pca_data$GROUP)) {
  colours <- if (group == "INH") c("#d6ebe4", "#228536") else if (group == "MIM") c("#f3dfee", "#57338A") else c("#d9d9d9", "#343434")
  take <- pca_data$GROUP == group
  pca_data$fill[take] <- scales::col_numeric(colours, domain = c(1, 6))(pca_data$TIME[take])
}
pca_plot <- ggplot(pca_data, aes(PC1, PC2)) +
  stat_ellipse(aes(group = reGROUP, colour = reGROUP), type = "norm", level = .99, linetype = 2, linewidth = .35, show.legend = FALSE) +
  geom_point(aes(fill = fill, shape = GROUP), size = 2, colour = "grey35", stroke = .35) +
  scale_fill_identity() + scale_colour_manual(values = c(INH = "#228536", CTRL = "black", MIM = "#57338A"), guide = "none") +
  scale_shape_manual(values = c(INH = 21, MIM = 21, inhCTRL = 22, mimCTRL = 23)) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_checkpoint$variance_percent[1]),
       y = sprintf("PC2 (%.1f%%)", pca_checkpoint$variance_percent[2]), fill = NULL, shape = NULL) + theme_figure()
gene_input <- read_input("genes")
gene_summary <- subset(gene_input$summaries, method == "VST_batch_corrected")
gene_summary$gene <- factor(gene_summary$gene, levels = c("EDA2R", "PTCHD4", "DUSP6", "EREG"))
gene_plot <- ggplot(gene_summary, aes(TIME, mean, colour = GROUP3, group = GROUP3)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = GROUP3), alpha = .12, colour = NA) +
  geom_line(linewidth = .25) + geom_point(size = 1.3) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .10, linewidth = .2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = .6) +
  facet_wrap(~gene, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = arm_colours) + scale_fill_manual(values = arm_colours, guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  labs(x = NULL, y = "mRNA (VST)", colour = NULL) + theme_figure()

# Jaccard values retain a common off-diagonal scale across treatment arms.
jaccard <- read_input("overlap")$original_data
recomputed_jaccard <- list()
for (arm in c("INH", "MIM")) {
  tables <- de$results[if (arm == "INH") 1:6 else 7:12]
  universe <- Reduce(intersect, lapply(tables, function(z)
    rownames(z)[is.finite(z$padj) & is.finite(z$log2FoldChange)]))
  significant <- vapply(tables, function(z) z[universe, "padj"] < .05, logical(length(universe)))
  for (t1 in 1:6) for (t2 in 1:6) {
    n_union <- sum(significant[, t1] | significant[, t2])
    value <- if (n_union) sum(significant[, t1] & significant[, t2]) / n_union else NA_real_
    recomputed_jaccard[[length(recomputed_jaccard) + 1L]] <- data.frame(arm, t1, t2, value, n_union)
  }
}
jaccard_check <- do.call(rbind, recomputed_jaccard)
stopifnot(identical(is.na(jaccard_check$value), is.na(jaccard$jaccard)),
          max(abs(jaccard_check$value - jaccard$jaccard), na.rm = TRUE) < 1e-14,
          all(jaccard_check$n_union == jaccard$union))
upper <- max(jaccard$jaccard[jaccard$t1 != jaccard$t2], na.rm = TRUE)
jaccard$colour <- "#D0D0D0"
for (arm in c("INH", "MIM")) {
  take <- jaccard$arm == arm & jaccard$t1 != jaccard$t2
  endpoint <- if (arm == "INH") "#248636" else "#923d8e"
  jaccard$colour[take] <- scales::col_numeric(c("white", endpoint), domain = c(0, upper))(jaccard$jaccard[take])
}
jaccard$area <- ifelse(jaccard$t1 == jaccard$t2, upper, jaccard$jaccard)
jaccard$text <- ifelse(jaccard$t1 == jaccard$t2, "", sprintf("%.2f", jaccard$jaccard))
jaccard_plot <- ggplot(jaccard, aes(t1, t2)) +
  geom_point(aes(size = area, fill = colour), shape = 21, stroke = 0) +
  geom_text(aes(label = text), size = 2.2) +
  scale_fill_identity() + scale_size_area(max_size = 10, limits = c(0, upper), guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  scale_y_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  facet_wrap(~arm, nrow = 1) + coord_fixed() + labs(x = NULL, y = NULL) + theme_figure()

# Preserved Cytoscape coordinates and node mappings define the network panels.
networks <- read_input("networks")
network_plot <- function(d, title) {
  n <- d$nodes
  e <- d$edges
  ggplot() +
    geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey70", linewidth = .12, alpha = .6) +
    ggforce::geom_circle(data = n, aes(x0 = display_x, y0 = display_y, r = diameter / 2, fill = fill, colour = border), n = 48, size = .08) +
    scale_fill_identity() + scale_colour_identity() + coord_fixed() + theme_void() + labs(title = title) +
    theme(plot.title = element_text(size = 9, hjust = .5))
}
network_panels <- lapply(names(networks$prepared), function(nm)
  network_plot(networks$prepared[[nm]], sub("MIR_", "", nm)))

# Expression fold changes for the epigenetic-regulator panel.
alluvial <- readRDS("data/main/alluvial/tf_n_o_results.rds")
regulator_matrix <- alluvial$o_log2fc
regulators <- as.data.frame(as.table(regulator_matrix))
names(regulators) <- c("gene", "condition", "log2FC")
regulator_plot <- ggplot(regulators, aes(condition, gene, fill = log2FC)) +
  geom_tile(colour = "grey25", linewidth = .2) +
  scale_fill_gradient2(low = "#262064", mid = "white", high = "#791B2D", limits = c(-1, 1), oob = scales::squish, name = "log2 fold change") +
  scale_y_discrete(limits = rev(rownames(regulator_matrix))) +
  labs(x = NULL, y = NULL) + coord_fixed() + theme_figure(9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())

# Signed-Wald transcriptome percentiles for all tested miR-29 targets.
target_data <- core$data$activity_results$target_ecdf
target_data$time <- as.integer(sub(".*_", "", target_data$condition))
target_data$rank_percentile <- 100 * target_data$percentile
target_data$display <- ifelse(target_data$group == "miR-29 targets", target_data$arm, "Background")
target_data$x <- target_data$time + ifelse(target_data$group == "miR-29 targets", .16, -.16)
target_plot <- ggplot(target_data, aes(x, rank_percentile, group = interaction(time, group), fill = display)) +
  geom_boxplot(width = .27, outlier.shape = NA, linewidth = .25) +
  geom_point(data = subset(target_data, group == "miR-29 targets"), aes(colour = arm),
    position = position_jitter(width = .07, seed = 29030908), size = .16, alpha = .45) +
  facet_wrap(~arm, nrow = 1) +
  scale_fill_manual(values = c(INH = "#C6E2C9", MIM = "#E2CCE9", Background = "#CCCCCC")) +
  scale_colour_manual(values = arm_colours) +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  labs(x = NULL, y = "Relative transcriptome rank (%)") + theme_figure() + theme(legend.position = "none")
target_tests <- subset(core$data$activity_results$target_camera_tests, assumed_correlation == .01)
write.csv(target_tests, "validation/target_CAMERA_statistics.csv", row.names = FALSE)
