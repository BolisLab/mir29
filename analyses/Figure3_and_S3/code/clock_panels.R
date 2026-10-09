# Final clock summaries and their matched-control inferential results.
main_clocks <- read_input("main_clocks")
clock_tests <- read_input("clock_stars")
mobile_clocks <- read_input("mobile_clocks")

mean_sem <- function(samples, value, groups) {
  do.call(rbind, lapply(split(samples, interaction(samples[groups], drop = TRUE)), function(d) {
    v <- d[[value]]
    cbind(d[1, groups, drop = FALSE], data.frame(n = length(v), mean = mean(v), sem = sd(v) / sqrt(length(v))))
  }))
}
quadratic_fit <- function(d, group_name, x_name = "x") {
  do.call(rbind, lapply(split(d, d[[group_name]]), function(z) {
    x <- z[[x_name]]
    model <- lm(z$mean ~ poly(x, 2, raw = TRUE))
    grid_x <- seq(min(x), max(x), length.out = 201)
    data.frame(group = z[[group_name]][1], x = grid_x,
               fitted = as.numeric(predict(model, newdata = data.frame(x = grid_x))))
  }))
}
plot_clock <- function(d, fit, title, y_label) {
  ggplot(d, aes(x, mean, colour = group, group = group)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = group), alpha = .12, colour = NA) +
    geom_line(aes(linetype = group), linewidth = .25) +
    geom_line(data = fit, aes(x, fitted, colour = group, linetype = group), inherit.aes = FALSE, linewidth = .85) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .12, linewidth = .25, linetype = "solid") +
    geom_point(size = 1.9, shape = 16) +
    scale_colour_manual(values = arm_colours) + scale_fill_manual(values = arm_colours, guide = "none") +
    scale_linetype_manual(values = arm_lines) +
    scale_x_continuous(breaks = time_positions, labels = paste0("T", 1:6)) +
    labs(x = NULL, y = y_label, title = title, colour = NULL, linetype = NULL) + theme_figure(10)
}
main_clock_panels <- list()
for (clock in c("RNAAgeCalc", "tAge")) {
  d <- if (clock == "RNAAgeCalc") main_clocks$rna_summary else main_clocks$summary_T1
  d$group <- d$grp
  d$x <- d$xdisp
  fit <- quadratic_fit(d, "group")
  main_clock_panels[[clock]] <- plot_clock(d, fit,
    if (clock == "tAge") "tAge (YuGene)" else clock,
    if (clock == "tAge") "Delta age from baseline (years)" else "Transcriptomic age (years)")
}

# Supplementary curves retain native predictions and separate control arms.
# Saved fitted y values are reused, so the main-style x spacing cannot refit them.
supplementary_clock_panels <- list()
for (clock in c("tAge", "RNAAgeCalc", "PASTA")) {
  p <- mobile_clocks$plots[[clock]]
  d <- as.data.frame(p$data)
  if (!"hours" %in% names(d)) d$hours <- d$TIME_H
  if (!"sem" %in% names(d)) d$sem <- d$SEM
  layer <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomSmooth"), logical(1)))
  fit <- p$layers[[layer]]$data
  if (clock == "tAge") {
    baseline <- d[d$TIME == 1 & d$GROUP %in% c("inhCTRL", "mimCTRL"), ]
    offset <- setNames(baseline$mean, baseline$arm)
    group_arm <- setNames(as.character(d$arm), d$GROUP)
    d$mean <- d$mean - offset[as.character(d$arm)]
    fit$y <- fit$y - offset[group_arm[fit$GROUP]]
    stopifnot(max(abs(d$mean[d$TIME == 1 & d$GROUP %in% c("inhCTRL", "mimCTRL")])) < 1e-10)
  }
  d$group <- as.character(d$GROUP)
  d$x <- map_time(d$hours)
  f <- data.frame(group = fit$GROUP, x = map_time(fit$x), fitted = fit$y)
  supplementary_clock_panels[[clock]] <- plot_clock(d, f,
    if (clock == "tAge") "tAge (YuGene)" else clock,
    if (clock == "tAge") "Delta age from baseline (years)" else "Transcriptomic age (years)")
}
write.csv(clock_tests$table, "validation/clock_timepoint_statistics.csv", row.names = FALSE)
write.csv(main_clocks$qvalues, "validation/clock_global_statistics.csv", row.names = FALSE)
