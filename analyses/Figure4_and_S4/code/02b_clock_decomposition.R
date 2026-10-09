suppressPackageStartupMessages(library(data.table))
bundle <- input("clocks/auc_inputs.rds")
meta <- as.data.table(bundle$metadata)
is_auc <- identical(bundle$estimand, "auc")
trap_mean <- function(y, hours) {
  stopifnot(length(y) == length(hours), all(diff(hours) > 0))
  sum(diff(hours) * (head(y, -1) + tail(y, -1)) / 2) / diff(range(hours))
}
trap_weights <- function(hours) {
  stopifnot(length(hours) >= 2L, all(diff(hours) > 0))
  interval <- diff(hours)
  (c(interval, 0) + c(0, interval)) / 2
}
order_clocks <- c("SkinBlood", "RetroAge", "CausAge", "DamAge", "Horvath2013")
n_clocks <- length(order_clocks)
clock_labels <- c(SkinBlood = "Skin & Blood", RetroAge = "RetroAge", CausAge = "CausAge",
                  DamAge = "DamAge", AdaptAge = "AdaptAge", Horvath2013 = "Horvath pan-tissue")
route_order <- c("down_negative_hyper", "down_positive_hypo", "up_positive_hyper", "up_negative_hypo")
get_route <- function(w, d) {
  fifelse(w < 0 & d > 0, "down_negative_hyper",
          fifelse(w > 0 & d < 0, "down_positive_hypo",
                   fifelse(w > 0 & d > 0, "up_positive_hyper",
                            fifelse(w < 0 & d < 0, "up_negative_hypo", "zero"))))
}
all_probes <- all_samples <- audits <- list()
for (clock_name in order_clocks) {
  model <- bundle$models[[clock_name]]
  x <- model$x
  w <- model$w[rownames(x)]
  stopifnot(all(is.finite(x)), !anyNA(w), !anyDuplicated(rownames(x)))
  pairs <- merge(meta[as.character(GROUP) == "INH", .(batch = as.character(BATCH), tp, treated = sample)],
                 meta[as.character(GROUP) == "CTRL", .(batch = as.character(BATCH), tp, control = sample)],
                 by = c("batch", "tp"))
  setorder(pairs, batch, tp)
  if (is_auc) {
    stopifnot(nrow(pairs) == 18L, all(table(pairs$batch) == 6L))
    pairs[, hours := meta$hours[match(treated, meta$sample)]]
    stopifnot(identical(pairs$hours, meta$hours[match(pairs$control, meta$sample)]),
              all(pairs$hours == c(72, 96, 168, 192, 264, 288)[pairs$tp]))
  } else {
    stopifnot(nrow(pairs) == 6L, all(table(pairs$batch) == 2L))
  }
  differences <- x[, pairs$treated, drop = FALSE] - x[, pairs$control, drop = FALSE]
  probe_rows <- list()
  for (batch_id in unique(pairs$batch)) {
    take <- pairs$batch == batch_id
    if (is_auc) {
      time_weights <- trap_weights(pairs$hours[take])
      delta <- as.numeric(differences[, take, drop = FALSE] %*% time_weights) / sum(time_weights)
    } else {
      delta <- rowMeans(differences[, take, drop = FALSE])
    }
    contribution <- w * delta
    probe_rows[[batch_id]] <- data.table(clock = clock_name, batch = batch_id,
      probe = rownames(x), weight = as.numeric(w), delta_beta = delta,
      contribution = as.numeric(contribution), route = get_route(w, delta))
  }
  all_probes[[clock_name]] <- rbindlist(probe_rows)
  linear <- as.numeric(model$intercept + w %*% x)
  predicted <- if (clock_name %in% c("SkinBlood", "Horvath2013")) ifelse(linear < 0, 21 * exp(linear) - 1, 21 * linear + 20) else linear
  samples <- data.table(clock = clock_name, sample = colnames(x), linear_predictor = linear,
                        predicted_score = predicted)
  all_samples[[clock_name]] <- merge(samples, meta[, .(sample, GROUP, BATCH, tp)], by = "sample")
  err <- max(abs(colSums(differences * w) -
                 (linear[match(pairs$treated, colnames(x))] - linear[match(pairs$control, colnames(x))])))
  stopifnot(err < 1e-8)
  if (is_auc) {
    # Independent trapezoid calculation from the reconstructed clock contrasts.
    auc_error <- max(vapply(unique(pairs$batch), function(batch_id) {
      pp <- pairs[batch == batch_id]
      expected <- trap_mean(linear[match(pp$treated, colnames(x))] -
                              linear[match(pp$control, colnames(x))], pp$hours)
      abs(expected - all_probes[[clock_name]][batch == batch_id, sum(contribution)])
    }, numeric(1)))
    stopifnot(auc_error < 1e-8)
  }
  audits[[clock_name]] <- data.table(clock = clock_name, n_present = length(w), n_required = model$required,
    n_positive_weights = sum(w > 0), n_negative_weights = sum(w < 0), n_zero_weights = sum(w == 0),
    missing_absolute_weight_fraction = model$missing_weight_fraction,
    max_delta_reconstruction_error = err, source_scale = model$scale)
  if (is_auc) audits[[clock_name]][, max_auc_reconstruction_error := auc_error]
}
probes <- rbindlist(all_probes)
if (is_auc) {
  duration <- diff(range(meta$hours))
  probes[, `:=`(auc_beta_hours = delta_beta * duration,
                 auc_contribution_hours = contribution * duration)]
}
sample_scores <- rbindlist(all_samples)
audit <- rbindlist(audits)
by_run <- probes[, .(contribution = sum(contribution), absolute_contribution = sum(abs(contribution)),
                     n_probes = .N), by = .(clock, batch, route)]
# Aggregate each probe within run BEFORE assigning its directional route:
# T5/T6 mean by default, or signed actual-time AUC / duration in AUC mode.
# Average absolute route totals over the three runs; normalize only afterwards.
grid <- CJ(clock = order_clocks, batch = unique(probes$batch), route = route_order)
by_run <- merge(grid, by_run, by = c("clock", "batch", "route"), all.x = TRUE)
for (column in c("contribution", "absolute_contribution", "n_probes"))
  set(by_run, which(is.na(by_run[[column]])), column, 0)
components <- by_run[, .(mean_contribution = mean(contribution),
                         mean_absolute = mean(absolute_contribution),
                         mean_n_probes = mean(n_probes)), by = .(clock, route)]
denom <- probes[, .(absolute_total = sum(abs(contribution)), net = sum(contribution)),
                by = .(clock, batch)][, .(denominator = mean(absolute_total)), by = clock]
components <- merge(components, denom, by = "clock")
components[, percent := 100 * mean_contribution / denominator]
stopifnot(all(denom$denominator > 0))
checks <- components[, .(sum_absolute_percent = sum(abs(percent))), by = clock]
stopifnot(max(abs(checks$sum_absolute_percent - 100)) < 1e-8)
net <- components[, .(net_percent = sum(percent),
                       down_percent = -sum(percent[percent < 0]),
                       up_percent = sum(percent[percent > 0]),
                       net_linear_shift = sum(mean_contribution)), by = clock]
net <- merge(net, audit, by = "clock")
net[, y := length(order_clocks) + 1 - match(clock, order_clocks)]
segments <- merge(components, net[, .(clock, y)], by = "clock")
segments[, route := factor(route, levels = route_order)]
setorder(segments, clock, route)
segments[, side := ifelse(percent < 0, "down", "up")]
segments[, edge_end := cumsum(percent), by = .(clock, side)]
segments[, edge_start := edge_end - percent]
segments[, `:=`(xmin = pmin(edge_start, edge_end), xmax = pmax(edge_start, edge_end))]
segments[, label := ifelse(abs(percent) >= 6, sprintf("%.0f%%", abs(percent)), "")]
segments[, label_color := ifelse(route %in% c("down_negative_hyper", "up_positive_hyper"), "white", "#263747")]
run_net <- probes[, .(net_linear_shift = sum(contribution)), by = .(clock, batch)]
run_net <- merge(run_net, denom, by = "clock")
run_net[, net_percent_common_denominator := 100 * net_linear_shift / denominator]
weight_sign <- probes[, .(contribution = sum(contribution)),
                       by = .(clock, batch, coefficient_sign = ifelse(weight > 0, "positive", "negative"))]
weight_sign <- weight_sign[, .(mean_contribution = mean(contribution)), by = .(clock, coefficient_sign)]
weight_sign <- merge(weight_sign, denom, by = "clock")
weight_sign[, net_percent := 100 * mean_contribution / denominator]

# Check observed age/score differences independently of route aggregation.
t <- sample_scores[as.character(GROUP) == "INH", .(clock, batch = as.character(BATCH), tp, treated = predicted_score)]
c <- sample_scores[as.character(GROUP) == "CTRL", .(clock, batch = as.character(BATCH), tp, control = predicted_score)]
observed <- merge(t, c, by = c("clock", "batch", "tp"))
observed[, difference := treated - control]
if (is_auc) {
  observed[, hours := meta$hours[match(paste(batch, tp), paste(as.character(meta$BATCH), meta$tp))]]
  setorder(observed, clock, batch, hours)
  observed_run <- observed[, .(mean_difference = trap_mean(difference, hours)), by = .(clock, batch)]
  observed_run[, auc_score_hours := mean_difference * duration]
} else {
  observed_run <- observed[, .(mean_difference = mean(difference)), by = .(clock, batch)]
}
observed_mean <- observed_run[, .(mean_clock_difference = mean(mean_difference),
                                sem_clock_difference = sd(mean_difference)/sqrt(.N)), by = clock]
net <- merge(net, observed_mean, by = "clock")

colors <- c(down_negative_hyper = "#286B9E", down_positive_hypo = "#B3D3E7",
            up_positive_hyper = "#C05842", up_negative_hypo = "#F0C4B6")
legend_text <- c(down_negative_hyper = "Hypermethylation; negative weight",
                 down_positive_hypo = "Hypomethylation; positive weight",
                 up_positive_hyper = "Hypermethylation; positive weight",
                 up_negative_hypo = "Hypomethylation; negative weight")
y_labels <- setNames(paste0(clock_labels[net$clock], "\n", net$n_present, "/", net$n_required, " features"), net$y)
plot <- ggplot() +
  annotate("rect", xmin = -100, xmax = 0, ymin = .5, ymax = n_clocks + .5, fill = "#F5F9FC") +
  annotate("rect", xmin = 0, xmax = 100, ymin = .5, ymax = n_clocks + .5, fill = "#FDF7F5") +
  geom_vline(xintercept = c(-100, -50, 50, 100), color = "#E1E6E9", linewidth = .3) +
  geom_vline(xintercept = 0, color = "#65727A", linewidth = .5) +
  geom_rect(data = segments, aes(xmin = xmin, xmax = xmax, ymin = y - .21, ymax = y + .21, fill = route),
            color = "white", linewidth = .35) +
  geom_text(data = segments, aes(x = (xmin + xmax)/2, y = y, label = label, color = label_color),
            size = 3.0, family = "Helvetica") +
  scale_color_identity() +
  geom_segment(data = net, aes(x = net_percent, xend = net_percent, y = y - .22, yend = y - .33),
               linewidth = .35, color = "#222222") +
  geom_point(data = net, aes(x = net_percent, y = y - .34), shape = 23, size = 2.6,
             fill = "#222222", color = "white", stroke = .35) +
  annotate("text", x = -50, y = n_clocks + .78, label = "Clock-decreasing", color = colors[1],
            size = 3.6, fontface = "bold", family = "Helvetica") +
  annotate("text", x = 50, y = n_clocks + .78, label = "Clock-increasing", color = colors[3],
            size = 3.6, fontface = "bold", family = "Helvetica") +
  scale_fill_manual(values = colors, breaks = route_order, labels = legend_text[route_order], name = NULL) +
  scale_x_continuous(limits = c(-102, 102), breaks = c(-100, -50, 0, 50, 100),
                     labels = c("-100", "-50", "0", "+50", "+100"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_len(n_clocks), labels = y_labels[as.character(seq_len(n_clocks))],
                     limits = c(.48, n_clocks + 1.02), expand = c(0, 0)) +
  labs(x = "Share of total absolute weighted contribution (%)", y = NULL,
       title = if (is_auc) "Components of the INH-induced clock shift: AUC" else "Components of the INH-induced clock shift",
       subtitle = if (is_auc) "T1-T6 signed AUC (72-288 h); own matched controls; n = 3 biological experiments" else "T5-T6 averaged within each biological experiment; own matched controls; n = 3",
       caption = if (is_auc) "Actual-time trapezoids at each CpG; AUC divided by duration. Diamonds: net contribution.\nShares describe integrated score composition, not years or evidence of rejuvenation." else "Diamonds: net shift on the same scale. Available features only; no significance filtering.\nShares describe score composition, not the magnitude of age change or rejuvenation.") +
  guides(fill = guide_legend(ncol = 2, byrow = FALSE)) +
  theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "#23333F", size = 10, lineheight = 1.0),
        axis.text.x = element_text(color = "#394B57"),
        axis.title.x = element_text(margin = margin(t = 8)),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, color = "#52626D", margin = margin(b = 8)),
        plot.caption = element_text(hjust = 0, size = 8, color = "#52626D", margin = margin(t = 8)),
        legend.position = "bottom", legend.text = element_text(size = 8.5),
        legend.key.size = unit(.38, "cm"), legend.spacing.x = unit(.12, "cm"),
        plot.margin = margin(10, 12, 8, 10))

save_panel(plot, "4k_clock_component_AUC", 8.5, 5.2)
save_table(segments, "4k_plotted_components")
save_table(net, "4k_clock_component_summary")
save_table(probes, "4k_probe_AUC_contributions")
save_table(by_run, "4k_experiment_components")
expected <- data.table::fread("data/clocks/expected_auc_components.tsv")
comparison <- merge(segments, expected, by = c("clock", "route"), suffixes = c("_new", "_old"))
assert_close(comparison$percent_new, comparison$percent_old, "Five-clock AUC route shares")
