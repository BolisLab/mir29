# Figure 4o: current 123-CpG external reference, raw-beta scoring and matched inference.
# The displayed global q tests the arithmetic mean matched displacement over
# T1-T6, not a temporal slope, an interaction, or a duration-weighted AUC.
read_fidelity123 <- function(name) {
  d <- read.delim(file.path("data/fidelity123", paste0(name, ".tsv")),
                 check.names = FALSE, stringsAsFactors = FALSE,
                 na.strings = c("NA", "\\N"))
  d[, !names(d) %in% c("", "\\N"), drop = FALSE]
}
reference123 <- read_fidelity123("HGF123_reference")
stopifnot(nrow(reference123) == 123L, !anyDuplicated(reference123$probe),
          all(is.finite(reference123$reference)))
ref <- setNames(reference123$reference, reference123$probe)

# Require all 123 probes, exactly once, in every scored sample. No imputation,
# recentering, batch correction or redefinition of the centroid is performed.
score_fidelity123 <- function(name, samples) {
  b <- read_fidelity123(name)
  stopifnot(setequal(b$probe, names(ref)), setequal(b$sample, samples),
            !anyDuplicated(b[c("probe", "sample")]), all(is.finite(b$beta)),
            all(table(b$sample) == 123L), all(b$beta >= 0 & b$beta <= 1))
  b$absolute_distance <- abs(b$beta - ref[b$probe])
  b %>% group_by(sample) %>%
    summarise(reference_distance_beta = mean(absolute_distance),
              n_reference_CpGs = n(), .groups = "drop")
}
one_sample <- function(v) {
  stopifnot(length(v) == 3L, all(is.finite(v)))
  test <- t.test(v, mu = 0, alternative = "two.sided")
  data.frame(n_experiments = length(v), mean_pp = mean(v),
             sem_pp = sd(v) / sqrt(length(v)), ci_low_pp = test$conf.int[1],
             ci_high_pp = test$conf.int[2], t = unname(test$statistic),
             df = unname(test$parameter), p = test$p.value)
}
test_run_summaries <- function(d) {
  d %>% group_by(endpoint, arm) %>% group_modify(~one_sample(.x$value)) %>%
    ungroup() %>% group_by(endpoint) %>%
    mutate(q = p.adjust(p, "BH"),
           test = "Two-sided one-sample t-test of matched run summaries against zero",
           q_family = "BH across INH and MIM, separately for this endpoint") %>% ungroup()
}
validate_reference_tests <- function(observed, expected, label) {
  if ("metric" %in% names(expected))
    expected <- expected[expected$metric == "reference_distance", , drop = FALSE]
  idx <- match(paste(observed$arm, observed$endpoint), paste(expected$arm, expected$endpoint))
  stopifnot(!anyNA(idx))
  for (nm in c("mean_pp", "sem_pp", "ci_low_pp", "ci_high_pp", "p"))
    assert_close(observed[[nm]], expected[[nm]][idx], paste(label, nm), 1e-10)
  assert_close(observed$q, expected$q_two_arms[idx], paste(label, "BH q"), 1e-10)
}
q_caption <- function(d, prefix = "Global q") {
  paste0(prefix, ": INH = ", formatC(d$q[match("INH", d$arm)], digits = 3, format = "f"),
         "; MIM = ", formatC(d$q[match("MIM", d$arm)], digits = 3, format = "f"))
}

fidelity_samples <- read_fidelity123("HGF_samples")
fidelity_samples <- fidelity_samples %>%
  left_join(score_fidelity123("HGF123_beta", fidelity_samples$sample), by = "sample")
stopifnot(nrow(fidelity_samples) == 72L, !anyNA(fidelity_samples$reference_distance_beta))
controls <- fidelity_samples %>% filter(GROUP == "CTRL") %>%
  select(EXPERIMENT, BATCH, TIME, control_sample = sample,
         control_distance_beta = reference_distance_beta)
stopifnot(!anyDuplicated(controls[c("EXPERIMENT", "BATCH", "TIME")]))
fidelity_delta <- fidelity_samples %>% filter(GROUP != "CTRL") %>%
  left_join(controls, by = c("EXPERIMENT", "BATCH", "TIME")) %>%
  mutate(arm = as.character(GROUP),
         delta_pp = 100 * (reference_distance_beta - control_distance_beta),
         cycle = ceiling(TIME / 2), hours = time_hours)
stopifnot(nrow(fidelity_delta) == 36L, !anyNA(fidelity_delta$delta_pp),
          !anyDuplicated(fidelity_delta[c("arm", "BATCH", "TIME")]))
cycle_exp <- fidelity_delta %>% group_by(arm, BATCH, cycle) %>%
  summarise(value = mean(delta_pp), n_timepoints = n(), .groups = "drop") %>%
  mutate(hours = c(84, 180, 276)[cycle])
stopifnot(all(cycle_exp$n_timepoints == 2))
cycle_summary <- cycle_exp %>% group_by(arm, cycle, hours) %>%
  summarise(mean = mean(value), sem = sd(value) / sqrt(n()), n = n(), .groups = "drop")
stopifnot(all(cycle_summary$n == 3))
full_means <- fidelity_delta %>% group_by(arm, BATCH) %>%
  summarise(value = mean(delta_pp), n_timepoints = n(), .groups = "drop") %>%
  mutate(endpoint = "all_T1_T6")
late_means <- fidelity_delta %>% filter(cycle == 3) %>% group_by(arm, BATCH) %>%
  summarise(value = mean(delta_pp), n_timepoints = n(), .groups = "drop") %>%
  mutate(endpoint = "late_T5_T6")
stopifnot(all(full_means$n_timepoints == 6), all(late_means$n_timepoints == 2))
run_means <- bind_rows(full_means, late_means)
primary_tests <- test_run_summaries(run_means)
global_tests <- filter(primary_tests, endpoint == "all_T1_T6")
late_tests <- filter(primary_tests, endpoint == "late_T5_T6")
validate_reference_tests(primary_tests, read_fidelity123("HGF_tests"), "123-CpG HGF")
cycle_average <- cycle_exp %>% group_by(arm, BATCH) %>%
  summarise(value = mean(value), .groups = "drop")
idx <- match(paste(full_means$arm, full_means$BATCH), paste(cycle_average$arm, cycle_average$BATCH))
assert_close(full_means$value, cycle_average$value[idx], "HGF full mean equals equal cycle mean")

# Quadratic interpolation and pointwise SEM are descriptive, not inferential.
cycle_fit <- cycle_summary %>% group_by(arm) %>% group_modify(~ {
  model <- lm(mean ~ poly(hours, 2, raw = TRUE), .x)
  grid <- seq(min(.x$hours), max(.x$hours), length.out = 300)
  data.frame(hours = grid, mean = predict(model, data.frame(hours = grid)),
             sem = approx(.x$hours, .x$sem, grid)$y)
}) %>% ungroup()
p_center <- ggplot(cycle_summary, aes(hours, mean, color = arm, group = arm)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey65", linewidth = .35) +
  geom_ribbon(data = cycle_fit, aes(ymin = mean - sem, ymax = mean + sem, fill = arm),
              color = NA, alpha = .14) +
  geom_line(data = cycle_fit, linewidth = .85) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 8, linewidth = .4) +
  geom_point(size = 3) + scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_x_continuous(breaks = c(84, 180, 276), labels = c("T1–T2", "T3–T4", "T5–T6")) +
  labs(x = NULL, y = "Deviation vs matched CTRL (beta percentage points)",
       title = "Deviation from fibroblast reference (123 CpGs)",
       subtitle = q_caption(global_tests)) + theme_figure()
save_panel(p_center, "4o_center_fibroblast_deviation", 5.1, 3.8)

late <- fidelity_delta %>% filter(cycle == 3) %>%
  group_by(arm) %>% mutate(batch_index = match(BATCH, sort(unique(BATCH))),
                         x = match(arm, c("INH", "MIM")) + .19 * (batch_index - 2) +
                           ifelse(TIME == 5, -.035, .035)) %>% ungroup()
late_plot_means <- late %>% group_by(arm, BATCH, batch_index) %>%
  summarise(value = mean(delta_pp), .groups = "drop") %>%
  mutate(x = match(arm, c("INH", "MIM")) + .19 * (batch_index - 2))
late_group <- late_plot_means %>% group_by(arm) %>%
  summarise(value = mean(value), .groups = "drop") %>% mutate(x = match(arm, c("INH", "MIM")))
p_right <- ggplot(late, aes(x, delta_pp, color = arm)) +
  geom_hline(yintercept = 0, color = "grey65", linetype = 2) +
  geom_line(aes(group = interaction(arm, BATCH)), linewidth = .4, alpha = .4) +
  geom_point(aes(shape = factor(TIME)), size = 2.5, alpha = .8) +
  geom_point(data = late_plot_means, aes(y = value, fill = arm), shape = 23, size = 4) +
  geom_segment(data = late_group, aes(x = x - .32, xend = x + .32, y = value, yend = value)) +
  scale_shape_manual(values = c("5" = 1, "6" = 0), labels = c("T5", "T6"), name = NULL) +
  scale_color_manual(values = COLORS, guide = "none") +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("INH", "MIM"), limits = c(.55, 2.45)) +
  labs(x = NULL, y = "Deviation vs matched CTRL (beta percentage points)",
       title = "Late response: T5 + T6", subtitle = q_caption(late_tests, "Late q")) + theme_figure()
save_panel(p_right, "4o_right_late_response", 4.2, 3.8)
save_panel(plot_grid(p_center, p_right, nrow = 1, rel_widths = c(5.1, 4.2)),
           "4o_center_and_right", 9.3, 3.8)
save_table(fidelity_samples, "4o_sample_reference_distances")
save_table(fidelity_delta, "4o_matched_deltas")
save_table(cycle_exp, "4o_cycle_means_per_experiment")
save_table(cycle_summary, "4o_cycle_mean_SEM")
save_table(run_means, "4o_run_means_all_and_late")
save_table(primary_tests, "4o_123CpG_primary_tests")
save_table(global_tests, "4o_global_tests")
save_table(late_tests, "4o_late_tests")
save_table(reference123, "4o_123CpG_reference_annotation")

# Retain the existing alternative trajectory summaries as diagnostics only;
# neither AUC is used for the current central-panel q.
normalized_auc <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / diff(range(x))
}
cycle_auc <- cycle_exp %>% group_by(arm, BATCH) %>%
  summarise(value = normalized_auc(hours, value), .groups = "drop") %>%
  mutate(endpoint = "cycle_mean_AUC")
full_auc <- fidelity_delta %>% group_by(arm, BATCH) %>%
  summarise(value = normalized_auc(hours, delta_pp), .groups = "drop") %>%
  mutate(endpoint = "full_trajectory_AUC")
cycle_tests <- test_run_summaries(cycle_auc)
full_tests <- test_run_summaries(full_auc)
validate_reference_tests(bind_rows(cycle_tests, full_tests),
                         read_fidelity123("HGF_tests"), "123-CpG HGF AUC diagnostics")
save_table(cycle_tests, "4o_cycle_AUC_tests")
save_table(full_tests, "4o_full_trajectory_AUC_sensitivity")

# Reconstruct selection from independent-cohort audit tables. The final bound
# is q <= 0.20 for the sex-adjusted fibroblast contrast; other filters are fixed.
source("code/upstream/reference_selection.R")
selection <- input("fidelity/selection_checkpoint.rds")
culture <- input("fidelity/culture_checkpoint.rds")
selected123 <- fibroblast_associated_sites(selection$identity_probe_audit)
stopifnot(length(selected123) == 123L, setequal(selected123, names(ref)),
          all(selected123 %in% culture$probes),
          all(culture$probes %in% selection$probes))
assert_close(unname(ref), unname(selection$centroid[names(ref)]), "123-CpG external centroid")
selected_audit <- as.data.frame(selection$identity_probe_audit)
selected_audit <- selected_audit[match(names(ref), selected_audit$probe), ]
# Historical membership flags refer to q <= 0.05; do not export them as the
# membership rule for the adopted q <= 0.20 reference.
selected_audit <- selected_audit[, !grepl("^identity_candidate_", names(selected_audit)), drop = FALSE]
selected_audit$reference_beta <- unname(ref)
save_table(selected_audit, "4o_123CpG_selection_audit")

# Independent donor-derived hTERT support: same fixed probes and centroid,
# matched within arm, experiment and collection. These are supplementary-data
# results, not a replacement or relabelling of an official S4 figure panel.
htert_reference_samples <- read_fidelity123("hTERT_samples") %>%
  left_join(score_fidelity123("hTERT123_beta", read_fidelity123("hTERT_samples")$sample), by = "sample")
stopifnot(nrow(htert_reference_samples) == 36L, !anyNA(htert_reference_samples$reference_distance_beta))
hc <- htert_reference_samples %>% filter(GROUP == "CTRL") %>%
  select(EXPERIMENT, BATCH, TP2, control_sample = sample,
         control_distance_beta = reference_distance_beta)
stopifnot(!anyDuplicated(hc[c("EXPERIMENT", "BATCH", "TP2")]))
htert_reference_delta <- htert_reference_samples %>% filter(GROUP != "CTRL") %>%
  left_join(hc, by = c("EXPERIMENT", "BATCH", "TP2")) %>%
  mutate(arm = GROUP, round = as.integer(TP2),
         delta_pp = 100 * (reference_distance_beta - control_distance_beta))
stopifnot(nrow(htert_reference_delta) == 18L, !anyNA(htert_reference_delta$delta_pp),
          !anyDuplicated(htert_reference_delta[c("arm", "BATCH", "round")]))
htert_reference_runs <- bind_rows(
  htert_reference_delta %>% group_by(arm, BATCH) %>%
    summarise(value = mean(delta_pp), n_timepoints = n(), .groups = "drop") %>%
    mutate(endpoint = "mean_three_rounds_pp"),
  htert_reference_delta %>% filter(round == 3) %>%
    transmute(arm, BATCH, value = delta_pp, n_timepoints = 1L, endpoint = "third_cycle_96h"))
stopifnot(all(htert_reference_runs$n_timepoints[htert_reference_runs$endpoint == "mean_three_rounds_pp"] == 3L))
htert_reference_tests <- test_run_summaries(htert_reference_runs)
validate_reference_tests(htert_reference_tests, read_fidelity123("hTERT_tests"), "123-CpG hTERT")
htert_reference_summary <- htert_reference_delta %>% group_by(arm, round) %>%
  summarise(mean_pp = mean(delta_pp), sem_pp = sd(delta_pp) / sqrt(n()),
            n_experiments = n(), .groups = "drop")
stopifnot(all(htert_reference_summary$n_experiments == 3L))
save_table(htert_reference_samples, "hTERT_123CpG_sample_reference_distances")
save_table(htert_reference_delta, "hTERT_123CpG_matched_deltas")
save_table(htert_reference_runs, "hTERT_123CpG_run_means_all_and_late")
save_table(htert_reference_summary, "hTERT_123CpG_cycle_mean_SEM")
save_table(htert_reference_tests, "hTERT_123CpG_reference_tests")

# HPF transfer remains descriptive (one matched pair per timepoint).
save_table(read_fidelity123("HPF_descriptive"), "HPF_123CpG_reference_distance_descriptive")
