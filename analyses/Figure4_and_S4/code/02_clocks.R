# Shared epigenetic clocks: frozen calibrated inputs; original scoring recipes.
clock_bundle <- input("clocks/primary_clocks.rds")
clock_scores <- as.data.frame(clock_bundle$scores)
clock_labels <- c(SkinBlood = "Skin & Blood DNAm age (years)", RetroAge = "RetroAge score",
                  CausAge = "CausAge score", DamAge = "Delta DamAge score", AdaptAge = "Delta AdaptAge score")
for (clock_name in names(clock_labels)) {
  d <- subset(clock_scores, clock == clock_name)
  d$score <- d$display_score
  d$series <- ifelse(d$GROUP == "CTRL", paste0(tolower(d$EXPERIMENT), "CTRL"), as.character(d$GROUP))
  panel_id <- c(SkinBlood = "S4c", RetroAge = "S4e", CausAge = "S4d_CausAge",
                DamAge = "S4d_DamAge", AdaptAge = "S4d_AdaptAge")[[clock_name]]
  save_panel(trajectory(d, clock_labels[clock_name], title = clock_name), panel_id, 5.3, 4.7)
  if (clock_name %in% c("SkinBlood", "RetroAge")) {
    d$series <- as.character(d$GROUP)
    save_panel(trajectory(d, clock_labels[clock_name], title = clock_name),
               ifelse(clock_name == "SkinBlood", "4i_right_SkinBlood", "4j_right_RetroAge"), 5.3, 4.7)
  }
}
save_panel(plot_grid(plotlist = panels[c("S4d_CausAge", "S4d_DamAge", "S4d_AdaptAge")], nrow = 1),
           "S4d_causal_clocks", 12, 4.5)
save_table(clock_scores, "epigenetic_clock_scores")
save_table(clock_bundle$coefficient_audit, "epigenetic_clock_coverage")

# Available-CpG linear predictors exactly preserve the historical missingness rule.
for (nm in names(clock_bundle$coefficients)) {
  co <- as.data.frame(clock_bundle$coefficients[[nm]])
  ids <- intersect(co$term, rownames(clock_bundle$small_beta_input))
  x <- clock_bundle$small_beta_input[ids, , drop = FALSE]
  w <- co$estimate[match(ids, co$term)]
  predicted <- as.numeric(co$estimate[co$term == "(Intercept)"] + crossprod(w, x))
  expected <- subset(clock_scores, clock == nm)
  assert_close(predicted, expected$score[match(colnames(x), expected$sample)], paste(nm, "linear predictor"))
}

# S4b: subtract the same experiment/timepoint control BEFORE averaging.
metrics <- input("primary/hgf_compartment_metrics_adjusted_verified_annotation.rds")
delta_panels <- list()
for (metric_name in c("entropy", "FAE")) {
  d <- subset(metrics, comp == "All_CpGs")
  d$value <- d[[metric_name]]
  controls <- d %>% filter(GROUP == "CTRL") %>% select(EXPERIMENT, BATCH, TP2, control = value)
  z <- d %>% filter(GROUP != "CTRL") %>% left_join(controls, by = c("EXPERIMENT", "BATCH", "TP2")) %>%
    mutate(score = (value - control) * ifelse(metric_name == "FAE", 100, 1), series = as.character(GROUP))
  stopifnot(nrow(z) == 36, !anyNA(z$score))
  delta_panels[[metric_name]] <- trajectory(z, ifelse(metric_name == "FAE", "Delta FAE (pp)", "Delta entropy"),
                                           title = metric_name) + geom_hline(yintercept = 0, linetype = 2, color = "grey60")
  save_table(z, paste0("S4b_", metric_name, "_matched_deltas"))
}
save_panel(plot_grid(plotlist = delta_panels, nrow = 1), "S4b_entropy_FAE_deltas", 9.5, 4.5)

# S4a: model-estimated own-T1 contrasts, retaining separate CTRL arms.
fit <- input("primary/global_mean_mixed_model.rds")
grid <- as.data.frame(fit$emm_full)
contrasts <- list()
for (group in unique(grid$mGROUP)) for (tp in unique(grid$TIME_factor)) {
  w <- numeric(nrow(grid))
  w[grid$mGROUP == group & grid$TIME_factor == tp] <- 1
  w[grid$mGROUP == group & grid$TIME_factor == 72] <-
    w[grid$mGROUP == group & grid$TIME_factor == 72] - 1
  contrasts[[paste(group, tp, sep = "__")]] <- w
}
est <- as.data.frame(summary(emmeans::contrast(fit$emm_full, contrasts), infer = TRUE, level = .90, adjust = "none"))
key <- strsplit(as.character(est$contrast), "__", fixed = TRUE)
est$group <- vapply(key, `[`, character(1), 1)
est$hours <- as.numeric(vapply(key, `[`, character(1), 2))
est$series <- c(inh_CTRL = "inhCTRL", inh_INH = "INH", mim_CTRL = "mimCTRL", mim_MIM = "MIM")[est$group]
est$x <- ORDINAL_X[match(est$hours, HOURS)]
save_table(est, "S4a_global_methylation_own_T1_contrasts")
p <- ggplot(est, aes(x, 1000 * estimate, color = series, group = series)) +
  geom_hline(yintercept = 0, color = "grey65", linetype = 2) +
  geom_ribbon(aes(ymin = 1000 * lower.CL, ymax = 1000 * upper.CL, fill = series), alpha = .12, color = NA) +
  geom_smooth(aes(linetype = series), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, linewidth = .8) +
  geom_errorbar(aes(ymin = 1000 * lower.CL, ymax = 1000 * upper.CL), width = .14, linewidth = .4) +
  geom_point(aes(shape = series), size = 2.6) + scale_color_manual(values = CS) +
  scale_fill_manual(values = CS, guide = "none") +
  scale_linetype_manual(values = c(INH = 1, MIM = 1, inhCTRL = 2, mimCTRL = 2), guide = "none") +
  scale_shape_manual(values = c(INH = 16, MIM = 16, inhCTRL = 1, mimCTRL = 1), guide = "none") +
  scale_x_continuous(breaks = ORDINAL_X, labels = paste0("T", 1:6)) +
  labs(x = NULL, y = "Delta mean beta from own T1 (x 1000)") + theme_figure()
save_panel(p, "S4a_global_methylation_separate_controls", 5.5, 4.6)
