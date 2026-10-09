# Requested main-panel statistics only: d, f, h, i-right, and l.
# Run from the reconstruction bundle root. Cached scores and figures are unchanged.
suppressPackageStartupMessages({
  library(lme4)
  library(pbkrtest)
  library(emmeans)
})

# Planned hypotheses, fixed before examining the resulting p-values:
# d/l: treatment-by-time interaction (5 df), allowing a baseline group offset.
# h/i: average treatment shift controlling categorical time (1 df), as in the
#      transcriptomic clocks. This is NOT a test of curve-shape equality.
# f: two-sided treatment contrasts against one pooled CTRL group at T6.
# Primary p-values use Kenward-Roger small-sample inference. Historical-style
# ML likelihood-ratio p-values are saved separately, never substituted post hoc.
# BH families: d=2, f=4, h=2, i-right=2, l=4. An all-14 BH audit is also saved.

meta <- metadata()
meta$bio_rep <- stringr::str_extract(meta$SHORT_ID, "P\\d+|exp\\d+")
stopifnot(!anyNA(meta$bio_rep), !anyDuplicated(meta$sample))
rep_map <- unique(meta[, c("bio_rep", "BATCH", "EXPERIMENT")])
stopifnot(nrow(rep_map) == 6L, !anyDuplicated(rep_map$bio_rep),
          !anyDuplicated(rep_map$BATCH))

models <- inputs <- diagnostics <- checks <- primary <- audit <- list()
fit_one <- function(formula, dd, key, reml = TRUE) {
  messages <- character()
  run_fit <- function(optimizer) withCallingHandlers(
      lmer(formula, data = dd, REML = reml,
           control = lmerControl(optimizer = optimizer,
             optCtrl = if (optimizer == "bobyqa") list(maxfun = 200000) else list(maxeval = 200000))),
      warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleMessage") })
  fit <- run_fit("bobyqa")
  converged <- function(f) all(f@optinfo$conv$opt == 0) &&
    all(grepl("boundary.*singular", f@optinfo$conv$lme4$messages))
  retried <- !converged(fit)
  if (retried) {
    messages <- c(messages, paste("Retrying", key, "with nloptwrap after convergence warning."))
    fit <- run_fit("nloptwrap")
  }
  conv <- fit@optinfo$conv$lme4$messages
  nonboundary <- conv[!grepl("boundary.*singular", conv)]
  xx <- model.matrix(nobars(formula), data = dd)
  stopifnot(qr(xx)$rank == ncol(xx), length(nonboundary) == 0L,
            all(fit@optinfo$conv$opt == 0), is.finite(sigma(fit)), sigma(fit) > 0)
  vv <- as.data.frame(VarCorr(fit))
  diagnostics[[key]] <<- data.frame(model = key, REML = reml, n = nrow(dd),
    n_runs = nlevels(dd$run), fixed_rank = qr(xx)$rank, singular = isSingular(fit),
    optimizer = fit@optinfo$optimizer, retried = retried,
    run_variance = vv$vcov[vv$grp == "run"], residual_variance = sigma(fit)^2,
    messages = paste(unique(c(messages, conv)), collapse = " | "))
  models[[key]] <<- fit
  fit
}

prepare <- function(z, field, key, arm) {
  dd <- z[z$EXPERIMENT == tolower(arm), , drop = FALSE]
  dd$score <- dd[[field]]
  dd$condition <- factor(ifelse(dd$GROUP == arm, "treated", "control"),
                         levels = c("control", "treated"))
  dd$run <- droplevels(factor(dd$BATCH))
  dd$time <- factor(dd$TP2, levels = 1:6)
  dd <- dd[order(dd$run, dd$time, dd$condition), ]
  counts <- with(dd, table(run, time, condition))
  stopifnot(nrow(dd) == 36L, nlevels(dd$run) == 3L, all(counts == 1L),
            all(is.finite(dd$score)), !anyDuplicated(dd$sample),
            all(dd$GROUP %in% c("CTRL", arm)))
  # Scale only for numerical conditioning; linear rescaling leaves tests invariant.
  dd$score_center <- mean(dd$score)
  dd$score_scale <- sd(dd$score)
  dd$y <- (dd$score - dd$score_center) / dd$score_scale
  inputs[[key]] <<- dd[, c("sample", "BATCH", "GROUP", "EXPERIMENT", "TP2",
                           "score", "score_center", "score_scale", "y", "condition", "run", "time")]
  checks[[key]] <<- data.frame(endpoint_arm = key, n = nrow(dd), n_runs = nlevels(dd$run),
    n_times = 6L, n_control = sum(dd$condition == "control"),
    n_treated = sum(dd$condition == "treated"), matched_cells_complete = all(counts == 1L))
  dd
}

curve_test <- function(z, field, panel, endpoint, hypothesis) {
  for (arm in c("INH", "MIM")) {
    key <- paste(endpoint, arm, sep = "__")
    message("Testing ", key)
    dd <- prepare(z, field, key, arm)
    null <- y ~ time + (1 | run)
    shift <- y ~ condition + time + (1 | run)
    interaction <- y ~ condition * time + (1 | run)
    full_formula <- if (hypothesis == "average_shift") shift else interaction
    reduced_formula <- if (hypothesis == "average_shift") null else shift
    full <- fit_one(full_formula, dd, paste0(key, "__full_REML"))
    reduced <- fit_one(reduced_formula, dd, paste0(key, "__reduced_REML"))
    kr <- KRmodcomp(full, reduced)
    ks <- kr$stats
    models[[paste0(key, "__KR")]] <<- kr
    p <- unname(ks$p.value)
    expected_df <- if (hypothesis == "average_shift") 1L else 5L
    stopifnot(is.finite(p), p >= 0, p <= 1, ks$ndf == expected_df)

    # Compatibility audit with the former ML-LRT inference convention.
    full_ml <- fit_one(full_formula, dd, paste0(key, "__full_ML"), FALSE)
    reduced_ml <- fit_one(reduced_formula, dd, paste0(key, "__reduced_ML"), FALSE)
    lr <- anova(reduced_ml, full_ml, refit = FALSE)
    stopifnot(lr$Df[2] == expected_df, lr$Chisq[2] >= -1e-8)
    shift_value <- if (hypothesis == "average_shift")
      unname(fixef(full)["conditiontreated"]) * dd$score_scale[1] else NA_real_

    # An arbitrary condition-specific baseline shift must leave the interaction
    # hypothesis invariant. This validates use of the cached delta-beta data in l.
    invariant <- NA_real_
    if (hypothesis == "time_interaction") {
      X <- model.matrix(nobars(full_formula), dd)
      X0 <- model.matrix(nobars(reduced_formula), dd)
      shifted <- dd$y + ifelse(dd$condition == "treated", .731, -.426)
      invariant <- max(abs(qr.resid(qr(X), dd$y) - qr.resid(qr(X), shifted)),
                       abs(qr.resid(qr(X0), dd$y) - qr.resid(qr(X0), shifted)))
      stopifnot(invariant < 1e-10)
    }
    primary[[key]] <<- data.frame(panel = panel, endpoint = endpoint, comparison = paste0(arm, " vs own CTRL"),
      arm = arm, hypothesis = hypothesis, test = "Kenward-Roger F", n_treated = 18L,
      n_control = 18L, n_runs = 3L, n_timepoints = 6L, estimate = shift_value,
      statistic = unname(ks$Fstat), numerator_df = unname(ks$ndf), denominator_df = unname(ks$ddf),
      p = p, singular = isSingular(full) || isSingular(reduced),
      input_state = if (panel == "d") "raw collapsed beta, sample means" else
        if (panel == "l") "adjusted noncollapsed beta, fixed CTRL-selected probes, condition T1 centered" else
          "frozen clock predictions from adjusted beta")
    audit[[key]] <<- data.frame(endpoint = endpoint, arm = arm, hypothesis = hypothesis,
      ML_LRT_chisq = lr$Chisq[2], ML_LRT_df = lr$Df[2], ML_LRT_p = lr$`Pr(>Chisq)`[2],
      baseline_shift_residual_invariance_error = invariant)
  }
}

global <- readRDS("data/primary/global_raw_collapsed_metrics.rds")
curve_test(global, "mean_beta", "d", "Global methylation", "time_interaction")
clocks <- readRDS("data/primary/clock_predictions.rds")
curve_test(clocks, "SkinBlood", "h", "Skin & Blood", "average_shift")
curve_test(clocks, "RetroAge", "i-right", "RetroAge", "average_shift")

movers <- readRDS("data/primary/ctrl_mover_trajectories.rds")$replicate_means
movers <- left_join(movers, rep_map, by = "bio_rep")
key_meta <- paste(meta$BATCH, meta$GROUP, meta$TP2, sep = "__")
key_movers <- paste(movers$BATCH, movers$GROUP, movers$TP2, sep = "__")
stopifnot(!anyDuplicated(key_meta), !anyNA(match(key_movers, key_meta)))
movers$sample <- meta$sample[match(key_movers, key_meta)]
for (direction in c("Hypo", "Hyper"))
  curve_test(movers[movers$direction == direction, ], "rep_mean", "l",
             paste(direction, "CTRL movers"), "time_interaction")

# T6: all six CTRL samples are a single fixed-effect group. Run intercepts retain
# the known pairing of each treatment with three of those controls; this avoids
# pretending partially paired observations are independent Welch-test samples.
metrics <- readRDS("data/primary/hgf_compartment_metrics_adjusted_verified_annotation.rds")
t6 <- subset(metrics, comp == "All_CpGs" & TP2 == 6)
stopifnot(nrow(t6) == 12L, all(table(t6$GROUP) == c(6L, 3L, 3L)))
for (field in c("entropy", "FAE")) {
  key <- paste0("f__", field)
  dd <- t6
  dd$score <- dd[[field]]
  dd$group <- factor(dd$GROUP, levels = c("CTRL", "INH", "MIM"))
  dd$run <- droplevels(factor(dd$BATCH))
  dd$y <- as.numeric(scale(dd$score))
  scale_sd <- sd(dd$score)
  stopifnot(all(table(dd$run) == 2L), !anyDuplicated(dd$sample), all(is.finite(dd$y)))
  fit <- fit_one(y ~ group + (1 | run), dd, paste0(key, "__REML"))
  inputs[[key]] <- dd[, c("sample", "BATCH", "GROUP", "EXPERIMENT", "TP2", "score", "y", "group", "run")]
  em <- emmeans(fit, ~ group, lmer.df = "kenward-roger")
  contrasts <- contrast(em, method = list("INH vs pooled CTRL" = c(-1, 1, 0),
                                          "MIM vs pooled CTRL" = c(-1, 0, 1)), adjust = "none")
  cc <- as.data.frame(summary(contrasts, infer = c(TRUE, TRUE), adjust = "none"))
  models[[paste0(key, "__contrasts")]] <- contrasts
  for (i in seq_len(nrow(cc))) {
    arm <- c("INH", "MIM")[i]
    primary[[paste(key, arm, sep = "__")]] <- data.frame(panel = "f", endpoint = field,
      comparison = cc$contrast[i], arm = arm, hypothesis = "T6_mean_difference",
      test = "Two-sided mixed-model t, Kenward-Roger df", n_treated = 3L, n_control = 6L,
      n_runs = 6L, n_timepoints = 1L, estimate = cc$estimate[i] * scale_sd,
      statistic = cc$t.ratio[i], numerator_df = 1, denominator_df = cc$df[i],
      p = cc$p.value[i], singular = isSingular(fit),
      input_state = "adjusted noncollapsed beta, sample entropy/FAE, all CTRL pooled")
  }
}

out <- bind_rows(primary)
out <- out[order(match(out$panel, c("d", "f", "h", "i-right", "l")),
                 match(out$endpoint, c("Global methylation", "entropy", "FAE", "Skin & Blood", "RetroAge", "Hypo CTRL movers", "Hyper CTRL movers")),
                 match(out$arm, c("INH", "MIM"))), ]
out$q_BH_panel <- ave(out$p, out$panel, FUN = function(p) p.adjust(p, "BH"))
out$BH_family_n <- ave(out$p, out$panel, FUN = length)
out$q_BH_all14_audit <- p.adjust(out$p, "BH")
stars <- function(p) ifelse(p < .0001, "****", ifelse(p < .001, "***",
                   ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
out$stars_q_BH_panel <- stars(out$q_BH_panel)
stopifnot(nrow(out) == 14L, all(is.finite(out$p)), all(out$q_BH_panel >= out$p - 1e-14),
          all(out$BH_family_n == c(2, 2, rep(4, 4), 2, 2, 2, 2, rep(4, 4))))
# Verify the nested-model hypotheses against explicitly selected coefficients,
# independently of the null-model formula construction above.
verification <- bind_rows(lapply(which(out$hypothesis != "T6_mean_difference"), function(i) {
  key <- paste(out$endpoint[i], out$arm[i], sep = "__")
  fit <- models[[paste0(key, "__full_REML")]]
  cn <- names(fixef(fit))
  selected <- if (out$hypothesis[i] == "average_shift") which(cn == "conditiontreated") else
    grep("^conditiontreated:time", cn)
  stopifnot(length(selected) == out$numerator_df[i])
  L <- diag(length(cn))[selected, , drop = FALSE]
  explicit_p <- KRmodcomp(fit, L)$stats$p.value
  err <- abs(log(explicit_p) - log(out$p[i]))
  stopifnot(err < 1e-6)
  data.frame(endpoint = out$endpoint[i], arm = out$arm[i], explicit_contrast_p = explicit_p,
             absolute_log_p_difference = err)
}))
save_table(out, "main_panel_statistics_primary")
save_table(verification, "main_panel_statistics_hypothesis_validation")
save_table(bind_rows(audit), "main_panel_statistics_ML_LRT_audit")
save_table(bind_rows(diagnostics), "main_panel_statistics_diagnostics")
save_table(bind_rows(checks), "main_panel_statistics_input_checks")
saveRDS(list(primary = out, ML_audit = bind_rows(audit), models = models, inputs = inputs,
             checks = checks, diagnostics = diagnostics, hypothesis_validation = verification),
        "validation/main_panel_statistics.rds")
capture.output(sessionInfo(), file = "validation/main_panel_statistics_sessionInfo.txt")
print(as.data.frame(out[, c("panel", "endpoint", "arm", "p", "q_BH_panel", "stars_q_BH_panel", "singular")]),
      row.names = FALSE, digits = 6)
