# Main Figure 4 only: compare own-arm and pooled-control inference.
# Run from the bundle root, after 15_main_panel_statistics.R.
# Frozen scores, beta inputs, probe selection, clock predictions and figures
# remain unchanged. This is a control-definition sensitivity, not raw-vs-adjusted.
suppressPackageStartupMessages({library(lme4); library(pbkrtest); library(emmeans)})
previous <- readRDS("validation/main_panel_statistics.rds")

models <- model_inputs <- diagnostics <- validation <- new_rows <- list()
fit_one <- function(formula, dd, key) {
  messages <- character()
  do_fit <- function(optimizer) withCallingHandlers(
    lmer(formula, data = dd, REML = TRUE, control = lmerControl(optimizer = optimizer,
      optCtrl = if (optimizer == "bobyqa") list(maxfun = 200000) else list(maxeval = 200000))),
    warning = function(w) {messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")},
    message = function(w) {messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleMessage")})
  ok <- function(f) all(f@optinfo$conv$opt == 0) &&
    all(grepl("boundary.*singular", f@optinfo$conv$lme4$messages))
  fit <- do_fit("bobyqa")
  retry <- !ok(fit)
  if (retry) fit <- do_fit("nloptwrap")
  X <- model.matrix(nobars(formula), dd)
  stopifnot(ok(fit), qr(X)$rank == ncol(X), is.finite(sigma(fit)), sigma(fit) > 0)
  v <- as.data.frame(VarCorr(fit))
  models[[key]] <<- fit
  diagnostics[[key]] <<- data.frame(model = key, n = nrow(dd), runs = nlevels(dd$run),
    singular = isSingular(fit), run_variance = v$vcov[v$grp == "run"],
    residual_variance = sigma(fit)^2, optimizer = fit@optinfo$optimizer, retried = retry,
    messages = paste(unique(messages), collapse = " | "))
  fit
}

endpoints <- unique(previous$primary$endpoint[previous$primary$panel != "f"])
for (endpoint in endpoints) {
  all_data <- bind_rows(previous$inputs[[paste0(endpoint, "__INH")]],
                        previous$inputs[[paste0(endpoint, "__MIM")]])
  stopifnot(nrow(all_data) == 72L, !anyDuplicated(all_data$sample))
  for (arm in c("INH", "MIM")) {
    key <- paste(endpoint, arm, "pooled", sep = "__")
    message("Testing ", key)
    reference <- previous$primary[previous$primary$endpoint == endpoint & previous$primary$arm == arm, ]
    stopifnot(nrow(reference) == 1L)
    dd <- all_data[all_data$GROUP %in% c("CTRL", arm), ]
    # The other treatment is excluded: the only additional observations relative
    # to the own-control test are the 18 controls from the other experiment arm.
    dd$condition <- factor(ifelse(dd$GROUP == arm, "treated", "control"),
                           levels = c("control", "treated"))
    dd$run <- droplevels(factor(dd$BATCH))
    dd$time <- factor(dd$TP2, levels = 1:6)
    dd <- dd[order(dd$run, dd$time, dd$condition), ]
    dd$score_center <- mean(dd$score)
    dd$score_scale <- sd(dd$score)
    dd$y <- (dd$score - dd$score_center) / dd$score_scale
    counts <- with(dd, table(run, time, condition))
    stopifnot(nrow(dd) == 54L, nlevels(dd$run) == 6L,
              all(table(dd$condition) == c(control = 36L, treated = 18L)),
              all(counts[, , "control"] == 1L),
              all(rowSums(counts[, , "treated"]) %in% c(0L, 6L)),
              sum(rowSums(counts[, , "treated"]) == 6L) == 3L,
              all(is.finite(dd$score)), !anyDuplicated(dd$sample))
    model_inputs[[key]] <- dd
    full_formula <- if (reference$hypothesis == "time_interaction")
      y ~ condition * time + (1 | run) else y ~ condition + time + (1 | run)
    reduced_formula <- if (reference$hypothesis == "time_interaction")
      y ~ condition + time + (1 | run) else y ~ time + (1 | run)
    full <- fit_one(full_formula, dd, paste0(key, "__full"))
    reduced <- fit_one(reduced_formula, dd, paste0(key, "__reduced"))
    kr <- KRmodcomp(full, reduced)
    ks <- kr$stats
    cn <- names(fixef(full))
    selected <- if (reference$hypothesis == "time_interaction")
      grep("^conditiontreated:time", cn) else which(cn == "conditiontreated")
    L <- diag(length(cn))[selected, , drop = FALSE]
    explicit_p <- KRmodcomp(full, L)$stats$p.value
    stopifnot(ks$ndf == reference$numerator_df, length(selected) == ks$ndf,
              is.finite(ks$p.value), abs(log(explicit_p) - log(ks$p.value)) < 1e-6)
    validation[[key]] <- data.frame(test = key, explicit_contrast_p = explicit_p,
      nested_model_p = ks$p.value, log_p_difference = abs(log(explicit_p) - log(ks$p.value)))
    row <- reference
    row$comparison <- paste0(arm, " vs pooled CTRL")
    row$n_control <- 36L
    row$n_runs <- 6L
    row$estimate <- if (reference$hypothesis == "average_shift")
      unname(fixef(full)["conditiontreated"]) * dd$score_scale[1] else NA_real_
    row$statistic <- unname(ks$Fstat)
    row$numerator_df <- unname(ks$ndf)
    row$denominator_df <- unname(ks$ddf)
    row$p <- unname(ks$p.value)
    row$singular <- isSingular(full) || isSingular(reduced)
    row$control_definition <- "pooled"
    new_rows[[key]] <- row
  }
}

# f was already pooled. Add own-arm T6 contrasts so all endpoints can be read
# in the same paired comparison table. The original pooled f results are reused.
for (endpoint in c("entropy", "FAE")) {
  all_data <- previous$inputs[[paste0("f__", endpoint)]]
  for (arm in c("INH", "MIM")) {
    key <- paste(endpoint, arm, "own", sep = "__")
    message("Testing ", key)
    reference <- previous$primary[previous$primary$endpoint == endpoint & previous$primary$arm == arm, ]
    dd <- all_data[all_data$EXPERIMENT == tolower(arm), ]
    dd$group <- factor(dd$GROUP, levels = c("CTRL", arm))
    dd$run <- droplevels(factor(dd$BATCH))
    dd$y <- as.numeric(scale(dd$score))
    stopifnot(nrow(dd) == 6L, nlevels(dd$run) == 3L, all(table(dd$run, dd$group) == 1L))
    fit <- fit_one(y ~ group + (1 | run), dd, key)
    model_inputs[[key]] <- dd
    emm <- emmeans(fit, ~ group, lmer.df = "kenward-roger")
    ct <- contrast(emm, method = list("treated vs own CTRL" = c(-1, 1)), adjust = "none")
    cc <- as.data.frame(summary(ct, infer = c(TRUE, TRUE), adjust = "none"))
    row <- reference
    row$comparison <- paste0(arm, " vs own CTRL")
    row$n_control <- 3L
    row$n_runs <- 3L
    row$estimate <- cc$estimate * sd(dd$score)
    row$statistic <- cc$t.ratio
    row$denominator_df <- cc$df
    row$p <- cc$p.value
    row$singular <- isSingular(fit)
    row$input_state <- "adjusted noncollapsed beta, sample entropy/FAE, own-arm CTRL"
    row$control_definition <- "own"
    new_rows[[key]] <- row
  }
}

old <- previous$primary
old$control_definition <- ifelse(old$panel == "f", "pooled", "own")
out <- bind_rows(old, bind_rows(new_rows))
out <- out[order(match(out$endpoint, unique(previous$primary$endpoint)),
                 match(out$arm, c("INH", "MIM")), match(out$control_definition, c("own", "pooled"))), ]
families <- interaction(out$panel, out$control_definition, drop = TRUE)
out$q_BH_panel <- ave(out$p, families, FUN = function(x) p.adjust(x, "BH"))
out$BH_family_n <- ave(out$p, families, FUN = length)
out$q_BH_all14_audit <- ave(out$p, out$control_definition, FUN = function(x) p.adjust(x, "BH"))
out$stars_q_BH_panel <- ifelse(out$q_BH_panel < .0001, "****", ifelse(out$q_BH_panel < .001, "***",
  ifelse(out$q_BH_panel < .01, "**", ifelse(out$q_BH_panel < .05, "*", "ns"))))
out$batch_model <- "random intercept for experimental run/BATCH"
stopifnot(nrow(out) == 28L, all(is.finite(out$p)), all(out$p >= 0 & out$p <= 1),
          all(out$q_BH_panel >= out$p - 1e-14), all(table(out$control_definition) == 14L))
# No previously delivered P or panel-wise q is silently changed.
kk <- function(x) paste(x$endpoint, x$arm, x$control_definition, sep = "__")
ix <- match(kk(old), kk(out))
stopifnot(!anyNA(ix), identical(out$p[ix], old$p),
          max(abs(out$q_BH_panel[ix] - old$q_BH_panel)) < 1e-14)

wide <- out %>% select(panel, endpoint, arm, hypothesis, control_definition,
                      p, q_BH_panel, stars_q_BH_panel, estimate, n_control, singular) %>%
  pivot_wider(names_from = control_definition,
              values_from = c(p, q_BH_panel, stars_q_BH_panel, estimate, n_control, singular))
save_table(out, "main_statistics_own_vs_pooled_long")
save_table(wide, "main_statistics_own_vs_pooled_comparison")
save_table(bind_rows(diagnostics), "main_statistics_own_vs_pooled_diagnostics")
save_table(bind_rows(validation), "main_statistics_own_vs_pooled_validation")
# Descriptive check of the control arms, without new significance tests.
controls <- bind_rows(lapply(endpoints, function(endpoint) {
  z <- bind_rows(previous$inputs[[paste0(endpoint, "__INH")]],
                 previous$inputs[[paste0(endpoint, "__MIM")]])
  z$endpoint <- endpoint
  z[z$GROUP == "CTRL", ]
})) %>% group_by(endpoint, EXPERIMENT, TP2) %>%
  summarise(n = n(), mean = mean(score), sd = sd(score), .groups = "drop")
save_table(controls, "main_statistics_control_arm_descriptives")
saveRDS(list(results = out, comparison = wide, models = models, inputs = model_inputs,
             diagnostics = diagnostics, validation = validation), "validation/main_statistics_own_vs_pooled.rds")
capture.output(sessionInfo(), file = "validation/main_statistics_own_vs_pooled_sessionInfo.txt")
print(as.data.frame(wide[, c("panel", "endpoint", "arm", "p_own", "q_BH_panel_own", "p_pooled", "q_BH_panel_pooled")]),
      row.names = FALSE, digits = 6)
