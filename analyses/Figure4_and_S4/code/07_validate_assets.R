# Independent scoring and selection checks using bundled compact source inputs.
source("code/upstream/reference_selection.R")
s <- input("fidelity/selection_checkpoint.rds")
culture <- input("fidelity/culture_checkpoint.rds")
a <- input("fidelity/initial_stability_audit.rds")
adopted_reference <- read.delim("data/fidelity123/HGF123_reference.tsv", check.names = FALSE,
                              na.strings = c("NA", "\\N"))
selected123 <- fibroblast_associated_sites(s$identity_probe_audit)
stopifnot(setequal(stable_reference_sites(a), s$probes),
          setequal(culture_stable_sites(culture$audit), culture$probes),
          length(selected123) == 123L,
          setequal(selected123, adopted_reference$probe),
          all(selected123 %in% culture$probes))
assert_close(s$centroid[adopted_reference$probe], adopted_reference$reference,
             "123-CpG centroid equals independent external reference")
cat("Three external-reference selection filters; adopted 123 CpGs\tPASS\t0\n",
    file = "validation/numeric_checks.tsv", append = TRUE)

model <- input("models/PASTA_linear_predictor.rds")
pasta <- input("clocks/htert_PASTA_inputs.rds")
stopifnot(setequal(names(model$weights), rownames(pasta$rank_matrix)))
pasta_prediction <- as.numeric(model$intercept + model$weights[rownames(pasta$rank_matrix)] %*% pasta$rank_matrix) * model$scale
assert_close(pasta_prediction, pasta$prediction, "hTERT PASTA frozen predictor")

horvath <- input("models/horvath_calibrated_input.rds")
models <- input("clocks/auc_inputs.rds")$models
for (nm in c("SkinBlood", "Horvath2013")) {
  model <- models[[nm]]
  x <- horvath$beta[names(model$w), , drop = FALSE]
  if (anyNA(x)) {
    missing <- which(is.na(x), arr.ind = TRUE)
    x[missing] <- rowMeans(x, na.rm = TRUE)[missing[, 1]]
  }
  predictor <- as.numeric(model$intercept + model$w %*% x)
  age <- ifelse(predictor < 0, 21 * exp(predictor) - 1, 21 * predictor + 20)
  expected <- if (nm == "SkinBlood") {
    scores <- input("primary/clock_predictions.rds")
    scores$SkinBlood[match(colnames(x), scores$sample)]
  } else horvath$age[colnames(x)]
  assert_close(age, expected, paste(nm, "all-sample calibrated predictor"))
}
retro <- input("primary/retroage_model_and_input.rds")
weights <- setNames(retro$coefficients$estimate, retro$coefficients$term)
predicted <- as.numeric(weights["(Intercept)"] + weights[rownames(retro$input)] %*% retro$input)
assert_close(predicted, retro$score[colnames(retro$input)], "RetroAge all-sample predictor")

# Primary statistics: own-control values keep the historical test families.
current <- readRDS("validation/main_statistics_own_vs_pooled.rds")$results
expected <- input("primary/main_statistics_own_vs_pooled.rds")$results
key <- function(x) paste(x$endpoint, x$arm, x$control_definition)
assert_close(current$p, expected$p[match(key(current), key(expected))], "Main-panel p values", 1e-5)
current$final_panel <- c(d = "4d / S4a", f = "4f", h = "4i-right / S4c",
                        `i-right` = "4j-right / S4e", l = "4l")[current$panel]
save_table(current, "final_main_statistics_all_control_definitions")
save_table(subset(current, control_definition == "own"), "final_main_statistics_own_controls")
