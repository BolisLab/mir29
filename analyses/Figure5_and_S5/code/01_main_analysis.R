# Build final panel tables from validated intermediates. Preserve the historical estimands, cohorts
# and score definitions. Do not recompute VST or resampling.
message("Preparing panels a-c")
a_cache <- read_data(paste0("miR29_Horvath_repeated_resampling/", "MIR29_CLOCK_RESULTS_Horvath_AGE_Acceleration_Pearson_Spearman_v1.rds"))
a_all <- do.call(rbind, lapply(a_cache, `[[`, "correlations"))
a_sample_n <- sum(vapply(a_cache, function(x) nrow(x$annotation), integer(1)))
rm(a_cache)
invisible(gc())
a_pan <- a_all[a_all$miRNA == "pan_miR29", ]
a_pan <- a_pan[order(a_pan$Cancer), ]
a_cancers <- a_pan$Cancer
check("a: 30 cancer types", length(a_cancers) == 30L)
a_metrics <- c("AGE", "Horvath", "Acceleration")
a_bubbles <- do.call(rbind, lapply(a_metrics, function(metric) {
  data.frame(Cancer = a_pan$Cancer, Metric = metric, rho = a_pan[[paste0("Spearman_", metric, "_Fisher")]])
}))
a_summary <- data.frame(Metric = a_metrics, rho = vapply(a_metrics, function(m) {
  fisher_mean(a_bubbles$rho[a_bubbles$Metric == m])
}, numeric(1)))
check("a: displayed summary correlations", identical(unname(round(a_summary$rho, 2)), c(0.09, 0.26, 0.23)),
  paste(round(a_summary$rho, 6), collapse = ", "))
redundant <- c("hsa-mir-29b-1", "hsa-mir-29b-2", "miR29AC")
a_distribution <- a_all[!a_all$miRNA %in% c(redundant, "pan_miR29") & is.finite(a_all$Spearman_Horvath_Fisher),
  c("Cancer", "miRNA", "N", "Spearman_Horvath_Fisher")]
names(a_distribution)[4] <- "rho"
a_rank <- do.call(rbind, lapply(split(a_all, a_all$miRNA), function(x) {
  r <- x$Spearman_Horvath_Fisher
  data.frame(miRNA = x$miRNA[1], K = sum(is.finite(r)), rho = fisher_mean(r))
}))
a_rank <- a_rank[!a_rank$miRNA %in% redundant & is.finite(a_rank$rho), ]
a_rank <- a_rank[order(a_rank$rho), ]
a_rank$Rank <- seq_len(nrow(a_rank))
check("a: no OV in final analysis", !"OV" %in% a_cancers)
mir_labels <- c(`hsa-mir-29a` = "miR-29a", miR29b_combined = "miR-29b", `hsa-mir-29c` = "miR-29c", pan_miR29 = "pan-miR-29")
a_rank$Feature <- unname(mir_labels[a_rank$miRNA])
a_distribution$Feature <- unname(mir_labels[a_distribution$miRNA])

b_final <- read_data("CCLE_miR29_DNAmAge_FINAL/CCLE_miR29_Age_Horvath_Skin_with_acceleration_N20_FINAL.rds")
b_tissues <- sort(b_final$settings$EligibleTissues)
b_bubbles <- subset(b_final$PerTissue, miRNA == "pan_miR29" & Phenotype %in% c("Age", "Horvath", "Skin"))
b_bubbles <- data.frame(Tissue = b_bubbles$Tissue, Metric = b_bubbles$Phenotype, rho = b_bubbles$rho_Fisher_resampled)
b_summary <- data.frame(Metric = c("Age", "Horvath", "Skin"), rho = vapply(c("Age", "Horvath", "Skin"),
  function(m) fisher_mean(b_bubbles$rho[b_bubbles$Metric == m]), numeric(1)))
b_rank <- b_final$Rankings$SkinAcceleration
b_rank <- b_rank[order(b_rank$rho_Fisher_equalTissue), ]
b_rank$rho <- b_rank$rho_Fisher_equalTissue
b_rank$Rank <- seq_len(nrow(b_rank))
b_labels <- c(MIMAT0000086 = "miR-29a", MIMAT0000100 = "miR-29b", MIMAT0000681 = "miR-29c", pan_miR29 = "pan-miR-29")
b_rank$Feature <- unname(b_labels[b_rank$miRNA])
check("b: 11 lineages and MIN_N=15", length(b_tissues) == 11L && b_final$settings$MIN_N == 15)
check("b: displayed summary correlations", identical(unname(round(b_summary$rho, 2)), c(0.03, 0.27, 0.41)))
check("b: pan-miR29 top Skin&Blood acceleration rank", tail(b_rank$miRNA, 1) == "pan_miR29")

c_final <- read_data("FIG6C_TCGA_RNA_FINAL/FIG6C_TCGA_RNA_Horvath_FINAL.rds")
c_rank <- c_final$meta
recalculated_cutoff <- mean(c_rank$rho) - 3 * sd(c_rank$rho)
check("c: exact negative-tail threshold", equal_numeric(recalculated_cutoff, c_final$threshold$cutoff))
c_185 <- c_rank[c_rank$rho < recalculated_cutoff, ]
check("c: 185 selected genes", nrow(c_185) == 185L)
check("c: selected identifiers equal saved selection", identical(c_185$FeatureID, c_final$agingRepressed_3SD$FeatureID))
c_matrix <- c_final$rhoCancerMatrix[c_rank$FeatureID, , drop = FALSE]
c_rebuilt <- tanh(rowMeans(atanh(pmin(pmax(c_matrix, -0.999999), 0.999999))))
check("c: meta-correlations rebuilt from per-cancer matrix", equal_numeric(c_rebuilt, c_rank$rho))
rm(c_matrix, c_rebuilt)
writeLines(c_185$FeatureID, file.path(BUNDLE, "data", "panel_c_185_ensembl_ids.txt"))
writeLines(unique(na.omit(c_185$SYMBOL)), file.path(BUNDLE, "data", "panel_c_185_gene_symbols.txt"))

message("Preparing paired expression panels e-f")
ef_final <- read_data("FIG6E_FINAL16_miR29_ERC_ge5pairs/FIG6E_FINAL16_miR29_ERC_ge5pairs.rds")
prepared <- setNames(lapply(CANCERS16, function(ca) read_data(file.path("FIG6E_TCGA_miR29_ERC/cache",
  paste0(ca, "_FIG6E_v1.rds")))), CANCERS16)
for (ca in CANCERS16) {
  x <- prepared[[ca]]
  check(paste0("MRC genes: ", ca), identical(names(attr(x$erc, "ERC_ID")), MRC_GENES))
}
ef_env <- new.env(parent = globalenv())
ef_env$MASTER_CANCERS16 <- CANCERS16
ef_env$ORGAN_MAP <- ef_final$settings$ORGAN_MAP
load_helpers(c("one_sample_per_patient", "run_histotype_meta", "build_organ_effects", "run_organ_meta",
  "weighted_kde2d"), ef_env, "paired_expression_and_density.R")
pair_tables <- lapply(CANCERS16, function(ca) {
  x <- prepared[[ca]]
  mT <- ef_env$one_sample_per_patient(x$mir, "01")
  mN <- ef_env$one_sample_per_patient(x$mir, "11")
  rT <- ef_env$one_sample_per_patient(x$erc, "01")
  rN <- ef_env$one_sample_per_patient(x$erc, "11")
  patients <- sort(Reduce(intersect, list(mT$Patient, mN$Patient, rT$Patient, rN$Patient)))
  d <- data.frame(Cancer = ca, Patient = patients, Pan29_Normal = mN$Pan29_Z[match(patients, mN$Patient)],
    Pan29_Tumor = mT$Pan29_Z[match(patients, mT$Patient)], ERC_Normal = rN$ERC_Z[match(patients,
      rN$Patient)], ERC_Tumor = rT$ERC_Z[match(patients, rT$Patient)])
  d$Delta_Pan29 <- d$Pan29_Tumor - d$Pan29_Normal
  d$Delta_ERC <- d$ERC_Tumor - d$ERC_Normal
  saved <- ef_final$paired_patient_level[ef_final$paired_patient_level$Cancer == ca, ]
  saved <- saved[match(d$Patient, saved$Patient), ]
  check(paste0("e/f exact paired scores: ", ca), equal_numeric(as.matrix(d[, 3:8]), as.matrix(saved[,
    3:8])))
  d
})
ef_pairs <- do.call(rbind, pair_tables)
check("e/f: 580 strict pairs", nrow(ef_pairs) == 580L)
e_effects <- do.call(rbind, lapply(CANCERS16, function(ca) do.call(rbind, lapply(c("pan-miR-29", "ERC"),
  function(readout) {
    x <- ef_pairs[ef_pairs$Cancer == ca, if (readout == "ERC")
      "Delta_ERC" else "Delta_Pan29"]
    tt <- t.test(x, mu = 0)
    data.frame(Cancer = ca, Readout = readout, Effect = mean(x), Lower = tt$conf.int[1], Upper = tt$conf.int[2],
      SE = sd(x)/sqrt(length(x)), N_pairs = length(x), P = tt$p.value)
  }))))
e_effects$FDR <- ave(e_effects$P, e_effects$Readout, FUN = function(p) p.adjust(p, "BH"))
for (col in c("Effect", "Lower", "Upper", "SE", "N_pairs", "P", "FDR")) check(paste0("e exact ", col),
  equal_numeric(e_effects[[col]], ef_final$paired[[col]]))
e_meta <- do.call(rbind, lapply(c("pan-miR-29", "ERC"), function(r) {
  h <- ef_env$run_histotype_meta(e_effects, r)
  o <- ef_env$run_organ_meta(e_effects, r)$meta
  do.call(rbind, lapply(c("Histotype", "Organ"), function(level) {
    m <- if (level == "Histotype")
      h else o
    data.frame(Readout = r, Level = level, K = m$k, Effect = as.numeric(m$b), CI_low = m$ci.lb, CI_high = m$ci.ub,
      P = m$pval)
  }))
}))
for (col in c("Effect", "CI_low", "CI_high", "P")) check(paste0("e exact meta ", col), equal_numeric(e_meta[[col]],
  subset(ef_final$meta_summary, Analysis == "Paired")[[col]]))
f_means <- merge(subset(e_effects, Readout == "pan-miR-29", c(Cancer, Effect)), subset(e_effects, Readout ==
  "ERC", c(Cancer, Effect)), by = "Cancer", suffixes = c("_miRNA", "_MRC"))
f_test <- cor.test(f_means$Effect_miRNA, f_means$Effect_MRC, method = "pearson")
check("f Pearson coefficient", abs(unname(f_test$estimate) - (-0.8063488)) < 1e-07)
check("f Pearson p-value", abs(f_test$p.value - 0.0001613243) < 1e-10)
f_weights <- 1/as.numeric(table(ef_pairs$Cancer)[ef_pairs$Cancer])
f_kde <- ef_env$weighted_kde2d(ef_pairs$Delta_Pan29, ef_pairs$Delta_ERC, f_weights, n = 200)
check("f equal cancer KDE weight", length(unique(round(tapply(f_weights, ef_pairs$Cancer, sum), 12))) ==
  1L)

message("Rebuilding both final age-adjusted fidelity forests")
h_env <- new.env(parent = globalenv())
load_helpers(c("tcga_sample_id_local", "tcga_patient_id_local", "tcga_type_code_local", "prepare_exact_clock_overlap",
  "build_analysis_data", "run_spline_residualized_spearman"), h_env, "fidelity_mirna.R")
load_helpers(c("tcga_sample_id_mrc", "tcga_patient_id_mrc", "tcga_type_code_mrc", "get_age_table_mrc",
  "get_mrc_table", "get_fidelity_table_mrc", "build_mrc_fidelity_data", "run_mrc_fidelity_cancer"),
  h_env, "fidelity_mrc.R")
h_env$allPred450K <- read_data("rstudio-export-7/allPred450K.RDS")
h_env$FIG6E_FINAL <- list(perCancer = prepared)
h_env$TCGA_miR29 <- setNames(vector("list", 15), CANCERS15)
h_env$TCGA_fidelity_main15 <- h_env$TCGA_miR29
reference_sets <- list()
reference_qc <- list()
old_cache <- "miR29_TCGA_epigenetic_fidelity_MAIN14/cache_final_autosomal_fidelity_MAIN14_N10_v1"
new_cache <- "miR29_TCGA_epigenetic_fidelity_FINAL16FLOW/cache_final_autosomal_fidelity_16flow_ge5_v2"
for (ca in CANCERS15) {
  dir <- if (ca == "CHOL")
    new_cache else old_cache
  h_env$TCGA_miR29[[ca]] <- read_data(file.path(dir, paste0("miR29_", ca, ".rds")))
  if (ca == "KICH") {
    x <- read_data(file.path(old_cache, "kidney3_shared_KIRC_KIRP_reference.rds"))
    fid <- list(status = x$status, scores = x$scores$KICH, stable_CpGs = x$stable_CpGs, reference = x$reference)
    reference_qc[[ca]] <- data.frame(Cancer = ca, Reference = "Shared KIRC+KIRP normal", N_CpGs = length(x$stable_CpGs))
  } else {
    file <- if (ca %in% c("KIRC", "KIRP"))
      paste0("sensitivity_", ca, "_own_reference.rds") else if (ca == "CHOL")
      "fidelity_own_CHOL_ge5pairs.rds" else paste0("fidelity_own_", ca, ".rds")
    fid <- read_data(file.path(dir, file))
    reference_qc[[ca]] <- data.frame(Cancer = ca, Reference = "Own adjacent normal", N_CpGs = length(fid$stable_CpGs))
  }
  check(paste0("Fidelity cache status: ", ca), identical(fid$status, "OK") && identical(h_env$TCGA_miR29[[ca]]$status,
    "OK"))
  reference_sets[[ca]] <- fid$stable_CpGs
  h_env$TCGA_fidelity_main15[[ca]] <- fid[c("status", "scores")]
}
h_left <- do.call(rbind, lapply(CANCERS15, h_env$run_spline_residualized_spearman))
h_right <- do.call(rbind, lapply(CANCERS15, h_env$run_mrc_fidelity_cancer))
check("h left: exact historical cohorts", identical(as.integer(h_left$N), as.integer(EXPECTED_H_LEFT_N[h_left$Cancer])))
check("h right: exact recovered cohorts", identical(as.integer(h_right$N), as.integer(EXPECTED_H_RIGHT_N[h_right$Cancer])))
h_organ_map <- setNames(CANCERS15, CANCERS15)
h_organ_map[c("KICH", "KIRC", "KIRP")] <- "Kidney"
h_organ_map[c("LUAD", "LUSC")] <- "Lung"
fidelity_meta <- function(d) {
  hist <- metafor::rma(yi = d$Fisher_z, sei = d$Fisher_SE, method = "REML", test = "knha")
  organs <- do.call(rbind, lapply(split(d, h_organ_map[d$Cancer]), function(x) {
    w <- 1/x$Fisher_SE^2
    data.frame(z = sum(w * x$Fisher_z)/sum(w), se = sqrt(1/sum(w)))
  }))
  organ <- metafor::rma(yi = organs$z, sei = organs$se, method = "REML", test = "knha")
  do.call(rbind, lapply(c("Histotype", "Organ"), function(level) {
    m <- if (level == "Histotype")
      hist else organ
    data.frame(Level = level, K = m$k, rho = tanh(as.numeric(m$b)), CI_low = tanh(m$ci.lb), CI_high = tanh(m$ci.ub),
      P = m$pval)
  }))
}
hl_meta <- fidelity_meta(h_left)
hr_meta <- fidelity_meta(h_right)
check("h: 15 histotypes / 12 organs", identical(as.integer(hl_meta$K), c(15L, 12L)) && identical(as.integer(hr_meta$K),
  c(15L, 12L)))
check("h left: approved historical pooled rho", abs(hl_meta$rho[1] - (-0.1451)) < 1e-04)
check("h left: approved historical pooled p", abs(hl_meta$P[1] - 0.0188) < 1e-04)

PANELS <- list(a = list(bubbles = a_bubbles, summary = a_summary, distribution = a_distribution, ranking = a_rank,
  cancers = a_cancers, N = a_sample_n), b = list(bubbles = b_bubbles, summary = b_summary, ranking = b_rank,
  tissues = b_tissues), c = list(ranking = c_rank, selected = c_185, threshold = recalculated_cutoff,
  targets = c_final$targets, network = c_final$oldCoreNetwork), e = list(effects = e_effects, meta = e_meta),
  f = list(means = f_means, pairs = ef_pairs, kde = f_kde, pearson = f_test), g = list(reference_qc = do.call(rbind,
    reference_qc), CpGs = reference_sets), h = list(left = h_left, right = h_right, left_meta = hl_meta,
    right_meta = hr_meta))
saveRDS(PANELS, file.path(BUNDLE, "data", "Figure5_panel_data.rds"))
jsonlite::write_json(CHECKS, file.path(BUNDLE, "validation", "numerical_checks.json"), pretty = TRUE,
  auto_unbox = TRUE)
writeLines(vapply(CHECKS, function(x) paste(if (x$pass) "PASS" else "FAIL", x$check, x$detail, sep = " | "),
  character(1)), file.path(BUNDLE, "validation", "numerical_checks.txt"))
writeLines(INPUT_FILES, file.path(BUNDLE, "provenance", "input_paths.txt"))
writeLines(capture.output(print(h_left), print(h_right), print(hl_meta), print(hr_meta)), file.path(BUNDLE,
  "validation", "fidelity_results.txt"))
