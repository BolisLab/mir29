# Rebuild only the analyses retained in the final Figure S5.
EXPECTED_S5 <- readRDS("data/expected_supplementary.rds")
EXPECTED_EXT <- readRDS("data/expected_clock_extensions.rds")
EXPECTED_BOX <- readRDS("data/expected_compact.rds")

# ---- S5a: tumor cohort accounting
annotations <- lapply(names(a_cache <- read_data("miR29_Horvath_repeated_resampling/MIR29_CLOCK_RESULTS_Horvath_AGE_Acceleration_Pearson_Spearman_v1.rds")),
  function(ca) {
    z <- a_cache[[ca]]$annotation
    z$Cancer <- ca
    z
  })
cohort <- do.call(rbind, annotations)
tcga_counts <- do.call(rbind, lapply(a_cancers, function(ca) {
  z <- cohort[cohort$Cancer == ca, ]
  data.frame(Cancer = ca, N = nrow(z), Primary = sum(z$TypeCode == "01"), Metastatic = sum(z$TypeCode ==
    "06"), UniquePatients = length(unique(z$PatientID)))
}))
check("S5a: 8428 unique tumor patients", nrow(cohort) == 8428L && length(unique(cohort$PatientID)) ==
  8428L && all(tcga_counts$N == tcga_counts$UniquePatients))
save_table(tcga_counts, "S5a_TCGA_cohort_counts")

# ---- S5b: TCGA equal-cancer Fisher meta-ranks
tcga_ranks <- lapply(c("AGE", "Horvath", "Acceleration"), function(metric) {
  z <- do.call(rbind, lapply(split(a_all, a_all$miRNA), function(d) {
    r <- d[[paste0("Spearman_", metric, "_Fisher")]]
    data.frame(miRNA = d$miRNA[1], K = sum(is.finite(r)), rho = fisher_mean(r))
  }))
  z <- z[!z$miRNA %in% redundant & is.finite(z$rho), ]
  z <- z[order(z$rho, z$miRNA), ]
  z$Rank <- seq_len(nrow(z))
  z$Feature <- unname(mir_labels[z$miRNA])
  z
})
names(tcga_ranks) <- c("AGE", "Horvath", "Acceleration")
for (metric in names(tcga_ranks)) {
  z <- tcga_ranks[[metric]]
  old <- EXPECTED_EXT$TCGA$ranks[[metric]]
  check(paste("S5b TCGA rank", metric), identical(z$miRNA, old$miRNA) && equal_numeric(z$rho, old$rho))
  save_table(z, paste0("TCGA_rank_", metric))
}

# ---- S5c/e: leave-one-stratum-out, using the original equal-stratum estimand
leave_one_stratum <- function(unit, rho, n) {
  full <- fisher_mean(rho)
  do.call(rbind, lapply(c("None", unit), function(excluded) {
    keep <- if (excluded == "None")
      rep(TRUE, length(unit)) else unit != excluded
    data.frame(Excluded = excluded, K = sum(keep), N = sum(n[keep]), rho = fisher_mean(rho[keep]),
      Full_rho = full)
  }))
}
pan_a <- a_all[a_all$miRNA == "pan_miR29", ]
pan_a <- pan_a[match(a_cancers, pan_a$Cancer), ]
loo_tcga <- leave_one_stratum(pan_a$Cancer, pan_a$Spearman_Horvath_Fisher, pan_a$N)
pan_b <- subset(b_final$PerTissue, miRNA == "pan_miR29" & Phenotype == "Skin")
pan_b <- pan_b[match(b_final$settings$EligibleTissues, pan_b$Tissue), ]
loo_ccle <- leave_one_stratum(pan_b$Tissue, pan_b$rho_Fisher_resampled, pan_b$N)
check("S5c: exact 31 TCGA LOO estimates", equal_numeric(loo_tcga$rho, EXPECTED_EXT$TCGA$loo$rho))
check("S5e: exact 12 CCLE LOO estimates", equal_numeric(loo_ccle$rho, EXPECTED_EXT$CCLE$loo$rho))
save_table(loo_tcga, "S5c_TCGA_leave_one_out")
save_table(loo_ccle, "S5e_CCLE_leave_one_out")

# ---- S5d/f: eligible cell lines and four CCLE ranks
ccle_counts <- data.frame(Tissue = b_tissues, N = vapply(b_tissues, function(ti) sum(b_final$Meta$TISSUE ==
  ti), integer(1)))
check("S5d: 429 lines in 11 eligible lineages", sum(ccle_counts$N) == 429L && nrow(ccle_counts) == 11L)
ccle_ranks <- lapply(c("Age", "Skin", "Horvath", "HorvathAcceleration"), function(metric) {
  z <- b_final$Rankings[[metric]]
  z <- data.frame(miRNA = z$miRNA, K = z$K, rho = z$rho_Fisher_equalTissue)
  z <- z[order(z$rho, z$miRNA), ]
  z$Rank <- seq_len(nrow(z))
  z$Feature <- unname(b_labels[z$miRNA])
  z
})
names(ccle_ranks) <- c("Age", "Skin", "Horvath", "HorvathAcceleration")
for (metric in names(ccle_ranks)) {
  old <- EXPECTED_EXT$CCLE$ranks[[metric]]
  z <- ccle_ranks[[metric]]
  check(paste("S5f CCLE rank", metric), identical(z$miRNA, old$miRNA) && equal_numeric(z$rho, old$rho))
  save_table(z, paste0("S5f_CCLE_rank_", metric))
}
save_table(ccle_counts, "S5d_CCLE_cohort_counts")

# ---- S5g: full RNA distributions, not the old rank-sampled display
summarize_rna_boxes <- function() {
  stats <- list()
  outliers <- list()
  highlights <- list()
  for (ca in a_cancers) {
    z <- c_final$perCancer[c_final$perCancer$Cancer == ca & is.finite(c_final$perCancer$rho_Fisher_resampled),
      ]
    z$Gene <- c_final$meta$SYMBOL[match(z$FeatureID, c_final$meta$FeatureID)]
    r <- z$rho_Fisher_resampled
    q <- quantile(r, c(0.25, 0.5, 0.75), type = 7, names = FALSE)
    is_out <- r < q[1] - 1.5 * (q[3] - q[1]) | r > q[3] + 1.5 * (q[3] - q[1])
    whiskers <- range(c(q[1], q[3], r[!is_out]))
    stats[[ca]] <- data.frame(Unit = ca, N = unique(z$N), N_features = nrow(z), Q25 = q[1], Median = q[2],
      Q75 = q[3], Low = whiskers[1], High = whiskers[2], N_outliers = sum(is_out))
    outliers[[ca]] <- data.frame(Unit = ca, ID = z$FeatureID[is_out], rho = r[is_out])
    ix <- which(z$Gene %in% c(MRC_GENES, "MIR29B2CHG", "LINC-PINT"))
    highlights[[ca]] <- data.frame(Unit = ca, ID = z$FeatureID[ix], rho = r[ix], Gene = z$Gene[ix],
      Feature = ifelse(z$Gene[ix] %in% MRC_GENES, "MRC genes", z$Gene[ix]))
  }
  list(summary = do.call(rbind, stats), outliers = do.call(rbind, outliers), highlights = do.call(rbind,
    highlights))
}
rna_boxes <- summarize_rna_boxes()
check("S5g: all 30 full RNA box statistics", equal_numeric(as.matrix(rna_boxes$summary[, -1]), as.matrix(EXPECTED_BOX$TCGA_RNA$summary[,
  -1])))
check("S5g: all 33264 outliers and 210 highlights", nrow(rna_boxes$outliers) == 33264L && nrow(rna_boxes$highlights) ==
  210L)
save_table(rna_boxes$summary, "S5g_RNA_box_statistics")
save_table(rna_boxes$highlights, "S5g_RNA_highlighted_genes")

# ---- S5h/i: strict paired cohorts and Pearson leave-one-cancer-out
pair_counts <- as.data.frame(table(ef_pairs$Cancer), stringsAsFactors = FALSE)
names(pair_counts) <- c("Cancer", "N_pairs")
loo_paired <- do.call(rbind, lapply(c("None", CANCERS16), function(excluded) {
  z <- if (excluded == "None")
    f_means else f_means[f_means$Cancer != excluded, ]
  ct <- cor.test(z$Effect_miRNA, z$Effect_MRC, method = "pearson")
  data.frame(Excluded = excluded, K = nrow(z), rho = unname(ct$estimate), p = ct$p.value, Full_rho = unname(f_test$estimate))
}))
old <- subset(EXPECTED_S5$robustness, Method == "pearson")
check("S5i: exact Pearson LOO estimates", equal_numeric(loo_paired$rho, old$rho))
save_table(pair_counts, "S5h_pair_counts")
save_table(loo_paired, "S5i_paired_Pearson_leave_one_out")

# ---- S5j: selected-reference QC, not an independent validation set
reference_metrics <- readRDS("data/reference_selection_metrics.rds")
reference_qc <- list()
reference_distributions <- list()
for (ca in CANCERS15) {
  if (ca == "KICH") {
    rel <- file.path(old_cache, "kidney3_shared_KIRC_KIRP_reference.rds")
    x <- read_data(rel)
    qc <- data.frame(Cancer = ca, N_normal_reference = x$N_normal, N_stable_CpGs = length(x$stable_CpGs),
      N_probes_after_NA = length(x$common_CpGs), Reference = "Shared KIRC+KIRP")
  } else {
    fname <- if (ca %in% c("KIRC", "KIRP"))
      paste0("sensitivity_", ca, "_own_reference.rds") else if (ca == "CHOL")
      "fidelity_own_CHOL_ge5pairs.rds" else paste0("fidelity_own_", ca, ".rds")
    rel <- file.path(if (ca == "CHOL")
      new_cache else old_cache, fname)
    x <- read_data(rel)
    qc <- x$qc[, c("Cancer", "N_normal_reference", "N_stable_CpGs", "N_probes_after_NA")]
    qc$Reference <- "Own adjacent normal"
  }
  ref <- reference_metrics[[rel]]
  stopifnot(!is.null(ref))
  for (metric in c("normal_sd_beta", "age_cor", "age_delta_beta_decade")) {
    v <- ref[[metric]]
    if (is.null(names(v))) {
      stopifnot(length(v) == length(ref$stable_CpGs))
      names(v) <- ref$stable_CpGs
    }
    v <- abs(v[x$stable_CpGs])
    stopifnot(all(is.finite(v)))
    cutoff <- switch(metric, normal_sd_beta = 0.05, age_cor = 0.25, age_delta_beta_decade = 0.01)
    check(paste("S5j original cutoff", ca, metric), all(v <= cutoff + 1e-12))
    q <- quantile(v, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), names = FALSE)
    reference_distributions[[paste(ca, metric)]] <- data.frame(Cancer = ca, Metric = metric, Min = q[1],
      P05 = q[2], Q25 = q[3], Median = q[4], Q75 = q[5], P95 = q[6], Max = q[7], N = length(v))
  }
  reference_qc[[ca]] <- qc
}
reference_qc <- do.call(rbind, reference_qc)
reference_distributions <- do.call(rbind, reference_distributions)
check("S5j all reference counts", equal_numeric(reference_qc$N_stable_CpGs, EXPECTED_S5$reference$qc$N_stable_CpGs))
check("S5j all distribution summaries", equal_numeric(as.matrix(reference_distributions[, -c(1, 2)]),
  as.matrix(EXPECTED_S5$reference$distributions[, -c(1, 2)])))
save_table(reference_qc, "S5j_reference_counts")
save_table(reference_distributions, "S5j_reference_distributions")
rm(reference_metrics)

# ---- S5k: retain the historical three-model FDR family
residual_spearman <- function(d, x, y, extra = NULL) {
  terms <- paste(c("ns(AGE,df=3)", extra), collapse = " + ")
  d <- d[complete.cases(d[, c(x, y, "AGE", extra)]), ]
  rx <- residuals(lm(as.formula(paste(x, "~", terms)), d))
  ry <- residuals(lm(as.formula(paste(y, "~", terms)), d))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  se <- 1/sqrt(nrow(d) - 6 - length(extra))
  data.frame(N = nrow(d), rho = rho, low = tanh(z - 1.96 * se), high = tanh(z + 1.96 * se), p = ct$p.value,
    se = se)
}
fidelity_age <- do.call(rbind, lapply(CANCERS15, function(ca) {
  d <- h_env$build_analysis_data(ca)
  d$AgeAccel <- residuals(lm(Horvath ~ AGE, d))
  save_table(d, paste0("fidelity_miRNA_cohort_", ca))
  save_table(h_env$build_mrc_fidelity_data(ca), paste0("fidelity_MRC_cohort_", ca))
  rbind(data.frame(Cancer = ca, Contrast = "Age acceleration ~ fidelity", residual_spearman(d, "AgeAccel",
    "FidelityLoss")), data.frame(Cancer = ca, Contrast = "miR-29 ~ fidelity | age", residual_spearman(d,
    "pan_miR29", "FidelityLoss")), data.frame(Cancer = ca, Contrast = "miR-29 ~ fidelity | age + acceleration",
    residual_spearman(d, "pan_miR29", "FidelityLoss", "AgeAccel")))
}))
fidelity_age$q <- p.adjust(fidelity_age$p, "BH")
fidelity_age_meta <- do.call(rbind, lapply(split(fidelity_age, fidelity_age$Contrast), function(d) {
  fit <- metafor::rma(yi = atanh(pmin(pmax(d$rho, -0.999999), 0.999999)), sei = d$se, method = "REML",
    test = "knha")
  data.frame(Contrast = d$Contrast[1], K = fit$k, rho = tanh(as.numeric(fit$b)), low = tanh(fit$ci.lb),
    high = tanh(fit$ci.ub), p = fit$pval, I2 = fit$I2)
}))
fidelity_age_meta$q <- p.adjust(fidelity_age_meta$p, "BH")
check("S5k all three pooled models unchanged", equal_numeric(as.matrix(fidelity_age_meta[, -1]), as.matrix(EXPECTED_S5$fidelity$age_meta[,
  -1])))
save_table(fidelity_age, "fidelity_age_all_45_tests")
save_table(fidelity_age_meta, "fidelity_age_all_three_meta_tests")
save_table(h_left, "Figure5h_miRNA_fidelity")
save_table(h_right, "Figure5h_MRC_fidelity")
save_table(hl_meta, "Figure5h_miRNA_pooled")
save_table(hr_meta, "Figure5h_MRC_pooled")
save_table(e_effects, "Figure5e_paired_effects")
save_table(e_meta, "Figure5e_pooled_effects")
save_table(ef_pairs, "Figure5ef_paired_patients")
save_table(c_185, "Figure5c_185_genes")
save_table(c_rank, "Figure5c_gene_ranking")
save_table(a_summary, "Figure5a_summary")
save_table(b_summary, "Figure5b_summary")
