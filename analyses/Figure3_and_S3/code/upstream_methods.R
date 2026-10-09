# Optional bulk-expression normalization from the bundled counts.
# The standard render uses the verified matrices to avoid repeating DESeq2 fits.
rebuild_pca_vst <- function(input) {
  metadata <- input$metadata[colnames(input$counts), , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(input$counts)), metadata, design = ~1)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  vst <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  protected <- model.matrix(~0 + reGROUP + TIME, metadata)
  adjusted <- limma::removeBatchEffect(vst, batch = metadata$EXPERIMENT, design = protected)
  adjusted <- limma::removeBatchEffect(adjusted, batch = metadata$BATCH, design = protected)
  list(vst = vst, corrected = adjusted)
}

# Historical per-arm DE model; no count-level subtraction is performed.
rebuild_arm_de <- function(input, arm) {
  metadata <- input$metadata[input$metadata$GROUP %in% c(arm, if (arm == "INH") "inhCTRL" else "mimCTRL"), ]
  metadata$SHORT_GROUP <- factor(metadata$SHORT_GROUP)
  metadata$BATCH <- factor(metadata$BATCH)
  counts <- round(as.matrix(input$counts[, rownames(metadata), drop = FALSE]))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, metadata, design = ~0 + SHORT_GROUP + BATCH)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  lapply(1:6, function(time) {
    treated <- unique(metadata$SHORT_GROUP[metadata$TIME == time & metadata$GROUP == arm])
    control <- unique(metadata$SHORT_GROUP[metadata$TIME == time & metadata$GROUP != arm])
    stopifnot(length(treated) == 1L, length(control) == 1L)
    as.data.frame(DESeq2::results(dds, contrast = c("SHORT_GROUP", as.character(treated), as.character(control)), alpha = .001))
  })
}

# Overall shift versus own-arm controls; repeated experimental runs share an intercept.
fit_global_shift <- function(samples) {
  samples$group <- factor(samples$group, levels = c("control", "treated"))
  samples$time <- factor(samples$time)
  null <- lme4::lmer(score ~ time + (1 | batch), samples, REML = FALSE)
  shift <- lme4::lmer(score ~ group + time + (1 | batch), samples, REML = FALSE)
  comparison <- anova(null, shift)
  list(p = comparison$`Pr(>Chisq)`[2], estimate = lme4::fixef(shift)["grouptreated"],
       singular = lme4::isSingular(null) || lme4::isSingular(shift))
}

# Timepoint comparisons use the complete treatment-by-time model.
fit_timepoint_tests <- function(samples, df_method = "kenward-roger") {
  samples$group <- factor(samples$group, levels = c("control", "treated"))
  samples$time <- factor(samples$time)
  model <- lme4::lmer(score ~ group * time + (1 | batch), samples, REML = FALSE)
  means <- emmeans::emmeans(model, ~group | time, lmer.df = df_method)
  tests <- as.data.frame(emmeans::contrast(means, list("treated - control" = c(-1, 1)), adjust = "none"))
  tests$p_holm6 <- p.adjust(tests$p.value, "holm", n = 6)
  tests
}
