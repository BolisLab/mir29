# Pure selection functions. Inputs are external cohort-level probe audits.
# The current replay applies the adopted q <= 0.20 bound without re-searching
# thresholds. Bound-based stability is not an equivalence test.
stable_reference_sites <- function(audit) {
  keep <- with(audit, detection_fraction >= .90 & site_abs <= .05 &
                 age_abs <= .001 & sex_abs <= .05)
  audit$cg[which(keep)]
}
culture_stable_sites <- function(audit) {
  keep <- with(audit, detection_fraction >= .90 &
    is.finite(delta_beta_full_range_unadjusted) &
    is.finite(delta_beta_full_range_batch_adjusted) &
    abs(delta_beta_full_range_unadjusted) <= .05 &
    abs(delta_beta_full_range_batch_adjusted) <= .05)
  audit$probe[which(keep)]
}
fibroblast_associated_sites <- function(audit, q_threshold = .20) {
  keep <- with(audit, q_gender_adjusted <= q_threshold & detection_gse31848 >= .90 &
    sign(delta_fibroblast_vs_nonfibro_unadjusted) ==
      sign(delta_fibroblast_vs_nonfibro_gender_adjusted) &
    abs(delta_fibroblast_vs_nonfibro_unadjusted) >= .10 &
    abs(delta_fibroblast_vs_nonfibro_gender_adjusted) >= .10 &
    abs(delta_fibroblast_vs_immune_external) >= .10 &
    sign(delta_fibroblast_vs_nonfibro_gender_adjusted) ==
      sign(delta_fibroblast_vs_immune_external))
  audit$probe[which(keep)]
}
# The discovery effect gives equal weight to epithelial and SMC/perivascular
# comparators, not to their unequal sample counts. Repeated arrays are averaged
# by cell line before this function. Both methylation directions are eligible.
fit_identity_contrast <- function(beta, metadata, adjust_gender = TRUE) {
  metadata$lineage_group <- factor(metadata$lineage_group,
    levels = c("Epithelial", "SMC_perivascular", "Fibroblast"))
  design <- model.matrix(if (adjust_gender) ~ 0 + lineage_group + gender else
                          ~ 0 + lineage_group, metadata)
  contrast <- limma::makeContrasts(
    lineage_groupFibroblast - .5 * lineage_groupEpithelial - .5 * lineage_groupSMC_perivascular,
    levels = design)
  limma::eBayes(limma::contrasts.fit(limma::lmFit(beta, design), contrast), robust = TRUE)
}
