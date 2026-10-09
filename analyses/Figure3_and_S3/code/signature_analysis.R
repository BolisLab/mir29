# Reconstruct the accepted partial correlations from pre-sex-correction VST.
definitions <- read_input("definitions")
curve_checkpoint <- read_input("curves")
simple_checkpoint <- read_input("simple_results")
authentic <- read_input("authentic")
fibroblast_vst <- read_input("unadjusted_fibroblast_vst")

partial_pearson <- function(expression, metadata) {
  stopifnot(identical(colnames(expression), rownames(metadata)))
  metadata$SEX <- factor(metadata$SEX)
  design <- model.matrix(~ SEX, metadata)
  residual_expression <- t(qr.resid(qr(design), t(expression)))
  residual_age <- qr.resid(qr(design), metadata$AGE)
  coefficient <- suppressWarnings(drop(cor(t(residual_expression), residual_age)))
  centered_norm <- function(m) sqrt(rowSums((m - rowMeans(m))^2))
  original_norm <- centered_norm(expression)
  informative <- original_norm > 1e-12 * pmax(1, sqrt(rowSums(expression^2))) &
    centered_norm(residual_expression) > 1e-12 * pmax(1, original_norm)
  coefficient[!informative] <- NA_real_
  setNames(coefficient, rownames(expression))
}
hf_ids <- rownames(definitions$partial_correlations)
partial_coefficients <- vapply(definitions$tissues, function(tissue) {
  z <- fibroblast_vst[[tissue]]
  partial_pearson(z$expression[hf_ids, , drop = FALSE], z$metadata)
}, numeric(length(hf_ids)))
rownames(partial_coefficients) <- hf_ids
stopifnot(identical(is.finite(partial_coefficients), is.finite(definitions$partial_correlations)),
          max(abs(partial_coefficients - definitions$partial_correlations), na.rm = TRUE) < 1e-10)

select_tails <- function(coefficient, symbols) {
  thresholds <- quantile(coefficient, c(.01, .99), na.rm = TRUE, names = FALSE)
  eligible <- is.finite(coefficient) & !symbols %in% c("LINC-PINT", "MIR29B2CHG")
  list(UP = eligible & coefficient > thresholds[2] & coefficient > 0,
       DN = eligible & coefficient < thresholds[1] & coefficient < 0,
       thresholds = thresholds)
}
aggregate_hf <- function(tissues) {
  r <- partial_coefficients[, tissues, drop = FALSE]
  n_valid <- rowSums(is.finite(r))
  r[!is.finite(r)] <- NA_real_
  score <- tanh(rowMeans(atanh(pmax(pmin(r, 1 - 1e-6), -1 + 1e-6)), na.rm = TRUE))
  score[n_valid < 5 | !is.finite(score)] <- NA_real_
  symbols <- definitions$run_definitions$full$table$symbol
  tails <- select_tails(score, symbols)
  data.frame(gene = hf_ids, symbol = symbols, r = score,
             UP = tails$UP, DN = tails$DN, n_valid = n_valid)
}
hf_curve <- aggregate_hf(definitions$tissues)
hgf_curve <- curve_checkpoint$source_curves$HGF
gingiva <- fibroblast_vst$Gingiva
hgf_recomputed <- partial_pearson(gingiva$expression[hgf_curve$gene, , drop = FALSE], gingiva$metadata)
stopifnot(max(abs(hgf_recomputed - hgf_curve$r), na.rm = TRUE) < 1e-10)
hgf_tails <- select_tails(hgf_recomputed, hgf_curve$symbol)
stopifnot(identical(unname(hgf_tails$UP), hgf_curve$UP), identical(unname(hgf_tails$DN), hgf_curve$DN))
signature_sets <- definitions$run_definitions$full$signatures
for (direction in c("UP", "DN")) {
  id <- if (direction == "UP") "new_panfib_up_R1" else "new_panfib_dn_R1"
  stopifnot(setequal(unique(hf_curve$symbol[hf_curve[[direction]]]), signature_sets[[id]]))
}
stopifnot(identical(signature_sets, support$signature_sets))

# Preserve the exact historical Monte Carlo NES, and use the accepted
# multilevel probabilities with a separately declared BH family.
main_statistics <- support$main_statistics
main_statistics$letter <- final_letter(main_statistics$signature)
main_statistics$label <- final_label(main_statistics$signature)
main_statistics <- main_statistics[match(signature_key$id, main_statistics$signature), ]
stopifnot(max(abs(main_statistics$NES - c(-1.917397, 1.900731, -2.040113, 2.061773, -1.821168, 1.763397))) < 1e-6,
          identical(as.integer(main_statistics$n_defined), c(354L, 355L, 149L, 149L, 125L, 125L)))
stopifnot(identical(as.integer(main_statistics$n_matched), c(311L, 351L, 149L, 148L, 124L, 125L)))
full_simple <- simple_checkpoint$runs$full$all_results
full_ml <- multilevel$variants$full
stopifnot(max(abs(main_statistics$NES - simple_checkpoint$runs$full$six$NES[
  match(main_statistics$signature, simple_checkpoint$runs$full$six$signature)])) < 1e-12)

# Verify every leave-one-tissue-out signature before plotting its cached GSEA.
for (tissue in definitions$tissues) {
  variant <- paste0("omit_", tissue)
  recomputed <- aggregate_hf(setdiff(definitions$tissues, tissue))
  expected <- definitions$run_definitions[[variant]]
  stopifnot(max(abs(recomputed$r - expected$score), na.rm = TRUE) < 1e-12,
            identical(unname(recomputed$UP), unname(expected$table$UP)),
            identical(unname(recomputed$DN), unname(expected$table$DN)))
}
loto_statistics <- multilevel$six_summary
loto_statistics$letter <- final_letter(loto_statistics$signature)
loto_statistics$label <- final_label(loto_statistics$signature)
write.csv(main_statistics, "validation/main_signature_statistics.csv", row.names = FALSE)
write.csv(loto_statistics, "validation/leave_one_tissue_out_statistics.csv", row.names = FALSE)
membership <- do.call(rbind, lapply(signature_key$id, function(id)
  data.frame(letter = final_letter(id), signature = final_label(id), gene = signature_sets[[id]])))
write.csv(membership, "validation/final_signature_membership.csv", row.names = FALSE)
