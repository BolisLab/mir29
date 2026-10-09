# Optional upstream recipes. Rendering starts from the frozen predictions.
# These functions are not called merely by sourcing this file.

predict_rnaagecalc <- function(counts, metadata, gene_annotation) {
  genes <- intersect(rownames(counts), rownames(gene_annotation)[
    is.finite(gene_annotation$length) & gene_annotation$length > 0])
  fpkm <- RNAAgeCalc::count2FPKM(counts[genes, , drop = FALSE],
    genelength = gene_annotation[genes, "length"], idtype = "ENSEMBL")
  design <- model.matrix(~ 0 + reGROUP + TIME, metadata)
  adjusted <- limma::removeBatchEffect(as.matrix(fpkm), batch = metadata$EXPERIMENT, design = design)
  adjusted <- limma::removeBatchEffect(adjusted, batch = metadata$BATCH, design = design)
  adjusted[adjusted < 0] <- 0
  RNAAgeCalc::predict_age(adjusted, exprtype = "FPKM", stype = "caucasian",
                         idtype = "ENSEMBL", signature = NULL, tissue = "all")
}
predict_pasta <- function(counts) {
  # Official model-gene filtering and rank normalization are applied ONCE.
  ranked <- pasta::filtering_age_model_genes_and_rank_norm(counts)
  pasta::predicting_age_score(t(ranked), model_type = "Pasta")
}
prepare_tage <- function(counts, metadata) {
  stopifnot(as.character(packageVersion("tAge")) == "1.1.0")
  eset <- tAge::make_ExpressionSet(counts, metadata, verbose = FALSE)
  processed <- tAge::tAge_preprocessing(eset, species = "human", gene_mapping_type = "Ensembl",
    control_group_column = "is_control", control_group_label = "CTRL", verbose = FALSE)
  centered <- processed$yugene
  for (arm in c("INH", "MIM")) {
    ix <- which(metadata$arm == arm)
    current <- tAge::control_subtraction(processed$yugene[, ix], "is_control", "CTRL", verbose = FALSE)
    Biobase::exprs(centered)[, ix] <- Biobase::exprs(current)
  }
  # Write this matrix and run code/upstream/replay_htert_tage.py in the locked
  # historical Python environment. Display centering/scaling is separate.
  Biobase::exprs(centered)
}
