# Exact statistical recipes for optional Monte Carlo replay.
# The normal notebook render validates the saved results without rerunning them.
prepare_collection <- function(sets) {
  pairs <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE))
  pairs <- unique(pairs[!is.na(pairs$term) & !is.na(pairs$gene), ])
  split(as.character(pairs$gene), as.character(pairs$term))
}
replay_simple <- function(variant = "full") {
  stopifnot(as.character(packageVersion("fgsea")) == "1.24.0")
  pathways <- prepare_collection(c(definitions$background, definitions$run_definitions[[variant]]$signatures))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(10403L)
  result <- as.data.frame(fgsea::fgseaSimple(pathways, definitions$ranked_stats,
    nperm = 10000L, minSize = 20L, maxSize = Inf, gseaParam = 1, scoreType = "std",
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)))
  stopifnot(nrow(result) == 4721L)
  result$q_global_BH <- p.adjust(result$pval, "BH", n = 4721L)
  result
}
build_multilevel_union <- function() {
  ids <- definitions$map$signature
  background <- prepare_collection(definitions$background)
  background <- background[vapply(background, function(s) length(intersect(s, names(definitions$ranked_stats))) >= 20L, logical(1))]
  fixed <- definitions$fixed_HGF_GTEx_signatures[ids[3:6]]
  varying <- key <- list()
  for (variant in names(definitions$run_definitions)) {
    aliases <- paste0("LOTO_ML_", variant, "_", c("UP", "DN"))
    varying <- c(varying, setNames(definitions$run_definitions[[variant]]$signatures[ids[1:2]], aliases))
    key[[variant]] <- setNames(c(aliases, ids[3:6]), ids)
  }
  list(pathways = prepare_collection(c(background, fixed, varying)), key = key, background = names(background))
}
replay_multilevel <- function() {
  stopifnot(as.character(packageVersion("fgsea")) == "1.24.0", .Platform$OS.type == "unix")
  union <- build_multilevel_union()
  stopifnot(length(union$pathways) == 4735L)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(29030910L)
  fgsea::fgseaMultilevel(union$pathways, definitions$ranked_stats,
    sampleSize = 501L, nPermSimple = 10000L, eps = 0, minSize = 20L,
    maxSize = Inf, gseaParam = 1, scoreType = "std",
    BPPARAM = BiocParallel::MulticoreParam(workers = 4L, progressbar = FALSE))
}

# Recompute BH from the saved raw multilevel probabilities, per tested family.
union_specification <- build_multilevel_union()
stopifnot(identical(union_specification$key, multilevel$variant_key))
for (variant in names(multilevel$variants)) {
  z <- multilevel$variants[[variant]]
  stopifnot(nrow(z$all_results) == 4721L,
    max(abs(p.adjust(z$all_results$pval, "BH", n = 4721L) - z$all_results$q_global_BH)) < 1e-14,
    max(abs(p.adjust(z$six$pval, "BH", n = 6L) - z$six$q_six_BH)) < 1e-14)
  expected <- simple_checkpoint$runs[[variant]]$all_results
  matched <- match(z$all_results$pathway, expected$pathway)
  stopifnot(max(abs(z$all_results$ES - expected$ES[matched])) < 1e-12,
    max(abs(z$all_results$NES_simple - expected$NES[matched])) < 1e-12)
}
