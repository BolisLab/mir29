# All paths and packages required by the processed-data reconstruction.
BUNDLE <- normalizePath(".")
DATA_ROOT <- file.path(BUNDLE, "data/source")
ROOT <- BUNDLE
PDF_DIR <- file.path(BUNDLE, "figures/reconstructed")
for (p in c(PDF_DIR, file.path(BUNDLE, c("data", "validation", "previews", "provenance", "code/helpers")))) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
packages <- c("ggplot2", "ggrepel", "cowplot", "metafor", "splines", "jsonlite")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({
  library(ggplot2)
  library(splines)
})
CANCERS16 <- c("BLCA", "BRCA", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD",
  "LUSC", "PRAD", "STAD", "THCA", "UCEC")
CANCERS15 <- setdiff(CANCERS16, "STAD")
MRC_GENES <- c("DNMT3A", "DNMT3B", "TET1", "TET3", "TDG")
EXPECTED_H_LEFT_N <- c(BLCA = 409, BRCA = 773, CHOL = 36, COAD = 293, ESCA = 184, HNSC = 522, KICH = 66,
  KIRC = 318, KIRP = 272, LIHC = 371, LUAD = 437, LUSC = 357, PRAD = 494, THCA = 506, UCEC = 427)
EXPECTED_H_RIGHT_N <- c(BLCA = 406, BRCA = 781, CHOL = 35, COAD = 292, ESCA = 162, HNSC = 501, KICH = 65,
  KIRC = 318, KIRP = 271, LIHC = 370, LUAD = 436, LUSC = 361, PRAD = 497, THCA = 504, UCEC = 427)
INPUT_FILES <- character()
read_data <- function(relative) {
  path <- file.path(DATA_ROOT, relative)
  stopifnot(file.exists(path))
  INPUT_FILES <<- unique(c(INPUT_FILES, path))
  readRDS(path)
}
CHECKS <- list()
check <- function(label, condition, detail = "") {
  ok <- isTRUE(condition)
  CHECKS[[length(CHECKS) + 1L]] <<- list(check = label, pass = ok, detail = as.character(detail))
  message(if (ok)
    "PASS " else "FAIL ", label, if (nzchar(detail))
    paste0(": ", detail))
  if (!ok)
    stop("Validation failed: ", label)
  invisible(ok)
}
equal_numeric <- function(a, b, tolerance = 1e-08) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance, check.attributes = FALSE))
}
fisher_mean <- function(r) {
  r <- r[is.finite(r)]
  if (!length(r))
    return(NA_real_)
  tanh(mean(atanh(pmin(pmax(r, -0.999999), 0.999999))))
}


cowplot::set_null_device("agg")
load_helpers <- function(required, envir, output_name) {
  source(file.path(BUNDLE, "code/helpers", output_name), local = envir)
  stopifnot(all(vapply(required, exists, logical(1), envir = envir, inherits = FALSE)))
}
save_table <- function(x, name) {
  write.table(as.data.frame(x), file.path("tables", paste0(name, ".tsv")), sep = "\t", row.names = FALSE,
    quote = FALSE, na = "NA")
}
