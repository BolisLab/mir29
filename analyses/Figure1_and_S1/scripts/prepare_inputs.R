#!/usr/bin/env Rscript

# Prepare the compact, self-contained inputs used by Figure1.Rmd.
#
# The script reads the unstranded-count column from the original STAR
# ReadsPerGene.out.tab files, applies only the sample exclusions used in the
# manuscript analysis, and stores compact integer matrices inside this project.
# The generated and copied files in data/ are the only inputs required to knit
# Figure1.Rmd. Sample-level RNA-seq inputs are rebuilt from STAR files; the
# cross-tissue aging, clock, and mature-miRNA objects are copied verbatim because
# they are the analysis-ready outputs distributed with the manuscript workflow.

options(stringsAsFactors = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required to prepare the inputs.")
}

command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
script_path <- sub("^--file=", "", script_argument[1])
script_dir <- dirname(normalizePath(script_path))
project_dir <- dirname(script_dir)
data_dir <- file.path(project_dir, "data")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Source locations can be overridden without editing this script.
ipsc_dir <- path.expand(Sys.getenv(
  "MIR29_IPSC_SOURCE",
  unset = "external_inputs/IPSC/GSE165176"
))
zga_dir <- path.expand(Sys.getenv(
  "MIR29_ZGA_SOURCE",
  unset = "external_inputs/ZGA"
))
aging_dir <- path.expand(Sys.getenv(
  "MIR29_AGING_SOURCE",
  unset = "external_inputs"
))

if (!dir.exists(ipsc_dir)) stop("iPSC source directory not found: ", ipsc_dir)
if (!dir.exists(zga_dir)) stop("ZGA source directory not found: ", zga_dir)
if (!dir.exists(aging_dir)) stop("Aging source directory not found: ", aging_dir)


# Read selected STAR files -----------------------------------------------------

read_star_directory <- function(directory, sample_ids) {
  star_files <- sort(list.files(
    directory,
    pattern = "ReadsPerGene\\.out\\.tab$",
    full.names = TRUE
  ))

  file_ids <- sub("\\.ReadsPerGene\\.out\\.tab$", "", basename(star_files))
  names(star_files) <- file_ids

  missing_ids <- setdiff(sample_ids, file_ids)
  if (length(missing_ids) > 0) {
    stop(
      "Missing STAR files in ", directory, ": ",
      paste(missing_ids, collapse = ", ")
    )
  }

  selected_files <- star_files[sample_ids]

  first_table <- data.table::fread(
    selected_files[1],
    skip = 4,
    header = FALSE,
    select = c(1, 2),
    showProgress = FALSE
  )

  gene_ids_raw <- first_table[[1]]
  gene_ids <- sub("\\..*$", "", gene_ids_raw)
  keep_genes <- grepl("^ENSG", gene_ids)

  count_matrix <- matrix(
    0L,
    nrow = sum(keep_genes),
    ncol = length(selected_files),
    dimnames = list(gene_ids[keep_genes], sample_ids)
  )

  count_matrix[, 1] <- as.integer(first_table[[2]][keep_genes])

  if (length(selected_files) > 1) {
    for (column_index in 2:length(selected_files)) {
      current_table <- data.table::fread(
        selected_files[column_index],
        skip = 4,
        header = FALSE,
        select = c(1, 2),
        showProgress = FALSE
      )

      current_gene_ids <- current_table[[1]]
      if (!identical(current_gene_ids, gene_ids_raw)) {
        stop("Gene order differs in: ", selected_files[column_index])
      }

      count_matrix[, column_index] <- as.integer(current_table[[2]][keep_genes])
    }
  }

  count_matrix
}


# iPSC reprogramming -----------------------------------------------------------

ipsc_annotation_file <- file.path(ipsc_dir, "ipsc_annot.txt")
ipsc_methylation_file <- file.path(ipsc_dir, "methylation.tsv")

ipsc_annotation <- read.delim(
  ipsc_annotation_file,
  row.names = 1,
  check.names = FALSE
)

ipsc_keep <- ipsc_annotation$TP %in% c(11, 13, 15) &
  grepl("Exp1", ipsc_annotation$EXPERIMENT)

ipsc_annotation <- ipsc_annotation[ipsc_keep, , drop = FALSE]
ipsc_counts <- read_star_directory(ipsc_dir, rownames(ipsc_annotation))

ipsc_methylation <- read.delim(
  ipsc_methylation_file,
  row.names = 1,
  check.names = FALSE
)

stopifnot(ncol(ipsc_counts) == 36L)
stopifnot(identical(colnames(ipsc_counts), rownames(ipsc_annotation)))
stopifnot(!anyNA(ipsc_counts))

ipsc_inputs <- list(
  counts = ipsc_counts,
  annotation = ipsc_annotation,
  methylation = ipsc_methylation,
  provenance = list(
    accession = "GSE165176",
    count_type = "STAR unstranded gene counts",
    sample_filter = "TP in {11, 13, 15}; Exp1; all six donors",
    prepared = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
)

saveRDS(
  ipsc_inputs,
  file.path(data_dir, "figure1_ipsc_inputs.rds"),
  compress = "xz"
)


# Zygotic development ----------------------------------------------------------

zga_annotation_file <- file.path(zga_dir, "zga_annot.txt")
zga_annotation <- read.delim(
  zga_annotation_file,
  row.names = 1,
  check.names = FALSE
)

zga_accessions <- c(
  "GSE101571",
  "GSE157834",
  "GSE71318",
  "GSE65481",
  "GSE133854"
)

zga_annotation <- zga_annotation[
  zga_annotation$DUPLICATED != 1 & zga_annotation$GEO %in% zga_accessions,
  ,
  drop = FALSE
]

zga_counts_by_accession <- vector("list", length(zga_accessions))
names(zga_counts_by_accession) <- zga_accessions

for (accession in zga_accessions) {
  accession_samples <- rownames(zga_annotation)[zga_annotation$GEO == accession]
  accession_dir <- file.path(zga_dir, accession)

  zga_counts_by_accession[[accession]] <- read_star_directory(
    accession_dir,
    accession_samples
  )
}

common_genes <- Reduce(intersect, lapply(zga_counts_by_accession, rownames))

zga_counts_by_accession <- lapply(
  zga_counts_by_accession,
  function(counts) counts[common_genes, , drop = FALSE]
)

zga_counts <- do.call(cbind, zga_counts_by_accession)
zga_counts <- zga_counts[, rownames(zga_annotation), drop = FALSE]

stopifnot(ncol(zga_counts) == 113L)
stopifnot(identical(colnames(zga_counts), rownames(zga_annotation)))
stopifnot(!anyNA(zga_counts))

zga_inputs <- list(
  counts = zga_counts,
  annotation = zga_annotation,
  provenance = list(
    accessions = zga_accessions,
    count_type = "STAR unstranded gene counts",
    sample_filter = "DUPLICATED != 1",
    samples_per_accession = table(zga_annotation$GEO),
    prepared = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
)

saveRDS(
  zga_inputs,
  file.path(data_dir, "figure1_zga_inputs.rds"),
  compress = "xz"
)


# Human-readable copies of the metadata ---------------------------------------

data.table::fwrite(
  data.table::as.data.table(ipsc_annotation, keep.rownames = "sample_id"),
  file.path(data_dir, "figure1_ipsc_samples.tsv"),
  sep = "\t"
)

data.table::fwrite(
  data.table::as.data.table(ipsc_methylation, keep.rownames = "sample_id"),
  file.path(data_dir, "figure1_ipsc_methylation.tsv"),
  sep = "\t"
)

data.table::fwrite(
  data.table::as.data.table(zga_annotation, keep.rownames = "sample_id"),
  file.path(data_dir, "figure1_zga_samples.tsv"),
  sep = "\t"
)


# Cross-tissue aging, clocks, and mature miRNAs -------------------------------

# These objects contain precomputed cross-tissue coefficients or compact source
# tables. In particular, coefficients.age.RDS is used downstream of the
# normalization/adjustment and age-constrained Monte Carlo procedure summarized
# in panel S1c; the sample-level controlled-access GTEx ages are not redistributed.
aging_files <- c(
  "hs_gene_annotations.RDS",
  "coefficients.age.RDS",
  "coefficients.age.sub7.RDS",
  "coefficients.clock.RDS",
  "epiclock_preds_GTEX.RDS",
  "clock_labels_classes.csv",
  "coefficients.mirnanew.RDS",
  "41467_2024_52943_MOESM6_ESM.xlsx",
  "GSE62776_all_smallRNA_reprogramming.counts_tracking.txt"
)

aging_sources <- file.path(aging_dir, aging_files)
missing_aging_files <- aging_sources[!file.exists(aging_sources)]

if (length(missing_aging_files) > 0) {
  stop(
    "Missing aging source files:\n",
    paste(missing_aging_files, collapse = "\n")
  )
}

copied <- file.copy(
  from = aging_sources,
  to = file.path(data_dir, aging_files),
  overwrite = TRUE,
  copy.date = TRUE
)

if (!all(copied)) {
  stop("One or more aging inputs could not be copied into: ", data_dir)
}

cat("Prepared Figure 1 inputs in:", data_dir, "\n")
cat("iPSC:", nrow(ipsc_counts), "genes x", ncol(ipsc_counts), "samples\n")
cat("ZGA: ", nrow(zga_counts), "genes x", ncol(zga_counts), "samples\n")
cat("Aging/clock/miRNA inputs copied:", length(aging_files), "files\n")
