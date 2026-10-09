#!/usr/bin/env Rscript
# Usage: Rscript run_hpf_dnam.R INPUT_DIRECTORY FRESH_OUTPUT_DIRECTORY
args_all <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args_all[grepl("^--file=",args_all)][1])
source(file.path(dirname(normalizePath(script)),"hpf_dnam_pipeline.R"))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: Rscript run_hpf_dnam.R INPUT_DIRECTORY FRESH_OUTPUT_DIRECTORY")
run_hpf_dnam(args[1],args[2])
