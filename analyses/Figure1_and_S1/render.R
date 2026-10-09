#!/usr/bin/env Rscript

# Render the complete Figure 1 notebook from the project directory.

options(stringsAsFactors = FALSE)

command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
script_path <- sub("^--file=", "", script_argument[1])
project_dir <- dirname(normalizePath(script_path))

setwd(project_dir)

rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc was not found. Open the project in RStudio or set RSTUDIO_PANDOC.")
}

rmarkdown::render(
  input = "Figure1.Rmd",
  output_file = "Figure1.html",
  envir = new.env(parent = globalenv()),
  clean = TRUE
)
