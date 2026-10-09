#!/usr/bin/env Rscript

# Render the portable Figure 2 main analytical notebook.

options(stringsAsFactors = FALSE)

command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
if (length(script_argument)) {
  project_dir <- dirname(normalizePath(sub("^--file=", "", script_argument[1])))
  setwd(project_dir)
}

rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc was not found. Open the project in RStudio or set RSTUDIO_PANDOC.")
}

rmarkdown::render(
  input = "Figure2_main.Rmd",
  output_file = "Figure2_main.html",
  envir = new.env(parent = globalenv()),
  clean = TRUE
)
