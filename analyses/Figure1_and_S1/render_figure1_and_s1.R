#!/usr/bin/env Rscript

# Render the complete Figure 1 + Supplementary Figure S1 notebook.

options(stringsAsFactors = FALSE)

command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
script_path <- sub("^--file=", "", script_argument[1])
project_dir <- dirname(normalizePath(script_path))
setwd(project_dir)

# Cairo does not always discover macOS fonts correctly. Keep the notebook's
# generic sans family, but provide a glyph-complete fallback for PDF exports.
# Respect an explicitly configured Fontconfig environment on other systems.
if (identical(Sys.info()[["sysname"]], "Darwin") &&
    !nzchar(Sys.getenv("FONTCONFIG_FILE"))) {
  font_config <- tempfile("figure1-fonts-", fileext = ".conf")
  font_cache <- file.path(tempdir(), "figure1-font-cache")
  dir.create(font_cache, showWarnings = FALSE, recursive = TRUE)
  writeLines(c(
    '<?xml version="1.0"?>',
    '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">',
    '<fontconfig>',
    '  <dir>/System/Library/Fonts</dir>',
    '  <dir>/Library/Fonts</dir>',
    sprintf('  <cachedir>%s</cachedir>', font_cache),
    '  <alias><family>sans-serif</family><prefer><family>Arial</family></prefer></alias>',
    '  <alias><family>sans</family><prefer><family>Arial</family></prefer></alias>',
    '</fontconfig>'
  ), font_config)
  Sys.setenv(FONTCONFIG_FILE = font_config)
}

rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc was not found. Open the project in RStudio or set RSTUDIO_PANDOC.")
}

rmarkdown::render(
  input = "Figure1_and_S1.Rmd",
  output_file = "Figure1_and_S1.html",
  envir = new.env(parent = globalenv()),
  clean = TRUE
)
