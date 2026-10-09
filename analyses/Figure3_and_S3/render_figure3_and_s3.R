# Open the RStudio project, then source this file or press Knit.
stopifnot(file.exists("Figure3_and_S3.Rmd"))
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install the rmarkdown package.")
}
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}
if (!rmarkdown::pandoc_available()) {
  stop("Open this project in RStudio or install Pandoc.")
}
rmarkdown::render(
  "Figure3_and_S3.Rmd",
  output_file = "Figure3_and_S3.html",
  envir = new.env(parent = globalenv()),
  clean = TRUE
)
source("code/validate_document.R")
