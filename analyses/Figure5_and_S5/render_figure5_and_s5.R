# Run from the bundle root, or open the RStudio project and source this file.
stopifnot(file.exists("Figure5_and_S5.Rmd"))
if (!requireNamespace("rmarkdown",quietly=TRUE))
  stop("Install rmarkdown before rendering.")
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc))
  Sys.setenv(RSTUDIO_PANDOC=rstudio_pandoc)
if (!rmarkdown::pandoc_available())
  stop("Open this project in RStudio or install Pandoc.")
rmarkdown::render("Figure5_and_S5.Rmd",output_file="Figure5_and_S5.html",
  envir=new.env(parent=globalenv()),clean=TRUE)
source("code/validate_document.R")
source("code/write_bundle_manifest.R")
