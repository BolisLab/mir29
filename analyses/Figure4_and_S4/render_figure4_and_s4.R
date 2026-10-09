# Open the RStudio project, then source this file or press Knit.
stopifnot(file.exists("Figure4_and_S4.Rmd"))
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install rmarkdown before rendering.")
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc))
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
if (!rmarkdown::pandoc_available()) stop("Open the project in RStudio or install Pandoc.")
rmarkdown::render("Figure4_and_S4.Rmd", output_file = "Figure4_and_S4.html",
                  envir = new.env(parent = globalenv()), clean = TRUE)
source("code/validate_document.R")
