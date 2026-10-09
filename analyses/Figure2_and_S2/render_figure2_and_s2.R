# Run in RStudio with the project folder as working directory.
# No commandArgs() or terminal launch is required.
stopifnot(file.exists("Figure2_and_S2.Rmd"))
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install the rmarkdown package.")
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
if (!rmarkdown::pandoc_available() && dir.exists(rstudio_pandoc))
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
if (!rmarkdown::pandoc_available()) stop("Open this project in RStudio (or install Pandoc).")
rmarkdown::render("Figure2_and_S2.Rmd", output_file = "Figure2_and_S2.html",
                  envir = new.env(parent = globalenv()), clean = TRUE)
