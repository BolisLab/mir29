# Run from the bundle project root. No HTML generation or upstream downloads.
stopifnot(file.exists("code/load_all.R"))
source("code/load_all.R")
message("Panel PDFs/PNGs: figures/reconstructed/; statistics: tables/; checks: validation/.")
