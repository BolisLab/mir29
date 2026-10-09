# Execute from the project root. Original final PDFs and checkpoints are unchanged.
source("code/load_all.R")
source("code/export_panels.R")
export_all_panels()
message("Validated panels exported to figures/reconstructed/.")
