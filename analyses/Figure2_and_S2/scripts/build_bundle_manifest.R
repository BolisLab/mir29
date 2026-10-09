# Rebuild a checksum inventory of the portable Figure 2 bundle.

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "Figure2_main.Rmd"))) {
  stop("Set the working directory to the Figure2_main project folder.")
}

all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
all_files <- all_files[!grepl("(^|/)(cache|\\.Rproj\\.user)(/|$)", all_files)]
all_files <- all_files[basename(all_files) != "bundle_manifest.tsv"]

manifest <- data.frame(
  relative_path = substring(all_files, nchar(root) + 2L),
  bytes = unname(file.info(all_files)$size),
  md5 = unname(tools::md5sum(all_files)),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$relative_path), ]

write.table(
  manifest,
  file.path(root, "bundle_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote bundle_manifest.tsv with ", nrow(manifest), " files.")
