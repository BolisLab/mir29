root <- normalizePath(getwd(), mustWork = TRUE)
exclude <- c(
  "Figure3_and_S3.html", "bundle_manifest.tsv",
  "logs/render.log", "logs/manifest_sessionInfo.txt"
)
files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
files <- files[file.info(files)$isdir %in% FALSE]
relative <- substring(files, nchar(root) + 2L)
keep <- !relative %in% exclude & !grepl("^cache/", relative)
files <- files[keep]
relative <- relative[keep]
info <- file.info(files)
manifest <- data.frame(
  relative_path = relative,
  bytes = unname(info$size),
  md5 = unname(tools::md5sum(files)),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$relative_path), ]
utils::write.table(
  manifest,
  file = file.path(root, "bundle_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
dir.create(file.path(root, "logs"), showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(root, "logs", "manifest_sessionInfo.txt"))
message("Wrote bundle_manifest.tsv with ", nrow(manifest), " files.")

