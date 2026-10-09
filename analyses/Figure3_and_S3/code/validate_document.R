# Post-render checks complement the numerical assertions run in load_all.R.
document <- readLines("Figure3_and_S3.Rmd", warn = FALSE)
markdown <- readLines("Figure3_and_S3.md", warn = FALSE)
stopifnot(sum(grepl("^## [a-m] —", document)) == 23L,
          !any(grepl("PLURYPROVA|questo per inhibitor|4 maggio", markdown)),
          all(file.exists(c("Figure3_and_S3.html", "validation/main_signature_statistics.csv",
                            "validation/final_signature_membership.csv"))))
stats <- read.csv("validation/main_signature_statistics.csv")
stopifnot(identical(stats$letter, LETTERS[2:7]),
  max(abs(stats$NES - c(-1.917397, 1.900731, -2.040113, 2.061773, -1.821168, 1.763397))) < 1e-6,
  all(stats$n_defined == c(354, 355, 149, 149, 125, 125)),
  all(stats$n_matched == c(311, 351, 149, 148, 124, 125)),
  all(stats$q_multilevel_global_4721 > 0))
report <- data.frame(check = c("13 main and 10 supplementary sections", "No rejected working comments",
  "Final B–G nomenclature", "Six accepted historical NES", "Final defined and matched sizes",
  "Nonzero multilevel FDR", "Completed HTML render"), status = "PASS")
write.csv(report, "validation/document_checks.csv", row.names = FALSE)
if (requireNamespace("digest", quietly = TRUE)) {
  files <- c("reference/Figure3_main_final.pdf", "reference/FigureS3_final.pdf", "Figure3_and_S3.Rmd", "Figure3_and_S3.html")
  hashes <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), character(1))
  stopifnot(hashes[1] == "1ed7fb706175b7bb41b3243cc69271223a943c09059130f8685b9e9afb91a7c8",
            hashes[2] == "4a070dccfbde1735d13c0ed2f905d421611d9a20321a5b5ba6d8566c505b663b")
  writeLines(paste(hashes, files, sep = "  "), "validation/final_sha256.txt")
}
print(report, row.names = FALSE)
