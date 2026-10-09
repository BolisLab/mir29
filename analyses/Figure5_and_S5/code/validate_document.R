# Structural validation complements numerical checks and visual review.
required <- c("Figure5_and_S5.Rmd", "Figure5_and_S5.md", "Figure5_and_S5.html")
stopifnot(all(file.exists(required)), all(file.info(required)$size > 10000))
html <- paste(readLines(required[3], warn = FALSE), collapse = "\n")
md <- paste(readLines(required[2], warn = FALSE), collapse = "\n")
stopifnot(grepl("Supplementary Figure S5", html, fixed = TRUE), grepl("data:image/png;base64", html,
  fixed = TRUE), grepl("code-folding", html, fixed = TRUE), !grepl("PLURYPROVA|questo per inhibitor|Error in |Execution halted",
  md))
expected <- c("Figure5a_TCGA_miRNA_age", "Figure5b_CCLE_miRNA_age", "Figure5c_RNA_age_ranking", "Figure5e_paired_expression",
  "Figure5f_paired_scatter_Pearson", "Figure5h_fidelity_associations", "S5a_TCGA_cohort", "S5b_TCGA_age_and_acceleration",
  "S5c_TCGA_leave_one_out", "S5d_CCLE_cohort", "S5e_CCLE_leave_one_out", "S5f_CCLE_four_ranks", "S5g_TCGA_RNA_distributions",
  "S5h_paired_cohort", "S5i_paired_leave_one_out", "S5j_reference_QC", "S5k_acceleration_fidelity")
for (ext in c(".pdf", ".png")) stopifnot(all(file.exists(file.path("figures/reconstructed", paste0(expected,
  ext)))))
for (name in expected) stopifnot(grepl(paste0(name, ".png"), md, fixed = TRUE))
current_main <- readRDS("data/Figure5_panel_data.rds")
expected_main <- readRDS("data/expected_main.rds")
checkpoint_checks <- vapply(names(current_main), function(panel) {
  isTRUE(all.equal(current_main[[panel]], expected_main[[panel]],
    tolerance=1e-8, check.attributes=FALSE))
}, logical(1))
stopifnot(all(checkpoint_checks))
writeLines(paste("PASS: complete main checkpoint",names(checkpoint_checks)),
  "validation/main_checkpoint_comparison.txt")
checks <- jsonlite::read_json("validation/numerical_checks.json", simplifyVector = TRUE)
stopifnot(all(checks$pass), nrow(checks) >= 139L)
writeLines(c("PASS: Rmd, expanded Markdown and HTML exist.", "PASS: embedded figure images and code-folding controls.",
  "PASS: all 17 analytical exports and their document images.", paste("PASS:", nrow(checks), "numerical checks."),
  "Two manual artwork panels are explicitly identified and retained.", paste("Validated:", format(Sys.time(),
    tz = "Europe/Rome"))), "validation/document_validation.txt")
