# Structural checks supplement, but do not replace, visual review.
required <- c("Figure4_and_S4.Rmd", "Figure4_and_S4.md", "Figure4_and_S4.html")
stopifnot(all(file.exists(required)), all(file.info(required)$size > 10000))
html <- paste(readLines("Figure4_and_S4.html", warn = FALSE), collapse = "\n")
md <- paste(readLines("Figure4_and_S4.md", warn = FALSE), collapse = "\n")
stopifnot(grepl("Supplementary Figure S4", html, fixed = TRUE),
          grepl("data:image/png;base64", html, fixed = TRUE),
          grepl("code-folding", html, fixed = TRUE),
          !grepl("PLURYPROVA|questo per inhibitor|Error in |Execution halted", md))
expected <- c("4a_PCA", "4b_temporal_DMP_counts", "4c_temporal_DMP_summary", "4d_global_methylation",
              "4e_beta_distributions", "4f_entropy_FAE", "4g_GTEx_context", "4h_compartment_heatmaps",
              "4i_left_compartment_burden", "4i_right_SkinBlood", "4j_left_repeat_heatmap", "4j_right_RetroAge",
              "4k_clock_component_AUC", "4l_CTRL_temporal_movers", "4m_temporal_slope_scatter",
              "4n_methylation_walk", "4o_center_and_right", "S4a_global_methylation_separate_controls",
              "S4b_entropy_FAE_deltas", "S4c", "S4d_causal_clocks", "S4e", "S4f_trajectory_outcomes",
              "S4g_CHIR_responsive_loci", "S4h_DKK2_aligned_locus", "S4i_DKK2_paired_responses",
              "S4i_DKK2_interaction_forest", "S4j_hTERT", "S4k_hTERT_transcriptomic_clocks",
              "S4m_HPF_retention", "S4n_HPF_residual_effect")
stopifnot(all(file.exists(file.path("figures/reconstructed", paste0(expected, ".pdf")))),
          all(file.exists(file.path("figures/reconstructed", paste0(expected, ".png")))))
checks <- read.delim("validation/numeric_checks.tsv")
stopifnot(all(trimws(checks$status) == "PASS"), nrow(checks) >= 20)
writeLines(c("PASS: Rmd, expanded Markdown and self-contained HTML exist.",
             "PASS: all 31 retained analytical panel/subpanel outputs exist.",
             "PASS: numerical validation checks passed.",
             "Editorial diagrams are separately identified and preserved from source PDFs.",
             paste("Validated:", format(Sys.time(), tz = "Europe/Rome"))),
           "validation/document_validation.txt")
