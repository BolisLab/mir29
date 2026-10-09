# Current curve/endpoint labels, applied after all statistical replay modules.
# This changes reconstruction annotations only: no scores, fitted curves,
# comparison families, reference-distance analyses or Figure 4l are changed.
# Clock curves display global q values, never individual-timepoint stars.

current_labels <- read.delim("CURRENT_STATISTICAL_LABELS.tsv", check.names = FALSE,
                             stringsAsFactors = FALSE)
main_stats <- read.delim("tables/main_statistics_own_vs_pooled_long.tsv",
                         check.names = FALSE, stringsAsFactors = FALSE)
epigenetic_stats <- read.delim("tables/S4_epigenetic_clock_global_statistics.tsv",
                               check.names = FALSE, stringsAsFactors = FALSE)
htert_stats <- read.delim("tables/S4k_clock_global_statistics.tsv",
                         check.names = FALSE, stringsAsFactors = FALSE)
primary_reference_stats <- read.delim("tables/4o_123CpG_primary_tests.tsv",
                                      check.names = FALSE, stringsAsFactors = FALSE)
htert_reference_stats <- read.delim("tables/hTERT_123CpG_reference_tests.tsv",
                                    check.names = FALSE, stringsAsFactors = FALSE)

# Match by endpoint, treatment and comparator, not by a historical panel letter.
# Old statistical tables call the current Skin & Blood panel h and RetroAge
# i-right. CURRENT_STATISTICAL_LABELS.tsv supplies the final artwork mapping.
handled_panels <- c("4d", "S4a", "4f", "4f legend matched comparator",
                    "4i right", "S4c", "4j right", "S4e", "S4d", "S4k",
                    "4o centre", "4o right", "hTERT reference-distance support")
annotation_validation <- list()
for (i in which(current_labels$panel %in% handled_panels)) {
  label <- current_labels[i, , drop = FALSE]
  if (label$panel %in% c("4o centre", "4o right", "hTERT reference-distance support")) {
    stopifnot(label$controls == "matched", label$endpoint == "123-CpG reference distance")
    if (label$panel == "hTERT reference-distance support") {
      endpoint <- switch(label$hypothesis,
        "Final-cycle 96-h displacement" = "third_cycle_96h",
        "Run-mean displacement over three cycles" = "mean_three_rounds_pp",
        stop("Unknown hTERT reference-distance hypothesis: ", label$hypothesis))
      source_stats <- htert_reference_stats
      source_name <- "hTERT_123CpG_reference_tests"
    } else {
      endpoint <- if (label$panel == "4o centre") "all_T1_T6" else "late_T5_T6"
      expected_hypothesis <- if (label$panel == "4o centre")
        "Run-mean displacement T1-T6" else "Run-mean displacement T5-T6"
      stopifnot(label$hypothesis == expected_hypothesis)
      source_stats <- primary_reference_stats
      source_name <- "4o_123CpG_primary_tests"
    }
    found <- source_stats[source_stats$endpoint == endpoint &
                            source_stats$arm == label$treatment, , drop = FALSE]
    stopifnot(nrow(found) == 1L, found$n_experiments == 3L, found$df == 2)
    replay_p <- found$p
    replay_q <- found$q
  } else if (label$panel == "S4k") {
    found <- htert_stats[htert_stats$clock == label$endpoint &
                          htert_stats$arm == label$treatment, , drop = FALSE]
    stopifnot(nrow(found) == 1L, label$controls == "matched")
    replay_p <- found$p_global_shift
    replay_q <- found$q_global_BH4
    source_name <- "S4k_clock_global_statistics"
  } else if (label$panel == "S4d") {
    found <- epigenetic_stats[epigenetic_stats$clock == label$endpoint &
                               epigenetic_stats$arm == label$treatment, , drop = FALSE]
    stopifnot(nrow(found) == 1L, label$controls == "matched")
    replay_p <- found$p_average_shift
    replay_q <- found$q_BH_two_arms
    source_name <- "S4_epigenetic_clock_global_statistics"
  } else {
    comparator <- if (label$controls == "matched") "own" else "pooled"
    found <- main_stats[main_stats$endpoint == label$endpoint &
                          main_stats$arm == label$treatment &
                          main_stats$control_definition == comparator, , drop = FALSE]
    stopifnot(nrow(found) == 1L, found$hypothesis == label$hypothesis)
    replay_p <- found$p
    replay_q <- found$q_BH_panel
    source_name <- "main_statistics_own_vs_pooled_long"
  }
  stopifnot(all(is.finite(c(replay_p, replay_q, label$P_value, label$q_value))),
            all(c(replay_p, replay_q, label$P_value, label$q_value) > 0))
  key <- paste(label$panel, label$endpoint, label$treatment, label$controls)
  # Log-scale comparison also detects errors in very small P/q values.
  assert_close(log(c(replay_p, replay_q)), log(c(label$P_value, label$q_value)),
               paste("Current statistical annotation:", key), tolerance = 1e-5)
  annotation_validation[[length(annotation_validation) + 1L]] <- data.frame(
    panel = label$panel, endpoint = label$endpoint, treatment = label$treatment,
    controls = label$controls, hypothesis = label$hypothesis,
    annotation_p = label$P_value, replay_p = replay_p,
    annotation_q = label$q_value, replay_q = replay_q,
    source_table = source_name, stringsAsFactors = FALSE)
}
save_table(bind_rows(annotation_validation), "current_annotation_validation")

format_annotation_q <- function(q) {
  if (!is.finite(q) || q < 0 || q > 1) stop("Invalid q value for annotation.")
  if (q < .001) {
    sub("e([+-])0+", "e\\1", formatC(q, format = "e", digits = 1))
  } else {
    format(signif(q, 3), scientific = FALSE, trim = TRUE)
  }
}
annotation_q <- function(panel, endpoint, arm, comparator) {
  z <- current_labels[current_labels$panel == panel &
                        current_labels$endpoint == endpoint &
                        current_labels$treatment == arm &
                        current_labels$controls == comparator, , drop = FALSE]
  stopifnot(nrow(z) == 1L)
  z$q_value
}
q_pair <- function(panel, endpoint, comparator) {
  paste0("INH = ", format_annotation_q(annotation_q(panel, endpoint, "INH", comparator)),
         "   MIM = ", format_annotation_q(annotation_q(panel, endpoint, "MIM", comparator)))
}
add_global_q <- function(plot, panel, endpoint, comparator,
                         hypothesis = "Global mean-shift") {
  plot + labs(subtitle = paste0(hypothesis, " q vs ", comparator, " CTRL\n",
                                q_pair(panel, endpoint, comparator))) +
    theme(plot.subtitle = element_text(size = 8.5, lineheight = 1.1,
                                       margin = margin(b = 7)),
          plot.margin = margin(6, 8, 6, 6))
}

# Main pooled curves and supplementary matched curves retain the same plotted
# data and uncertainty bands. Only subtitles are added outside the data area.
curve_spec <- data.frame(
  plot = c("4d_global_methylation", "S4a_global_methylation_separate_controls",
           "4i_right_SkinBlood", "S4c", "4j_right_RetroAge", "S4e"),
  panel = c("4d", "S4a", "4i right", "S4c", "4j right", "S4e"),
  endpoint = c("Global methylation", "Global methylation", "Skin & Blood",
               "Skin & Blood", "RetroAge", "RetroAge"),
  comparator = c("pooled", "matched", "pooled", "matched", "pooled", "matched"),
  width = c(5.5, 5.5, 5.3, 5.3, 5.3, 5.3),
  height = c(3.2, 5.1, 5.2, 5.2, 5.2, 5.2), stringsAsFactors = FALSE)
for (i in seq_len(nrow(curve_spec))) {
  spec <- curve_spec[i, ]
  stopifnot(spec$plot %in% names(panels))
  hypothesis <- if (spec$endpoint == "Global methylation")
    "Global interaction" else "Global mean-shift"
  annotated <- add_global_q(panels[[spec$plot]], spec$panel, spec$endpoint,
                            spec$comparator, hypothesis)
  save_panel(annotated, spec$plot, spec$width, spec$height)
}

# Causal-clock panels have matched controls and their own two-arm BH families.
for (clock in c("CausAge", "DamAge", "AdaptAge")) {
  name <- paste0("S4d_", clock)
  stopifnot(name %in% names(panels))
  save_panel(add_global_q(panels[[name]], "S4d", clock, "matched"), name, 5.3, 5.2)
}
save_panel(plot_grid(plotlist = panels[c("S4d_CausAge", "S4d_DamAge", "S4d_AdaptAge")],
                     nrow = 1), "S4d_causal_clocks", 13.5, 5.1)

# hTERT RNA clocks: display inputs are unchanged. The inference remains the
# established ML-LRT global shift, BH across two shifts plus two interactions
# within each clock. Rebuild individual panels to annotate the combined view.
htert_display <- as.data.frame(input("clocks/htert_clocks.rds")$display)
htert_display$TP2 <- htert_display$TIME
htert_display$series <- as.character(htert_display$GROUP)
htert_annotated <- lapply(c("RNAAgeCalc", "tAge", "PASTA"), function(clock) {
  plot <- trajectory(htert_display[htert_display$clock == clock, , drop = FALSE],
                     if (clock == "tAge") "Relative tAge (model-scaled units)" else
                       paste(clock, "score"), xpos = 1:3, title = clock)
  add_global_q(plot, "S4k", clock, "matched")
})
save_panel(plot_grid(plotlist = htert_annotated, nrow = 1),
           "S4k_hTERT_transcriptomic_clocks", 13.5, 5.1)

# Figure 4f is a T6 endpoint, not a global curve. Put q values in each facet
# heading, with matched-control values in the caption as in the figure legend.
stopifnot("4f_entropy_FAE" %in% names(panels))
endpoint_labels <- setNames(vapply(c("entropy", "FAE"), function(endpoint) {
  paste0(if (endpoint == "entropy") "Entropy" else "FAE", "\n",
         "q(INH) = ", format_annotation_q(annotation_q("4f", endpoint, "INH", "pooled")),
         "\nq(MIM) = ", format_annotation_q(annotation_q("4f", endpoint, "MIM", "pooled")))
}, character(1)), c("entropy", "FAE"))
matched_caption <- paste0(
  "Matched CTRL q: entropy ", q_pair("4f legend matched comparator", "entropy", "matched"),
  "\nFAE ", q_pair("4f legend matched comparator", "FAE", "matched"))
endpoint_plot <- panels[["4f_entropy_FAE"]] +
  facet_wrap(~name, scales = "free_y", nrow = 1,
             labeller = as_labeller(endpoint_labels)) +
  labs(title = "T6 vs pooled CTRL", caption = matched_caption) +
  theme(strip.text = element_text(size = 8.5, face = "plain", lineheight = 1.1),
        plot.caption = element_text(size = 7.5, hjust = 0, margin = margin(t = 7)),
        plot.margin = margin(6, 8, 6, 6))
save_panel(endpoint_plot, "4f_entropy_FAE", 5.5, 3.6)

# No loop above includes Figure 4l or any reference-distance panel. Their
# existing reconstructions are preserved by this annotation-only module.
