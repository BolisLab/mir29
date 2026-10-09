---
title: "MIR29 trades youthful plasticity for epigenetic fidelity"
subtitle: "Figure 5 and Supplementary Figure S5 — methods, data and executable reconstruction"
date: "1 October 2026"
output:
  html_document:
    self_contained: true
    toc: true
    toc_float: true
    toc_depth: 2
    code_folding: show
    code_download: true
    fig_caption: true
    highlight: tango
    css: styles.css
    keep_md: true
---

> **Repository snapshot — input bundle pending.** This chapter preserves the existing analytical narrative, code and rendered outputs. References below to bundled input data describe the full analytical bundle, not this lightweight repository snapshot. Required inputs are not yet fully distributed here, so a fresh clone cannot regenerate the analysis. See [module instructions](README.md).

# Scope and execution

This document accompanies the final **Figure 5 and Figure S5**. It reconstructs
the associations linking miR-29 expression, DNA methylation age, the methylome
remodelling core (MRC), and deviation from normal-tissue methylation references
in human cancers and cancer cell lines. Figure S5 documents cohort sizes,
alternative age readouts, leave-one-stratum-out sensitivity, and reference
selection. These observational associations do not establish causal direction.

The final Illustrator pages are preserved unchanged in `reference/`.
Every analytical plot below is regenerated from the bundled numerical inputs.
The Cytoscape network in Figure 5d and the conceptual schematic in Figure 5g
are explicitly identified as retained artwork. Page geometry, organ icons and
other manually arranged decorations are not analytical outputs.

Open `Figure5_main_plus_supplementary.Rproj` in RStudio and knit this file, or run:


```r
source("render_figure5_and_s5.R")
# Rebuild numerical summaries, statistical tests and plots without HTML:
source("run_all.R")
```

Setup executes the clean modules once. Every displayed code block is read
directly from those executable files. Supplementary sections contain analysis
and plotting code, not just images. Paths are relative to the bundle root;
knitting does not download data, install software, recompute genome-wide
preprocessing or require the original server.

## Reproduction boundary and provenance

The supported starting point is the **verified processed-data checkpoint**.
The bundle includes per-cancer resampled correlations and annotations,
CCLE results and metadata, RNA correlation matrices, per-sample paired
expression scores, tumor-level fidelity scores, and the selected-reference
CpG vectors and selection metrics. Independent historical outputs are also
included as numerical regression checks.

The reconstruction recalculates meta-ranks, the 185-gene selection, strict
patient matching, paired tests, Pearson associations, residual-Spearman
associations, random-effects summaries, leave-one-out results, FDR values,
boxplot statistics and reference-QC summaries. It preserves, rather than
refits, upstream VST normalization, the 100 resampling iterations, clock
predictions, and tumor-versus-reference RMS scores.

Raw count matrices and genome-wide tumor/normal beta matrices are **not**
bundled. Original server scripts are retained in
`provenance/upstream_archive/`; they are not the clean executable entry
point and may require their original server inputs. The matrix-level
reference algorithm is separately extracted in clean form below. A successful
render certifies processed-data reconstruction, not a new raw-data/server
run or recreation of Cytoscape's exact edge layout.

`provenance/input_sources.tsv` records original paths, transformations,
sizes and MD5 checksums. Unused VST matrices were removed from the portable
miRNA cache without changing correlations or annotations. Reference vectors
were collected into `data/reference_selection_metrics.rds` without changing
their values. Source-page checksums and artwork crop coordinates are recorded
in `provenance/reference_sources.json`.

## Analysis map and conventions

Here **P** denotes Horvath's pan-tissue clock and **S** denotes Skin & Blood.
The word "acceleration" must not be used to interchange these clocks.
Acceleration is the residual of a linear DNAmAge-versus-chronological-age
model fitted within the relevant cancer or cell-line lineage.

| Panels | Units / strata | Readout and aggregation |
|:--|:--|:--|
| 5a; S5a–c | 8,428 tumor patients; 30 cancers | Resampled Spearman correlations; equal-cancer Fisher aggregation |
| 5b; S5d–f | 429 cell lines; 11 lineages | Resampled Spearman correlations; equal-lineage Fisher aggregation |
| 5c; S5g | RNA features within 30 cancers | Equal-cancer Fisher rank; native within-cancer distributions |
| 5e–f; S5h–i | 580 strict tumor–normal pairs; 16 cancers | Paired score differences; cancer-mean Pearson association |
| 5g–h; S5j–k | 15 reference-compatible histotypes | RMS deviation in M values; age-spline residual correlations |

For correlations \(r_j\), the equal-stratum Fisher summary is
\[
\bar r=\tanh\!\left[\frac{1}{K}\sum_{j=1}^{K}\operatorname{atanh}(r_j)\right].
\]
Each eligible cancer or lineage has equal weight. This is not a correlation
obtained by pooling all patients, nor an inverse-variance random-effects
meta-analysis. Finite available entries are used for feature ranks; their
numbers of contributing strata are retained in the tables.

Random-effects pooling is used separately for the paired expression and
fidelity forests, using REML with Knapp–Hartung inference. The two
organ-aggregation procedures remain distinct: paired expression uses an
initial within-organ random-effects summary; fidelity combines histotypes
within organ by inverse-variance Fisher weighting before between-organ
pooling.

All tests are two-sided. Nominal \(p\), BH-adjusted \(q\), correlation
intervals, CpG-distribution quantiles and descriptive leave-one-out connectors
are distinguished explicitly. The 185-gene three-SD threshold is an
operational tail selection, **not an FDR threshold**.

## Shared setup and plotting functions

The modules stop if inputs are absent or a numerical check fails.
Historical names such as `ERC` do not change the biological definition:
MRC comprises **DNMT3A, DNMT3B, TET1, TET3 and TDG** throughout.


```r
# All paths and packages required by the processed-data reconstruction.
BUNDLE <- normalizePath(".")
DATA_ROOT <- file.path(BUNDLE, "data/source")
ROOT <- BUNDLE
PDF_DIR <- file.path(BUNDLE, "figures/reconstructed")
for (p in c(PDF_DIR, file.path(BUNDLE, c("data", "validation", "previews", "provenance", "code/helpers")))) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
packages <- c("ggplot2", "ggrepel", "cowplot", "metafor", "splines", "jsonlite")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({
  library(ggplot2)
  library(splines)
})
CANCERS16 <- c("BLCA", "BRCA", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD",
  "LUSC", "PRAD", "STAD", "THCA", "UCEC")
CANCERS15 <- setdiff(CANCERS16, "STAD")
MRC_GENES <- c("DNMT3A", "DNMT3B", "TET1", "TET3", "TDG")
EXPECTED_H_LEFT_N <- c(BLCA = 409, BRCA = 773, CHOL = 36, COAD = 293, ESCA = 184, HNSC = 522, KICH = 66,
  KIRC = 318, KIRP = 272, LIHC = 371, LUAD = 437, LUSC = 357, PRAD = 494, THCA = 506, UCEC = 427)
EXPECTED_H_RIGHT_N <- c(BLCA = 406, BRCA = 781, CHOL = 35, COAD = 292, ESCA = 162, HNSC = 501, KICH = 65,
  KIRC = 318, KIRP = 271, LIHC = 370, LUAD = 436, LUSC = 361, PRAD = 497, THCA = 504, UCEC = 427)
INPUT_FILES <- character()
read_data <- function(relative) {
  path <- file.path(DATA_ROOT, relative)
  stopifnot(file.exists(path))
  INPUT_FILES <<- unique(c(INPUT_FILES, path))
  readRDS(path)
}
CHECKS <- list()
check <- function(label, condition, detail = "") {
  ok <- isTRUE(condition)
  CHECKS[[length(CHECKS) + 1L]] <<- list(check = label, pass = ok, detail = as.character(detail))
  message(if (ok)
    "PASS " else "FAIL ", label, if (nzchar(detail))
    paste0(": ", detail))
  if (!ok)
    stop("Validation failed: ", label)
  invisible(ok)
}
equal_numeric <- function(a, b, tolerance = 1e-08) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance, check.attributes = FALSE))
}
fisher_mean <- function(r) {
  r <- r[is.finite(r)]
  if (!length(r))
    return(NA_real_)
  tanh(mean(atanh(pmin(pmax(r, -0.999999), 0.999999))))
}


cowplot::set_null_device("agg")
load_helpers <- function(required, envir, output_name) {
  source(file.path(BUNDLE, "code/helpers", output_name), local = envir)
  stopifnot(all(vapply(required, exists, logical(1), envir = envir, inherits = FALSE)))
}
save_table <- function(x, name) {
  write.table(as.data.frame(x), file.path("tables", paste0(name, ".tsv")), sep = "\t", row.names = FALSE,
    quote = FALSE, na = "NA")
}
```


```r
# Figure 5 drawing layer. All values are read from the validated panel tables.  Graphical assembly
# is regenerated, not inferred from pixels in Illustrator.
if (!exists("PANELS")) PANELS <- readRDS(file.path(BUNDLE, "data", "Figure5_panel_data.rds"))
COLORS <- c(`miR-29a` = "#E8D833", `miR-29b` = "#91C7D3", `miR-29c` = "#E74337", `pan-miR-29` = "#111111")
panel_theme <- function(size = 10) theme_classic(base_size = size, base_family = "sans") + theme(plot.title = element_text(size = size +
  1, face = "bold", margin = margin(b = 8)), axis.line = element_line(linewidth = 0.35, colour = "#4A4A4A"),
  axis.ticks = element_line(linewidth = 0.3, colour = "#4A4A4A"), plot.margin = margin(8, 10, 8, 10),
  legend.position = "bottom", legend.title = element_blank())
save_panel <- function(plot, name, width, height) {
  ggsave(file.path(PDF_DIR, paste0(name, ".pdf")), plot, width = width, height = height, device = cairo_pdf,
    bg = "white", symbolfamily = grDevices::cairoSymbolFont("Arial", usePUA = FALSE))
  ggsave(file.path(PDF_DIR, paste0(name, ".png")), plot, width = width, height = height, dpi = 180,
    bg = "white")
}
combine_panels <- function(plots, widths, title) {
  # Align panel regions exactly even when one graph has a bottom legend and another has top-axis
  # labels. Shared cancer rows must have identical height.
  grobs <- lapply(plots, ggplotGrob)
  above <- lapply(grobs, function(g) sum(g$heights[seq_len(g$layout$t[g$layout$name == "panel"] - 1)]))
  below <- lapply(grobs, function(g) sum(g$heights[seq.int(g$layout$b[g$layout$name == "panel"] + 1,
    length(g$heights))]))
  max_above <- do.call(grid::unit.pmax, above)
  max_below <- do.call(grid::unit.pmax, below)
  grobs <- lapply(seq_along(grobs), function(i) {
    g <- gtable::gtable_add_rows(grobs[[i]], max_above - above[[i]], pos = 0)
    gtable::gtable_add_rows(g, max_below - below[[i]], pos = -1)
  })
  body <- cowplot::plot_grid(plotlist = grobs, nrow = 1, rel_widths = widths)
  heading <- cowplot::ggdraw() + cowplot::draw_label(title, x = 0.015, hjust = 0, fontface = "bold",
    size = 14)
  cowplot::plot_grid(heading, body, ncol = 1, rel_heights = c(0.075, 1))
}
rank_plot <- function(d, ylab, title, ylim = NULL) {
  highlighted <- d[!is.na(d$Feature), ]
  background <- d[is.na(d$Feature), ]
  zerorank <- max(c(1, d$Rank[d$rho <= 0]))
  p <- ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#FAE9E9") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#E7EDF5") + geom_hline(yintercept = 0,
    linewidth = 0.35) + geom_vline(xintercept = zerorank, linetype = "dashed", linewidth = 0.3, color = "grey60") +
    geom_point(data = background, aes(Rank, rho, colour = rho > 0), size = 0.5) + scale_colour_manual(values = c(`FALSE` = "#426FAD",
    `TRUE` = "#D95C64"), guide = "none") + geom_point(data = highlighted, aes(Rank, rho, fill = Feature,
    shape = Feature), size = 2.5, colour = "#555555", stroke = 0.3) + scale_fill_manual(values = COLORS,
    guide = "none") + scale_shape_manual(values = c(`miR-29a` = 21, `miR-29b` = 21, `miR-29c` = 21,
    `pan-miR-29` = 24), guide = "none") + ggrepel::geom_text_repel(data = highlighted, aes(Rank,
    rho, label = Feature), size = 2.7, direction = "y", nudge_x = -0.22 * max(d$Rank), hjust = 1,
    seed = 29, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.18, segment.size = 0.25) +
    labs(x = "miRNA meta-rank", y = ylab, title = title) + panel_theme(9)
  if (!is.null(ylim))
    p <- p + coord_cartesian(ylim = ylim)
  p
}
bubble_plot <- function(d, summary, row_names, metric_order, metric_labels, labels = row_names, title = NULL,
  summary_label = "Pan-cancer") {
  n <- length(row_names)
  d$row <- n + 1 - match(d$Row, row_names)
  d$x <- match(d$Metric, metric_order)
  summary$x <- match(summary$Metric, metric_order)
  summary$row <- -1.15
  display <- rbind(d[, c("row", "x", "rho")], summary[, c("row", "x", "rho")])
  ggplot(display, aes(x, row)) + annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -1.6, ymax = -0.65,
    fill = "#F0F0F0", colour = "grey40", linewidth = 0.3) + geom_hline(yintercept = 0.05, linewidth = 0.25,
    linetype = "dashed", colour = "grey65") + geom_vline(xintercept = 1:3, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(size = abs(rho), fill = rho), shape = 21, colour = "grey45", stroke = 0.25) +
    scale_size_area(max_size = 5, limits = c(0, max(abs(display$rho)) + 0.01), guide = "none") +
    scale_fill_gradient2(low = "#426FAD", mid = "white", high = "#D65359", midpoint = 0, limits = c(-0.55,
      0.55), oob = scales::squish, guide = "none") + geom_text(data = summary, aes(label = sprintf("%.2f",
    rho)), nudge_y = -0.9, size = 2.8) + scale_x_continuous(breaks = 1:3, labels = metric_labels,
    position = "top", limits = c(0.5, 3.5), expand = c(0, 0)) + scale_y_continuous(breaks = c(-1.15,
    seq_len(n)), labels = c(summary_label, rev(labels)), limits = c(-2.7, n + 0.7), expand = c(0,
    0)) + labs(x = NULL, y = NULL, title = title) + panel_theme(9) + theme(axis.line = element_blank(),
    axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0), axis.text.y = element_text(size = 8.5),
    plot.margin = margin(8, 4, 8, 5))
}
```

# Figure 5 — main

## a — miR-29 and methylation age across human cancers

The analysis uses primary and metastatic TCGA tumors with matched miRNA,
chronological age and Horvath DNAmAge. One sample is retained per patient,
preferring primary tumor when both types are available. Adjacent normals
are excluded. There are **8,428 unique patients in 30 cancers**; ovarian
cancer is not part of the retained analysis.

The upstream workflow sums duplicated sample IDs at the raw-count level.
The two miR-29b loci are summed into a combined b feature; pan-miR-29 is
the sum of miR-29a, miR-29b-1, miR-29b-2 and miR-29c before DESeq2 variance
stabilization. The preserved implementation uses a blind VST with fallback
fitting when required. Within each cancer, 100 random 80% subsamples
produce correlations with age, DNAmAge and acceleration, averaged on the
Fisher-z scale. The upstream minimum cohort size is 25. The exact
processed output is replayed here.

Left: pan-miR-29 correlations within each cancer and their equal-cancer
summary. Center: native-miRNA DNAmAge correlation distributions with the
three family members highlighted. Redundant b-locus entries and the
auxiliary AC composite are excluded from this display. Right: the
cross-cancer rank, including pan-miR-29, contains 1,764 entries.
Points are ranked observations, not a fitted curve.


```r
# Build final panel tables from validated intermediates. Preserve the historical estimands, cohorts
# and score definitions. Do not recompute VST or resampling.
message("Preparing panels a-c")
a_cache <- read_data(paste0("miR29_Horvath_repeated_resampling/", "MIR29_CLOCK_RESULTS_Horvath_AGE_Acceleration_Pearson_Spearman_v1.rds"))
a_all <- do.call(rbind, lapply(a_cache, `[[`, "correlations"))
a_sample_n <- sum(vapply(a_cache, function(x) nrow(x$annotation), integer(1)))
rm(a_cache)
invisible(gc())
a_pan <- a_all[a_all$miRNA == "pan_miR29", ]
a_pan <- a_pan[order(a_pan$Cancer), ]
a_cancers <- a_pan$Cancer
check("a: 30 cancer types", length(a_cancers) == 30L)
a_metrics <- c("AGE", "Horvath", "Acceleration")
a_bubbles <- do.call(rbind, lapply(a_metrics, function(metric) {
  data.frame(Cancer = a_pan$Cancer, Metric = metric, rho = a_pan[[paste0("Spearman_", metric, "_Fisher")]])
}))
a_summary <- data.frame(Metric = a_metrics, rho = vapply(a_metrics, function(m) {
  fisher_mean(a_bubbles$rho[a_bubbles$Metric == m])
}, numeric(1)))
check("a: displayed summary correlations", identical(unname(round(a_summary$rho, 2)), c(0.09, 0.26, 0.23)),
  paste(round(a_summary$rho, 6), collapse = ", "))
redundant <- c("hsa-mir-29b-1", "hsa-mir-29b-2", "miR29AC")
a_distribution <- a_all[!a_all$miRNA %in% c(redundant, "pan_miR29") & is.finite(a_all$Spearman_Horvath_Fisher),
  c("Cancer", "miRNA", "N", "Spearman_Horvath_Fisher")]
names(a_distribution)[4] <- "rho"
a_rank <- do.call(rbind, lapply(split(a_all, a_all$miRNA), function(x) {
  r <- x$Spearman_Horvath_Fisher
  data.frame(miRNA = x$miRNA[1], K = sum(is.finite(r)), rho = fisher_mean(r))
}))
a_rank <- a_rank[!a_rank$miRNA %in% redundant & is.finite(a_rank$rho), ]
a_rank <- a_rank[order(a_rank$rho), ]
a_rank$Rank <- seq_len(nrow(a_rank))
check("a: no OV in final analysis", !"OV" %in% a_cancers)
mir_labels <- c(`hsa-mir-29a` = "miR-29a", miR29b_combined = "miR-29b", `hsa-mir-29c` = "miR-29c", pan_miR29 = "pan-miR-29")
a_rank$Feature <- unname(mir_labels[a_rank$miRNA])
a_distribution$Feature <- unname(mir_labels[a_distribution$miRNA])
```


```r
# a: identical cancer ordering across bubble and distribution plots.
a <- PANELS$a
ab <- a$bubbles
ab$Row <- ab$Cancer
pa_bubble <- bubble_plot(ab, a$summary, a$cancers, c("AGE", "Horvath", "Acceleration"), c("Age", "DNAmAge (P)",
  "Age acc. (P)"), title = "Σ(miR-29)")
ad <- a$distribution
ad$row <- length(a$cancers) + 1 - match(ad$Cancer, a$cancers)
pa_box <- ggplot(ad, aes(rho, row, group = Cancer)) + geom_vline(xintercept = 0, linetype = "dashed",
  colour = "grey60", linewidth = 0.3) + geom_boxplot(orientation = "y", width = 0.65, fill = "#626868",
  colour = "#444848", linewidth = 0.25, outlier.size = 0.18, outlier.alpha = 0.5) + geom_point(data = ad[!is.na(ad$Feature),
  ], aes(fill = Feature), shape = 21, size = 1.7, stroke = 0.25, colour = "#666666") + scale_fill_manual(values = COLORS,
  breaks = names(COLORS)[1:3]) + scale_y_continuous(limits = c(-2.7, length(a$cancers) + 0.7), expand = c(0,
  0)) + scale_x_continuous(breaks = c(-0.5, 0, 0.5)) + labs(x = "miRNA expression ~ DNAmAge [ρ]",
  y = NULL, title = "Within-cancer correlations") + panel_theme(9) + theme(axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.text = element_text(size = 8),
  legend.key.width = grid::unit(0.65, "cm"))
pa_rank <- rank_plot(a$ranking, "miRNA ~ DNAmAge [ρ]", "Cross-cancer meta-rank")
pa <- combine_panels(list(pa_bubble, pa_box, pa_rank), c(1.15, 1.6, 1.05), "a | TCGA miR-29 and epigenetic age")
save_panel(pa, "Figure5a_TCGA_miRNA_age", 9.3, 7.1)
```

<img src="figures/reconstructed/Figure5a_TCGA_miRNA_age.png" width="100%" style="display: block; margin: auto;" />


|Metric       |      rho|
|:------------|--------:|
|AGE          | 0.090054|
|Horvath      | 0.256610|
|Acceleration | 0.230591|

## b — miR-29 and age readouts in cancer cell lines

The CCLE checkpoint contains **429 eligible lines from 11 lineages**.
Eligibility requires complete age, Horvath and Skin & Blood values and
at least 15 lines per lineage; the historical filename contains "N20",
but the saved setting and final cohort use **MIN_N = 15**.
Chronological age is donor age, not culture duration.

Native miRNA values come from the processed NanoString input. Within each
lineage, the three miR-29 measurements are individually standardized and
averaged to form pan-miR-29. This differs from the TCGA raw-count sum.
The saved 100 × 80% repeated-subsampling analysis is aggregated equally
across lineages.

Bubbles show chronological age, Horvath pan-tissue age and Skin & Blood
age. The main rank uses **Skin & Blood acceleration**, not Horvath
acceleration. Figure S5f shows the other retained readouts.


```r
b_final <- read_data("CCLE_miR29_DNAmAge_FINAL/CCLE_miR29_Age_Horvath_Skin_with_acceleration_N20_FINAL.rds")
b_tissues <- sort(b_final$settings$EligibleTissues)
b_bubbles <- subset(b_final$PerTissue, miRNA == "pan_miR29" & Phenotype %in% c("Age", "Horvath", "Skin"))
b_bubbles <- data.frame(Tissue = b_bubbles$Tissue, Metric = b_bubbles$Phenotype, rho = b_bubbles$rho_Fisher_resampled)
b_summary <- data.frame(Metric = c("Age", "Horvath", "Skin"), rho = vapply(c("Age", "Horvath", "Skin"),
  function(m) fisher_mean(b_bubbles$rho[b_bubbles$Metric == m]), numeric(1)))
b_rank <- b_final$Rankings$SkinAcceleration
b_rank <- b_rank[order(b_rank$rho_Fisher_equalTissue), ]
b_rank$rho <- b_rank$rho_Fisher_equalTissue
b_rank$Rank <- seq_len(nrow(b_rank))
b_labels <- c(MIMAT0000086 = "miR-29a", MIMAT0000100 = "miR-29b", MIMAT0000681 = "miR-29c", pan_miR29 = "pan-miR-29")
b_rank$Feature <- unname(b_labels[b_rank$miRNA])
check("b: 11 lineages and MIN_N=15", length(b_tissues) == 11L && b_final$settings$MIN_N == 15)
check("b: displayed summary correlations", identical(unname(round(b_summary$rho, 2)), c(0.03, 0.27, 0.41)))
check("b: pan-miR29 top Skin&Blood acceleration rank", tail(b_rank$miRNA, 1) == "pan_miR29")
```


```r
# b: chronological age, Horvath age, Skin&Blood age; ranking uses Skin&Blood acceleration.
b <- PANELS$b
bb <- b$bubbles
bb$Row <- bb$Tissue
tissue_labels <- c("breast", "CNS", "blood", "intestine", "liver", "lung", "esophagus", "ovary", "pancreas",
  "skin", "stomach")
pb_bubble <- bubble_plot(bb, b$summary, b$tissues, c("Age", "Horvath", "Skin"), c("Age", "DNAmAge (P)",
  "DNAmAge (S)"), labels = tissue_labels, title = "Σ(miR-29)", summary_label = "Pan-lineage")
pb_rank <- rank_plot(b$ranking, "miRNA ~ age acceleration (S) [ρ]", "Skin & Blood acceleration")
pb <- combine_panels(list(pb_bubble, pb_rank), c(1.1, 1), "b | Cancer cell lines: 11 lineages")
save_panel(pb, "Figure5b_CCLE_miRNA_age", 6.8, 4.8)
```

<img src="figures/reconstructed/Figure5b_CCLE_miRNA_age.png" width="100%" style="display: block; margin: auto;" />


|Metric  |      rho|
|:-------|--------:|
|Age     | 0.034284|
|Horvath | 0.267813|
|Skin    | 0.408831|

## c — RNA age associations and the 185-gene lower tail

RNA counts were processed upstream with DESeq2 variance stabilization,
and nonvariable features removed within cancer. Repeated-subsample
Spearman associations with Horvath DNAmAge were computed within each
cancer. This section rebuilds the equal-cancer Fisher mean from the
saved complete 30-cancer correlation matrix, rather than reading the
final rank alone.

The 44,398 ranked features are ordered by their meta-correlation.
The selected set is defined by
\[
r_g < \operatorname{mean}(r)-3\,\operatorname{SD}(r).
\]
Exactly **185 Ensembl features** satisfy this rule. Their IDs are checked
against the original selection and exported with symbols and the full
rank. Association with this epigenetic-age readout is not proof of
repression during chronological aging or of a causal role.


```r
c_final <- read_data("FIG6C_TCGA_RNA_FINAL/FIG6C_TCGA_RNA_Horvath_FINAL.rds")
c_rank <- c_final$meta
recalculated_cutoff <- mean(c_rank$rho) - 3 * sd(c_rank$rho)
check("c: exact negative-tail threshold", equal_numeric(recalculated_cutoff, c_final$threshold$cutoff))
c_185 <- c_rank[c_rank$rho < recalculated_cutoff, ]
check("c: 185 selected genes", nrow(c_185) == 185L)
check("c: selected identifiers equal saved selection", identical(c_185$FeatureID, c_final$agingRepressed_3SD$FeatureID))
c_matrix <- c_final$rhoCancerMatrix[c_rank$FeatureID, , drop = FALSE]
c_rebuilt <- tanh(rowMeans(atanh(pmin(pmax(c_matrix, -0.999999), 0.999999))))
check("c: meta-correlations rebuilt from per-cancer matrix", equal_numeric(c_rebuilt, c_rank$rho))
rm(c_matrix, c_rebuilt)
writeLines(c_185$FeatureID, file.path(BUNDLE, "data", "panel_c_185_ensembl_ids.txt"))
writeLines(unique(na.omit(c_185$SYMBOL)), file.path(BUNDLE, "data", "panel_c_185_gene_symbols.txt"))
```


```r
# c: actual gene meta-ranks and the negative-tail selection.
c <- PANELS$c
cr <- c$ranking
show_genes <- c("MIR29B2CHG", "LINC-PINT", MRC_GENES)
ch <- cr[cr$SYMBOL %in% show_genes & !is.na(cr$SYMBOL), ]
ch <- ch[order(ch$SYMBOL, -abs(ch$rho)), ]
ch <- ch[!duplicated(ch$SYMBOL), ]
ch$Colour <- ifelse(ch$SYMBOL %in% MRC_GENES, "#D36E93", ifelse(ch$SYMBOL == "LINC-PINT", "#AE8EB7",
  "#E2A1A1"))
pc <- ggplot(cr, aes(MetaRank, rho)) + annotate("rect", xmin = 0, xmax = max(c$selected$MetaRank), ymin = -Inf,
  ymax = c$threshold, fill = "#F5E696", alpha = 0.7) + geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_line(linewidth = 0.4, colour = "#333333") + geom_point(data = ch, aes(fill = Colour), shape = 21,
  size = 2.6, stroke = 0.35, colour = "#555555") + scale_fill_identity() + ggrepel::geom_text_repel(data = ch,
  aes(label = SYMBOL), seed = 29, size = 3.2, fontface = "italic", nudge_x = ifelse(ch$rho < 0, 8500,
    -7500), direction = "y", min.segment.length = 0, segment.size = 0.25, box.padding = 0.35, max.overlaps = Inf) +
  annotate("text", x = 16000, y = -0.25, label = "185 genes\nbelow mean - 3 SD", size = 3.4, hjust = 0) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000, 40000), labels = c("0", "10k", "20k", "30k",
    "40k")) + labs(x = "Gene meta-rank", y = "Gene expression ~ DNAmAge [meta-ρ]", title = "c | Pan-cancer transcriptome ranking") +
  panel_theme(10)
save_panel(pc, "Figure5c_RNA_age_ranking", 5.1, 4.6)
```

<img src="figures/reconstructed/Figure5c_RNA_age_ranking.png" width="100%" style="display: block; margin: auto;" />


| Ranked_features| Cancers|    Cutoff| Selected_features|
|---------------:|-------:|---------:|-----------------:|
|           44398|      30| -0.132566|               185|

The complete annotated selection is in `tables/Figure5c_185_genes.tsv`;
exact network input IDs are in `data/panel_c_185_ensembl_ids.txt`.

## d — network context and the methylome remodelling core

The approved Cytoscape network was assembled from the selected 185-gene
set. The final diagram highlights the five MRC genes: DNMT3A, DNMT3B,
TET1, TET3 and TDG. The analytical input list is reproducible in panel c.
Exact interaction retrieval, edge filtering and Cytoscape layout are
not certified by this R reconstruction. The vector artwork is therefore
preserved directly, without substituting an inferred network.


```r
MRC_GENES <- c('DNMT3A','DNMT3B','TET1','TET3','TDG')
network_input <- readLines('data/panel_c_185_ensembl_ids.txt')
stopifnot(length(network_input)==185L)
# Approved artwork: reference/Figure5d_Cytoscape_artwork.pdf
```

<img src="reference/Figure5d_Cytoscape_artwork.png" width="100%" style="display: block; margin: auto;" />

## e — paired tumor–normal expression differences

Within each cancer, primary tumor (01) and adjacent normal (11) samples
are matched by patient. A patient is retained only if both tissue states
are available for **both miRNA and RNA**. This yields 580 strict pairs
from 16 histotypes, not the 8,428 tumor-only patients in panel a.

Upstream pan-miR-29 VST values are standardized over the combined
tumor-plus-normal universe within cancer. For MRC, each of the five
gene VST measurements is standardized, their sample-wise mean is
calculated, and the composite is standardized again. The cached name
`ERC_Z` refers to this MRC score.

The effect is the mean paired difference
\(\Delta=Z_{\mathrm{tumor}}-Z_{\mathrm{normal}}\). A one-sample t test of
these patient differences is equivalent to a paired t test. Intervals
are 95% t intervals. BH correction is applied across 16 cancers
**separately for each readout**. The effect is in the underlying
standardized-score units, not Cohen's \(d_z\).

Bottom diamonds show REML/Knapp–Hartung summaries over 16 histotypes
and 13 organs. For the organ analysis, kidney and lung histotypes are
first pooled within organ using the preserved random-effects procedure.


```r
one_sample_per_patient <- function(D, sample_type) {
  X <- D[D$Type == sample_type, , drop = FALSE]
  if (nrow(X) == 0) {
    return(X)
  }
  X <- X[order(X$Patient, X$Sample), , drop = FALSE]
  X <- X[!duplicated(X$Patient), , drop = FALSE]
  X
}

run_histotype_meta <- function(D, readout) {
  X <- D[D$Readout == readout, , drop = FALSE]
  X <- X[match(MASTER_CANCERS16, X$Cancer), , drop = FALSE]
  m <- metafor::rma(yi = Effect, sei = SE, data = X, method = "REML", test = "knha", slab = Cancer)
  stopifnot(m$k == 16L)
  m
}

build_organ_effects <- function(D, readout) {
  X <- D[D$Readout == readout, , drop = FALSE]
  X$Organ <- unname(ORGAN_MAP[X$Cancer])
  if (anyNA(X$Organ)) {
    stop("Missing organ mapping.")
  }
  SPLIT <- split(X, X$Organ)
  O <- do.call(rbind, lapply(names(SPLIT), function(org) {
    Z <- SPLIT[[org]]
    if (nrow(Z) == 1) {
      return(data.frame(Organ = org, Histotypes = Z$Cancer, N_histotypes = 1L, Effect = Z$Effect,
        SE = Z$SE, CI_low = Z$Effect - 1.96 * Z$SE, CI_high = Z$Effect + 1.96 * Z$SE, tau2_within = 0,
        Q_p_within = NA_real_, stringsAsFactors = FALSE))
    }
    m_org <- metafor::rma(yi = Effect, sei = SE, data = Z, method = "REML", test = "z")
    data.frame(Organ = org, Histotypes = paste(Z$Cancer, collapse = "/"), N_histotypes = nrow(Z),
      Effect = as.numeric(m_org$b[1]), SE = m_org$se, CI_low = m_org$ci.lb, CI_high = m_org$ci.ub,
      tau2_within = m_org$tau2, Q_p_within = m_org$QEp, stringsAsFactors = FALSE)
  }))
  rownames(O) <- NULL
  O
}

run_organ_meta <- function(D, readout) {
  O <- build_organ_effects(D, readout)
  m <- metafor::rma(yi = Effect, sei = SE, data = O, method = "REML", test = "knha", slab = Organ)
  stopifnot(m$k == 13L)
  list(effects = O, meta = m)
}

weighted_kde2d <- function(x, y, weights, n = 200, h_x = NULL, h_y = NULL, expand = 0.1) {
  ok <- is.finite(x) & is.finite(y) & is.finite(weights) & weights > 0
  x <- x[ok]
  y <- y[ok]
  weights <- weights[ok]
  weights <- weights/sum(weights)
  if (is.null(h_x)) {
    h_x <- 1.1 * stats::bw.nrd0(x)
  }
  if (is.null(h_y)) {
    h_y <- 1.1 * stats::bw.nrd0(y)
  }
  if (!is.finite(h_x) || h_x <= 0) {
    h_x <- sd(x) * 0.3
  }
  if (!is.finite(h_y) || h_y <= 0) {
    h_y <- sd(y) * 0.3
  }
  xr <- range(x)
  yr <- range(y)
  xpad <- diff(xr) * expand
  ypad <- diff(yr) * expand
  gx <- seq(xr[1] - xpad, xr[2] + xpad, length.out = n)
  gy <- seq(yr[1] - ypad, yr[2] + ypad, length.out = n)
  Z <- matrix(0, nrow = length(gx), ncol = length(gy))
  norm_const <- 1/(2 * pi * h_x * h_y)
  for (i in seq_along(x)) {
    kx <- exp(-0.5 * ((gx - x[i])/h_x)^2)
    ky <- exp(-0.5 * ((gy - y[i])/h_y)^2)
    Z <- Z + weights[i] * outer(kx, ky)
  }
  Z <- Z * norm_const
  GRID <- expand.grid(x = gx, y = gy)
  GRID$density <- as.vector(Z)
  GRID$density_scaled <- GRID$density/max(GRID$density, na.rm = TRUE)
  list(grid = GRID, h_x = h_x, h_y = h_y)
}
```


```r
message("Preparing paired expression panels e-f")
ef_final <- read_data("FIG6E_FINAL16_miR29_ERC_ge5pairs/FIG6E_FINAL16_miR29_ERC_ge5pairs.rds")
prepared <- setNames(lapply(CANCERS16, function(ca) read_data(file.path("FIG6E_TCGA_miR29_ERC/cache",
  paste0(ca, "_FIG6E_v1.rds")))), CANCERS16)
for (ca in CANCERS16) {
  x <- prepared[[ca]]
  check(paste0("MRC genes: ", ca), identical(names(attr(x$erc, "ERC_ID")), MRC_GENES))
}
ef_env <- new.env(parent = globalenv())
ef_env$MASTER_CANCERS16 <- CANCERS16
ef_env$ORGAN_MAP <- ef_final$settings$ORGAN_MAP
load_helpers(c("one_sample_per_patient", "run_histotype_meta", "build_organ_effects", "run_organ_meta",
  "weighted_kde2d"), ef_env, "paired_expression_and_density.R")
pair_tables <- lapply(CANCERS16, function(ca) {
  x <- prepared[[ca]]
  mT <- ef_env$one_sample_per_patient(x$mir, "01")
  mN <- ef_env$one_sample_per_patient(x$mir, "11")
  rT <- ef_env$one_sample_per_patient(x$erc, "01")
  rN <- ef_env$one_sample_per_patient(x$erc, "11")
  patients <- sort(Reduce(intersect, list(mT$Patient, mN$Patient, rT$Patient, rN$Patient)))
  d <- data.frame(Cancer = ca, Patient = patients, Pan29_Normal = mN$Pan29_Z[match(patients, mN$Patient)],
    Pan29_Tumor = mT$Pan29_Z[match(patients, mT$Patient)], ERC_Normal = rN$ERC_Z[match(patients,
      rN$Patient)], ERC_Tumor = rT$ERC_Z[match(patients, rT$Patient)])
  d$Delta_Pan29 <- d$Pan29_Tumor - d$Pan29_Normal
  d$Delta_ERC <- d$ERC_Tumor - d$ERC_Normal
  saved <- ef_final$paired_patient_level[ef_final$paired_patient_level$Cancer == ca, ]
  saved <- saved[match(d$Patient, saved$Patient), ]
  check(paste0("e/f exact paired scores: ", ca), equal_numeric(as.matrix(d[, 3:8]), as.matrix(saved[,
    3:8])))
  d
})
ef_pairs <- do.call(rbind, pair_tables)
check("e/f: 580 strict pairs", nrow(ef_pairs) == 580L)
e_effects <- do.call(rbind, lapply(CANCERS16, function(ca) do.call(rbind, lapply(c("pan-miR-29", "ERC"),
  function(readout) {
    x <- ef_pairs[ef_pairs$Cancer == ca, if (readout == "ERC")
      "Delta_ERC" else "Delta_Pan29"]
    tt <- t.test(x, mu = 0)
    data.frame(Cancer = ca, Readout = readout, Effect = mean(x), Lower = tt$conf.int[1], Upper = tt$conf.int[2],
      SE = sd(x)/sqrt(length(x)), N_pairs = length(x), P = tt$p.value)
  }))))
e_effects$FDR <- ave(e_effects$P, e_effects$Readout, FUN = function(p) p.adjust(p, "BH"))
for (col in c("Effect", "Lower", "Upper", "SE", "N_pairs", "P", "FDR")) check(paste0("e exact ", col),
  equal_numeric(e_effects[[col]], ef_final$paired[[col]]))
e_meta <- do.call(rbind, lapply(c("pan-miR-29", "ERC"), function(r) {
  h <- ef_env$run_histotype_meta(e_effects, r)
  o <- ef_env$run_organ_meta(e_effects, r)$meta
  do.call(rbind, lapply(c("Histotype", "Organ"), function(level) {
    m <- if (level == "Histotype")
      h else o
    data.frame(Readout = r, Level = level, K = m$k, Effect = as.numeric(m$b), CI_low = m$ci.lb, CI_high = m$ci.ub,
      P = m$pval)
  }))
}))
for (col in c("Effect", "CI_low", "CI_high", "P")) check(paste0("e exact meta ", col), equal_numeric(e_meta[[col]],
  subset(ef_final$meta_summary, Analysis == "Paired")[[col]]))
f_means <- merge(subset(e_effects, Readout == "pan-miR-29", c(Cancer, Effect)), subset(e_effects, Readout ==
  "ERC", c(Cancer, Effect)), by = "Cancer", suffixes = c("_miRNA", "_MRC"))
f_test <- cor.test(f_means$Effect_miRNA, f_means$Effect_MRC, method = "pearson")
check("f Pearson coefficient", abs(unname(f_test$estimate) - (-0.8063488)) < 1e-07)
check("f Pearson p-value", abs(f_test$p.value - 0.0001613243) < 1e-10)
f_weights <- 1/as.numeric(table(ef_pairs$Cancer)[ef_pairs$Cancer])
f_kde <- ef_env$weighted_kde2d(ef_pairs$Delta_Pan29, ef_pairs$Delta_ERC, f_weights, n = 200)
check("f equal cancer KDE weight", length(unique(round(tapply(f_weights, ef_pairs$Cancer, sum), 12))) ==
  1L)
```


```r
# Generic compact forest. Both the uncertainty intervals and pooling are unchanged.
forest_plot <- function(d, meta, title, xlab, limits, show_labels = TRUE, expression_panel = FALSE) {
  rows <- if (expression_panel)
    CANCERS16 else CANCERS15
  d$row <- length(rows) + 1 - match(d$Cancer, rows)
  d$alpha <- if (expression_panel)
    ifelse(d$FDR < 0.05, 1, 0.35) else ifelse(d$CI_low * d$CI_high > 0, 1, 0.4)
  d$colour <- ifelse(d$Effect >= 0, "#BD3D4B", "#355EA6")
  meta$row <- c(-0.5, -1.6)
  diamonds <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    data.frame(x = c(meta$CI_low[i], meta$Effect[i], meta$CI_high[i], meta$Effect[i]), y = meta$row[i] +
      c(0, 0.25, 0, -0.25), group = i, colour = ifelse(meta$Effect[i] >= 0, "#BD3D4B", "#355EA6"))
  }))
  p <- ggplot() + annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F2F8FB") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FCF5F6") + geom_vline(xintercept = 0,
    linetype = "dashed", colour = "grey65", linewidth = 0.35) + geom_segment(data = d, aes(x = CI_low,
    xend = CI_high, y = row, yend = row), colour = "#777777", linewidth = 0.4) + geom_point(data = d,
    aes(Effect, row, colour = colour, alpha = alpha), size = 1.8) + geom_hline(yintercept = 0.35,
    linewidth = 0.35) + geom_polygon(data = diamonds, aes(x, y, group = group, fill = colour), linewidth = 0.25) +
    scale_colour_identity() + scale_fill_identity() + scale_alpha_identity() + scale_y_continuous(breaks = c(-1.6,
    -0.5, seq_along(rows)), labels = c("Organ", "Histotype", rev(rows)), limits = c(-2.3, length(rows) +
    0.7), expand = c(0, 0)) + coord_cartesian(xlim = limits) + labs(x = xlab, y = NULL, title = title) +
    panel_theme(10) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 11, face = "plain"))
  if (!show_labels)
    p <- p + theme(axis.text.y = element_blank())
  p
}
ed <- PANELS$e$effects
em <- PANELS$e$meta
ed$CI_low <- ed$Lower
ed$CI_high <- ed$Upper
pe1 <- forest_plot(subset(ed, Readout == "pan-miR-29"), subset(em, Readout == "pan-miR-29"), "Σ(miR-29)",
  "Paired difference (SD units)", c(-3.2, 2.1), TRUE, TRUE)
pe2 <- forest_plot(subset(ed, Readout == "ERC"), subset(em, Readout == "ERC"), "MRC", "Paired difference (SD units)",
  c(-2, 4.1), FALSE, TRUE)
pe <- combine_panels(list(pe1, pe2), c(1, 1), "e | Matched tumor-normal expression")
save_panel(pe, "Figure5e_paired_expression", 6.7, 5.4)
```

<img src="figures/reconstructed/Figure5e_paired_expression.png" width="100%" style="display: block; margin: auto;" />


|Readout    |Level     |  K|   Effect|   CI_low|  CI_high|       P|
|:----------|:---------|--:|--------:|--------:|--------:|-------:|
|pan-miR-29 |Histotype | 16| -0.72143| -1.22867| -0.21420| 0.00841|
|pan-miR-29 |Organ     | 13| -0.77559| -1.29910| -0.25208| 0.00725|
|MRC        |Histotype | 16|  0.99749|  0.51349|  1.48149| 0.00052|
|MRC        |Organ     | 13|  1.03521|  0.49487|  1.57555| 0.00129|

Full per-cancer effects, nominal p values, q values and patient-level
differences are exported in the `Figure5e_*` and
`Figure5ef_paired_patients.tsv` tables.

## f — association of paired miR-29 and MRC changes

Each labeled point is the mean paired effect for one of the 16 cancers,
using the same patients as panel e. The statistic is **Pearson correlation
across the 16 cancer means**, not a patient-level or Spearman correlation.

The background is a two-dimensional density of all 580 patient pairs.
Each pair has weight \(1/n_c\), so each cancer contributes the same total
weight despite unequal numbers of pairs. Contours show density, not a
confidence region. The line is a descriptive ordinary least-squares fit
to the cancer means. Figure S5i examines omission of each cancer.
The calculation and density helper are displayed in panel e.


```r
# f: original weighted KDE; Pearson coefficient and p-value are kept together.
f <- PANELS$f
pf <- ggplot() + geom_contour_filled(data = f$kde$grid, aes(x, y, z = density_scaled), breaks = c(0.05,
  0.15, 0.3, 0.5, 0.75, Inf), alpha = 0.85) + scale_fill_manual(values = c("#F4F8F4", "#E1EEDB", "#C9DFCF",
  "#ACCED0", "#86B8C2"), guide = "none") + geom_smooth(data = f$means, aes(Effect_miRNA, Effect_MRC),
  method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, colour = "#555555") + geom_point(data = f$means,
  aes(Effect_miRNA, Effect_MRC), colour = "#D84147", size = 1.9) + ggrepel::geom_text_repel(data = f$means,
  aes(Effect_miRNA, Effect_MRC, label = Cancer), seed = 29, size = 3, segment.size = 0.25, min.segment.length = 0,
  box.padding = 0.23, max.overlaps = Inf) + annotate("text", x = 1.8, y = 4.2, hjust = 1, label = sprintf("Pearson r = %.2f\np = %.4f",
  unname(f$pearson$estimate), f$pearson$p.value), size = 3.4) + coord_cartesian(xlim = c(-3.5, 2.1),
  ylim = c(-2.3, 4.6)) + labs(x = "Σ(miR-29) Δ(Tumor - Normal)", y = "MRC difference (Tumor - Normal)",
  title = "f | Matched tumor-normal pairs") + panel_theme(10)
save_panel(pf, "Figure5f_paired_scatter_Pearson", 5.3, 5)
```

<img src="figures/reconstructed/Figure5f_paired_scatter_Pearson.png" width="100%" style="display: block; margin: auto;" />


| Cancers| Paired_patients| Pearson_r|   CI_low|  CI_high|        P|
|-------:|---------------:|---------:|--------:|--------:|--------:|
|      16|             580| -0.806349| -0.93023| -0.51749| 0.000161|

## g — stable normal references and fidelity loss

The reference is defined separately for each eligible tumor histotype
using adjacent-normal methylation. The upstream pipeline restricts to
autosomal probes, removes probes exceeding 5% missing values and applies
KNN imputation as implemented in the archived workflow. CpGs must satisfy
all three normal-sample conditions:

- beta standard deviation \(\leq0.05\);
- absolute Pearson correlation with chronological age \(\leq0.25\);
- absolute linear beta change per decade \(\leq0.01\).

These are **stability-selected, tissue-matched references**, not
tissue-specific marker CpGs. Passing cutoffs in the same normal samples
does not establish independent age invariance. CpG counts and normal
sample numbers appear in S5j.

KICH uses a shared normal-kidney reference from **160 KIRC plus 45 KIRP
normal samples**. KIRC and KIRP use their own normal references in the
final analysis. CHOL uses its recovered reference with nine normal
samples. STAD is not included in the final fidelity analysis.

Beta values are clipped to \([10^{-6},1-10^{-6}]\), then converted to
\(M=\log_2[\beta/(1-\beta)]\). Each reference M value is the median of
the normal-sample M values. For tumor \(s\),
\[
\mathrm{FidelityLoss}_s =
\sqrt{\frac{1}{P}\sum_{i=1}^{P}
 (M_{is}-\operatorname{median}_{n\in\mathrm{normal}}M_{in})^2}.
\]
Higher RMS deviation means greater reference deviation, or lower
operational epigenetic fidelity. It is not DNAmAge, a percentage,
a signed methylation change or the 43-CpG fibroblast score in Figure 4.
Cellular composition and other tumor-associated changes can contribute
to bulk-tissue reference deviation.

These extracted functions specify matrix-level selection and scoring.
They are **not executed during knitting**, because full normal/tumor
beta matrices are not included. The reconstruction starts from saved
sample scores and checks selected-reference metrics. The prerequisite
`BN` is already QC-filtered, imputed and autosomal; these functions
do not perform that preprocessing.


```r
# Optional matrix-level reference derivation; not run during document rendering.  BN must contain
# prefiltered, imputed autosomal beta values from normal samples.
MAX_NORMAL_SD_BETA <- 0.05
MAX_ABS_AGE_COR <- 0.25
MAX_AGE_DBETA_DECADE <- 0.01
MIN_STABLE_CPGS <- 5000
ESCAPE_BETA_FLOOR <- 0.05
DBETA_DIRECTIONAL <- 0.1
BLOCK_SIZE <- 50
beta_to_M <- function(x, eps = 1e-06) {
  x <- pmin(pmax(x, eps), 1 - eps)
  log2(x/(1 - x))
}

derive_stable_reference <- function(BN, normal_age) {
  if (ncol(BN) != length(normal_age)) {
    stop("Normal methylation matrix and age vector have different lengths.")
  }

  normal_sd_beta <- matrixStats::rowSds(BN)
  is_stable <- normal_sd_beta <= MAX_NORMAL_SD_BETA

  age_centered <- normal_age - mean(normal_age)
  age_ss <- sum(age_centered^2)
  if (!is.finite(age_ss) || age_ss <= 0)
    stop("Normal age has no usable variance.")

  age_num <- as.vector(BN %*% age_centered)
  age_delta_beta_decade <- 10 * age_num/age_ss

  denom_cor <- sqrt(((ncol(BN) - 1) * normal_sd_beta^2) * age_ss)

  age_cor <- age_num/denom_cor
  age_cor[!is.finite(age_cor)] <- 0

  is_age_independent <- abs(age_cor) <= MAX_ABS_AGE_COR & abs(age_delta_beta_decade) <= MAX_AGE_DBETA_DECADE

  stable <- is_stable & is_age_independent
  stable_cpgs <- rownames(BN)[stable]

  if (length(stable_cpgs) < MIN_STABLE_CPGS) {
    stop("Too few stable age-independent CpGs: ", length(stable_cpgs))
  }

  BNs <- BN[stable_cpgs, , drop = FALSE]

  ref_beta <- matrixStats::rowMedians(BNs)
  ref_beta_mad <- matrixStats::rowMads(BNs, center = ref_beta, constant = 1.4826)
  ref_M <- matrixStats::rowMedians(beta_to_M(BNs))

  names(ref_beta) <- stable_cpgs
  names(ref_beta_mad) <- stable_cpgs
  names(ref_M) <- stable_cpgs

  list(stable_CpGs = stable_cpgs, ref_beta = ref_beta, ref_beta_mad = ref_beta_mad, ref_M = ref_M,
    normal_sd_beta = normal_sd_beta[stable], age_cor = age_cor[stable], age_delta_beta_decade = age_delta_beta_decade[stable])
}


score_against_reference <- function(B, ann, reference, block_size = BLOCK_SIZE) {
  cpgs <- intersect(reference$stable_CpGs, rownames(B))
  if (length(cpgs) < MIN_STABLE_CPGS) {
    stop("Too few reference CpGs available in target matrix: ", length(cpgs))
  }

  B <- B[cpgs, , drop = FALSE]
  ref_beta <- reference$ref_beta[cpgs]
  ref_beta_mad <- reference$ref_beta_mad[cpgs]
  ref_M <- reference$ref_M[cpgs]
  escape_width <- pmax(3 * ref_beta_mad, ESCAPE_BETA_FLOOR)

  n <- ncol(B)
  RMS_M <- MeanAbsM <- MeanAbsBeta <- MeanSignedBeta <- EscapeFraction <- Hyper10 <- Hypo10 <- rep(NA_real_,
    n)

  starts <- seq(1, n, by = block_size)

  for (s in starts) {
    jj <- s:min(s + block_size - 1, n)
    BB <- B[, jj, drop = FALSE]

    dBeta <- sweep(BB, 1, ref_beta, "-")
    dM <- sweep(beta_to_M(BB), 1, ref_M, "-")

    RMS_M[jj] <- sqrt(colMeans(dM^2))
    MeanAbsM[jj] <- colMeans(abs(dM))
    MeanAbsBeta[jj] <- colMeans(abs(dBeta))
    MeanSignedBeta[jj] <- colMeans(dBeta)
    EscapeFraction[jj] <- colMeans(abs(dBeta) > escape_width)
    Hyper10[jj] <- colMeans(dBeta >= DBETA_DIRECTIONAL)
    Hypo10[jj] <- colMeans(dBeta <= -DBETA_DIRECTIONAL)

    rm(BB, dBeta, dM)
  }

  score_df <- data.frame(Stable_FidelityLoss_RMS_M = RMS_M, Stable_EscapeFraction = EscapeFraction,
    Stable_MeanAbsDeltaM = MeanAbsM, Stable_MeanAbsDeltaBeta = MeanAbsBeta, Stable_MeanSignedDeltaBeta = MeanSignedBeta,
    Stable_HyperFraction10 = Hyper10, Stable_HypoFraction10 = Hypo10, row.names = colnames(B))

  cbind(ann, score_df[ann$MethID, , drop = FALSE])
}
```

<img src="reference/Figure5g_reference_artwork.png" width="100%" style="display: block; margin: auto;" />

This schematic is preserved artwork. The numerical score and selection
rules, not the cartoon geometry, define the method.

## h — miR-29, MRC and age-adjusted reference deviation

The forests have different assay-overlap cohorts. The miRNA analysis
requires primary-tumor fidelity, pan-miR-29, age and a finite Horvath
prediction. Horvath availability is an **inclusion rule**, not an
adjustment covariate here. The MRC analysis requires MRC, fidelity and
age but not miRNA or Horvath. One primary sample per patient is retained.
Original sample counts are checked separately; the cohorts are not
silently replaced by their intersection.

The fidelity miRNA checkpoint uses DESeq2 size-factor-normalized counts:
the miR-29a, combined miR-29b and miR-29c abundances are summed, then
transformed with log2(x + 1). This is not the VST composite used for the
paired expression scores. The MRC scores retain their five-gene
standardization described in panel e.

Within cancer, each variable is regressed on a three-degree-of-freedom
natural spline of age; Spearman correlation is then calculated between
the residual vectors. This is correlation of residualized values, not
a rank-first partial correlation. The historical nominal p value is
the asymptotic `cor.test(..., exact=FALSE)` result. Fisher intervals use
\(SE_z=1/\sqrt{N-6}\), accounting approximately for the three spline
degrees of freedom. These are preserved approximation conventions,
not an exact finite-sample permutation test.

Pooling uses REML/Knapp–Hartung on Fisher-z correlations, with
**15-histotype and 12-organ** summaries. Unlike panel e, the first stage
for the fidelity organ summary is inverse-variance weighting within organ.


```r
tcga_sample_id_local <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    paste(z[1:4], collapse = "-")
  }, character(1))
}

tcga_patient_id_local <- function(x) {
  substr(tcga_sample_id_local(x), 1, 12)
}

tcga_type_code_local <- function(x) {
  sid <- tcga_sample_id_local(x)
  vapply(strsplit(sid, "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    substr(z[4], 1, 2)
  }, character(1))
}

prepare_exact_clock_overlap <- function(cancer) {
  P <- allPred450K[[cancer]]
  if (is.null(P)) {
    stop("No allPred450K object for ", cancer)
  }
  P <- as.data.frame(P)
  age_idx <- which(colnames(P) == "AGE")
  hor_idx <- which(colnames(P) == "horvath")
  if (length(age_idx) == 0) {
    stop(cancer, ": AGE column not found in allPred450K")
  }
  if (length(hor_idx) == 0) {
    stop(cancer, ": horvath column not found in allPred450K")
  }
  age_idx <- age_idx[1]
  hor_idx <- hor_idx[1]
  A <- data.frame(MethID = rownames(P), SampleID = tcga_sample_id_local(rownames(P)), PatientID = tcga_patient_id_local(rownames(P)),
    TypeCode = tcga_type_code_local(rownames(P)), AGE = suppressWarnings(as.numeric(P[, age_idx])),
    Horvath = suppressWarnings(as.numeric(P[, hor_idx])), stringsAsFactors = FALSE)
  A <- A[A$TypeCode == "01" & is.finite(A$AGE) & is.finite(A$Horvath) & !is.na(A$SampleID), , drop = FALSE]
  A <- A[!duplicated(A$SampleID), , drop = FALSE]
  A <- A[order(A$PatientID, A$SampleID), , drop = FALSE]
  A <- A[!duplicated(A$PatientID), , drop = FALSE]
  A
}

build_analysis_data <- function(cancer) {
  Fobj <- TCGA_fidelity_main15[[cancer]]
  Mobj <- TCGA_miR29[[cancer]]
  if (is.null(Fobj) || is.null(Fobj$status) || Fobj$status != "OK") {
    stop("Invalid fidelity object for ", cancer)
  }
  if (is.null(Mobj) || is.null(Mobj$status) || Mobj$status != "OK") {
    stop("Invalid miRNA object for ", cancer)
  }
  F <- Fobj$scores[Fobj$scores$SampleType == "PrimaryTumor", c("SampleID", "Stable_FidelityLoss_RMS_M"),
    drop = FALSE]
  names(F)[names(F) == "Stable_FidelityLoss_RMS_M"] <- "FidelityLoss"
  M <- Mobj$table[, c("SampleID", "pan_miR29"), drop = FALSE]
  A <- prepare_exact_clock_overlap(cancer)
  A <- A[, c("SampleID", "PatientID", "AGE", "Horvath"), drop = FALSE]
  D <- merge(F, M, by = "SampleID")
  D <- merge(D, A, by = "SampleID")
  D <- D[complete.cases(D[, c("pan_miR29", "FidelityLoss", "AGE", "Horvath")]), , drop = FALSE]
  D <- D[order(D$PatientID, D$SampleID), , drop = FALSE]
  D <- D[!duplicated(D$PatientID), , drop = FALSE]
  rownames(D) <- NULL
  D
}

run_spline_residualized_spearman <- function(cancer) {
  D <- build_analysis_data(cancer)
  N <- nrow(D)
  if (N <= 6) {
    stop(cancer, ": insufficient N for nonlinear residualized correlation.")
  }
  fit_pan <- lm(pan_miR29 ~ ns(AGE, df = 3), data = D)
  fit_fid <- lm(FidelityLoss ~ ns(AGE, df = 3), data = D)
  pan_resid <- residuals(fit_pan)
  fid_resid <- residuals(fit_fid)
  ct <- suppressWarnings(cor.test(pan_resid, fid_resid, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  rho_clip <- min(max(rho, -0.999999), 0.999999)
  z <- atanh(rho_clip)
  se_z <- 1/sqrt(N - 6)
  ci_low <- tanh(z - qnorm(0.975) * se_z)
  ci_high <- tanh(z + qnorm(0.975) * se_z)
  data.frame(Cancer = cancer, N = N, rho = rho, CI_low = ci_low, CI_high = ci_high, pvalue = ct$p.value,
    Fisher_z = z, Fisher_SE = se_z, stringsAsFactors = FALSE)
}
```


```r
tcga_sample_id_mrc <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    paste(z[1:4], collapse = "-")
  }, character(1))
}

tcga_patient_id_mrc <- function(x) {
  substr(tcga_sample_id_mrc(x), 1, 12)
}

tcga_type_code_mrc <- function(x) {
  sid <- tcga_sample_id_mrc(x)
  vapply(strsplit(sid, "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    substr(z[4], 1, 2)
  }, character(1))
}

get_age_table_mrc <- function(cancer) {
  P <- allPred450K[[cancer]]
  if (is.null(P)) {
    stop(cancer, ": allPred450K missing.")
  }
  P <- as.data.frame(P)
  age_idx <- which(colnames(P) == "AGE")
  if (length(age_idx) == 0) {
    stop(cancer, ": AGE column not found.")
  }
  age_idx <- age_idx[1]
  D <- data.frame(SampleID = tcga_sample_id_mrc(rownames(P)), PatientID = tcga_patient_id_mrc(rownames(P)),
    TypeCode = tcga_type_code_mrc(rownames(P)), AGE = suppressWarnings(as.numeric(P[, age_idx])),
    stringsAsFactors = FALSE)
  D <- D[D$TypeCode == "01" & is.finite(D$AGE) & !is.na(D$PatientID), , drop = FALSE]
  D <- D[order(D$PatientID, D$SampleID), , drop = FALSE]
  D <- D[!duplicated(D$PatientID), , drop = FALSE]
  D[, c("PatientID", "AGE"), drop = FALSE]
}

get_mrc_table <- function(cancer) {
  PC <- FIG6E_FINAL$perCancer[[cancer]]
  if (is.null(PC)) {
    stop(cancer, ": absent from FIG6E_FINAL$perCancer.")
  }
  X <- if (!is.null(PC$mrc)) {
    PC$mrc
  } else {
    PC$erc
  }
  if (is.null(X)) {
    stop(cancer, ": MRC/ERC table missing.")
  }
  X <- as.data.frame(X)
  value_candidates <- c("MRC_Z", "MRC_score", "MRC_Score", "MRC", "ERC_Z", "ERC_score", "ERC_Score",
    "ERC_raw")
  value_col <- value_candidates[value_candidates %in% names(X)][1]
  if (length(value_col) == 0 || is.na(value_col)) {
    stop(cancer, ": MRC/ERC score column not identified.")
  }
  if ("Patient" %in% names(X)) {
    PatientID <- as.character(X$Patient)
  } else if ("PatientID" %in% names(X)) {
    PatientID <- as.character(X$PatientID)
  } else if ("Sample" %in% names(X)) {
    PatientID <- substr(as.character(X$Sample), 1, 12)
  } else if ("SampleID" %in% names(X)) {
    PatientID <- substr(as.character(X$SampleID), 1, 12)
  } else {
    stop(cancer, ": no patient/sample identifier in MRC table.")
  }
  if ("Type" %in% names(X)) {
    TypeCode <- as.character(X$Type)
  } else if ("TypeCode" %in% names(X)) {
    TypeCode <- as.character(X$TypeCode)
  } else if ("Sample" %in% names(X)) {
    TypeCode <- tcga_type_code_mrc(X$Sample)
  } else if ("SampleID" %in% names(X)) {
    TypeCode <- tcga_type_code_mrc(X$SampleID)
  } else {
    stop(cancer, ": cannot identify tumor rows in MRC table.")
  }
  OUT <- data.frame(PatientID = PatientID, TypeCode = TypeCode, MRC = suppressWarnings(as.numeric(X[[value_col]])),
    stringsAsFactors = FALSE)
  OUT <- OUT[OUT$TypeCode == "01" & is.finite(OUT$MRC) & !is.na(OUT$PatientID), , drop = FALSE]
  OUT <- OUT[order(OUT$PatientID), , drop = FALSE]
  OUT <- OUT[!duplicated(OUT$PatientID), , drop = FALSE]
  OUT[, c("PatientID", "MRC"), drop = FALSE]
}

get_fidelity_table_mrc <- function(cancer) {
  Fobj <- TCGA_fidelity_main15[[cancer]]
  if (is.null(Fobj) || is.null(Fobj$status) || Fobj$status != "OK") {
    stop(cancer, ": invalid fidelity object.")
  }
  S <- as.data.frame(Fobj$scores)
  S <- S[S$SampleType == "PrimaryTumor", , drop = FALSE]
  if (!"Stable_FidelityLoss_RMS_M" %in% names(S)) {
    stop(cancer, ": Stable_FidelityLoss_RMS_M missing.")
  }
  if ("PatientID" %in% names(S)) {
    PatientID <- as.character(S$PatientID)
  } else {
    PatientID <- substr(as.character(S$SampleID), 1, 12)
  }
  OUT <- data.frame(PatientID = PatientID, FidelityLoss = suppressWarnings(as.numeric(S$Stable_FidelityLoss_RMS_M)),
    stringsAsFactors = FALSE)
  OUT <- OUT[is.finite(OUT$FidelityLoss) & !is.na(OUT$PatientID), , drop = FALSE]
  OUT <- OUT[!duplicated(OUT$PatientID), , drop = FALSE]
  OUT
}

build_mrc_fidelity_data <- function(cancer) {
  M <- get_mrc_table(cancer)
  F <- get_fidelity_table_mrc(cancer)
  A <- get_age_table_mrc(cancer)
  D <- merge(M, F, by = "PatientID")
  D <- merge(D, A, by = "PatientID")
  D <- D[complete.cases(D[, c("MRC", "FidelityLoss", "AGE")]), , drop = FALSE]
  rownames(D) <- NULL
  D
}

run_mrc_fidelity_cancer <- function(cancer) {
  D <- build_mrc_fidelity_data(cancer)
  N <- nrow(D)
  if (N <= 6) {
    stop(cancer, ": insufficient sample size.")
  }
  fit_mrc <- lm(MRC ~ ns(AGE, df = 3), data = D)
  fit_fid <- lm(FidelityLoss ~ ns(AGE, df = 3), data = D)
  mrc_res <- residuals(fit_mrc)
  fid_res <- residuals(fit_fid)
  ct <- suppressWarnings(cor.test(mrc_res, fid_res, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  rho_clip <- min(max(rho, -0.999999), 0.999999)
  fisher_z <- atanh(rho_clip)
  fisher_se <- 1/sqrt(N - 6)
  CI_low <- tanh(fisher_z - qnorm(0.975) * fisher_se)
  CI_high <- tanh(fisher_z + qnorm(0.975) * fisher_se)
  data.frame(Cancer = cancer, N = N, rho = rho, CI_low = CI_low, CI_high = CI_high, pvalue = ct$p.value,
    Fisher_z = fisher_z, Fisher_SE = fisher_se, stringsAsFactors = FALSE)
}
```


```r
message("Rebuilding both final age-adjusted fidelity forests")
h_env <- new.env(parent = globalenv())
load_helpers(c("tcga_sample_id_local", "tcga_patient_id_local", "tcga_type_code_local", "prepare_exact_clock_overlap",
  "build_analysis_data", "run_spline_residualized_spearman"), h_env, "fidelity_mirna.R")
load_helpers(c("tcga_sample_id_mrc", "tcga_patient_id_mrc", "tcga_type_code_mrc", "get_age_table_mrc",
  "get_mrc_table", "get_fidelity_table_mrc", "build_mrc_fidelity_data", "run_mrc_fidelity_cancer"),
  h_env, "fidelity_mrc.R")
h_env$allPred450K <- read_data("rstudio-export-7/allPred450K.RDS")
h_env$FIG6E_FINAL <- list(perCancer = prepared)
h_env$TCGA_miR29 <- setNames(vector("list", 15), CANCERS15)
h_env$TCGA_fidelity_main15 <- h_env$TCGA_miR29
reference_sets <- list()
reference_qc <- list()
old_cache <- "miR29_TCGA_epigenetic_fidelity_MAIN14/cache_final_autosomal_fidelity_MAIN14_N10_v1"
new_cache <- "miR29_TCGA_epigenetic_fidelity_FINAL16FLOW/cache_final_autosomal_fidelity_16flow_ge5_v2"
for (ca in CANCERS15) {
  dir <- if (ca == "CHOL")
    new_cache else old_cache
  h_env$TCGA_miR29[[ca]] <- read_data(file.path(dir, paste0("miR29_", ca, ".rds")))
  if (ca == "KICH") {
    x <- read_data(file.path(old_cache, "kidney3_shared_KIRC_KIRP_reference.rds"))
    fid <- list(status = x$status, scores = x$scores$KICH, stable_CpGs = x$stable_CpGs, reference = x$reference)
    reference_qc[[ca]] <- data.frame(Cancer = ca, Reference = "Shared KIRC+KIRP normal", N_CpGs = length(x$stable_CpGs))
  } else {
    file <- if (ca %in% c("KIRC", "KIRP"))
      paste0("sensitivity_", ca, "_own_reference.rds") else if (ca == "CHOL")
      "fidelity_own_CHOL_ge5pairs.rds" else paste0("fidelity_own_", ca, ".rds")
    fid <- read_data(file.path(dir, file))
    reference_qc[[ca]] <- data.frame(Cancer = ca, Reference = "Own adjacent normal", N_CpGs = length(fid$stable_CpGs))
  }
  check(paste0("Fidelity cache status: ", ca), identical(fid$status, "OK") && identical(h_env$TCGA_miR29[[ca]]$status,
    "OK"))
  reference_sets[[ca]] <- fid$stable_CpGs
  h_env$TCGA_fidelity_main15[[ca]] <- fid[c("status", "scores")]
}
h_left <- do.call(rbind, lapply(CANCERS15, h_env$run_spline_residualized_spearman))
h_right <- do.call(rbind, lapply(CANCERS15, h_env$run_mrc_fidelity_cancer))
check("h left: exact historical cohorts", identical(as.integer(h_left$N), as.integer(EXPECTED_H_LEFT_N[h_left$Cancer])))
check("h right: exact recovered cohorts", identical(as.integer(h_right$N), as.integer(EXPECTED_H_RIGHT_N[h_right$Cancer])))
h_organ_map <- setNames(CANCERS15, CANCERS15)
h_organ_map[c("KICH", "KIRC", "KIRP")] <- "Kidney"
h_organ_map[c("LUAD", "LUSC")] <- "Lung"
fidelity_meta <- function(d) {
  hist <- metafor::rma(yi = d$Fisher_z, sei = d$Fisher_SE, method = "REML", test = "knha")
  organs <- do.call(rbind, lapply(split(d, h_organ_map[d$Cancer]), function(x) {
    w <- 1/x$Fisher_SE^2
    data.frame(z = sum(w * x$Fisher_z)/sum(w), se = sqrt(1/sum(w)))
  }))
  organ <- metafor::rma(yi = organs$z, sei = organs$se, method = "REML", test = "knha")
  do.call(rbind, lapply(c("Histotype", "Organ"), function(level) {
    m <- if (level == "Histotype")
      hist else organ
    data.frame(Level = level, K = m$k, rho = tanh(as.numeric(m$b)), CI_low = tanh(m$ci.lb), CI_high = tanh(m$ci.ub),
      P = m$pval)
  }))
}
hl_meta <- fidelity_meta(h_left)
hr_meta <- fidelity_meta(h_right)
check("h: 15 histotypes / 12 organs", identical(as.integer(hl_meta$K), c(15L, 12L)) && identical(as.integer(hr_meta$K),
  c(15L, 12L)))
check("h left: approved historical pooled rho", abs(hl_meta$rho[1] - (-0.1451)) < 1e-04)
check("h left: approved historical pooled p", abs(hl_meta$P[1] - 0.0188) < 1e-04)
```


```r
# h: each forest retains its original (different) assay-availability cohort.
hl <- PANELS$h$left
hl$Effect <- hl$rho
hr <- PANELS$h$right
hr$Effect <- hr$rho
hm1 <- PANELS$h$left_meta
hm1$Effect <- hm1$rho
hm2 <- PANELS$h$right_meta
hm2$Effect <- hm2$rho
ph1 <- forest_plot(hl, hm1, "Σ(miR-29) ~ ΔM", "Age-adjusted Spearman's ρ", c(-0.55, 0.55))
ph2 <- forest_plot(hr, hm2, "MRC ~ ΔM", "Age-adjusted Spearman's ρ", c(-0.55, 0.85), FALSE)
ph <- combine_panels(list(ph1, ph2), c(1, 1), "h | Epigenetic fidelity loss")
save_panel(ph, "Figure5h_fidelity_associations", 6.7, 5.4)

saveRDS(list(a = pa, b = pb, c = pc, e = pe, f = pf, h = ph), file.path(BUNDLE, "data", "Figure5_ggplots.rds"))
```

<img src="figures/reconstructed/Figure5h_fidelity_associations.png" width="100%" style="display: block; margin: auto;" />


|Readout    |Level     |  K|       rho|    CI_low|   CI_high|        P|
|:----------|:---------|--:|---------:|---------:|---------:|--------:|
|pan-miR-29 |Histotype | 15| -0.145050| -0.258084| -0.028098| 0.018817|
|pan-miR-29 |Organ     | 12| -0.167357| -0.273702| -0.056970| 0.006800|
|MRC        |Histotype | 15|  0.178639|  0.084109|  0.269983| 0.001258|
|MRC        |Organ     | 12|  0.164717|  0.070284|  0.256223| 0.002853|



|Cancer | N_miRNA| N_MRC|
|:------|-------:|-----:|
|BLCA   |     409|   406|
|BRCA   |     773|   781|
|CHOL   |      36|    35|
|COAD   |     293|   292|
|ESCA   |     184|   162|
|HNSC   |     522|   501|
|KICH   |      66|    65|
|KIRC   |     318|   318|
|KIRP   |     272|   271|
|LIHC   |     371|   370|
|LUAD   |     437|   436|
|LUSC   |     357|   361|
|PRAD   |     494|   497|
|THCA   |     506|   504|
|UCEC   |     427|   427|

Per-cancer cohorts and nominal p values are saved in `tables/`.
The association does not itself establish direct maintenance of fidelity
by miR-29, or demonstrate lineage conversion.

# Figure S5 — supplementary

## a — TCGA cohort sizes and sample types

These are the exact tumor samples used for Figure 5a: 8,428 samples,
all from unique patients, comprising 8,084 primary and 344 metastatic
tumors without matched normals. BRCA has 773 primary plus 1 metastatic
sample; SKCM has 97 primary plus 343 metastatic samples. Other cancers
contribute primary tumors only. Slash labels mean **primary / metastatic**,
not samples / unique patients.


```r
# Rebuild only the analyses retained in the final Figure S5.
EXPECTED_S5 <- readRDS("data/expected_supplementary.rds")
EXPECTED_EXT <- readRDS("data/expected_clock_extensions.rds")
EXPECTED_BOX <- readRDS("data/expected_compact.rds")
```


```r
# Final Figure S5: retained panels only. Illustrator-only decorations are not redrawn.
stopifnot(requireNamespace("ggrastr", quietly = TRUE))
tissue_labels <- c(breast = "Breast", central_nervous_system = "CNS", haematopoietic_and_lymphoid_tissue = "Blood / lymphoid",
  large_intestine = "Intestine", liver = "Liver", lung = "Lung", oesophagus = "Esophagus", ovary = "Ovary",
  pancreas = "Pancreas", skin = "Skin", stomach = "Stomach")
supp_theme <- function(size = 10) panel_theme(size) + theme(plot.title = element_text(size = size, face = "plain"))
```


```r
# ---- S5a: tumor cohort accounting
annotations <- lapply(names(a_cache <- read_data("miR29_Horvath_repeated_resampling/MIR29_CLOCK_RESULTS_Horvath_AGE_Acceleration_Pearson_Spearman_v1.rds")),
  function(ca) {
    z <- a_cache[[ca]]$annotation
    z$Cancer <- ca
    z
  })
cohort <- do.call(rbind, annotations)
tcga_counts <- do.call(rbind, lapply(a_cancers, function(ca) {
  z <- cohort[cohort$Cancer == ca, ]
  data.frame(Cancer = ca, N = nrow(z), Primary = sum(z$TypeCode == "01"), Metastatic = sum(z$TypeCode ==
    "06"), UniquePatients = length(unique(z$PatientID)))
}))
check("S5a: 8428 unique tumor patients", nrow(cohort) == 8428L && length(unique(cohort$PatientID)) ==
  8428L && all(tcga_counts$N == tcga_counts$UniquePatients))
save_table(tcga_counts, "S5a_TCGA_cohort_counts")
```


```r
# ---- S5a
d <- rbind(data.frame(Cancer = tcga_counts$Cancer, Type = "Primary", N = tcga_counts$Primary), data.frame(Cancer = tcga_counts$Cancer,
  Type = "Metastatic", N = tcga_counts$Metastatic))
d$Cancer <- factor(d$Cancer, levels = rev(a_cancers))
d$Type <- factor(d$Type, levels = c("Primary", "Metastatic"))
ann <- tcga_counts
ann$Cancer <- factor(ann$Cancer, levels = rev(a_cancers))
ann$Label <- ifelse(ann$Metastatic > 0, paste0(ann$Primary, " / ", ann$Metastatic), ann$Primary)
p <- ggplot(d, aes(N, Cancer, fill = Type)) + geom_col(width = 0.7) + geom_text(data = ann, inherit.aes = FALSE,
  mapping = aes(x = N, y = Cancer, label = Label), hjust = -0.12, size = 2.6, colour = "#477CB4") +
  scale_fill_manual(values = c(Primary = "#4C7FB8", Metastatic = "#DE983C")) + scale_x_continuous(expand = expansion(mult = c(0,
  0.28))) + labs(title = "TCGA: miRNA and DNAmAge", x = "Patients", y = NULL) + supp_theme(9)
save_panel(p, "S5a_TCGA_cohort", 5, 7)
```

<img src="figures/reconstructed/S5a_TCGA_cohort.png" width="100%" style="display: block; margin: auto;" />

|Cancer |   N| Primary| Metastatic| UniquePatients|
|:------|---:|-------:|----------:|--------------:|
|ACC    |  80|      80|          0|             80|
|BLCA   | 409|     409|          0|            409|
|BRCA   | 774|     773|          1|            774|
|CESC   | 307|     307|          0|            307|
|CHOL   |  36|      36|          0|             36|
|COAD   | 293|     293|          0|            293|
|DLBC   |  47|      47|          0|             47|
|ESCA   | 184|     184|          0|            184|
|HNSC   | 522|     522|          0|            522|
|KICH   |  66|      66|          0|             66|
|KIRC   | 318|     318|          0|            318|
|KIRP   | 272|     272|          0|            272|
|LGG    | 511|     511|          0|            511|
|LIHC   | 371|     371|          0|            371|
|LUAD   | 437|     437|          0|            437|
|LUSC   | 357|     357|          0|            357|
|MESO   |  87|      87|          0|             87|
|PAAD   | 178|     178|          0|            178|
|PCPG   | 179|     179|          0|            179|
|PRAD   | 494|     494|          0|            494|
|READ   |  97|      97|          0|             97|
|SARC   | 259|     259|          0|            259|
|SKCM   | 440|      97|        343|            440|
|STAD   | 383|     383|          0|            383|
|TGCT   | 134|     134|          0|            134|
|THCA   | 506|     506|          0|            506|
|THYM   | 123|     123|          0|            123|
|UCEC   | 427|     427|          0|            427|
|UCS    |  57|      57|          0|             57|
|UVM    |  80|      80|          0|             80|

## b — TCGA chronological-age and acceleration meta-ranks

The same per-cancer results and feature rules as Figure 5a are applied
to chronological age (left) and Horvath pan-tissue acceleration (right).
Each panel has 1,764 entries. All marks are individual ranked
observations, without an interpolating or smoothing curve.
Colors indicate correlation sign; the vertical line marks its zero crossing.


```r
# ---- S5b: TCGA equal-cancer Fisher meta-ranks
tcga_ranks <- lapply(c("AGE", "Horvath", "Acceleration"), function(metric) {
  z <- do.call(rbind, lapply(split(a_all, a_all$miRNA), function(d) {
    r <- d[[paste0("Spearman_", metric, "_Fisher")]]
    data.frame(miRNA = d$miRNA[1], K = sum(is.finite(r)), rho = fisher_mean(r))
  }))
  z <- z[!z$miRNA %in% redundant & is.finite(z$rho), ]
  z <- z[order(z$rho, z$miRNA), ]
  z$Rank <- seq_len(nrow(z))
  z$Feature <- unname(mir_labels[z$miRNA])
  z
})
names(tcga_ranks) <- c("AGE", "Horvath", "Acceleration")
for (metric in names(tcga_ranks)) {
  z <- tcga_ranks[[metric]]
  old <- EXPECTED_EXT$TCGA$ranks[[metric]]
  check(paste("S5b TCGA rank", metric), identical(z$miRNA, old$miRNA) && equal_numeric(z$rho, old$rho))
  save_table(z, paste0("TCGA_rank_", metric))
}
```


```r
# ---- S5b
tcga_rank_plot <- function(d, ylabel) {
  h <- d[!is.na(d$Feature), ]
  background <- d[is.na(d$Feature), ]
  h$Label <- c(`miR-29a` = "29a", `miR-29b` = "29b", `miR-29c` = "29c", `pan-miR-29` = "Σ(miR-29)")[h$Feature]
  zero <- max(d$Rank[d$rho <= 0])
  ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#EAEDF7") + annotate("rect",
    xmin = zero, xmax = Inf, ymin = 0, ymax = Inf, fill = "#FAF0F0") + geom_hline(yintercept = 0,
    linewidth = 0.4) + geom_vline(xintercept = zero, linetype = "dashed", linewidth = 0.3) + geom_point(data = background,
    aes(Rank, rho, colour = rho > 0, size = 0.35 + 4 * abs(rho)), alpha = 0.88) + scale_colour_manual(values = c(`FALSE` = "#446F9E",
    `TRUE` = "#D65242"), guide = "none") + scale_size_identity() + geom_point(data = h, aes(Rank,
    rho, fill = Feature, shape = Feature), size = 3, stroke = 0.3) + scale_fill_manual(values = COLORS,
    guide = "none") + scale_shape_manual(values = c(`miR-29a` = 21, `miR-29b` = 21, `miR-29c` = 21,
    `pan-miR-29` = 24), guide = "none") + ggrepel::geom_text_repel(data = h, aes(Rank, rho, label = Label),
    nudge_x = -320, direction = "y", hjust = 1, size = 3.2, seed = 29, min.segment.length = 0, box.padding = 0.4,
    max.overlaps = Inf) + annotate("text", x = 2040, y = -0.214, label = "N = 8,428", hjust = 1,
    size = 3.3, fontface = "italic") + scale_x_continuous(limits = c(-300, 2100), breaks = c(0, 500,
    1000, 1500), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.225, 0.265), breaks = seq(-0.2,
    0.2, 0.1), expand = c(0, 0)) + labs(title = "Cross-cancer meta-rank", x = "miRNA meta-rank",
    y = ylabel) + supp_theme(10) + theme(aspect.ratio = 2.53, panel.border = element_rect(fill = NA,
    colour = "grey25", linewidth = 0.35), axis.line = element_blank())
}
p <- cowplot::plot_grid(tcga_rank_plot(tcga_ranks$AGE, "miRNA expression ~ chronological age [ρ]"),
  tcga_rank_plot(tcga_ranks$Acceleration, "miRNA expression ~ age acceleration (P) [ρ]"), nrow = 1,
  align = "hv", axis = "tblr")
save_panel(p, "S5b_TCGA_age_and_acceleration", 7.2, 8.5)
```

<img src="figures/reconstructed/S5b_TCGA_age_and_acceleration.png" width="100%" style="display: block; margin: auto;" />

## c — TCGA leave-one-cancer-out sensitivity

The pan-miR-29 versus Horvath DNAmAge summary is recomputed after omitting
each of 30 cancers. The full estimate is approximately 0.25661.
This examines dependence of the equal-cancer summary on one cancer;
it does not resample patients or provide external validation.

Horizontal segments connect omission estimates to the full estimate:
they are **not confidence intervals**. The same calculation is used for S5e.


```r
# ---- S5c/e: leave-one-stratum-out, using the original equal-stratum estimand
leave_one_stratum <- function(unit, rho, n) {
  full <- fisher_mean(rho)
  do.call(rbind, lapply(c("None", unit), function(excluded) {
    keep <- if (excluded == "None")
      rep(TRUE, length(unit)) else unit != excluded
    data.frame(Excluded = excluded, K = sum(keep), N = sum(n[keep]), rho = fisher_mean(rho[keep]),
      Full_rho = full)
  }))
}
pan_a <- a_all[a_all$miRNA == "pan_miR29", ]
pan_a <- pan_a[match(a_cancers, pan_a$Cancer), ]
loo_tcga <- leave_one_stratum(pan_a$Cancer, pan_a$Spearman_Horvath_Fisher, pan_a$N)
pan_b <- subset(b_final$PerTissue, miRNA == "pan_miR29" & Phenotype == "Skin")
pan_b <- pan_b[match(b_final$settings$EligibleTissues, pan_b$Tissue), ]
loo_ccle <- leave_one_stratum(pan_b$Tissue, pan_b$rho_Fisher_resampled, pan_b$N)
check("S5c: exact 31 TCGA LOO estimates", equal_numeric(loo_tcga$rho, EXPECTED_EXT$TCGA$loo$rho))
check("S5e: exact 12 CCLE LOO estimates", equal_numeric(loo_ccle$rho, EXPECTED_EXT$CCLE$loo$rho))
save_table(loo_tcga, "S5c_TCGA_leave_one_out")
save_table(loo_ccle, "S5e_CCLE_leave_one_out")
```


```r
# ---- S5c/e/i: horizontal segments show movement from the full estimate, not uncertainty
loo_plot_final <- function(d, title, labels = NULL, colour = "#396896", xlabel = "Equal-stratum Fisher mean of Spearman ρ",
  full_label = "All cancers") {
  if (is.null(labels))
    labels <- setNames(d$Excluded[-1], d$Excluded[-1])
  yl <- c(None = full_label, setNames(paste("Omit", labels[d$Excluded[-1]]), d$Excluded[-1]))
  d$Row <- factor(d$Excluded, levels = rev(d$Excluded))
  ggplot(d, aes(rho, Row)) + geom_vline(xintercept = d$Full_rho[1], linetype = "dashed", colour = "grey50",
    linewidth = 0.4) + geom_segment(aes(x = Full_rho, xend = rho, yend = Row), colour = colour, linewidth = 0.45,
    alpha = 0.7) + geom_point(aes(shape = Excluded == "None"), colour = colour, fill = "white", size = 2.5) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 23), guide = "none") + scale_y_discrete(labels = yl) +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) + labs(title = title, x = xlabel,
    y = NULL) + supp_theme(9)
}
p <- loo_plot_final(loo_tcga, "Leave-one-cancer-out: Σ(miR-29) ~ DNAmAge (P)")
save_panel(p, "S5c_TCGA_leave_one_out", 6, 7.5)
```

<img src="figures/reconstructed/S5c_TCGA_leave_one_out.png" width="100%" style="display: block; margin: auto;" />

## d — eligible CCLE lineages

Counts come from the complete-case metadata used in the final analysis,
not all available CCLE lines. They sum to 429 across 11 lineages.
The eligibility threshold is 15, irrespective of the "N20" filename.


```r
# ---- S5d/f: eligible cell lines and four CCLE ranks
ccle_counts <- data.frame(Tissue = b_tissues, N = vapply(b_tissues, function(ti) sum(b_final$Meta$TISSUE ==
  ti), integer(1)))
check("S5d: 429 lines in 11 eligible lineages", sum(ccle_counts$N) == 429L && nrow(ccle_counts) == 11L)
ccle_ranks <- lapply(c("Age", "Skin", "Horvath", "HorvathAcceleration"), function(metric) {
  z <- b_final$Rankings[[metric]]
  z <- data.frame(miRNA = z$miRNA, K = z$K, rho = z$rho_Fisher_equalTissue)
  z <- z[order(z$rho, z$miRNA), ]
  z$Rank <- seq_len(nrow(z))
  z$Feature <- unname(b_labels[z$miRNA])
  z
})
names(ccle_ranks) <- c("Age", "Skin", "Horvath", "HorvathAcceleration")
for (metric in names(ccle_ranks)) {
  old <- EXPECTED_EXT$CCLE$ranks[[metric]]
  z <- ccle_ranks[[metric]]
  check(paste("S5f CCLE rank", metric), identical(z$miRNA, old$miRNA) && equal_numeric(z$rho, old$rho))
  save_table(z, paste0("S5f_CCLE_rank_", metric))
}
save_table(ccle_counts, "S5d_CCLE_cohort_counts")
```


```r
# ---- S5d
d <- ccle_counts
d$Tissue <- factor(d$Tissue, levels = b_tissues)
p <- ggplot(d, aes(Tissue, N)) + geom_col(fill = "#58A397", width = 0.7) + geom_text(aes(label = N),
  vjust = -0.4, size = 3, colour = "#438F88") + scale_x_discrete(labels = tissue_labels) + scale_y_continuous(expand = expansion(mult = c(0,
  0.12))) + labs(x = NULL, y = "Eligible cell lines") + supp_theme(10) + theme(axis.text.x = element_text(angle = 90,
  hjust = 1, vjust = 0.5))
save_panel(p, "S5d_CCLE_cohort", 6, 3.4)
```

<img src="figures/reconstructed/S5d_CCLE_cohort.png" width="100%" style="display: block; margin: auto;" />

|Tissue                             |   N|
|:----------------------------------|---:|
|breast                             |  44|
|central_nervous_system             |  27|
|haematopoietic_and_lymphoid_tissue |  96|
|large_intestine                    |  32|
|liver                              |  17|
|lung                               | 103|
|oesophagus                         |  15|
|ovary                              |  22|
|pancreas                           |  24|
|skin                               |  33|
|stomach                            |  16|

## e — CCLE leave-one-lineage-out sensitivity

This concerns **pan-miR-29 versus Skin & Blood DNAmAge**, not
acceleration. One eligible lineage is removed at a time from the
equal-lineage Fisher summary. The full estimate is approximately 0.40883.
The 12 estimates are in `tables/S5e_CCLE_leave_one_out.tsv`.
These are cell-line lineages, not individual tumors. Segment interpretation
is the same as S5c; its shared calculation is displayed there.


```r
# ---- S5e
p <- loo_plot_final(loo_ccle, "Leave-one-lineage-out: Σ(miR-29) ~ DNAmAge (S)", tissue_labels, "#45998C",
  full_label = "All lineages")
save_panel(p, "S5e_CCLE_leave_one_out", 6, 4.5)
```

<img src="figures/reconstructed/S5e_CCLE_leave_one_out.png" width="100%" style="display: block; margin: auto;" />

## f — four complementary CCLE age ranks

Top row: donor chronological age and Skin & Blood DNAmAge.
Bottom row: Horvath pan-tissue DNAmAge and its acceleration.
Each plot has 735 entries with the miR-29 members and pan-miR-29
highlighted, using equal-lineage Fisher aggregation as in Figure 5b.
Rank calculations are shown in S5d.

No fitted line is drawn through ranks. To keep Illustrator responsive,
only the dense background point layer is rasterized at 600 dpi.
Axes, highlights and text remain vector objects in the PDF.
All coordinates are retained: rasterization is not subsampling.


```r
# ---- S5f
ccle_rank_plot <- function(d, title) {
  h <- d[!is.na(d$Feature), ]
  h$Label <- c(`miR-29a` = "29a", `miR-29b` = "29b", `miR-29c` = "29c", `pan-miR-29` = "Σ(miR-29)")[h$Feature]
  ggplot(d, aes(Rank, rho)) + geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) + ggrastr::rasterise(geom_point(colour = "#B6C5CA",
    size = 0.6), dpi = 600, dev = "cairo") + geom_point(data = h, aes(fill = Feature, shape = Feature),
    size = 2.6, colour = "grey35", stroke = 0.25) + scale_fill_manual(values = COLORS, guide = "none") +
    scale_shape_manual(values = c(`miR-29a` = 21, `miR-29b` = 21, `miR-29c` = 21, `pan-miR-29` = 23),
      guide = "none") + ggrepel::geom_text_repel(data = h, aes(label = Label), seed = 29, nudge_x = -215,
    direction = "y", size = 2.7, min.segment.length = 0, max.overlaps = Inf) + scale_x_continuous(breaks = c(0,
    250, 500, 750)) + labs(title = title, x = "miRNA meta-rank", y = "Spearman ρ") + supp_theme(10)
}
titles <- c(Age = "miRNA ~ chronological age", Skin = "miRNA ~ DNAmAge (S)", Horvath = "miRNA ~ DNAmAge (P)",
  HorvathAcceleration = "miRNA ~ age acceleration (P)")
p <- cowplot::plot_grid(plotlist = lapply(names(titles), function(m) ccle_rank_plot(ccle_ranks[[m]],
  titles[[m]])), ncol = 2, align = "hv", axis = "tblr")
save_panel(p, "S5f_CCLE_four_ranks", 8, 7.4)
```

<img src="figures/reconstructed/S5f_CCLE_four_ranks.png" width="100%" style="display: block; margin: auto;" />

## g — within-cancer RNA correlation distributions

Horizontal boxplots summarize **all eligible RNA features** within 30
cancers, using resampled Spearman correlations with Horvath DNAmAge.
These are distributions over features, not uncertainty intervals for
a pooled correlation. There are 1,673,475 feature-by-cancer observations;
feature counts can differ across cancers.

Boxes show median and first/third quartiles; whiskers reach the most
extreme observations within 1.5 interquartile ranges. All 33,264
outliers are plotted. The five MRC genes, MIR29B2CHG and LINC-PINT
produce 210 highlighted gene–cancer pairs. MRC markers are yellow,
matching the final montage.

Only the gray outlier layer is rasterized at 600 dpi. Boxes, highlights,
labels and axes remain editable vectors, reducing Illustrator object
counts without dropping outliers or changing boxplot statistics.


```r
# ---- S5g: full RNA distributions, not the old rank-sampled display
summarize_rna_boxes <- function() {
  stats <- list()
  outliers <- list()
  highlights <- list()
  for (ca in a_cancers) {
    z <- c_final$perCancer[c_final$perCancer$Cancer == ca & is.finite(c_final$perCancer$rho_Fisher_resampled),
      ]
    z$Gene <- c_final$meta$SYMBOL[match(z$FeatureID, c_final$meta$FeatureID)]
    r <- z$rho_Fisher_resampled
    q <- quantile(r, c(0.25, 0.5, 0.75), type = 7, names = FALSE)
    is_out <- r < q[1] - 1.5 * (q[3] - q[1]) | r > q[3] + 1.5 * (q[3] - q[1])
    whiskers <- range(c(q[1], q[3], r[!is_out]))
    stats[[ca]] <- data.frame(Unit = ca, N = unique(z$N), N_features = nrow(z), Q25 = q[1], Median = q[2],
      Q75 = q[3], Low = whiskers[1], High = whiskers[2], N_outliers = sum(is_out))
    outliers[[ca]] <- data.frame(Unit = ca, ID = z$FeatureID[is_out], rho = r[is_out])
    ix <- which(z$Gene %in% c(MRC_GENES, "MIR29B2CHG", "LINC-PINT"))
    highlights[[ca]] <- data.frame(Unit = ca, ID = z$FeatureID[ix], rho = r[ix], Gene = z$Gene[ix],
      Feature = ifelse(z$Gene[ix] %in% MRC_GENES, "MRC genes", z$Gene[ix]))
  }
  list(summary = do.call(rbind, stats), outliers = do.call(rbind, outliers), highlights = do.call(rbind,
    highlights))
}
rna_boxes <- summarize_rna_boxes()
check("S5g: all 30 full RNA box statistics", equal_numeric(as.matrix(rna_boxes$summary[, -1]), as.matrix(EXPECTED_BOX$TCGA_RNA$summary[,
  -1])))
check("S5g: all 33264 outliers and 210 highlights", nrow(rna_boxes$outliers) == 33264L && nrow(rna_boxes$highlights) ==
  210L)
save_table(rna_boxes$summary, "S5g_RNA_box_statistics")
save_table(rna_boxes$highlights, "S5g_RNA_highlighted_genes")
```


```r
# ---- S5g
s <- rna_boxes$summary
o <- rna_boxes$outliers
h <- rna_boxes$highlights
for (obj in c("s", "o", "h")) {
  z <- get(obj)
  z$Row <- 31 - match(z$Unit, a_cancers)
  assign(obj, z)
}
h$Row <- h$Row + seq(-0.16, 0.16, length.out = 7)[match(h$Gene, c(MRC_GENES, "MIR29B2CHG", "LINC-PINT"))]
p <- ggplot(s, aes(x = Row, group = Unit)) + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
  linewidth = 0.35) + geom_boxplot(aes(ymin = Low, lower = Q25, middle = Median, upper = Q75, ymax = High),
  stat = "identity", orientation = "x", width = 0.64, fill = "#626868", colour = "#444848", linewidth = 0.25) +
  ggrastr::rasterise(geom_point(data = o, aes(y = rho), size = 0.16, alpha = 0.45, colour = "#6C7479"),
    dpi = 600, dev = "cairo") + geom_point(data = h, aes(x = Row, y = rho, fill = Feature), shape = 21,
  size = 1.9, colour = "grey35", stroke = 0.22) + scale_fill_manual(values = c(`MRC genes` = "#E8D833",
  MIR29B2CHG = "#E2A1A1", `LINC-PINT` = "#AE8EB7")) + scale_x_continuous(breaks = 1:30, labels = rev(a_cancers),
  limits = c(0.4, 30.6), expand = c(0, 0)) + scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1,
  1), expand = expansion(mult = 0.015)) + coord_flip() + labs(x = NULL, y = "Gene expression ~ DNAmAge (P) [ρ]",
  fill = NULL) + supp_theme(9) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
  legend.position = "top", legend.text = element_text(size = 8))
save_panel(p, "S5g_TCGA_RNA_distributions", 6.7, 7.7)
```

<img src="figures/reconstructed/S5g_TCGA_RNA_distributions.png" width="100%" style="display: block; margin: auto;" />

Box summaries and highlighted values are exported as
`S5g_RNA_box_statistics.tsv` and `S5g_RNA_highlighted_genes.tsv`.

## h — strict paired tumor–normal cohort

Counts refer to the patients used in both Figure 5e and 5f.
Each must have primary-tumor and adjacent-normal data for both miRNA
and RNA. The 16 cancer counts sum to 580. Reference-normal counts
in S5j are a different cohort and cannot replace these paired counts.


```r
# ---- S5h/i: strict paired cohorts and Pearson leave-one-cancer-out
pair_counts <- as.data.frame(table(ef_pairs$Cancer), stringsAsFactors = FALSE)
names(pair_counts) <- c("Cancer", "N_pairs")
loo_paired <- do.call(rbind, lapply(c("None", CANCERS16), function(excluded) {
  z <- if (excluded == "None")
    f_means else f_means[f_means$Cancer != excluded, ]
  ct <- cor.test(z$Effect_miRNA, z$Effect_MRC, method = "pearson")
  data.frame(Excluded = excluded, K = nrow(z), rho = unname(ct$estimate), p = ct$p.value, Full_rho = unname(f_test$estimate))
}))
old <- subset(EXPECTED_S5$robustness, Method == "pearson")
check("S5i: exact Pearson LOO estimates", equal_numeric(loo_paired$rho, old$rho))
save_table(pair_counts, "S5h_pair_counts")
save_table(loo_paired, "S5i_paired_Pearson_leave_one_out")
```


```r
# ---- S5h
d <- pair_counts
d$Cancer <- factor(d$Cancer, levels = rev(CANCERS16))
p <- ggplot(d, aes(N_pairs, Cancer)) + geom_col(fill = "#9E6CAF", width = 0.7) + geom_text(aes(label = N_pairs),
  hjust = -0.15, size = 3, colour = "#86589A") + scale_x_continuous(expand = expansion(mult = c(0,
  0.2))) + labs(title = "Matched tumor-normal pairs", x = "Paired patients", y = NULL) + supp_theme(10)
save_panel(p, "S5h_paired_cohort", 4.5, 5.2)
```

<img src="figures/reconstructed/S5h_paired_cohort.png" width="100%" style="display: block; margin: auto;" />

|Cancer | N_pairs|
|:------|-------:|
|BLCA   |      19|
|BRCA   |     103|
|CHOL   |       8|
|COAD   |       8|
|ESCA   |       8|
|HNSC   |      42|
|KICH   |      24|
|KIRC   |      71|
|KIRP   |      32|
|LIHC   |      49|
|LUAD   |      19|
|LUSC   |      38|
|PRAD   |      52|
|STAD   |      27|
|THCA   |      59|
|UCEC   |      21|

## i — paired-effect association after omitting one cancer

The Pearson association between cancer-mean paired pan-miR-29 and MRC
effects is recalculated after removing each cancer. The full estimate is
\(r=-0.80635\). Only Pearson is retained in the final S5. Each omission
leaves 15 cancer means; patients in remaining cancers are not rematched.

The calculation is shown in S5h. Segments indicate movement from the
full correlation, **not confidence intervals**. Nominal p values are
exported without treating this descriptive sensitivity display as a
new discovery family.


```r
# ---- S5i
p <- loo_plot_final(loo_paired, "Leave-one-cancer-out: ΔΣ(miR-29) ~ ΔMRC", colour = "#8139A4", xlabel = "Pearson r (paired cancer means)")
save_panel(p, "S5i_paired_leave_one_out", 5.8, 5.8)
```

<img src="figures/reconstructed/S5i_paired_leave_one_out.png" width="100%" style="display: block; margin: auto;" />

## j — normal reference size and retained-CpG properties

Left: retained CpGs and normal sample counts for each reference.
KICH's asterisk means **205 shared kidney normals (160 KIRC + 45 KIRP)**,
not 205 KICH donors.

Right: the retained CpGs' beta variability and **absolute** age slope.
The point is the median, thick interval the 25th–75th percentiles,
and thin interval the 5th–95th percentiles. These are CpG-distribution
ranges, not standard errors or confidence intervals. Dashed cutoffs:
SD(beta) = 0.05 and \(|\Delta\beta|/\mathrm{decade}=0.01\).

The absolute age-correlation cutoff of 0.25 remains in selection even
though its diagnostic panel is omitted from the final montage.
All three cutoffs are checked for every selected CpG. This documents
selection properties, not independent validation of age invariance.
Large retained sets are compatible with broadly stable loci rather
than uniquely tissue-specific markers.


```r
# ---- S5j: selected-reference QC, not an independent validation set
reference_metrics <- readRDS("data/reference_selection_metrics.rds")
reference_qc <- list()
reference_distributions <- list()
for (ca in CANCERS15) {
  if (ca == "KICH") {
    rel <- file.path(old_cache, "kidney3_shared_KIRC_KIRP_reference.rds")
    x <- read_data(rel)
    qc <- data.frame(Cancer = ca, N_normal_reference = x$N_normal, N_stable_CpGs = length(x$stable_CpGs),
      N_probes_after_NA = length(x$common_CpGs), Reference = "Shared KIRC+KIRP")
  } else {
    fname <- if (ca %in% c("KIRC", "KIRP"))
      paste0("sensitivity_", ca, "_own_reference.rds") else if (ca == "CHOL")
      "fidelity_own_CHOL_ge5pairs.rds" else paste0("fidelity_own_", ca, ".rds")
    rel <- file.path(if (ca == "CHOL")
      new_cache else old_cache, fname)
    x <- read_data(rel)
    qc <- x$qc[, c("Cancer", "N_normal_reference", "N_stable_CpGs", "N_probes_after_NA")]
    qc$Reference <- "Own adjacent normal"
  }
  ref <- reference_metrics[[rel]]
  stopifnot(!is.null(ref))
  for (metric in c("normal_sd_beta", "age_cor", "age_delta_beta_decade")) {
    v <- ref[[metric]]
    if (is.null(names(v))) {
      stopifnot(length(v) == length(ref$stable_CpGs))
      names(v) <- ref$stable_CpGs
    }
    v <- abs(v[x$stable_CpGs])
    stopifnot(all(is.finite(v)))
    cutoff <- switch(metric, normal_sd_beta = 0.05, age_cor = 0.25, age_delta_beta_decade = 0.01)
    check(paste("S5j original cutoff", ca, metric), all(v <= cutoff + 1e-12))
    q <- quantile(v, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), names = FALSE)
    reference_distributions[[paste(ca, metric)]] <- data.frame(Cancer = ca, Metric = metric, Min = q[1],
      P05 = q[2], Q25 = q[3], Median = q[4], Q75 = q[5], P95 = q[6], Max = q[7], N = length(v))
  }
  reference_qc[[ca]] <- qc
}
reference_qc <- do.call(rbind, reference_qc)
reference_distributions <- do.call(rbind, reference_distributions)
check("S5j all reference counts", equal_numeric(reference_qc$N_stable_CpGs, EXPECTED_S5$reference$qc$N_stable_CpGs))
check("S5j all distribution summaries", equal_numeric(as.matrix(reference_distributions[, -c(1, 2)]),
  as.matrix(EXPECTED_S5$reference$distributions[, -c(1, 2)])))
save_table(reference_qc, "S5j_reference_counts")
save_table(reference_distributions, "S5j_reference_distributions")
rm(reference_metrics)
```


```r
# ---- S5j
d <- reference_qc
d$Cancer <- factor(d$Cancer, levels = rev(CANCERS15))
d$Nlabel <- paste0("N = ", d$N_normal_reference, ifelse(d$Cancer == "KICH", "*", ""))
p0 <- ggplot(d, aes(N_stable_CpGs/1000, Cancer)) + geom_col(fill = "#5A91A3", width = 0.7) + geom_text(aes(label = sprintf("%.1fk",
  N_stable_CpGs/1000)), hjust = -0.15, size = 2.8) + geom_text(aes(x = 4, label = Nlabel), hjust = 0,
  colour = "white", size = 2.8) + scale_x_continuous(limits = c(0, 365), breaks = c(0, 100, 200, 300),
  expand = c(0, 0)) + labs(x = "Retained CpGs (thousands)", y = NULL, title = "Normal-tissue references") +
  supp_theme(10)
ref_summary_plot <- function(metric, cutoff, xlabel) {
  z <- reference_distributions[reference_distributions$Metric == metric, ]
  z$Cancer <- factor(z$Cancer, levels = rev(CANCERS15))
  ggplot(z, aes(Median, Cancer)) + geom_vline(xintercept = cutoff, colour = "#C84242", linetype = "dashed",
    linewidth = 0.4) + geom_segment(aes(x = P05, xend = P95, yend = Cancer), colour = "#AAC3CC",
    linewidth = 0.4) + geom_segment(aes(x = Q25, xend = Q75, yend = Cancer), colour = "#5AA38C",
    linewidth = 1.3) + geom_point(size = 1.3, colour = "#526C6B") + scale_x_continuous(limits = c(0,
    cutoff * 1.08)) + labs(x = xlabel, y = NULL, title = paste("Cutoff:", cutoff)) + supp_theme(10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
}
p <- cowplot::plot_grid(p0, ref_summary_plot("normal_sd_beta", 0.05, "SD(β)"), ref_summary_plot("age_delta_beta_decade",
  0.01, "|Δβ| / decade"), nrow = 1, rel_widths = c(1.85, 1, 1), align = "h", axis = "tb")
save_panel(p, "S5j_reference_QC", 10, 6)
```

<img src="figures/reconstructed/S5j_reference_QC.png" width="100%" style="display: block; margin: auto;" />

|Cancer | N_normal_reference| N_stable_CpGs| N_probes_after_NA|Reference           |
|:------|------------------:|-------------:|-----------------:|:-------------------|
|BLCA   |                 21|        164300|            373446|Own adjacent normal |
|BRCA   |                 97|        241182|            369818|Own adjacent normal |
|CHOL   |                  9|        145573|            375859|Own adjacent normal |
|COAD   |                 38|        182067|            387050|Own adjacent normal |
|ESCA   |                 16|        110803|            374010|Own adjacent normal |
|HNSC   |                 50|        216324|            386972|Own adjacent normal |
|KICH   |                205|        286113|            372864|Shared KIRC+KIRP    |
|KIRC   |                160|        293407|            378714|Own adjacent normal |
|KIRP   |                 45|        249720|            387716|Own adjacent normal |
|LIHC   |                 50|        222508|            361787|Own adjacent normal |
|LUAD   |                 32|        243492|            386267|Own adjacent normal |
|LUSC   |                 42|        259906|            376653|Own adjacent normal |
|PRAD   |                 50|        209510|            380698|Own adjacent normal |
|THCA   |                 56|        271883|            396538|Own adjacent normal |
|UCEC   |                 34|        207399|            374732|Own adjacent normal |

## k — age acceleration and reference deviation

The retained forest tests Horvath pan-tissue acceleration versus
fidelity loss. Acceleration is recalculated as the residual of Horvath
DNAmAge on age in each exact primary-tumor overlap cohort. Both
acceleration and fidelity are then residualized against a three-df
age spline before Spearman correlation. This is explicitly an
**age-adjusted** association.

The historical residual-correlation approximation and REML/Knapp–Hartung
pooling are preserved. The pooled result is approximately
\(\rho=-0.02078\), 95% CI \([-0.10881,0.06756]\), with
**p = q = 0.6220752**. This nonsignificant summary does not establish
equivalence or rule out cancer-specific relationships; heterogeneity
between cancers is substantial.

FDR retains the original analysis family: 45 cancer-level tests
(three models × 15 cancers), and separately three pooled tests.
Displaying only one forest does not remove the other models from
correction. The two miR-29 models are recalculated and tabulated for
provenance, not reintroduced as supplementary panels.


```r
# ---- S5k: retain the historical three-model FDR family
residual_spearman <- function(d, x, y, extra = NULL) {
  terms <- paste(c("ns(AGE,df=3)", extra), collapse = " + ")
  d <- d[complete.cases(d[, c(x, y, "AGE", extra)]), ]
  rx <- residuals(lm(as.formula(paste(x, "~", terms)), d))
  ry <- residuals(lm(as.formula(paste(y, "~", terms)), d))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  se <- 1/sqrt(nrow(d) - 6 - length(extra))
  data.frame(N = nrow(d), rho = rho, low = tanh(z - 1.96 * se), high = tanh(z + 1.96 * se), p = ct$p.value,
    se = se)
}
fidelity_age <- do.call(rbind, lapply(CANCERS15, function(ca) {
  d <- h_env$build_analysis_data(ca)
  d$AgeAccel <- residuals(lm(Horvath ~ AGE, d))
  save_table(d, paste0("fidelity_miRNA_cohort_", ca))
  save_table(h_env$build_mrc_fidelity_data(ca), paste0("fidelity_MRC_cohort_", ca))
  rbind(data.frame(Cancer = ca, Contrast = "Age acceleration ~ fidelity", residual_spearman(d, "AgeAccel",
    "FidelityLoss")), data.frame(Cancer = ca, Contrast = "miR-29 ~ fidelity | age", residual_spearman(d,
    "pan_miR29", "FidelityLoss")), data.frame(Cancer = ca, Contrast = "miR-29 ~ fidelity | age + acceleration",
    residual_spearman(d, "pan_miR29", "FidelityLoss", "AgeAccel")))
}))
fidelity_age$q <- p.adjust(fidelity_age$p, "BH")
fidelity_age_meta <- do.call(rbind, lapply(split(fidelity_age, fidelity_age$Contrast), function(d) {
  fit <- metafor::rma(yi = atanh(pmin(pmax(d$rho, -0.999999), 0.999999)), sei = d$se, method = "REML",
    test = "knha")
  data.frame(Contrast = d$Contrast[1], K = fit$k, rho = tanh(as.numeric(fit$b)), low = tanh(fit$ci.lb),
    high = tanh(fit$ci.ub), p = fit$pval, I2 = fit$I2)
}))
fidelity_age_meta$q <- p.adjust(fidelity_age_meta$p, "BH")
check("S5k all three pooled models unchanged", equal_numeric(as.matrix(fidelity_age_meta[, -1]), as.matrix(EXPECTED_S5$fidelity$age_meta[,
  -1])))
save_table(fidelity_age, "fidelity_age_all_45_tests")
save_table(fidelity_age_meta, "fidelity_age_all_three_meta_tests")
save_table(h_left, "Figure5h_miRNA_fidelity")
save_table(h_right, "Figure5h_MRC_fidelity")
save_table(hl_meta, "Figure5h_miRNA_pooled")
save_table(hr_meta, "Figure5h_MRC_pooled")
save_table(e_effects, "Figure5e_paired_effects")
save_table(e_meta, "Figure5e_pooled_effects")
save_table(ef_pairs, "Figure5ef_paired_patients")
save_table(c_185, "Figure5c_185_genes")
save_table(c_rank, "Figure5c_gene_ranking")
save_table(a_summary, "Figure5a_summary")
save_table(b_summary, "Figure5b_summary")
```


```r
# ---- S5k
d <- subset(fidelity_age, Contrast == "Age acceleration ~ fidelity")
m <- subset(fidelity_age_meta, Contrast == "Age acceleration ~ fidelity")
d$Row <- factor(d$Cancer, levels = rev(c("Pooled", CANCERS15)))
m$Row <- factor("Pooled", levels = levels(d$Row))
p <- ggplot(d, aes(rho, Row)) + annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#E7EDF5") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FAE9E9") + geom_vline(xintercept = 0,
  linetype = "dashed", colour = "grey55", linewidth = 0.35) + geom_segment(aes(x = low, xend = high,
  yend = Row), linewidth = 0.45, colour = "grey50") + geom_point(aes(colour = rho > 0, alpha = q <
  0.05), size = 2) + geom_segment(data = m, aes(x = low, xend = high, yend = Row), linewidth = 0.7) +
  geom_point(data = m, shape = 18, size = 3) + scale_colour_manual(values = c(`FALSE` = "#426FAD",
  `TRUE` = "#D95C64"), guide = "none") + scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1),
  guide = "none") + scale_y_discrete(labels = c(Pooled = "Pooled", setNames(CANCERS15, CANCERS15))) +
  labs(title = sprintf("Age acceleration (P) ~ fidelity loss\nPooled q = %.3f", m$q), x = "Age-adjusted Spearman ρ",
    y = NULL) + supp_theme(10)
save_panel(p, "S5k_acceleration_fidelity", 5.3, 6.1)
```

<img src="figures/reconstructed/S5k_acceleration_fidelity.png" width="100%" style="display: block; margin: auto;" />

|Contrast                                    |  K|       rho|       low|      high|        p|       I2|        q|
|:-------------------------------------------|--:|---------:|---------:|---------:|--------:|--------:|--------:|
|Age acceleration ~ fidelity                 | 15| -0.020784| -0.108809|  0.067564| 0.622075| 88.34567| 0.622075|
|miR-29 ~ fidelity &#124; age                | 15| -0.145050| -0.258084| -0.028098| 0.018817| 93.42797| 0.028225|
|miR-29 ~ fidelity &#124; age + acceleration | 15| -0.150101| -0.260681| -0.035635| 0.014007| 92.95205| 0.028225|

The per-cancer results for all 45 tests are in
`tables/fidelity_age_all_45_tests.tsv`.

# Validation, source pages and software

## Numerical validation

The reconstruction stops on a failed check. Checks cover cohort sizes,
feature selection, paired scores/tests, Fisher ranks, fidelity cohorts,
leave-one-out estimates, full RNA box statistics, reference cutoffs and
the original three-model fidelity FDR family. Floating-point comparisons
use explicit tolerances. Historical outputs serve as regression checks;
an unexecuted raw-data pipeline is not certified.


| Checks| Passed| Main_analytical_exports| Supplementary_analytical_exports| Retained_artwork_panels|
|------:|------:|-----------------------:|--------------------------------:|-----------------------:|
|    139|    139|                       6|                               11|                       2|


```r
# Run from the root of this bundle. No network or original workstation paths are required.
for (directory in c("tables", "validation", "figures/reconstructed")) dir.create(directory, recursive = TRUE,
  showWarnings = FALSE)
for (module in c("00_setup.R", "01_main_analysis.R", "02_main_plots.R", "03_supplementary_analysis.R",
  "04_supplementary_plots.R")) {
  message("Running ", module)
  source(file.path("code", module), local = TRUE)
}
jsonlite::write_json(CHECKS, "validation/numerical_checks.json", pretty = TRUE, auto_unbox = TRUE)
writeLines(vapply(CHECKS, function(x) paste(if (x$pass) "PASS" else "FAIL", x$check, x$detail, sep = " | "),
  character(1)), "validation/numerical_checks.txt")
writeLines(capture.output(sessionInfo()), "validation/sessionInfo.txt")
saveRDS(list(tcga_counts = tcga_counts, tcga_ranks = tcga_ranks, loo_tcga = loo_tcga, ccle_counts = ccle_counts,
  ccle_ranks = ccle_ranks, loo_ccle = loo_ccle, rna_boxes = rna_boxes, pair_counts = pair_counts, loo_paired = loo_paired,
  reference_qc = reference_qc, reference_distributions = reference_distributions, fidelity_age = fidelity_age,
  fidelity_age_meta = fidelity_age_meta), "data/FigureS5_panel_data.rds")
message("Completed ", length(CHECKS), " numerical checks.")
```

## Montage notes

Submitted PDF pages remain unchanged. Reconstructed plots and this
document use the verified analytical labels and values:

| Panel | Clarification |
|:--|:--|
| 5a | miRNA meta-rank, not mRNA; TCGA abbreviations LIHC and PCPG |
| S5e | Leave-one-lineage-out for CCLE, not leave-one-tumor-out |
| S5j | Absolute beta change per decade; ranges are CpG quantiles |
| S5k | Pooled q = 0.622, not the montage's 0.662; age adjustment explicit |

These are transcription/label clarifications, not changes to the
analysis or its significance.

## Approved main and supplementary layouts

These full pages are reference images of the supplied final PDFs,
not the analytical reconstructions displayed above.

<img src="reference/Figure5_main_final.png" width="100%" style="display: block; margin: auto;" />

<img src="reference/FigureS5_final.png" width="100%" style="display: block; margin: auto;" />

## Software and files

The tested environment is recorded below and in `validation/sessionInfo.txt`.
The bundle contains the RStudio project, source R Markdown, expanded
Markdown and HTML, executable modules, checkpoints, tables, PDFs/PNGs,
unchanged source pages and validation/provenance records.

HTML is self-contained for reading. Re-execution requires the full
bundle, not the HTML or downloaded Rmd alone. Reconstructed PDFs are in
`figures/reconstructed/`. The two dense supplementary layers identified
above use 600-dpi raster backgrounds; other elements remain vectors.


```
## R version 4.2.3 (2023-03-15)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
## 
## attached base packages:
## [1] splines   stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.5.1
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.10         highr_0.10          vipor_0.4.5         ggrastr_1.0.2      
##  [5] pillar_1.9.0        bslib_0.9.0         compiler_4.2.3      RColorBrewer_1.1-3 
##  [9] jquerylib_0.1.4     mathjaxr_1.6-0      tools_4.2.3         metadat_1.2-0      
## [13] digest_0.6.31       jsonlite_1.8.4      evaluate_0.20       lifecycle_1.0.3    
## [17] tibble_3.2.1        gtable_0.3.6        nlme_3.1-162        lattice_0.20-45    
## [21] mgcv_1.8-42         pkgconfig_2.0.3     rlang_1.1.4         Matrix_1.5-3       
## [25] cli_3.6.1           ggrepel_0.9.3       yaml_2.3.7          beeswarm_0.4.0     
## [29] xfun_0.41           fastmap_1.2.0       withr_2.5.0         metafor_4.4-0      
## [33] dplyr_1.1.4         knitr_1.42          systemfonts_1.2.3   generics_0.1.3     
## [37] vctrs_0.6.5         sass_0.4.9          isoband_0.2.7       grid_4.2.3         
## [41] tidyselect_1.2.0    cowplot_1.1.1       glue_1.6.2          R6_2.5.1           
## [45] textshaping_1.0.1   fansi_1.0.4         ggbeeswarm_0.7.2    rmarkdown_2.21     
## [49] farver_2.1.1        magrittr_2.0.3      scales_1.4.0        htmltools_0.5.9    
## [53] numDeriv_2016.8-1.1 ragg_1.2.5          labeling_0.4.2      utf8_1.2.3         
## [57] cachem_1.1.0        Cairo_1.6-2
```
