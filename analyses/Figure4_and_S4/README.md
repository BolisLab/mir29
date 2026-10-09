# Figure 4 and Supplementary Figure S4

## Methylome remodelling and somatic-reference distance

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure4_and_S4.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from methylation and predictor checkpoints and the 123-CpG reference extract. The reference-distance and pooled/matched statistics are recalculated by the source modules; raw IDAT preprocessing and genome-wide discovery are separate upstream steps.

**The external input bundle is not yet distributed with this repository.** Data-availability and redistribution review must be completed before promising an executable fresh clone. Text inside the archived reading chapter describing inputs as “bundled” refers to the original full analytical bundle. No statistical analyses were rerun during this repository staging step.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure4_and_s4.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Required input locations will be supplied with the data manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `Biobase`, `GenomicRanges`, `RNAAgeCalc`, `cowplot`, `data.table`, `dplyr`, `emmeans`, `ggplot2`, `grid`, `knitr`, `limma`, `lme4`, `matrixStats`, `mgcv`, `nlme`, `pasta`, `pbkrtest`, `rmarkdown`, `scales`, `scico`, `stats`, `stringr`, `tAge`, `tibble`, `tidyr`, `tools`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.

## Current statistical labels

See [CURRENT_STATISTICAL_LABELS.tsv](CURRENT_STATISTICAL_LABELS.tsv) for the approved panel-to-test mapping.

The active somatic reference contains **123 CpGs**, selected with lineage-comparison q ≤ 0.20. The HGF whole-series test is an arithmetic mean-displacement test, **not a treatment-by-time slope test**. The hTERT final-cycle and whole-series tests remain distinct. Figure 4l is unchanged.

