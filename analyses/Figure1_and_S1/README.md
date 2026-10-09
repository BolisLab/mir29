# Figure 1 and Supplementary Figure S1

## Ageing and rejuvenation

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure1_and_S1.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from selected STAR gene-count files and analysis-ready developmental, ageing and mature-miRNA checkpoints; it does not repeat all raw-read alignment or clock training.

**The external input bundle is not yet distributed with this repository.** Data-availability and redistribution review must be completed before promising an executable fresh clone. Text inside the archived reading chapter describing inputs as “bundled” refers to the original full analytical bundle. No statistical analyses were rerun during this repository staging step.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure1_and_s1.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Required input locations will be supplied with the data manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `ComplexHeatmap`, `DESeq2`, `SingleCellExperiment`, `circlize`, `cowplot`, `data.table`, `dendextend`, `dplyr`, `edgeR`, `forcats`, `gghalves`, `ggplot2`, `ggplotify`, `ggrastr`, `ggrepel`, `grDevices`, `grid`, `knitr`, `limma`, `org.Hs.eg.db`, `png`, `ragg`, `readxl`, `rmarkdown`, `scales`, `scater`, `scico`, `slingshot`, `svglite`, `tidyr`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.

