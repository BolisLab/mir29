# Figure 5 and Supplementary Figure S5

## Cancer and tumour-cell associations

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure5_and_S5.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from tumour and cell-line analytical checkpoints and cohort summaries; it does not repeat the complete upstream acquisition and processing of TCGA/CCLE data.

**Inputs: `MIR29_Figure5_S5_inputs_v1.zip`**, uploaded to the unpublished, restricted Zenodo draft with reserved DOI [10.5281/zenodo.23266285](https://doi.org/10.5281/zenodo.23266285). The DOI is not an active public download until publication; reviewer access is granted separately. See the [data/installation guide](../../data/README.md) and [file manifest](../../data/zenodo-manifest.json). Extract the archive into the repository root. No individual GTEx donor ages are distributed.

**Current execution boundary:** Core processed-input replay passed 139 assertions. Full plotting/rendering and upstream raw-data calculations are separate verification levels.

**Verification:** 139 numerical assertions passed; all seven main-panel numerical objects matched expected values at tolerance 1e-8. This does not certify a complete raw-data-to-figure rerun.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure5_and_s5.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Input locations and hashes are listed in the manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `cowplot`, `ggplot2`, `ggrastr`, `ggrepel`, `grDevices`, `grid`, `gtable`, `jsonlite`, `knitr`, `matrixStats`, `metafor`, `rmarkdown`, `scales`, `splines`, `stats`, `tools`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.

