# Figure 2 and Supplementary Figure S2

## Evolutionary conservation

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure2_and_S2.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from analysis-ready miRNA and methylation correlation objects, sequence catalogs and cached phylogeny; it does not rerun the nf-core small-RNA workflows or whole-genome alignments.

**Inputs: `MIR29_Figure2_S2_inputs_v1.zip`**, last confirmed in the unpublished, restricted Zenodo draft on 9 October 2026 with reserved DOI [10.5281/zenodo.23266285](https://doi.org/10.5281/zenodo.23266285). The DOI is not an active public download until publication; reviewer access is granted separately. See the [data/installation guide](../../data/README.md) and [file manifest](../../data/zenodo-manifest.json). Extract the archive into the repository root. No individual GTEx donor ages are distributed.

**Current execution boundary:** Processed inputs and required reference PNGs are included. Raw small-RNA processing and genome alignment are separate upstream workflows.

**Verification:** Species, 76-seed ranking, methylation and bootstrap/leave-one-out numerical checks passed. This does not certify a complete raw-data-to-figure rerun.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure2_and_s2.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Input locations and hashes are listed in the manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `ape`, `cowplot`, `dplyr`, `gghalves`, `ggplot2`, `ggrastr`, `ggtext`, `ggtree`, `grDevices`, `grid`, `knitr`, `purrr`, `readr`, `rmarkdown`, `rstudioapi`, `scales`, `stats`, `stringr`, `svglite`, `tibble`, `tidyr`, `tools`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.
