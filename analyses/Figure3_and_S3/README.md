# Figure 3 and Supplementary Figure S3

## Transcriptional responses to miR-29 perturbation

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure3_and_S3.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from growth counts, RNA-count/VST matrices, gene sets, activity and clock predictions, network coordinates and compact single-cell summaries; it does not rerun read alignment, single-cell integration or external clock training.

**Inputs: `MIR29_Figure3_S3_inputs_v1.zip`**, last confirmed in the unpublished, restricted Zenodo draft on 9 October 2026 with reserved DOI [10.5281/zenodo.23266285](https://doi.org/10.5281/zenodo.23266285). The DOI is not an active public download until publication; reviewer access is granted separately. See the [data/installation guide](../../data/README.md) and [file manifest](../../data/zenodo-manifest.json). Extract the archive into the repository root. No individual GTEx donor ages are distributed.

**Current execution boundary:** Processed numerical inputs are included. When `digest` is installed, the standard renderer automatically runs a post-render hash check requiring two original reference PDFs not in the package. Without those PDFs, the command will fail at this check even if HTML rendering has completed.

**Verification:** Dependency and serialized-input inspection complete; full chapter replay not run in package preparation. This does not certify a complete raw-data-to-figure rerun.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure3_and_s3.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Input locations and hashes are listed in the manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `BiocParallel`, `DESeq2`, `SummarizedExperiment`, `circlize`, `digest`, `emmeans`, `fgsea`, `ggalluvial`, `ggforce`, `ggplot2`, `grDevices`, `grid`, `gridExtra`, `knitr`, `limma`, `lme4`, `matrixStats`, `minpack.lm`, `png`, `rmarkdown`, `scales`, `tools`, `utils`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.

## Current statistical labels

See [CURRENT_STATISTICAL_LABELS.tsv](CURRENT_STATISTICAL_LABELS.tsv) for the approved panel-to-test mapping.

Existing reconstruction plots may retain older labels; consult that mapping for current global pooled/matched q values. Updating the plots is a remaining release check, not an analysis rerun claimed here. The source narrative uses the corrected donor age of 60 years.
