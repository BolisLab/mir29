# Figure 4 and Supplementary Figure S4

## Methylome remodelling and somatic-reference distance

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure4_and_S4.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from methylation and predictor checkpoints and the 123-CpG reference extract. The reference-distance and pooled/matched statistics are recalculated by the source modules; raw IDAT preprocessing and genome-wide discovery are separate upstream steps.

**Inputs: `MIR29_Figure4_S4_inputs_v1.zip`**, uploaded to the unpublished, restricted Zenodo draft with reserved DOI [10.5281/zenodo.23266285](https://doi.org/10.5281/zenodo.23266285). The DOI is not an active public download until publication; reviewer access is granted separately. See the [data/installation guide](../../data/README.md) and [file manifest](../../data/zenodo-manifest.json). Extract the archive into the repository root. No individual GTEx donor ages are distributed.

**Current execution boundary:** Full unmodified execution still reads withheld data/primary/gtex_sample_metrics.rds. Existing aggregate windows are included but do not automatically bypass this read.

**Verification:** 123-CpG selection, HGF/hTERT subset matrices and existing aggregate GTEx windows verified. This does not certify a complete raw-data-to-figure rerun.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure4_and_s4.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Input locations and hashes are listed in the manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `Biobase`, `GenomicRanges`, `RNAAgeCalc`, `cowplot`, `data.table`, `dplyr`, `emmeans`, `ggplot2`, `grid`, `knitr`, `limma`, `lme4`, `matrixStats`, `mgcv`, `nlme`, `pasta`, `pbkrtest`, `rmarkdown`, `scales`, `scico`, `stats`, `stringr`, `tAge`, `tibble`, `tidyr`, `tools`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.

## Current statistical labels

See [CURRENT_STATISTICAL_LABELS.tsv](CURRENT_STATISTICAL_LABELS.tsv) for the approved panel-to-test mapping.

The active somatic reference contains **123 CpGs**, selected with lineage-comparison q ≤ 0.20. The HGF whole-series test is an arithmetic mean-displacement test, **not a treatment-by-time slope test**. The hTERT final-cycle and whole-series tests remain distinct. Figure 4l is unchanged.

