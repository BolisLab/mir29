# Figure 1 and Supplementary Figure S1

## Ageing and rejuvenation

[Read the analytical chapter](READABLE.md) · [R Markdown source](Figure1_and_S1.Rmd)

## What is included

This private-review snapshot includes the analytical source, the existing rendered chapter, plot assets and selected aggregate reference/statistical tables. Most plots are extracted losslessly from the original chapter; SVG plots larger than 500 kB are displayed as 2,400-pixel PNG previews for efficient browsing. This presentation conversion does not recompute results. It does not include large input matrices, participant-level metadata, original figure PDFs or internal working archives.

## Reproduction boundary

Starts from compact iPSC/development count matrices and analysis-ready ageing and mature-miRNA checkpoints. Reconstruction of those matrices from STAR gene-count files is documented but not executed by the routine notebook; raw-read alignment and clock training are separate.

**Inputs: `MIR29_Figure1_S1_inputs_v1.zip`**, uploaded to the unpublished, restricted Zenodo draft with reserved DOI [10.5281/zenodo.23266285](https://doi.org/10.5281/zenodo.23266285). The DOI is not an active public download until publication; reviewer access is granted separately. See the [data/installation guide](../../data/README.md) and [file manifest](../../data/zenodo-manifest.json). Extract the archive into the repository root. No individual GTEx donor ages are distributed.

**Current execution boundary:** Full notebook initialization still reads withheld data/epiclock_preds_GTEX.RDS. Authorised local input is required; the aggregate S1 helper can run without it.

**Verification:** Aggregate S1 helper: assertions and six panel constructions passed. This does not certify a complete raw-data-to-figure rerun.

## Intended execution, after input installation

From this module directory, in an R environment with the recorded dependencies:

```sh
Rscript render_figure1_and_s1.R
```

This regenerates the chapter from its analytical checkpoints; it does not repeat all upstream raw-data processing. Inspect the source and the module-specific scripts before execution. Input locations and hashes are listed in the manifest; do not substitute a different input or reference silently.

## Dependencies

Directly detected R dependencies (not a validated lockfile): `ComplexHeatmap`, `DESeq2`, `SingleCellExperiment`, `circlize`, `cowplot`, `data.table`, `dendextend`, `dplyr`, `edgeR`, `forcats`, `gghalves`, `ggplot2`, `ggplotify`, `ggrastr`, `ggrepel`, `grDevices`, `grid`, `knitr`, `limma`, `org.Hs.eg.db`, `png`, `ragg`, `readxl`, `rmarkdown`, `scales`, `scater`, `scico`, `slingshot`, `svglite`, `tidyr`.

Pandoc, rmarkdown and knitr are required for notebook rendering. Where present, `RECORDED_R_SESSION.txt` records the original analytical environment, not a newly tested installation. Optional Python predictor replay may require its separate recorded environment.
