# Data access and analytical starting points

The analyses use five figure-specific input packages. This repository separates
code and documentation from large numerical checkpoints and from data subject
to access restrictions. A readable notebook is not, by itself, evidence that a
fresh checkout can execute every analysis.

The five packages were last confirmed in an **unpublished, restricted Zenodo
draft on 9 October 2026**, with reserved DOI
[10.5281/zenodo.23266285](https://doi.org/10.5281/zenodo.23266285).
The DOI becomes active on publication; no public download was available at that
check. The deposit's current publication status has not been rechecked for this update.
Editors and reviewers require separately granted draft-preview access. A secret
preview link must never be committed here. This deposit contains data, not code.

The RNA-seq and DNA-methylation data generated in this study have also been
deposited in GEO. Accessions are withheld from this public-facing repository
until release; reviewer links are provided separately to the journal.

The [input catalogue](input-catalogue.json) records module-level access and
reproduction boundaries. The [file manifest](zenodo-manifest.json) lists every
archive member, its size, purpose and SHA-256. The
[checksum file](zenodo-SHA256SUMS.txt) is an unchanged copy of the deposit's
`SHA256SUMS.txt`; its filenames refer to the original files downloaded from Zenodo.

## Packages and installation

| Package | Size (MB) | RDS objects |
|---|---:|---:|
| `MIR29_Figure1_S1_inputs_v1.zip` | 83.8 | 8 |
| `MIR29_Figure2_S2_inputs_v1.zip` | 35.1 | 11 |
| `MIR29_Figure3_S3_inputs_v1.zip` | 265.9 | 25 |
| `MIR29_Figure4_S4_inputs_v1.zip` | 184.4 | 36 |
| `MIR29_Figure5_S5_inputs_v1.zip` | 344.2 | 56 |

The deposit also contains `README.txt`, `MANIFEST.json` and `SHA256SUMS.txt`.
Download the desired packages and these three documentation files. Verify the
checksums in the download directory (with all files present):

```sh
shasum -a 256 -c SHA256SUMS.txt
```

Extract each ZIP into the **root of this repository**, not inside a figure's
`data/` folder. The archive already supplies `analyses/FigureN_and_SN/data/...`
and any missing reference-image paths. No source code is overwritten. Compare
hashes before replacing existing input files.

Version v1 was matched to analysis code at commit
`ada517267dce6d2665710d79f3d73b54f0cba50b`; subsequent documentation-only changes
do not change that numerical compatibility. Input bytes were preserved and
verified against the source files; package installation is not a full raw-data
reproduction test.

Modules are located under `analyses/Figure1_and_S1` through
`analyses/Figure5_and_S5`. Input paths in the catalogue are relative to their
module directory, not this top-level `data/` directory.

## What can be reproduced

| Module | Supported analytical starting point | Execution boundary |
|---|---|---|
| Figure 1 + S1 | Compact iPSC/development count matrices, metadata and cross-tissue summary coefficients | The current full notebook reads withheld `epiclock_preds_GTEX.RDS`; authorised local GTEx inputs are required. The aggregate S1 helper runs without it. |
| Figure 2 + S2 | Species/seed coefficients, phylogeny, methylation inputs and sequence references | Numerical input set packaged, including required reference PNGs; raw small-RNA processing and genome searches are separate. |
| Figure 3 + S3 | Expression/signature, clock, activity, growth and single-cell checkpoints | Numerical inputs packaged. If `digest` is installed, the renderer automatically checks two original reference PDFs not in the deposit; this check fails without them, even after successful HTML rendering. |
| Figure 4 + S4 | Frozen methylation/expression inputs and the current 123-CpG reference | Full unmodified replay reads withheld `gtex_sample_metrics.rds`. Included aggregate windows do not automatically bypass that read. |
| Figure 5 + S5 | TCGA/CCLE correlations, paired scores and reference distances | Core numerical replay passed 139 assertions; upstream normalization, predictors and resampling remain frozen inputs. |

The checkpoint workflows recalculate documented panel-level statistics and
plots. They do **not** automatically repeat raw read/IDAT processing, original
normalization, all resampling, clock training, or manual figure assembly.
Upstream recipes and their additional inputs are identified separately.

## Access-controlled data

Exact GTEx donor ages and other controlled-access participant metadata must not
be committed to this repository or uploaded as release assets, including to a
private repository. Their use requires the relevant dbGaP authorization and an
appropriate computing environment. Derived sample-level tables can retain the
same restriction: a precomputed object is not necessarily an aggregate object.

The data packages exclude individual GTEx donor ages, including the two
sample-level clock objects identified above. GTEx inputs in the packages are
aggregate results, not individual donor-age tables. Other public GEO, TCGA and
CCLE/DepMap inputs can retain individual sample annotations; the packages are
not universally aggregate-only. Source-specific attribution and data-use terms
remain applicable. Study-level dbGaP authorization does not grant public
redistribution rights for controlled GTEx metadata.

## Sequence and small-RNA resources

The Figure 2 sequence catalogues cover 266 precursor records, their 266 mature
miRNAs, the 257 clustered precursor subset, and 60 genomic fragments around
cg10501210. The catalogues retain assembly, chromosome/scaffold, coordinates,
sequence and locus annotations where available. Counts describe records, not
unique sequences or necessarily distinct species.

The small-RNA workflow distinguishes eight raw-read reprocessed datasets from
three sources used as published processed matrices. The recorded upstream
version is nf-core/smrnaseq 2.4.1. Genome references, exact adapter settings,
read-length settings and source accessions are documented per species; the
historical settings should not be replaced with assumed uniform defaults.

## Release checklist

1. Finalize attribution, source-specific terms and release licences.
2. Resolve or retain explicit documentation of the GTEx/PDF execution boundaries.
3. Test full chapter rendering in the documented environments.
4. Activate the deposited DOI and public downloads only after author approval.
5. Report checkpoint reproducibility separately from raw-data reproducibility.

No automatic download or authentication is performed by this catalogue.
