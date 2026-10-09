# Data access and analytical starting points

The analyses use five figure-specific input packages. This repository separates
code and documentation from large numerical checkpoints and from data subject
to access restrictions. A readable notebook is not, by itself, evidence that a
fresh checkout can execute every analysis.

The machine-readable [input catalogue](input-catalogue.json) records the current
availability and reproduction boundary of each package. Download locations will
be added after the distribution files have been reviewed and deposited; no
public data release is claimed at this stage.

Modules are located under `analyses/Figure1_and_S1` through
`analyses/Figure5_and_S5`. Input paths in the catalogue are relative to their
module directory, not this top-level `data/` directory.

## What can be reproduced

| Module | Supported analytical starting point | Preparation still required for a distributable release |
|---|---|---|
| Figure 1 + S1 | Compact iPSC/development count matrices, aligned metadata and cross-tissue summary coefficients | Package the recovered inputs; replace restricted GTEx sample-level clock/age input with a permitted aggregate checkpoint, or document authorized local execution |
| Figure 2 + S2 | Species coefficient tables, seed references, phylogenetic annotations, methylation coefficients and sequence catalogues | Package the recovered runtime inputs; verify that current corrected species inputs and final sequence catalogues remain authoritative |
| Figure 3 + S3 | Verified expression/signature, clock, activity, growth and single-cell plotting checkpoints | Package the recovered numerical inputs and preserve the current global pooled/matched statistical annotations |
| Figure 4 + S4 | Frozen methylation and expression checkpoints, including the current 123-CpG reference | Deposit an approved checkpoint package; preserve the distinction between mean-displacement and late-endpoint tests |
| Figure 5 + S5 | Cancer/CCLE correlation, paired-expression and reference-distance checkpoints | Review and package the recovered processed inputs; document which upstream calculations are frozen |

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

Only approved aggregate results should be distributed here. Study-level
dbGaP authorization does not grant public redistribution rights. Human
sample-level records from other resources also require a source-specific
sharing check before release.

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

1. Review redistribution rights and remove restricted sample-level information.
2. Freeze the input package for each module and record its checksum.
3. Add a persistent, version-specific download location to the catalogue.
4. Test each documented entry point from a clean directory with those inputs.
5. Report checkpoint reproducibility separately from raw-data reproducibility.

No automatic download or authentication is performed by this catalogue.
