# Software environments

The analysis modules were developed in R with chapter-specific dependencies.
Their recorded session information is evidence of the original environment;
it is not replaced with an untested lockfile generated from current package versions.

The documentation build requires Python 3 (standard library) and Quarto.
It has been prepared using Quarto 1.2.335. The documentation renderer does not
execute the analysis code. Analytical R dependencies and optional Python clock
environments will be pinned after reconstruction of the matching input bundles.

Before the first public release:

1. Recover session/package versions for each chapter.
2. Separate R/Bioconductor, optional Python predictors and Nextflow/container workflows.
3. Create and test appropriate environment locks or container references.
4. Record successful execution checks, hardware requirements and analysis runtimes.

Do not treat documentation-render success as evidence that these steps have been completed.
