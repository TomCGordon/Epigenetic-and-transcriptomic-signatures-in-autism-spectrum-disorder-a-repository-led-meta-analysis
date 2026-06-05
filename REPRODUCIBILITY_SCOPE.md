# Reproducibility Scope

This code package is intended to make the public-data analyses transparent and reproducible in R.

## What Is Fully Scripted

- Public repository file staging or download where public source URLs were available.
- Gene-level methylation and expression summary generation from public matrices, count files, signal files, or WGBS/RRBS promoter summaries.
- Dataset-level Hedges' g effect-size calculation.
- Primary and sensitivity meta-analysis model fitting with `metafor`.
- FDR correction and modified Knapp-Hartung interval classification.
- Cross-omic convergence and pathway enrichment analyses.
- Main figure generation.

## What Requires External Data Availability

Some public repository files are large and are not bundled with the code. Users should either allow scripts to download them from public repositories or provide a local source folder through the environment variables listed in `docs/RUN_ORDER.md`.

Controlled-access datasets are not analysed by this code package.

## Important Boundaries

- Single-dataset results are descriptive and are not treated as pooled meta-analyses.
- GSE36315 is retained only as a custom-annotated sensitivity analysis.
- Methylation and expression data generally come from different participants; convergence is tissue-matched summary-level convergence.
- Where the repository-available source was a processed matrix, analyses begin from that processed matrix rather than from raw IDAT or raw intensity files.
