# Epigenetic and transcriptomic signatures in autism spectrum disorder: a repository-led meta-analysis

This repository contains the code, documentation and supporting materials for a repository-led meta-analysis of epigenetic and transcriptomic signatures associated with autism spectrum disorder.

The study synthesises publicly available DNA methylation and gene-expression datasets across blood, brain, placenta and lymphoblastoid cell-line samples. It examines the reproducibility of molecular differences between autistic and non-autistic samples, as well as convergence across tissues, omic layers, genes and biological pathways.

The repository is being developed alongside the associated manuscript. Additional processing, analysis, validation and reporting materials will be added as they are prepared for public release.

## Repository contents

### [`Dataset_Download_Scripts/`](./Dataset_Download_Scripts/)

Contains the scripts used to retrieve and organise the publicly available datasets and reference resources included in the meta-analysis.

Equivalent implementations are provided in both Python and R. The scripts support:

- downloading datasets from NCBI GEO, ArrayExpress and institutional repositories
- downloading and preparing genomic reference resources
- checksum verification
- automatic archive extraction
- branch-specific or dataset-specific downloads
- reproducible organisation of raw data and annotation files

Detailed installation instructions, command-line options, dataset accession numbers and output structures are provided in the folder’s [README](./Dataset_Download_Scripts/README.md).

## Study components

The complete repository will include materials relating to:

- dataset identification and acquisition
- DNA methylation data processing
- gene-expression data processing
- dataset-level statistical analyses
- within-omic meta-analyses
- cross-omic and cross-tissue comparisons
- pathway and gene-set enrichment analyses
- sensitivity and robustness analyses
- manuscript figures, tables and supplementary materials

These components will be linked from this README as they are added.

## Reproducibility

The raw biological datasets are not stored directly in this repository because of their size and because they are already available through public data repositories. Instead, the provided scripts retrieve the original files from their source repositories and organise them into the directory structure used by the analysis pipeline.

## Licence

The code in this repository is released under the [MIT Licence](./LICENSE).
