# R Package Roles

The scripts use standard R and Bioconductor packages commonly used in public omics reanalysis and meta-analysis workflows.

## Repository Retrieval And File Handling

- `GEOquery`: GEO metadata, series matrix and supplementary file retrieval.
- `data.table`, `readr`, `readxl`, `openxlsx`, `xml2`: large table, workbook and XML parsing.

## Expression Processing

- `limma`: microarray-oriented expression data handling and annotation-compatible workflows.
- `edgeR`: count-based RNA-seq processing, including TMM/logCPM handling where count data were available.
- `AnnotationDbi`, `org.Hs.eg.db`, `huex10sttranscriptcluster.db`: gene identifier and platform annotation mapping.

## Methylation And Genomic Annotation

- `GEOquery`: methylation repository source retrieval.
- `data.table`: large methylation table processing.
- `Biostrings`, `Rbowtie`: sequence-based probe remapping for the labelled GSE36315 sensitivity workflow.

## Meta-Analysis

- `metafor`: production DerSimonian-Laird random-effects meta-analysis and modified Knapp-Hartung interval estimation.

## Pathway Enrichment And Convergence

- `clusterProfiler`, `ReactomePA`, `msigdbr`, `fgsea`, `gprofiler2`: over-representation, Reactome/MSigDB and ranked enrichment analyses.

## Figures

- `ggplot2`, `dplyr`, `tidyr`, `forcats`, `stringr`, `scales`, `viridis`, `ggrepel`, `ggalluvial`, `igraph`, `ggraph`, `tidygraph`, `patchwork`: plotting, network summaries, alluvial plots and figure layout.

For final manuscript references, use `citation("package")` in R for the packages reported in the Methods and supplement.

