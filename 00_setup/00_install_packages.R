#!/usr/bin/env Rscript

# Install the R packages used by the public-data analysis scripts.
# This script is intentionally explicit so the GitHub repository remains easy
# to reproduce on a clean machine.

cran_packages <- c(
  "data.table", "readr", "readxl", "xml2", "openxlsx",
  "metafor", "ggplot2", "dplyr", "tidyr", "forcats",
  "stringr", "scales", "viridis", "ggrepel", "ggalluvial",
  "igraph", "ggraph", "tidygraph", "patchwork"
)

bioc_packages <- c(
  "GEOquery", "limma", "edgeR", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "ReactomePA", "msigdbr", "fgsea",
  "gprofiler2", "Biostrings", "Rbowtie", "huex10sttranscriptcluster.db"
)

install_if_missing <- function(pkgs, repos = getOption("repos")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = repos)
}

install_if_missing(cran_packages, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package installation/check complete.")

