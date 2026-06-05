#!/usr/bin/env Rscript

# Check that the main reproducibility packages are available and record versions.

packages <- c(
  "data.table", "readr", "readxl", "xml2", "openxlsx",
  "GEOquery", "limma", "edgeR", "AnnotationDbi", "org.Hs.eg.db",
  "metafor", "clusterProfiler", "ReactomePA", "msigdbr", "fgsea",
  "gprofiler2", "Biostrings", "Rbowtie", "huex10sttranscriptcluster.db",
  "ggplot2", "dplyr", "tidyr", "forcats", "stringr", "scales",
  "viridis", "ggrepel", "ggalluvial", "igraph", "ggraph",
  "tidygraph", "patchwork"
)

versions <- data.frame(
  package = packages,
  available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)

dir.create("manifests", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(versions, file.path("manifests", "R_package_versions_check.csv"), row.names = FALSE)

print(versions)
if (any(!versions$available)) {
  stop("Some required packages are missing. See manifests/R_package_versions_check.csv")
}

message("Environment check passed.")

