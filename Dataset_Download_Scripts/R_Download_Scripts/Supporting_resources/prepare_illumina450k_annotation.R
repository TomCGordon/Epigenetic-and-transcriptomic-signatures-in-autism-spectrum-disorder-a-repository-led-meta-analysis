#!/usr/bin/env Rscript

# Reconstruct the compact Illumina 450K annotation table used by the
# methylation workflows from the public Bioconductor annotation package.

args <- commandArgs(trailingOnly = TRUE)
value_for <- function(prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else sub(prefix, "", hit[[1]], fixed = TRUE)
}

output_root <- value_for(
  "--output-root=",
  file.path(getwd(), "downloaded_public_inputs")
)
dry_run <- "--dry-run" %in% args
destination <- file.path(
  output_root,
  "03_Required_Annotations_and_Metadata",
  "data_raw_annotation",
  "illumina450k_annotation_core.csv"
)

if (dry_run) {
  message("[dry-run] construct Illumina 450K annotation -> ", destination)
  quit(save = "no", status = 0)
}

required <- c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Install the required Bioconductor package(s) before running this script: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})
data(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  package = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)
annotation <- minfi::getAnnotation(
  IlluminaHumanMethylation450kanno.ilmn12.hg19
)

required_columns <- c(
  "UCSC_RefGene_Name",
  "UCSC_RefGene_Group",
  "Relation_to_Island",
  "chr",
  "pos"
)
missing_columns <- setdiff(required_columns, colnames(annotation))
if (length(missing_columns)) {
  stop(
    "Required annotation columns were absent: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

output <- data.frame(
  Name = rownames(annotation),
  annotation[, required_columns, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(output, destination, row.names = FALSE, na = "")
message("Wrote ", nrow(output), " annotation rows to ", destination)
