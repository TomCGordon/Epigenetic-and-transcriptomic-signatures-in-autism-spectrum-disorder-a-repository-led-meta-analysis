#!/usr/bin/env Rscript

# Download the public reference resources used by the analysis workflows.
args <- commandArgs(trailingOnly = TRUE)
value_for <- function(prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else sub(prefix, "", hit[[1]], fixed = TRUE)
}
output_root <- value_for("--output-root=", file.path(getwd(), "downloaded_public_inputs"))
overwrite <- "--overwrite" %in% args
dry_run <- "--dry-run" %in% args
retries <- as.integer(value_for("--retries=", "3"))
timeout <- as.integer(value_for("--timeout=", "7200"))
normalise_path <- function(path) gsub("[\\\\/]+", .Platform$file.sep, path)

resources <- data.frame(
  remote_url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL16nnn/GPL16686/soft/GPL16686_family.soft.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13388/soft/GPL13388_family.soft.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15207/soft/GPL15207_family.soft.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15314/soft/GPL15314_family.soft.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15314/soft/GPL15314_family.soft.gz",
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36315/matrix/GSE36315_series_matrix.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1708/annot/GPL1708.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL3nnn/GPL3427/annot/GPL3427.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6883/annot/GPL6883.annot.gz"
  ),
  destination = c(
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL10558.annot.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL16686_family.soft.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL570.annot.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL6244.annot.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL6480.annot.gz",
    "03_Required_Annotations_and_Metadata/branch_specific_gene_nomenclature/blood_expression/hgnc_complete_set.txt",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL13388_family.soft.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL15207_family.soft.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL15314_family.soft.gz",
    "03_Required_Annotations_and_Metadata/GSE36315_custom_annotation_source/GPL15314_family.soft.gz",
    "03_Required_Annotations_and_Metadata/branch_specific_gene_nomenclature/brain_expression/hgnc_complete_set.txt",
    "03_Required_Annotations_and_Metadata/GSE36315_custom_annotation_source/gencode.v19.lncRNA_transcripts.fa.gz",
    "03_Required_Annotations_and_Metadata/GSE36315_custom_annotation_source/gencode.v19.pc_transcripts.fa.gz",
    "03_Required_Annotations_and_Metadata/GSE36315_custom_annotation_source/GSE36315_series_matrix.txt.gz",
    "03_Required_Annotations_and_Metadata/data_raw_annotation/hg18_refGene.txt.gz",
    "03_Required_Annotations_and_Metadata/data_raw_annotation/hg19_refGene.txt.gz",
    "03_Required_Annotations_and_Metadata/data_raw_annotation/hg38_refGene.txt.gz",
    "03_Required_Annotations_and_Metadata/data_raw_annotation/UCSC_hg38_refGene.txt.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL13158.annot.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL1708.annot.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL3427.annot.gz",
    "03_Required_Annotations_and_Metadata/expression_platform_annotations/GPL6883.annot.gz"
  ),
  stringsAsFactors = FALSE
)

if (is.na(retries) || retries < 1L) stop("--retries must be a positive integer.", call. = FALSE)
if (is.na(timeout) || timeout < 1L) stop("--timeout must be a positive integer.", call. = FALSE)
for (i in seq_len(nrow(resources))) {
  row <- resources[i, ]
  destination <- file.path(output_root, normalise_path(row$destination))
  if (dry_run) {
    message("[dry-run] ", row$remote_url, " -> ", destination)
    next
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination) && file.info(destination)$size > 0 && !overwrite) next
  temporary <- paste0(destination, ".part")
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, old_timeout))
  ok <- FALSE
  for (attempt in seq_len(retries)) {
    ok <- tryCatch({
      utils::download.file(row$remote_url, temporary, mode = "wb", method = "libcurl")
      file.exists(temporary) && file.info(temporary)$size > 0
    }, error = function(error) FALSE)
    if (ok) break
    if (file.exists(temporary)) unlink(temporary)
    if (attempt < retries) Sys.sleep(min(30, 2^attempt))
  }
  options(timeout = old_timeout)
  if (!ok) stop("Download failed for ", row$remote_url, call. = FALSE)
  if (file.exists(destination)) unlink(destination)
  if (!file.rename(temporary, destination)) stop("Could not stage ", destination, call. = FALSE)
}
if (!dry_run) message("Reference resources downloaded successfully.")
