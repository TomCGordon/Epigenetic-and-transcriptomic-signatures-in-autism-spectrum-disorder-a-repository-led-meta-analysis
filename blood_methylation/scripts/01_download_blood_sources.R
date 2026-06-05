#!/usr/bin/env Rscript

# Download the public repository files used by the blood methylation workflow.
# The large GSE140730 WGBS CpG reports are handled by
# 03_process_GSE140730_wgbs_from_geo_cpg_reports.R because they are streamed
# sample-by-sample rather than downloaded as one monolithic archive.

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(Sys.getenv("ASD_REPO_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
raw_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_RAW_DIR", unset = file.path(package_root, "data_raw")),
                         winslash = "/", mustWork = FALSE)
qc_dir <- normalizePath(file.path(package_root, "qc"), winslash = "/", mustWork = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

geo_series_bucket <- function(accession) {
  n <- as.integer(sub("^GSE", "", accession))
  paste0("GSE", floor(n / 1000), "nnn")
}

download_if_needed <- function(url, dest, force = FALSE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && !force) {
    return(data.table(url = url, file = dest, status = "already_present",
                      bytes = file.info(dest)$size))
  }
  message("Downloading ", url)
  ok <- tryCatch({
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) {
    message("Download failed: ", conditionMessage(e))
    FALSE
  })
  data.table(url = url, file = dest, status = if (ok) "downloaded" else "failed",
             bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_)
}

force_download <- toupper(Sys.getenv("FORCE_DOWNLOAD", "FALSE")) == "TRUE"

series_accessions <- c("GSE109905", "GSE113967", "GSE108785", "GSE27044")
series_results <- rbindlist(lapply(series_accessions, function(acc) {
  bucket <- geo_series_bucket(acc)
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
                 bucket, acc, acc)
  dest <- file.path(raw_dir, acc, paste0(acc, "_series_matrix.txt.gz"))
  download_if_needed(url, dest, force_download)
}), fill = TRUE)

# GSE83424 is most reliably represented by the GEO MINiML family archive, which
# contains the family XML plus sample table files used by the broad gene-level
# workflow.
miniml_accessions <- c("GSE83424", "GSE140730")
miniml_results <- rbindlist(lapply(miniml_accessions, function(acc) {
  bucket <- geo_series_bucket(acc)
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/miniml/%s_family.xml.tgz",
                 bucket, acc, acc)
  dest <- file.path(raw_dir, acc, paste0(acc, "_family.xml.tgz"))
  res <- download_if_needed(url, dest, force_download)
  if (file.exists(dest)) {
    out_dir <- file.path(raw_dir, acc)
    utils::untar(dest, exdir = out_dir)
    if (acc == "GSE83424") {
      sample_dir <- file.path(out_dir, "sample_tables")
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      tbls <- list.files(out_dir, pattern = "^GSM.*-tbl-1\\.txt$", full.names = TRUE)
      if (length(tbls)) file.copy(tbls, sample_dir, overwrite = TRUE)
    }
  }
  res
}), fill = TRUE)

annotation_src <- file.path(repo_root, "illumina450k_annotation_core.csv")
annotation_dest <- file.path(raw_dir, "annotation", "illumina450k_annotation_core.csv")
dir.create(dirname(annotation_dest), recursive = TRUE, showWarnings = FALSE)
annotation_status <- if (file.exists(annotation_src)) {
  file.copy(annotation_src, annotation_dest, overwrite = TRUE)
  "copied_from_repository_snapshot"
} else {
  "missing_repository_annotation"
}

manifest <- rbindlist(list(series_results, miniml_results), fill = TRUE)
manifest[, checked_at := format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")]
manifest <- rbind(
  manifest,
  data.table(url = "repository_snapshot:illumina450k_annotation_core.csv",
             file = annotation_dest,
             status = annotation_status,
             bytes = if (file.exists(annotation_dest)) file.info(annotation_dest)$size else NA_real_,
             checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  fill = TRUE
)

data.table::fwrite(manifest, file.path(qc_dir, "blood_source_download_manifest.csv"))

if (any(manifest$status == "failed") || any(manifest$status == "missing_repository_annotation")) {
  stop("One or more required blood methylation source files could not be retrieved. See qc/blood_source_download_manifest.csv",
       call. = FALSE)
}

message("Blood methylation source retrieval complete.")
