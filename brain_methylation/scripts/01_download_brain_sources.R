#!/usr/bin/env Rscript

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
source(file.path(script_dir, "lib", "brain_methylation_functions.R"))

raw_dir <- file.path(package_root, "data_raw")
array_dir <- file.path(raw_dir, "arrays")
wgbs_dir <- file.path(raw_dir, "WGBS")
annotation_dir <- file.path(raw_dir, "annotation")
qc_dir <- file.path(package_root, "qc")
dir_create(array_dir)
dir_create(wgbs_dir)
dir_create(annotation_dir)
dir_create(qc_dir)

geo_series_url <- function(gse) {
  prefix <- sub("([0-9]{3})[0-9]+$", "\\1", sub("^GSE", "", gse))
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE", prefix, "nnn/", gse, "/matrix/", gse, "_series_matrix.txt.gz")
}

normalise_public_url <- function(url) {
  ifelse(is.na(url) | !nzchar(url), url, sub("^ftp://", "https://", url))
}

optional_cache_file <- function(env_var, file_name) {
  root <- Sys.getenv(env_var, unset = "")
  if (!nzchar(root)) return(NA_character_)
  file.path(root, file_name)
}

array_sources <- data.table(
  accession = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608"),
  source_kind = c("series_matrix", "series_matrix", "series_matrix", "processed_matrix", "series_matrix", "processed_matrix", "series_matrix"),
  local_fallback = NA_character_,
  dest_name = c(
    "GSE53162_series_matrix.txt.gz",
    "GSE53924_series_matrix.txt.gz",
    "GSE80017_series_matrix.txt.gz",
    "GSE131706_Matrix_processed.csv.gz",
    "GSE242427_series_matrix.txt.gz",
    "GSE278285_MatrixProcessed_Avg_Beta.txt.gz",
    "GSE38608_series_matrix.txt.gz"
  )
)
array_sources[, local_fallback := vapply(dest_name, optional_cache_file, character(1), env_var = "BRAIN_ARRAY_SOURCE_CACHE")]
array_sources[, url := fifelse(source_kind == "series_matrix", geo_series_url(accession), NA_character_)]
array_sources[accession == "GSE131706", url := "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131706/suppl/GSE131706_Matrix_processed.csv.gz"]
array_sources[accession == "GSE278285", url := "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278285/suppl/GSE278285_MatrixProcessed_Avg_Beta.txt.gz"]

manifest <- rbindlist(lapply(seq_len(nrow(array_sources)), function(i) {
  row <- array_sources[i]
  download_or_copy(row$url, row$local_fallback, file.path(array_dir, row$dest_name))
}), fill = TRUE)
manifest[, accession := array_sources$accession]
manifest[, omic_source := "array_or_HM27"]

# Stage companion series-matrix metadata for processed-matrix-only routes.
metadata_sources <- data.table(
  accession = c("GSE131706", "GSE278285"),
  dest_name = c("GSE131706_series_matrix.txt.gz", "GSE278285_series_matrix.txt.gz"),
  url = geo_series_url(c("GSE131706", "GSE278285"))
)
metadata_sources[, local_fallback := vapply(dest_name, optional_cache_file, character(1), env_var = "BRAIN_ARRAY_SOURCE_CACHE")]
metadata_manifest <- rbindlist(lapply(seq_len(nrow(metadata_sources)), function(i) {
  row <- metadata_sources[i]
  download_or_copy(row$url, row$local_fallback, file.path(array_dir, row$dest_name))
}), fill = TRUE)
metadata_manifest[, accession := metadata_sources$accession]
metadata_manifest[, omic_source := "series_matrix_metadata"]
manifest <- rbind(manifest, metadata_manifest, fill = TRUE)

annotation_candidates <- c(Sys.getenv("ILLUMINA450K_ANNOTATION", unset = ""),
                           file.path(annotation_dir, "illumina450k_annotation_core.csv"))
annotation_source <- annotation_candidates[file.exists(annotation_candidates)][1]
if (is.na(annotation_source)) {
  stop("Could not find illumina450k_annotation_core.csv. Provide it at data_raw/annotation/ or set ILLUMINA450K_ANNOTATION to a CSV containing Name, UCSC_RefGene_Name and UCSC_RefGene_Group.")
}
annotation_dest <- file.path(annotation_dir, "illumina450k_annotation_core.csv")
if (!file.exists(annotation_dest)) file.copy(annotation_source, annotation_dest, overwrite = TRUE)
manifest <- rbind(manifest, data.table(dest = annotation_dest, status = "copied_annotation",
                                       source_used = annotation_source, bytes = file.info(annotation_dest)$size,
                                       accession = "annotation", omic_source = "annotation"), fill = TRUE)

build_wgbs_sample_manifest <- function(gse) {
  g <- getGEO(gse, GSEMatrix = FALSE, AnnotGPL = FALSE, getGPL = FALSE)
  rows <- lapply(GSMList(g), function(gsm) {
    meta <- Meta(gsm)
    chars <- meta$characteristics_ch1 %||% character()
    title <- meta$title %||% NA_character_
    source_name <- meta$source_name_ch1 %||% NA_character_
    group <- infer_group(c(title, source_name, chars))
    include <- FALSE
    reason <- ""
    if (gse == "GSE109875") {
      include <- group %in% c("ASD", "Control")
      reason <- if (include) "BA9 ASD/control WGBS sample" else "not ASD/control"
    } else if (gse == "GSE81541") {
      text <- tolower(paste(title, source_name, chars, collapse = " "))
      include <- grepl("brain_idioaut|brain_control", text) && group %in% c("ASD", "Control") &&
        !grepl("braindfba|syndrom|cell|neuron|oligodendrocyte", text)
      reason <- if (include) "Idiopathic autism/control brain WGBS sample" else "excluded non-primary GSE81541 route"
    }
    data.table(
      accession = gse,
      sample_id = meta$geo_accession %||% NA_character_,
      title = title,
      source_name = source_name,
      group = group,
      include = include,
      exclusion_reason = reason,
      supplementary_file = meta$supplementary_file_1 %||% NA_character_
    )
  })
  rbindlist(rows, fill = TRUE)
}

wgbs_manifest <- rbind(build_wgbs_sample_manifest("GSE109875"), build_wgbs_sample_manifest("GSE81541"), fill = TRUE)
fwrite(wgbs_manifest, file.path(qc_dir, "brain_WGBS_GEO_sample_manifest.csv"))

wgbs_source_cache <- Sys.getenv("BRAIN_WGBS_SOURCE_CACHE", unset = "")
wgbs_stage_manifest <- rbindlist(lapply(seq_len(nrow(wgbs_manifest[include == TRUE])), function(i) {
  row <- wgbs_manifest[include == TRUE][i]
  short_name <- if (row$accession == "GSE109875") {
    paste0(row$sample_id, ".cpg_report.txt.gz")
  } else {
    paste0(row$sample_id, ".bed.gz")
  }
  dest <- file.path(wgbs_dir, row$accession, short_name)
  fallback <- if (nzchar(wgbs_source_cache)) file.path(wgbs_source_cache, row$accession, basename(row$supplementary_file)) else NA_character_
  ans <- download_or_copy(normalise_public_url(row$supplementary_file), fallback, dest)
  ans[, `:=`(sample_id = row$sample_id, group = row$group, original_file_name = basename(row$supplementary_file))]
  ans
}), fill = TRUE)
if (nrow(wgbs_stage_manifest)) {
  wgbs_stage_manifest[, accession := wgbs_manifest[include == TRUE]$accession]
  wgbs_stage_manifest[, omic_source := "WGBS_public_CpG_or_BED_file"]
}
manifest <- rbind(manifest, wgbs_stage_manifest, fill = TRUE)
fwrite(wgbs_stage_manifest, file.path(qc_dir, "brain_WGBS_staged_source_manifest.csv"))

fwrite(manifest, file.path(qc_dir, "brain_source_download_manifest.csv"))
message("Brain methylation source files staged.")
