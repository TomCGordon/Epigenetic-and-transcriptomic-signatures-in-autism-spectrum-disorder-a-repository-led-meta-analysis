#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})
options(timeout = max(3600, getOption("timeout")))

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "lib", "placenta_lcl_methylation_functions.R"))

raw_dir <- file.path(package_root, "data_raw")
placenta_dir <- file.path(raw_dir, "placenta")
lcl_dir <- file.path(raw_dir, "lcl")
annotation_dir <- file.path(raw_dir, "annotation")
qc_dir <- file.path(package_root, "qc")
dir_create(placenta_dir)
dir_create(lcl_dir)
dir_create(annotation_dir)
dir_create(qc_dir)

geo_series_url <- function(gse) {
  number <- as.integer(sub("^GSE", "", gse))
  prefix <- floor(number / 1000)
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE", prefix, "nnn/", gse, "/matrix/", gse, "_series_matrix.txt.gz")
}

optional_source_archive_file <- function(env_var, file_name) {
  root <- Sys.getenv(env_var, unset = "")
  if (!nzchar(root)) return(NA_character_)
  file.path(root, file_name)
}

get_gsm_table <- function(gse) {
  g <- getGEO(gse, GSEMatrix = FALSE, AnnotGPL = FALSE, getGPL = FALSE)
  rows <- lapply(GSMList(g), function(gsm) {
    meta <- Meta(gsm)
    chars <- paste(meta[["characteristics_ch1"]] %||% character(), collapse = " | ")
    data.table(
      accession = gse,
      sample_id = meta[["geo_accession"]][1] %||% NA_character_,
      title = meta[["title"]][1] %||% NA_character_,
      source_name = meta[["source_name_ch1"]][1] %||% NA_character_,
      characteristics = chars,
      supplementary_file = meta[["supplementary_file_1"]][1] %||% NA_character_
    )
  })
  rbindlist(rows, fill = TRUE)
}

gse178203 <- get_gsm_table("GSE178203")
gse178203[, group := fifelse(grepl("genotype:\\s*ASD", characteristics), "ASD",
                             fifelse(grepl("genotype:\\s*TD", characteristics), "Control", NA_character_))]
gse178203[, include := group %in% c("ASD", "Control") & grepl("Placenta", characteristics, ignore.case = TRUE)]
gse178203[, build := "hg38"]
gse178203[, source_file_url := normalise_public_url(supplementary_file)]
gse178203[, staged_file := file.path(placenta_dir, "GSE178203", paste0(sample_id, ".CpG_report.txt.gz"))]

gse67615 <- get_gsm_table("GSE67615")
gse67615[, group := fifelse(grepl("diagnosis:\\s*autism", characteristics, ignore.case = TRUE), "ASD",
                            fifelse(grepl("diagnosis:\\s*typical", characteristics, ignore.case = TRUE), "Control", NA_character_))]
gse67615[, include := group %in% c("ASD", "Control")]
gse67615[, build := "hg19"]
gse67615[, source_file_url := normalise_public_url(supplementary_file)]
gse67615[, staged_file := file.path(placenta_dir, "GSE67615", basename(source_file_url))]

placenta_manifest <- rbind(gse178203, gse67615, fill = TRUE)
write_csv(placenta_manifest, file.path(qc_dir, "placenta_GEO_sample_manifest.csv"))

lcl_matrix_sources <- data.table(
  accession = c("GSE34099", "GSE99935"),
  dest_name = c("GSE34099_series_matrix.txt.gz",
                "GSE99935_Matrix_Transposed_normalized_MeDIP_data_bkgd_subtracted.txt.gz"),
  url = c(geo_series_url("GSE34099"),
          "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99935/suppl/GSE99935_Matrix_Transposed_normalized_MeDIP_data_bkgd_subtracted.txt.gz"),
  local_fallback = c(
    optional_source_archive_file("LCL_SOURCE_ARCHIVE", "GSE34099_series_matrix.txt.gz"),
    optional_source_archive_file("LCL_SOURCE_ARCHIVE", "GSE99935_Matrix_Transposed_normalized_MeDIP_data_bkgd_subtracted.txt.gz")
  )
)
lcl_manifest <- rbindlist(lapply(seq_len(nrow(lcl_matrix_sources)), function(i) {
  row <- lcl_matrix_sources[i]
  ans <- download_or_copy(row$url, row$local_fallback, file.path(lcl_dir, row$dest_name))
  ans[, accession := row$accession]
  ans[, omic_source := "LCL_public_matrix"]
  ans
}), fill = TRUE)
write_csv(lcl_manifest, file.path(qc_dir, "lcl_source_matrix_manifest.csv"))

gse34099_meta <- get_gsm_table("GSE34099")
gse34099_meta[, group := fifelse(grepl("disease state:\\s*Autism", characteristics), "ASD",
                                 fifelse(grepl("disease state:\\s*Control", characteristics), "Control", NA_character_))]
gse34099_meta[, include := group %in% c("ASD", "Control")]
write_csv(gse34099_meta, file.path(qc_dir, "GSE34099_GEO_sample_manifest.csv"))

gse99935_meta <- get_gsm_table("GSE99935")
gse99935_meta[, group := fifelse(grepl("disease status:\\s*Autistic", characteristics), "ASD",
                                 fifelse(grepl("disease status:\\s*Control", characteristics), "Control", NA_character_))]
gse99935_meta[, assay_component := fifelse(grepl("_MeDIP$", title), "MeDIP",
                                           fifelse(grepl("_Input$", title), "Input", NA_character_))]
gse99935_meta[, subject_code := sub("^[AC]", "", sub("_(MeDIP|Input)$", "", title))]
gse99935_meta[, include := group %in% c("ASD", "Control") & assay_component == "MeDIP"]
write_csv(gse99935_meta, file.path(qc_dir, "GSE99935_GEO_sample_manifest.csv"))

refgene_sources <- data.table(
  build = c("hg18", "hg19", "hg38"),
  url = c(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
  ),
  local_fallback = c(
    optional_source_archive_file("REFGENE_SOURCE_ARCHIVE", "hg18_refGene.txt.gz"),
    optional_source_archive_file("REFGENE_SOURCE_ARCHIVE", "hg19_refGene.txt.gz"),
    optional_source_archive_file("REFGENE_SOURCE_ARCHIVE", "hg38_refGene.txt.gz")
  )
)
refgene_manifest <- rbindlist(lapply(seq_len(nrow(refgene_sources)), function(i) {
  row <- refgene_sources[i]
  dest <- file.path(annotation_dir, paste0(row$build, "_refGene.txt.gz"))
  ans <- download_or_copy(row$url, row$local_fallback, dest)
  ans[, `:=`(accession = row$build, omic_source = "refGene_annotation")]
  ans
}), fill = TRUE)

write_csv(refgene_manifest, file.path(qc_dir, "refGene_source_manifest.csv"))
message("Placenta/LCL source manifests and public matrix files staged.")
