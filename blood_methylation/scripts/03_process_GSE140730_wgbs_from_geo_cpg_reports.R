#!/usr/bin/env Rscript

# Build per-sample promoter-gene methylation summaries for GSE140730 directly
# from public GEO Bismark CpG reports. This is intentionally resource-heavy:
# the CpG reports are large, so the script processes samples one at a time and
# writes reusable per-sample outputs.

suppressPackageStartupMessages({
  library(data.table)
  library(xml2)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
raw_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_RAW_DIR", unset = file.path(package_root, "data_raw")),
                         winslash = "/", mustWork = FALSE)
processed_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_PROCESSED_DIR", unset = file.path(package_root, "data_processed")),
                               winslash = "/", mustWork = FALSE)
qc_dir <- normalizePath(file.path(package_root, "qc"), winslash = "/", mustWork = FALSE)
wgbs_dir <- file.path(processed_dir, "GSE140730_wgbs")
per_sample_dir <- file.path(wgbs_dir, "per_sample_outputs")
dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

family_xml <- file.path(raw_dir, "GSE140730", "GSE140730_family.xml")
coord_path <- file.path(processed_dir, "annotation", "GSE140730_hg38_promoter_coordinates_from_UCSC_refGene.csv")
if (!file.exists(family_xml)) stop("Missing GSE140730 family XML: ", family_xml, "\nRun 01_download_blood_sources.R first.", call. = FALSE)
if (!file.exists(coord_path)) stop("Missing WGBS promoter coordinates: ", coord_path, "\nRun 02_build_wgbs_promoter_coordinates.R first.", call. = FALSE)

clean_text <- function(x) trimws(gsub("[\r\n\t]+", " ", x))

parse_gse140730_xml <- function(path) {
  doc <- xml2::read_xml(path)
  samples <- xml2::xml_find_all(doc, "//*[local-name()='Sample']")
  rows <- lapply(samples, function(s) {
    chars <- xml2::xml_find_all(s, ".//*[local-name()='Characteristics']")
    char_dt <- data.table(
      tag = xml2::xml_attr(chars, "tag"),
      value = clean_text(xml2::xml_text(chars))
    )
    val <- function(target_tag) {
      out <- char_dt[tag == target_tag, value]
      if (length(out)) out[[1]] else NA_character_
    }
    supp <- clean_text(xml2::xml_text(xml2::xml_find_all(s, ".//Supplementary-Data")))
    data.table(
      sample_id = xml2::xml_attr(s, "iid"),
      title = clean_text(xml2::xml_text(xml2::xml_find_first(s, "./*[local-name()='Title']"))),
      sample_set = val("sample set"),
      tissue = val("tissue"),
      diagnosis = val("diagnosis"),
      sex = val("Sex"),
      study = val("study"),
      report_url = clean_text(xml2::xml_text(xml2::xml_find_all(s, ".//*[local-name()='Supplementary-Data']")))[1]
    )
  })
  rbindlist(rows, fill = TRUE)
}

samples <- parse_gse140730_xml(family_xml)
samples[, phenotype_group := fifelse(diagnosis == "ASD", "ASD",
                                     fifelse(diagnosis %in% c("TD", "Control"), "TD", "Exclude"))]
samples <- samples[phenotype_group %in% c("ASD", "TD") & grepl("cord", tissue, ignore.case = TRUE)]
sample_limit <- as.integer(Sys.getenv("BLOOD_WGBS_SAMPLE_LIMIT", "0"))
if (is.finite(sample_limit) && sample_limit > 0) samples <- samples[seq_len(min(.N, sample_limit))]
fwrite(samples, file.path(qc_dir, "GSE140730_wgbs_sample_manifest.csv"))

promoters <- fread(coord_path, showProgress = FALSE)
promoters <- promoters[, .(
  gene = gene_symbol,
  chromosome,
  promoter_start = as.integer(promoter_start),
  promoter_end = as.integer(promoter_end)
)]
promoters <- unique(promoters[!is.na(gene) & !is.na(chromosome) & is.finite(promoter_start) & is.finite(promoter_end)])
setkey(promoters, chromosome, promoter_start, promoter_end)

local_report_path <- function(sample_id, url) {
  fname <- basename(url)
  candidates <- c(
    file.path(raw_dir, "GSE140730", "cpg_reports", fname),
    file.path(raw_dir, "GSE140730", fname),
    file.path(getwd(), fname)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) && !is.na(hit)) hit else NA_character_
}

download_report_if_requested <- function(url) {
  out <- file.path(raw_dir, "GSE140730", "cpg_reports", basename(url))
  if (file.exists(out)) return(out)
  if (toupper(Sys.getenv("BLOOD_WGBS_DOWNLOAD_REPORTS", "FALSE")) != "TRUE") return(NA_character_)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  message("Downloading WGBS CpG report: ", basename(url))
  utils::download.file(url, out, mode = "wb", quiet = FALSE)
  out
}

read_cpg_report <- function(sample_id, url) {
  local <- local_report_path(sample_id, url)
  if (is.na(local)) local <- download_report_if_requested(url)
  source_ref <- if (!is.na(local) && file.exists(local)) local else url
  message("Reading CpG report for ", sample_id, " from ", ifelse(grepl("^https?://", source_ref), "remote GEO URL", "local file"))
  x <- fread(source_ref, header = FALSE, select = 1:5, showProgress = FALSE,
             col.names = c("chromosome", "position", "strand", "methylated_count", "unmethylated_count"))
  x[, `:=`(
    position = as.integer(position),
    methylated_count = as.numeric(methylated_count),
    unmethylated_count = as.numeric(unmethylated_count)
  )]
  x <- x[chromosome %chin% unique(promoters$chromosome)]
  x[, depth := methylated_count + unmethylated_count]
  x <- x[is.finite(position) & is.finite(depth) & depth > 0]
  x[, beta := methylated_count / depth]
  x[, `:=`(start = position, end = position)]
  setkey(x, chromosome, start, end)
  list(data = x, source_ref = source_ref)
}

summarise_one_sample <- function(row) {
  sample_id <- row$sample_id
  out_file <- file.path(per_sample_dir, paste0(sample_id, "_promoter_long.csv"))
  if (file.exists(out_file) && toupper(Sys.getenv("FORCE_RERUN", "FALSE")) != "TRUE") {
    return(fread(out_file, showProgress = FALSE))
  }
  cpg <- read_cpg_report(sample_id, row$report_url)
  hits <- foverlaps(cpg$data, promoters, by.x = c("chromosome", "start", "end"),
                    by.y = c("chromosome", "promoter_start", "promoter_end"), nomatch = 0L)
  if (!nrow(hits)) {
    out <- data.table(sample_id = sample_id, phenotype_group = row$phenotype_group,
                      gene = character(), cpg_count = integer(),
                      methylated_count = numeric(), unmethylated_count = numeric(),
                      mean_methylation_beta = numeric(), finite_flag = logical(),
                      notes = character(), source_ref = character())
  } else {
    hits <- unique(hits[, .(gene, chromosome, position, methylated_count, unmethylated_count, beta)])
    out <- hits[, .(
      cpg_count = .N,
      methylated_count = sum(methylated_count, na.rm = TRUE),
      unmethylated_count = sum(unmethylated_count, na.rm = TRUE),
      mean_methylation_beta = mean(beta, na.rm = TRUE)
    ), by = gene]
    out[, `:=`(
      sample_id = sample_id,
      phenotype_group = row$phenotype_group,
      finite_flag = is.finite(mean_methylation_beta),
      notes = "",
      source_ref = cpg$source_ref
    )]
    setcolorder(out, c("sample_id", "phenotype_group", "gene", "cpg_count",
                       "methylated_count", "unmethylated_count", "mean_methylation_beta",
                       "finite_flag", "notes", "source_ref"))
  }
  fwrite(out, out_file)
  out
}

if (!nrow(samples)) stop("No ASD/TD cord-blood samples were parsed from GSE140730 family XML.", call. = FALSE)

pieces <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  message("Processing GSE140730 sample ", i, "/", nrow(samples), ": ", samples$sample_id[[i]])
  pieces[[i]] <- summarise_one_sample(samples[i])
  gc()
}

combined <- rbindlist(pieces, fill = TRUE)
combined_out <- file.path(wgbs_dir, "GSE140730_per_sample_promoter_long.csv")
fwrite(combined, combined_out)

qc <- data.table(
  metric = c("samples_processed", "ASD_samples", "TD_control_samples", "unique_genes_with_any_WGBS_value", "combined_output"),
  value = c(nrow(samples), sum(samples$phenotype_group == "ASD"), sum(samples$phenotype_group == "TD"),
            uniqueN(combined$gene), combined_out)
)
fwrite(qc, file.path(qc_dir, "GSE140730_wgbs_processing_QC.csv"))

message("GSE140730 WGBS promoter extraction complete: ", combined_out)
