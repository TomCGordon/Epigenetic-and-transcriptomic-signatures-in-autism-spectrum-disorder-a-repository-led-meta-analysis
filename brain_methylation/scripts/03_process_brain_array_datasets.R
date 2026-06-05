#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "lib", "brain_methylation_functions.R"))

array_dir <- file.path(package_root, "data_raw", "arrays")
annotation_dir <- file.path(package_root, "data_processed", "annotation")
processed_dir <- file.path(package_root, "data_processed", "arrays")
qc_dir <- file.path(package_root, "qc")
dir_create(processed_dir)
dir_create(qc_dir)

probe_map <- fread(file.path(annotation_dir, "brain_450k_promoter_probe_gene_map.csv"))
universe <- fread(file.path(annotation_dir, "brain_harmonised_20960_gene_universe.csv"))$gene

dataset_info <- data.table(
  accession = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608"),
  matrix_file = file.path(array_dir, c(
    "GSE53162_series_matrix.txt.gz",
    "GSE53924_series_matrix.txt.gz",
    "GSE80017_series_matrix.txt.gz",
    "GSE131706_Matrix_processed.csv.gz",
    "GSE242427_series_matrix.txt.gz",
    "GSE278285_MatrixProcessed_Avg_Beta.txt.gz",
    "GSE38608_series_matrix.txt.gz"
  )),
  metadata_file = file.path(array_dir, c(
    "GSE53162_series_matrix.txt.gz",
    "GSE53924_series_matrix.txt.gz",
    "GSE80017_series_matrix.txt.gz",
    "GSE131706_series_matrix.txt.gz",
    "GSE242427_series_matrix.txt.gz",
    "GSE278285_series_matrix.txt.gz",
    "GSE38608_series_matrix.txt.gz"
  )),
  platform = c(rep("Illumina HumanMethylation450 BeadChip", 6), "Illumina HumanMethylation27 BeadChip"),
  assay_class = c(rep("450K", 6), "HM27"),
  brain_region = c(
    "Post-mortem prefrontal/temporal/cerebellum grouped route",
    "Post-mortem Brodmann area 10/24 cortex",
    "Post-mortem prefrontal cortex",
    "Post-mortem subventricular zone of lateral ventricles",
    "Post-mortem dorsal raphe",
    "Post-mortem cerebellum",
    "Post-mortem cerebellum and occipital cortex"
  ),
  broader_region_group = c("mixed", "cortex", "cortex", "subventricular zone", "dorsal raphe", "cerebellum", "mixed"),
  collapse_subject = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)
)

sample_id_from_title <- function(title, fallback) {
  id <- sub("_.*$", "", title)
  ifelse(is.na(id) | !nzchar(id), fallback, id)
}

subject_from_title <- function(title, sample_id) {
  x <- gsub("_rep[0-9]+$", "", title, ignore.case = TRUE)
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  out <- parts[length(parts)]
  ifelse(is.na(out) | !nzchar(out), sample_id, out)
}

parse_subregion <- function(accession, title, source_name, characteristics) {
  text <- tolower(paste(title, source_name, characteristics, collapse = " "))
  if (grepl("prefrontal", text)) return("prefrontal cortex")
  if (grepl("temporal", text)) return("temporal cortex")
  if (grepl("frontal cortex \\(ba10\\)|\\bba10\\b", text)) return("frontal cortex BA10")
  if (grepl("cingulate cortex \\(ba24\\)|\\bba24\\b", text)) return("cingulate cortex BA24")
  if (grepl("occipital", text)) return("occipital cortex")
  if (grepl("cerebell", text)) return("cerebellum")
  if (grepl("subventricular", text)) return("subventricular zone")
  if (grepl("dorsal raphe|raphe", text)) return("dorsal raphe")
  NA_character_
}

slugify <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

build_sample_meta <- function(accession, metadata_file, matrix_sample_cols) {
  meta <- read_geo_sample_metadata(metadata_file)
  for (nm in c("geo_accession", "title", "source_name_ch1", "characteristics_ch1", "description")) {
    if (!nm %in% names(meta)) meta[, (nm) := NA_character_]
  }
  meta[, sample_id := geo_accession %||% title]
  meta[, group := vapply(seq_len(.N), function(i) infer_group(c(title[i], source_name_ch1[i], characteristics_ch1[i], description[i])), character(1))]
  meta[, region := vapply(characteristics_ch1, extract_after_colon, character(1), pattern = "brain region|brodmann|tissue subtype")]
  if (accession %in% c("GSE131706", "GSE278285")) {
    # These processed matrices use submitter sample labels, not GSM accessions.
    meta[, sample_id := sample_id_from_title(title, geo_accession)]
    if (length(matrix_sample_cols) == nrow(meta)) {
      meta[, sample_id := matrix_sample_cols]
    }
  }
  if (accession == "GSE53924") {
    meta[, sample_id := geo_accession]
    meta[, subject_id := source_name_ch1]
  } else if (accession == "GSE38608") {
    meta[, sample_id := geo_accession]
    meta[, subject_id := vapply(seq_len(.N), function(i) subject_from_title(title[i], geo_accession[i]), character(1))]
  } else {
    meta[, subject_id := sample_id]
  }
  meta[, analysis_region := vapply(seq_len(.N), function(i) {
    parse_subregion(accession, title[i], source_name_ch1[i], characteristics_ch1[i])
  }, character(1))]
  meta[, .(sample_id, geo_accession, title, source_name_ch1, characteristics_ch1, group, region, subject_id, analysis_region)]
}

region_qc_outputs <- list()
region_summary_outputs <- list()
region_effect_outputs <- list()

process_dataset <- function(info) {
  accession <- info$accession
  message("Processing ", accession)
  if (accession %in% c("GSE131706", "GSE278285")) {
    beta_dt <- fread(info$matrix_file, showProgress = TRUE)
    setnames(beta_dt, names(beta_dt)[1], "ID_REF")
  } else {
    beta_dt <- read_geo_matrix_table(info$matrix_file)
  }
  sample_cols <- setdiff(names(beta_dt), "ID_REF")
  beta_dt <- normalise_numeric_matrix(beta_dt, sample_cols, blanks_as_zero = FALSE)
  sample_meta <- build_sample_meta(accession, info$metadata_file, sample_cols)
  sample_meta <- sample_meta[sample_id %in% sample_cols & group %in% c("ASD", "Control")]
  if (!nrow(sample_meta)) stop("No ASD/control samples retained for ", accession)
  beta_dt <- beta_dt[, c("ID_REF", sample_meta$sample_id), with = FALSE]
  gm <- build_gene_matrix(beta_dt, probe_map, sample_meta$sample_id)
  summary_dt <- summarise_gene_matrix(
    dataset_id = accession,
    gene_dt = gm$gene_matrix,
    sample_meta = sample_meta,
    feature_counts = gm$feature_counts,
    brain_region = info$brain_region,
    broader_region_group = info$broader_region_group,
    platform = info$platform,
    assay_class = info$assay_class,
    collapse_subject = info$collapse_subject
  )
  summary_dt <- merge(data.table(gene = universe), summary_dt, by = "gene", all.x = TRUE)
  summary_dt[is.na(accession), `:=`(
    accession = info$accession,
    brain_region = info$brain_region,
    broader_region_group = info$broader_region_group,
    platform = info$platform,
    assay_class = info$assay_class,
    finite_summary = "no",
    missingness_reason = "gene not represented by finite promoter values"
  )]
  effects <- hedges_effects(summary_dt)
  fwrite(sample_meta, file.path(processed_dir, paste0(accession, "_sample_classification.csv")))
  fwrite(summary_dt, file.path(processed_dir, paste0(accession, "_gene_summary_statistics.csv")))
  fwrite(effects, file.path(processed_dir, paste0(accession, "_effect_sizes.csv")))

  regions <- sort(unique(na.omit(sample_meta$analysis_region)))
  for (rg in regions) {
    region_meta <- sample_meta[analysis_region == rg]
    if (sum(region_meta$group == "ASD") < 2 || sum(region_meta$group == "Control") < 2) next
    sub_id <- paste0(accession, "_", slugify(rg))
    broader <- infer_region_group(rg)
    sub_summary <- summarise_gene_matrix(
      dataset_id = sub_id,
      gene_dt = gm$gene_matrix,
      sample_meta = region_meta,
      feature_counts = gm$feature_counts,
      brain_region = paste0("Region-specific subroute: ", rg),
      broader_region_group = broader,
      platform = info$platform,
      assay_class = info$assay_class,
      collapse_subject = info$collapse_subject
    )
    sub_summary <- merge(data.table(gene = universe), sub_summary, by = "gene", all.x = TRUE)
    sub_summary[is.na(accession), `:=`(
      accession = sub_id,
      brain_region = paste0("Region-specific subroute: ", rg),
      broader_region_group = broader,
      platform = info$platform,
      assay_class = info$assay_class,
      finite_summary = "no",
      missingness_reason = "gene not represented by finite promoter values"
    )]
    sub_effects <- hedges_effects(sub_summary)
    fwrite(sub_summary, file.path(processed_dir, paste0(sub_id, "_gene_summary_statistics.csv")))
    fwrite(sub_effects, file.path(processed_dir, paste0(sub_id, "_effect_sizes.csv")))
    region_summary_outputs[[sub_id]] <<- sub_summary
    region_effect_outputs[[sub_id]] <<- sub_effects
    region_qc_outputs[[sub_id]] <<- data.table(
      accession = sub_id,
      parent_accession = accession,
      region = rg,
      ASD_n_or_subjects = if (info$collapse_subject) uniqueN(region_meta[group == "ASD"]$subject_id) else sum(region_meta$group == "ASD"),
      control_n_or_subjects = if (info$collapse_subject) uniqueN(region_meta[group == "Control"]$subject_id) else sum(region_meta$group == "Control"),
      effect_size_rows = nrow(sub_effects),
      use = "region/subtissue sensitivity only"
    )
  }

  data.table(
    accession = accession,
    retained_samples = nrow(sample_meta),
    ASD_n_or_subjects = if (info$collapse_subject) uniqueN(sample_meta[group == "ASD"]$subject_id) else sum(sample_meta$group == "ASD"),
    control_n_or_subjects = if (info$collapse_subject) uniqueN(sample_meta[group == "Control"]$subject_id) else sum(sample_meta$group == "Control"),
    finite_gene_summaries = sum(summary_dt$finite_summary == "yes", na.rm = TRUE),
    effect_size_rows = nrow(effects),
    collapse_subject = info$collapse_subject,
    region_subroute_count = length(regions)
  )
}

qc <- rbindlist(lapply(seq_len(nrow(dataset_info)), function(i) process_dataset(dataset_info[i])), fill = TRUE)
fwrite(qc, file.path(qc_dir, "brain_array_dataset_processing_QC.csv"))

region_qc <- rbindlist(region_qc_outputs, fill = TRUE)
region_summaries <- rbindlist(region_summary_outputs, fill = TRUE)
region_effects <- rbindlist(region_effect_outputs, fill = TRUE)
if (nrow(region_qc)) fwrite(region_qc, file.path(qc_dir, "brain_array_region_subroute_QC.csv"))
if (nrow(region_summaries)) fwrite(region_summaries, file.path(processed_dir, "brain_array_region_subroute_gene_summary_statistics.csv"))
if (nrow(region_effects)) fwrite(region_effects, file.path(processed_dir, "brain_array_region_subroute_effect_sizes.csv"))

all_summary <- rbindlist(lapply(dataset_info$accession, function(id) {
  fread(file.path(processed_dir, paste0(id, "_gene_summary_statistics.csv")))
}), fill = TRUE)
all_effects <- rbindlist(lapply(dataset_info$accession, function(id) {
  fread(file.path(processed_dir, paste0(id, "_effect_sizes.csv")))
}), fill = TRUE)
fwrite(all_summary, file.path(processed_dir, "brain_array_all_dataset_gene_summary_statistics.csv"))
fwrite(all_effects, file.path(processed_dir, "brain_array_all_dataset_effect_sizes.csv"))
message("Brain array/HM27 datasets processed.")
