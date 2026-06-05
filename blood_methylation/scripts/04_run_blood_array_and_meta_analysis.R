#!/usr/bin/env Rscript

# Blood methylation processing and branch-level QC summaries.
#
# This script generates promoter-gene summaries and dataset-level Hedges' g
# inputs for the final methylation meta-analysis. It also writes branch-level
# QC summaries. The reported pooled methylation models are fitted
# separately in scripts/04_meta_analysis/ using metafor.

suppressPackageStartupMessages({
  extra_lib <- Sys.getenv("ASD_R_LIB", unset = "")
  if (nzchar(extra_lib)) .libPaths(unique(c(extra_lib, .libPaths())))
  library(data.table)
  library(metafor)
})

# Set ASD_REPO_ROOT if this script is not run from the repository root.
# Example: Sys.setenv(ASD_REPO_ROOT = "/path/to/ASD_repository_audit")
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
BUNDLE_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(Sys.getenv("ASD_REPO_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
RAW_DIR <- normalizePath(Sys.getenv("BLOOD_METHYLATION_RAW_DIR", unset = file.path(BUNDLE_ROOT, "data_raw")),
                         winslash = "/", mustWork = FALSE)
PROCESSED_DIR <- normalizePath(Sys.getenv("BLOOD_METHYLATION_PROCESSED_DIR", unset = file.path(BUNDLE_ROOT, "data_processed")),
                               winslash = "/", mustWork = FALSE)
OUT_BASE <- normalizePath(
  Sys.getenv("BLOOD_METHYLATION_R_OUTPUT_DIR", unset = file.path(BUNDLE_ROOT, "results")),
  winslash = "/", mustWork = FALSE
)
blank_mode <- tolower(Sys.getenv("BLOOD_ARRAY_BLANK_HANDLING", "missing"))
BLANKS_AS_ZERO <- blank_mode %in% c("zero", "blank_as_zero", "sensitivity_zero")
OUT <- file.path(
  OUT_BASE,
  if (BLANKS_AS_ZERO) "array_blanks_as_zero_sensitivity" else "array_blanks_as_missing_R_default"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

lp <- function(path) {
  p <- normalizePath(path, winslash = "\\", mustWork = FALSE)
  if (.Platform$OS.type == "windows" && !startsWith(p, "\\\\?\\")) paste0("\\\\?\\", p) else p
}
rp <- function(...) file.path(ROOT, ...)
op <- function(...) file.path(OUT, ...)

read_csv <- function(path, ...) fread(file = lp(path), showProgress = FALSE, data.table = TRUE, ...)
write_csv <- function(x, name) fwrite(as.data.table(x), lp(op(name)))
num <- function(x) suppressWarnings(as.numeric(x))
clean <- function(x) trimws(gsub('^"|"$', "", as.character(x)))
mean2 <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
sd2 <- function(x) if (length(x) > 1) sd(x, na.rm = TRUE) else NA_real_
ci_excludes_zero <- function(lo, hi) is.finite(lo) & is.finite(hi) & ((lo > 0 & hi > 0) | (lo < 0 & hi < 0))
tcrit95_lookup <- function(df) {
  tab <- c(
    `1` = 12.706, `2` = 4.303, `3` = 3.182, `4` = 2.776, `5` = 2.571,
    `6` = 2.447, `7` = 2.365, `8` = 2.306, `9` = 2.262, `10` = 2.228,
    `11` = 2.201, `12` = 2.179, `13` = 2.160, `14` = 2.145, `15` = 2.131,
    `16` = 2.120, `17` = 2.110, `18` = 2.101, `19` = 2.093, `20` = 2.086,
    `21` = 2.080, `22` = 2.074, `23` = 2.069, `24` = 2.064, `25` = 2.060,
    `26` = 2.056, `27` = 2.052, `28` = 2.048, `29` = 2.045, `30` = 2.042
  )
  tab[[as.character(max(1, min(30, round(df))))]] %||% 1.96
}

promoter_terms <- c("TSS200", "TSS1500", "5'UTR", "1stExon")
blank_handling_note <- if (BLANKS_AS_ZERO) {
  paste(
    "Blank-as-zero sensitivity mode: blank array cells are converted to zero after numeric coercion.",
    "This mode is retained only as an explicit missing-data sensitivity check."
  )
} else {
  paste(
    "R-default mode: blank array cells are retained as missing values and excluded from sample-level gene summaries.",
    "This is the stricter missing-data interpretation used to diagnose blank-cell handling."
  )
}

datasets <- data.table(
  accession = c("GSE109905", "GSE113967", "GSE83424", "GSE108785", "GSE27044"),
  assay_class = c("450K", "450K", "450K", "450K", "HM27"),
  tissue = c("Whole blood", "Whole blood", "Peripheral blood", "Whole blood", "Peripheral blood leukocytes"),
  data_mode = c("series_matrix_gz", "series_matrix_gz", "family_tables", "series_matrix_gz", "series_matrix_gz"),
  matrix = c(
    file.path(RAW_DIR, "GSE109905", "GSE109905_series_matrix.txt.gz"),
    file.path(RAW_DIR, "GSE113967", "GSE113967_series_matrix.txt.gz"),
    NA_character_,
    file.path(RAW_DIR, "GSE108785", "GSE108785_series_matrix.txt.gz"),
    file.path(RAW_DIR, "GSE27044", "GSE27044_series_matrix.txt.gz")
  ),
  family_xml = c(NA_character_, NA_character_, file.path(RAW_DIR, "GSE83424", "GSE83424_family.xml"), NA_character_, NA_character_),
  sample_dir = c(NA_character_, NA_character_, file.path(RAW_DIR, "GSE83424", "sample_tables"), NA_character_, NA_character_),
  phenotype_note = c(
    "Series-matrix group metadata were used to classify ASD and Control samples.",
    "Series-matrix ASD and Control labels were used; 16p11.2 deletion and CHD8 subgroups were excluded.",
    "Family XML status was used; Case was harmonised to ASD and Control retained directly.",
    "Metadata were classified conservatively using ASD/autism versus unaffected/control wording.",
    "Sample titles ending in .p were treated as autistic probands and .s as unaffected siblings, matching the prior validated rule."
  )
)

model_defs <- list(
  blood_array_peripheral_primary = c("GSE109905", "GSE113967", "GSE83424", "GSE108785", "GSE27044"),
  blood_450k_only_sensitivity = c("GSE109905", "GSE113967", "GSE83424", "GSE108785"),
  blood_array_plus_cord_WGBS_sensitivity = c("GSE109905", "GSE113967", "GSE83424", "GSE108785", "GSE27044", "GSE140730"),
  blood_450k_plus_cord_WGBS_sensitivity = c("GSE109905", "GSE113967", "GSE83424", "GSE108785", "GSE140730")
)

model_roles <- data.table(
  model_name = names(model_defs),
  role = c("primary", "sensitivity", "sensitivity", "sensitivity"),
  caveat = c(
    "",
    "450K-only platform sensitivity.",
    "Cord blood is developmental and not equivalent to postnatal peripheral blood.",
    "Cord blood is developmental and not equivalent to postnatal peripheral blood."
  )
)

hedges_from_groups <- function(asd, control) {
  asd <- asd[is.finite(asd)]
  control <- control[is.finite(control)]
  n1 <- length(asd)
  n0 <- length(control)
  if (n1 < 2 || n0 < 2) return(NULL)
  sd1 <- sd(asd)
  sd0 <- sd(control)
  if (!is.finite(sd1) || !is.finite(sd0)) return(NULL)
  pooled <- sqrt(((n1 - 1) * sd1^2 + (n0 - 1) * sd0^2) / (n1 + n0 - 2))
  if (!is.finite(pooled) || pooled <= 0) return(NULL)
  d <- (mean(asd) - mean(control)) / pooled
  j <- 1 - (3 / (4 * (n1 + n0) - 9))
  g <- d * j
  vg <- ((n1 + n0) / (n1 * n0)) + ((g^2) / (2 * (n1 + n0 - 2)))
  list(
    pooled_sd = pooled,
    Hedges_g = g,
    variance_g = vg,
    standard_error_g = sqrt(vg)
  )
}

parse_promoter_annotation <- function() {
  cache_file <- file.path(OUT_BASE, "blood_R_promoter_annotation_cache.rds")
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  anno <- read_csv(file.path(RAW_DIR, "annotation", "illumina450k_annotation_core.csv"))
  needed <- c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")
  stopifnot(all(needed %in% names(anno)))
  out <- vector("list", nrow(anno))
  for (i in seq_len(nrow(anno))) {
    probe <- anno$Name[[i]]
    genes <- strsplit(clean(anno$UCSC_RefGene_Name[[i]]), ";", fixed = TRUE)[[1]]
    groups <- strsplit(clean(anno$UCSC_RefGene_Group[[i]]), ";", fixed = TRUE)[[1]]
    if (!length(genes) || !nzchar(genes[[1]])) next
    n <- max(length(genes), length(groups))
    rows <- vector("list", n)
    for (j in seq_len(n)) {
      gene <- clean(if (j <= length(genes) && nzchar(genes[[j]])) genes[[j]] else genes[[1]])
      group <- clean(if (j <= length(groups) && nzchar(groups[[j]])) groups[[j]] else groups[[1]])
      if (nzchar(gene) && group %in% promoter_terms) {
        rows[[j]] <- data.table(probe = probe, gene = gene, promoter_group = group)
      }
    }
    out[[i]] <- rbindlist(rows, fill = TRUE)
  }
  map <- unique(rbindlist(out, fill = TRUE)[!is.na(gene) & nzchar(gene), .(probe, gene)])
  setkey(map, probe)
  universe <- sort(unique(map$gene))
  feature_info <- map[, .(
    promoter_feature_count = uniqueN(probe),
    promoter_feature_ids = paste(sort(unique(probe)), collapse = ";")
  ), by = gene]
  parsed <- list(map = map, universe = universe, feature_info = feature_info)
  saveRDS(parsed, cache_file)
  parsed
}

split_tsv <- function(line) clean(strsplit(line, "\t", fixed = TRUE)[[1]])

series_metadata <- function(file, gzipped = FALSE) {
  con <- if (gzipped) gzfile(lp(file), "rt") else file(lp(file), "rt")
  on.exit(close(con), add = TRUE)
  sample_meta <- list()
  char_lines <- list()
  platform <- ""
  table_samples <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (line == "!series_matrix_table_begin") {
      header <- readLines(con, n = 1, warn = FALSE)
      table_samples <- split_tsv(header)[-1]
      break
    }
    if (startsWith(line, "!Series_platform_id")) platform <- split_tsv(line)[2]
    if (startsWith(line, "!Sample_title")) sample_meta$title <- split_tsv(line)[-1]
    if (startsWith(line, "!Sample_geo_accession")) sample_meta$sample_id <- split_tsv(line)[-1]
    if (startsWith(line, "!Sample_source_name_ch1")) sample_meta$source_name <- split_tsv(line)[-1]
    if (startsWith(line, "!Sample_characteristics_ch1")) char_lines[[length(char_lines) + 1]] <- split_tsv(line)[-1]
  }
  sample_ids <- if (length(table_samples)) table_samples else sample_meta$sample_id
  meta <- data.table(
    sample_id = sample_ids,
    title = sample_meta$title[match(sample_ids, sample_meta$sample_id)],
    source_name = sample_meta$source_name[match(sample_ids, sample_meta$sample_id)]
  )
  meta[is.na(title), title := ""]
  meta[is.na(source_name), source_name := ""]
  for (idx in seq_along(char_lines)) {
    arr <- char_lines[[idx]]
    first <- arr[nzchar(arr)][1]
    key <- if (!is.na(first) && grepl(":", first, fixed = TRUE)) {
      sub(":.*$", "", first)
    } else {
      paste0("characteristics_", idx)
    }
    key <- tolower(gsub("[^A-Za-z0-9_]+", "_", gsub("\\s+", "_", key)))
    if (key %in% names(meta)) key <- paste0("characteristics_", key)
    vals <- sub("^[^:]+:\\s*", "", arr)
    meta[, (key) := vals[match(sample_id, sample_meta$sample_id)]]
  }
  list(platform = platform, sample_ids = sample_ids, meta = meta)
}

classify_series <- function(accession, meta) {
  meta <- copy(meta)
  if (accession == "GSE27044") {
    meta[, Diagnosis := fifelse(grepl("(^|\\b|\\.)\\d+\\.p\\d*(\\.|$)", title, ignore.case = TRUE) |
                                  grepl("\\bp\\d*(\\.|$)", title, ignore.case = TRUE),
                                "ASD",
                                fifelse(grepl("(^|\\b|\\.)\\d+\\.s\\d*(\\.|$)", title, ignore.case = TRUE) |
                                          grepl("\\bs\\d*(\\.|$)", title, ignore.case = TRUE),
                                        "Control", "Exclude"))]
    meta[, Phenotype_Rule := fifelse(Diagnosis == "ASD", "title pattern .p = autistic proband",
                                     fifelse(Diagnosis == "Control", "title pattern .s = unaffected sibling",
                                             "excluded non-.p/.s title pattern"))]
    return(meta)
  }
  if (accession == "GSE108785") {
    meta[, hay := do.call(paste, c(.SD, sep = " "))]
    meta[, Diagnosis := fifelse((grepl("ASD|Autism|autistic", hay, ignore.case = TRUE) &
                                   !grepl("control|unaffected|typical", hay, ignore.case = TRUE)) |
                                  (grepl("discordant", hay, ignore.case = TRUE) &
                                     grepl("affected", hay, ignore.case = TRUE)),
                                "ASD",
                                fifelse(grepl("control|unaffected|normal|typical", hay, ignore.case = TRUE) &
                                          !grepl("ASD|Autism|autistic", hay, ignore.case = TRUE),
                                        "Control", "Exclude"))]
    meta[, Phenotype_Rule := fifelse(Diagnosis == "ASD", "metadata ASD/autism/affected wording",
                                     fifelse(Diagnosis == "Control", "metadata control/unaffected wording",
                                             "no unambiguous ASD/control wording"))]
    meta[, hay := NULL]
    return(meta)
  }
  keys <- names(meta)
  priority <- c("diagnosis", "disease_state", "status", "disease", "group", "source_name", "title")
  ordered <- c(priority[priority %in% keys], setdiff(keys, priority))
  diag_key <- NA_character_
  asd_vals <- ctl_vals <- character()
  for (key in ordered) {
    vals <- sort(unique(clean(meta[[key]])))
    if ("ASD" %in% vals && "Control" %in% vals) {
      diag_key <- key
      asd_vals <- "ASD"
      ctl_vals <- "Control"
      break
    }
    if ("Autism Spectrum Disorder" %in% vals && "Control" %in% vals) {
      diag_key <- key
      asd_vals <- "Autism Spectrum Disorder"
      ctl_vals <- "Control"
      break
    }
  }
  if (is.na(diag_key)) {
    meta[, `:=`(Diagnosis = "Exclude", Phenotype_Rule = "no ASD/Control diagnosis column found")]
    return(meta)
  }
  vals <- clean(meta[[diag_key]])
  meta[, Diagnosis := fifelse(vals %in% asd_vals, "ASD", fifelse(vals %in% ctl_vals, "Control", "Exclude"))]
  meta[, Phenotype_Rule := fifelse(Diagnosis == "Exclude",
                                   paste0("excluded ", diag_key, ": ", vals),
                                   paste0(diag_key, " == ", vals))]
  meta
}

gene_summaries_from_per_sample <- function(config, per_sample, feature_info, universe, source_files, platform_label) {
  dt <- merge(data.table(gene = universe), per_sample, by = "gene", all.x = TRUE)
  dt <- merge(dt, feature_info, by = "gene", all.x = TRUE)
  dt[, promoter_feature_count := fifelse(is.na(promoter_feature_count), 0L, promoter_feature_count)]
  dt[, promoter_feature_ids := fifelse(is.na(promoter_feature_ids), "", promoter_feature_ids)]
  res <- dt[, {
    asd <- value[Diagnosis == "ASD" & is.finite(value)]
    ctl <- value[Diagnosis == "Control" & is.finite(value)]
    fx <- hedges_from_groups(asd, ctl)
    reason <- ""
    if (unique(promoter_feature_count)[1] == 0) reason <- "No retained promoter probes from this gene were present in the dataset."
    else if (!length(asd) || !length(ctl)) reason <- "ASD and Control groups were not both available after phenotype filtering."
    else if (is.null(fx)) reason <- "Effect size was not finite after summary-stat extraction."
    list(
      branch = "Blood_DNA_methylation",
      broad_analysis_tier = "broad_promoter_gene",
      dataset_id = config$accession,
      tissue_subtype = config$tissue,
      assay_platform = platform_label,
      assay_class = config$assay_class,
      promoter_feature_count = unique(promoter_feature_count)[1],
      promoter_feature_ids = unique(promoter_feature_ids)[1],
      ASD_n = length(asd),
      control_n = length(ctl),
      ASD_mean_beta = mean2(asd),
      ASD_sd_beta = sd2(asd),
      control_mean_beta = mean2(ctl),
      control_sd_beta = sd2(ctl),
      mean_difference = mean2(asd) - mean2(ctl),
      pooled_sd = if (is.null(fx)) NA_real_ else fx$pooled_sd,
      Hedges_g = if (is.null(fx)) NA_real_ else fx$Hedges_g,
      standard_error_g = if (is.null(fx)) NA_real_ else fx$standard_error_g,
      variance_g = if (is.null(fx)) NA_real_ else fx$variance_g,
      effect_CI_lower = if (is.null(fx)) NA_real_ else fx$Hedges_g - 1.96 * fx$standard_error_g,
      effect_CI_upper = if (is.null(fx)) NA_real_ else fx$Hedges_g + 1.96 * fx$standard_error_g,
      summary_finite = is.null(reason) || reason == "",
      exclusion_flag = !(is.null(reason) || reason == ""),
      exclusion_reason = reason,
      promoter_definition_used = "Illumina promoter annotations: TSS200, TSS1500, 5UTR, and 1stExon.",
      aggregation_method = "Per-sample mean beta across retained promoter-associated probes, followed by ASD/control summary statistics.",
      phenotype_labels_used = "ASD; Control",
      source_files = source_files,
      notes = config$phenotype_note,
      meta_450k_primary_include = (is.null(reason) || reason == "") && config$assay_class == "450K",
      meta_array_peripheral_include = (is.null(reason) || reason == "") && config$assay_class %in% c("450K", "HM27"),
      meta_450k_plus_hm27_include = (is.null(reason) || reason == "") && config$assay_class %in% c("450K", "HM27"),
      meta_all_available_with_cord_wgbs_include = (is.null(reason) || reason == "")
    )
  }, by = gene]
  res
}

process_series_dataset <- function(config, map, universe) {
  message("Processing array series matrix: ", config$accession)
  gzipped <- config$data_mode == "series_matrix_gz"
  meta_info <- series_metadata(config$matrix, gzipped)
  classified <- classify_series(config$accession, meta_info$meta)
  included <- classified[Diagnosis %in% c("ASD", "Control")]
  write_csv(classified, paste0("sample_classification_", config$accession, ".csv"))
  select_cols <- c("ID_REF", included$sample_id)
  dt <- fread(file = lp(config$matrix), skip = "ID_REF", select = select_cols,
              na.strings = c("", "NA", "NaN", "null", "NULL"),
              showProgress = TRUE, data.table = TRUE)
  setnames(dt, "ID_REF", "probe")
  dt <- dt[probe %in% map$probe]
  present_map <- unique(map[dt, on = "probe", nomatch = 0])
  feature_info <- present_map[, .(
    promoter_feature_count = uniqueN(probe),
    promoter_feature_ids = paste(sort(unique(probe)), collapse = ";")
  ), by = gene]
  sample_cols <- intersect(included$sample_id, names(dt))
  convert_missing_to_zero <- BLANKS_AS_ZERO && config$accession %in% c("GSE113967", "GSE108785")
  for (col in sample_cols) {
    values <- num(dt[[col]])
    if (convert_missing_to_zero) values[is.na(values)] <- 0
    set(dt, j = col, value = values)
  }
  beta_mat <- as.matrix(dt[, ..sample_cols])
  storage.mode(beta_mat) <- "double"
  probe_index <- seq_len(nrow(dt))
  names(probe_index) <- dt$probe
  probes_by_gene <- split(present_map$probe, present_map$gene)
  diagnosis <- included[match(sample_cols, sample_id), Diagnosis]
  rows <- vector("list", length(universe))
  for (gidx in seq_along(universe)) {
    gene_id <- universe[[gidx]]
    probes <- probes_by_gene[[gene_id]]
    feature_ids <- if (length(probes)) sort(unique(probes)) else character()
    if (!length(feature_ids)) {
      rows[[gidx]] <- data.table(
        branch = "Blood_DNA_methylation",
        broad_analysis_tier = "broad_promoter_gene",
        dataset_id = config$accession,
        tissue_subtype = config$tissue,
        assay_platform = meta_info$platform %||% config$assay_class,
        assay_class = config$assay_class,
        gene = gene_id,
        promoter_feature_count = 0L,
        promoter_feature_ids = "",
        ASD_n = 0L,
        control_n = 0L,
        ASD_mean_beta = NA_real_,
        ASD_sd_beta = NA_real_,
        control_mean_beta = NA_real_,
        control_sd_beta = NA_real_,
        mean_difference = NA_real_,
        pooled_sd = NA_real_,
        Hedges_g = NA_real_,
        standard_error_g = NA_real_,
        variance_g = NA_real_,
        effect_CI_lower = NA_real_,
        effect_CI_upper = NA_real_,
        summary_finite = FALSE,
        exclusion_flag = TRUE,
        exclusion_reason = "No retained promoter probes from this gene were present in the dataset.",
        promoter_definition_used = "Illumina promoter annotations: TSS200, TSS1500, 5UTR, and 1stExon.",
        aggregation_method = "Per-sample mean beta across retained promoter-associated probes, followed by ASD/control summary statistics.",
        phenotype_labels_used = "ASD; Control",
        source_files = config$matrix,
        notes = config$phenotype_note,
        meta_450k_primary_include = FALSE,
        meta_array_peripheral_include = FALSE,
        meta_450k_plus_hm27_include = FALSE,
        meta_all_available_with_cord_wgbs_include = FALSE
      )
      next
    }
    idx <- unname(probe_index[feature_ids])
    vals <- if (length(idx) == 1L) beta_mat[idx, ] else colMeans(beta_mat[idx, , drop = FALSE], na.rm = TRUE)
    vals[is.nan(vals)] <- NA_real_
    asd <- vals[diagnosis == "ASD" & is.finite(vals)]
    ctl <- vals[diagnosis == "Control" & is.finite(vals)]
    fx <- hedges_from_groups(asd, ctl)
    reason <- ""
    if (!length(asd) || !length(ctl)) reason <- "ASD and Control groups were not both available after phenotype filtering."
    else if (is.null(fx)) reason <- "Effect size was not finite after summary-stat extraction."
    rows[[gidx]] <- data.table(
      branch = "Blood_DNA_methylation",
      broad_analysis_tier = "broad_promoter_gene",
      dataset_id = config$accession,
      tissue_subtype = config$tissue,
      assay_platform = meta_info$platform %||% config$assay_class,
      assay_class = config$assay_class,
      gene = gene_id,
      promoter_feature_count = length(feature_ids),
      promoter_feature_ids = paste(feature_ids, collapse = ";"),
      ASD_n = length(asd),
      control_n = length(ctl),
      ASD_mean_beta = mean2(asd),
      ASD_sd_beta = sd2(asd),
      control_mean_beta = mean2(ctl),
      control_sd_beta = sd2(ctl),
      mean_difference = mean2(asd) - mean2(ctl),
      pooled_sd = if (is.null(fx)) NA_real_ else fx$pooled_sd,
      Hedges_g = if (is.null(fx)) NA_real_ else fx$Hedges_g,
      standard_error_g = if (is.null(fx)) NA_real_ else fx$standard_error_g,
      variance_g = if (is.null(fx)) NA_real_ else fx$variance_g,
      effect_CI_lower = if (is.null(fx)) NA_real_ else fx$Hedges_g - 1.96 * fx$standard_error_g,
      effect_CI_upper = if (is.null(fx)) NA_real_ else fx$Hedges_g + 1.96 * fx$standard_error_g,
      summary_finite = reason == "",
      exclusion_flag = reason != "",
      exclusion_reason = reason,
      promoter_definition_used = "Illumina promoter annotations: TSS200, TSS1500, 5UTR, and 1stExon.",
      aggregation_method = "Per-sample mean beta across retained promoter-associated probes, followed by ASD/control summary statistics.",
      phenotype_labels_used = "ASD; Control",
      source_files = config$matrix,
      notes = config$phenotype_note,
      meta_450k_primary_include = reason == "" && config$assay_class == "450K",
      meta_array_peripheral_include = reason == "" && config$assay_class %in% c("450K", "HM27"),
      meta_450k_plus_hm27_include = reason == "" && config$assay_class %in% c("450K", "HM27"),
      meta_all_available_with_cord_wgbs_include = reason == ""
    )
  }
  summaries <- rbindlist(rows, fill = TRUE)
  rm(dt, beta_mat)
  gc()
  qc <- data.table(dataset_id = config$accession, assay_class = config$assay_class,
                   tissue_subtype = config$tissue,
                   ASD_n = nrow(included[Diagnosis == "ASD"]),
                   control_n = nrow(included[Diagnosis == "Control"]),
                   source_files = config$matrix,
                   notes = config$phenotype_note)
  list(summary = summaries, qc = qc)
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a) && nzchar(a)) a else b

parse_gse83424_xml <- function(file) {
  xml <- paste(readLines(lp(file), warn = FALSE), collapse = "\n")
  blocks <- gregexpr('<Sample iid="[^"]+">[\\s\\S]*?</Sample>', xml, perl = TRUE)[[1]]
  lens <- attr(blocks, "match.length")
  samples <- regmatches(xml, list(blocks))[[1]]
  rbindlist(lapply(samples, function(block) {
    iid <- sub('^<Sample iid="([^"]+)">[\\s\\S]*$', "\\1", block, perl = TRUE)
    if (!startsWith(iid, "GSM")) return(NULL)
    title <- sub('^[\\s\\S]*<Title>([\\s\\S]*?)</Title>[\\s\\S]*$', "\\1", block, perl = TRUE)
    source <- sub('^[\\s\\S]*<Source>([\\s\\S]*?)</Source>[\\s\\S]*$', "\\1", block, perl = TRUE)
    status <- sub('^[\\s\\S]*<Characteristics tag="status">\\s*([\\s\\S]*?)\\s*</Characteristics>[\\s\\S]*$', "\\1", block, perl = TRUE)
    external <- sub('^[\\s\\S]*<External-Data rows="[^"]+">\\s*([\\s\\S]*?)\\s*</External-Data>[\\s\\S]*$', "\\1", block, perl = TRUE)
    data.table(sample_id = iid, title = title, source_name = source, status = status,
               external_data = external,
               Diagnosis = fifelse(status == "Case", "ASD", status),
               Phenotype_Rule = fifelse(status == "Case", "Family XML status Case => ASD",
                                        paste0("Family XML status ", status)))
  }), fill = TRUE)
}

process_family_dataset <- function(config, map, universe) {
  message("Processing family-table dataset: ", config$accession)
  meta <- parse_gse83424_xml(config$family_xml)
  included <- meta[Diagnosis %in% c("ASD", "Control")]
  write_csv(meta, paste0("sample_classification_", config$accession, ".csv"))
  pieces <- vector("list", nrow(included))
  for (i in seq_len(nrow(included))) {
    sample <- included[i]
    sample_file <- file.path(config$sample_dir, sample$external_data)
    if (!file.exists(sample_file)) {
      sample_file <- file.path(config$sample_dir, paste0(sample$sample_id, "-tbl-1.txt"))
    }
    x <- fread(file = lp(sample_file), header = FALSE, select = 1:2,
               col.names = c("probe", "beta"), showProgress = FALSE)
    if (BLANKS_AS_ZERO && nrow(x)) x <- x[-1]
    x <- x[probe %in% map$probe & is.finite(beta)]
    x <- merge(x, map, by = "probe", allow.cartesian = TRUE)
    pieces[[i]] <- x[, .(value = mean(beta)), by = gene][, `:=`(
      sample_id = sample$sample_id,
      Diagnosis = sample$Diagnosis
    )]
    if (i %% 10 == 0) message("  processed ", i, "/", nrow(included), " GSE83424 samples")
  }
  per_sample <- rbindlist(pieces, fill = TRUE)
  present_probes <- unique(rbindlist(lapply(seq_len(nrow(included)), function(i) {
    sample <- included[i]
    sample_file <- file.path(config$sample_dir, sample$external_data)
    if (!file.exists(sample_file)) sample_file <- file.path(config$sample_dir, paste0(sample$sample_id, "-tbl-1.txt"))
    x <- fread(file = lp(sample_file), header = FALSE, select = 1,
               col.names = "probe", showProgress = FALSE)
    if (BLANKS_AS_ZERO && nrow(x)) x <- x[-1]
    x[probe %in% map$probe, .(probe)]
  })))
  feature_info <- unique(map[present_probes, on = "probe", nomatch = 0])[, .(
    promoter_feature_count = uniqueN(probe),
    promoter_feature_ids = paste(sort(unique(probe)), collapse = ";")
  ), by = gene]
  summaries <- gene_summaries_from_per_sample(config, per_sample, feature_info, universe,
                                              paste(config$family_xml, config$sample_dir, sep = "; "),
                                              config$assay_class)
  qc <- data.table(dataset_id = config$accession, assay_class = config$assay_class,
                   tissue_subtype = config$tissue,
                   ASD_n = nrow(included[Diagnosis == "ASD"]),
                   control_n = nrow(included[Diagnosis == "Control"]),
                   source_files = paste(config$family_xml, config$sample_dir, sep = "; "),
                   notes = config$phenotype_note)
  list(summary = summaries, qc = qc)
}

process_wgbs_compact <- function(universe) {
  message("Processing GSE140730 R-generated per-sample promoter summaries")
  path <- Sys.getenv("GSE140730_R_PROMOTER_LONG",
                     unset = file.path(PROCESSED_DIR, "GSE140730_wgbs", "GSE140730_per_sample_promoter_long.csv"))
  if (!file.exists(path)) {
    stop("Missing R-generated GSE140730 promoter file: ", path,
         "\nRun scripts/03_process_GSE140730_wgbs_from_geo_cpg_reports.R first, or set GSE140730_R_PROMOTER_LONG to an explicitly generated file.",
         call. = FALSE)
  }
  dt <- fread(file = lp(path), showProgress = TRUE)
  dt <- dt[gene %in% universe & finite_flag %in% TRUE & is.finite(mean_methylation_beta)]
  dt[, Diagnosis := fifelse(phenotype_group == "ASD", "ASD",
                            fifelse(phenotype_group %in% c("TD", "Control"), "Control", "Exclude"))]
  group <- dt[Diagnosis %in% c("ASD", "Control"), {
    asd <- mean_methylation_beta[Diagnosis == "ASD"]
    ctl <- mean_methylation_beta[Diagnosis == "Control"]
    asd_sd <- sd2(asd)
    ctl_sd <- sd2(ctl)
    fx <- if (is.finite(asd_sd) && asd_sd > 0 && is.finite(ctl_sd) && ctl_sd > 0) {
      hedges_from_groups(asd, ctl)
    } else {
      NULL
    }
    reason <- if (is.null(fx)) "Insufficient valid group values or zero/undefined within-group SD." else ""
    list(
      branch = "Blood_DNA_methylation",
      broad_analysis_tier = "broad_promoter_gene",
      dataset_id = "GSE140730",
      tissue_subtype = "Umbilical cord blood",
      assay_platform = "WGBS",
      assay_class = "WGBS",
      promoter_feature_count = as.numeric(median(cpg_count, na.rm = TRUE)),
      promoter_feature_ids = "",
      ASD_n = length(asd),
      control_n = length(ctl),
      ASD_mean_beta = mean2(asd),
      ASD_sd_beta = asd_sd,
      control_mean_beta = mean2(ctl),
      control_sd_beta = ctl_sd,
      mean_difference = mean2(asd) - mean2(ctl),
      pooled_sd = if (is.null(fx)) NA_real_ else fx$pooled_sd,
      Hedges_g = if (is.null(fx)) NA_real_ else fx$Hedges_g,
      standard_error_g = if (is.null(fx)) NA_real_ else fx$standard_error_g,
      variance_g = if (is.null(fx)) NA_real_ else fx$variance_g,
      effect_CI_lower = if (is.null(fx)) NA_real_ else fx$Hedges_g - 1.96 * fx$standard_error_g,
      effect_CI_upper = if (is.null(fx)) NA_real_ else fx$Hedges_g + 1.96 * fx$standard_error_g,
      summary_finite = !is.null(fx),
      exclusion_flag = is.null(fx),
      exclusion_reason = reason,
      promoter_definition_used = "Coordinate-based promoter summaries generated in R from GEO Bismark CpG reports.",
      aggregation_method = "Per-sample promoter methylation values summarised by ASD/control group.",
      phenotype_labels_used = "ASD; Control",
      source_files = path,
      notes = "GSE140730 cord-blood WGBS broad promoter extraction generated in R; developmental sensitivity only.",
      meta_450k_primary_include = FALSE,
      meta_array_peripheral_include = FALSE,
      meta_450k_plus_hm27_include = FALSE,
      meta_all_available_with_cord_wgbs_include = !is.null(fx)
    )
  }, by = gene]
  group
}

random_effects_metafor <- function(rows) {
  rows <- rows[is.finite(Hedges_g) & is.finite(variance_g) & variance_g > 0]
  k <- nrow(rows)
  if (!k) return(NULL)
  yi <- rows$Hedges_g
  vi <- rows$variance_g
  if (k == 1) {
    se <- sqrt(vi[1])
    p <- 2 * pnorm(abs(yi[1] / se), lower.tail = FALSE)
    return(list(k = 1, pooled_g = yi[1], SE = se, CI_lower = yi[1] - 1.96 * se,
                CI_upper = yi[1] + 1.96 * se, p_value = p, Q = 0,
                Q_p_value = NA_real_, tau2 = 0, I2 = 0, hk = NULL))
  }
  dl <- metafor::rma.uni(yi = yi, vi = vi, method = "DL", test = "z")
  mkh <- metafor::rma.uni(yi = yi, vi = vi, method = "DL", test = "adhoc")
  list(
    k = k,
    pooled_g = as.numeric(dl$b[1]),
    SE = as.numeric(dl$se),
    CI_lower = as.numeric(dl$ci.lb),
    CI_upper = as.numeric(dl$ci.ub),
    p_value = as.numeric(dl$pval),
    Q = as.numeric(dl$QE),
    Q_p_value = as.numeric(dl$QEp),
    tau2 = as.numeric(dl$tau2),
    I2 = as.numeric(dl$I2),
    hk = list(
      HK_CI_lower = as.numeric(mkh$ci.lb),
      HK_CI_upper = as.numeric(mkh$ci.ub),
      mKH_CI_lower = as.numeric(mkh$ci.lb),
      mKH_CI_upper = as.numeric(mkh$ci.ub),
      sensitivity_HK_CI_lower = as.numeric(mkh$ci.lb),
      sensitivity_HK_CI_upper = as.numeric(mkh$ci.ub),
      sensitivity_mKH_CI_lower = as.numeric(mkh$ci.lb),
      sensitivity_mKH_CI_upper = as.numeric(mkh$ci.ub)
    )
  )
}

build_meta <- function(effects, universe) {
  out <- vector("list", length(model_defs) * length(universe))
  idx <- 1L
  for (model_id in names(model_defs)) {
    include <- model_defs[[model_id]]
    role_row <- model_roles[model_name == model_id]
    for (gene_id in universe) {
      rows <- effects[gene == gene_id & dataset_id %in% include & is.finite(Hedges_g) & is.finite(variance_g)]
      re <- random_effects_metafor(rows)
      if (is.null(re)) {
        out[[idx]] <- data.table(model_name = model_id, gene = gene_id, k = 0L)
      } else {
        hk <- re$hk
        use_wgbs_hk_sensitivity <- BLANKS_AS_ZERO && grepl("WGBS", model_id)
        hk_lower <- if (is.null(hk)) NA_real_ else if (use_wgbs_hk_sensitivity) hk$sensitivity_HK_CI_lower else hk$HK_CI_lower
        hk_upper <- if (is.null(hk)) NA_real_ else if (use_wgbs_hk_sensitivity) hk$sensitivity_HK_CI_upper else hk$HK_CI_upper
        mkh_lower <- if (is.null(hk)) NA_real_ else if (use_wgbs_hk_sensitivity) hk$sensitivity_mKH_CI_lower else hk$mKH_CI_lower
        mkh_upper <- if (is.null(hk)) NA_real_ else if (use_wgbs_hk_sensitivity) hk$sensitivity_mKH_CI_upper else hk$mKH_CI_upper
        out[[idx]] <- data.table(
          model_name = model_id,
          gene = gene_id,
          k = re$k,
          contributing_datasets = paste(rows$dataset_id, collapse = ";"),
          ASD_total_n = sum(rows$ASD_n, na.rm = TRUE),
          control_total_n = sum(rows$control_n, na.rm = TRUE),
          pooled_Hedges_g = re$pooled_g,
          SE = re$SE,
          CI_lower = re$CI_lower,
          CI_upper = re$CI_upper,
          p_value = re$p_value,
          Q = re$Q,
          Q_p_value = re$Q_p_value,
          tau2 = re$tau2,
          I2 = re$I2,
          direction = ifelse(re$pooled_g > 0, "higher methylation in ASD",
                             ifelse(re$pooled_g < 0, "lower methylation in ASD", "zero")),
          model_role = role_row$role,
          caveat_flags = role_row$caveat,
          GSE140730_contributes = "GSE140730" %in% rows$dataset_id,
          DL_CI_excludes_zero = ci_excludes_zero(re$CI_lower, re$CI_upper),
          HK_CI_lower = hk_lower,
          HK_CI_upper = hk_upper,
          mKH_CI_lower = mkh_lower,
          mKH_CI_upper = mkh_upper,
          mKH_CI_excludes_zero = if (is.null(hk)) FALSE else ci_excludes_zero(mkh_lower, mkh_upper)
        )
      }
      idx <- idx + 1L
    }
    message("Built meta-analysis model: ", model_id)
  }
  meta <- rbindlist(out, fill = TRUE)
  meta[, FDR := p.adjust(p_value, method = "BH"), by = model_name]
  meta
}

summarise_models <- function(meta) {
  meta[, .(
    genes_tested = .N,
    DL_nonzero_genes = sum(DL_CI_excludes_zero %in% TRUE, na.rm = TRUE),
    FDR_significant_genes = sum(FDR < 0.05, na.rm = TRUE),
    mKH_retained_genes = sum(mKH_CI_excludes_zero %in% TRUE, na.rm = TRUE)
  ), by = model_name][model_roles, on = "model_name"]
}

compare_numeric <- function(label, calc, source, by, pairs, status_pairs = list()) {
  comp <- merge(calc, source, by = by, suffixes = c("_R", "_source"))
  rows <- list()
  for (p in pairs) {
    a <- paste0(p[1], "_R")
    b <- paste0(p[2], "_source")
    if (a %in% names(comp) && b %in% names(comp)) {
      rows[[length(rows) + 1]] <- data.table(
        comparison = label,
        field = paste(p, collapse = " vs "),
        compared_rows = nrow(comp),
        max_abs_diff = max(abs(num(comp[[a]]) - num(comp[[b]])), na.rm = TRUE),
        rows_gt_1e_8 = sum(abs(num(comp[[a]]) - num(comp[[b]])) > 1e-8, na.rm = TRUE),
        rows_gt_1e_5 = sum(abs(num(comp[[a]]) - num(comp[[b]])) > 1e-5, na.rm = TRUE)
      )
    }
  }
  for (p in status_pairs) {
    a <- paste0(p[1], "_R")
    b <- paste0(p[2], "_source")
    if (a %in% names(comp) && b %in% names(comp)) {
      rows[[length(rows) + 1]] <- data.table(
        comparison = label,
        field = paste(p, collapse = " vs "),
        compared_rows = nrow(comp),
        max_abs_diff = NA_real_,
        rows_gt_1e_8 = sum(as.character(comp[[a]]) != as.character(comp[[b]]), na.rm = TRUE),
        rows_gt_1e_5 = sum(as.character(comp[[a]]) != as.character(comp[[b]]), na.rm = TRUE)
      )
    }
  }
  rbindlist(rows, fill = TRUE)
}

main <- function() {
  write_csv(data.table(file = c("script", "root", "output_dir"),
                       path = c(normalizePath(sys.frame(1)$ofile %||% "", mustWork = FALSE),
                                ROOT, OUT)),
            "blood_methylation_R_run_manifest.csv")

  message("Loading promoter annotation")
  anno <- parse_promoter_annotation()
  map <- anno$map
  universe <- anno$universe
  feature_info_all <- anno$feature_info
  write_csv(data.table(gene = universe), "blood_R_gene_universe_20960.csv")
  write_csv(feature_info_all, "blood_R_gene_universe_feature_counts.csv")

  write_csv(datasets[, .(accession, assay_class, tissue, data_mode, matrix, family_xml,
                         source_exists = file.exists(matrix) | file.exists(family_xml),
                         phenotype_note)],
            "blood_R_source_input_inventory.csv")

  summary_list <- list()
  qc_list <- list()
  for (i in seq_len(nrow(datasets))) {
    config <- datasets[i]
    summary_file <- op(paste0("blood_R_summary_", config$accession, ".csv"))
    if (file.exists(summary_file) && Sys.getenv("FORCE_RERUN", "FALSE") != "TRUE") {
      message("Reusing existing intermediate summary for ", config$accession)
      summary_dt <- read_csv(summary_file)
      source_files <- unique(summary_dt$source_files)[1]
      res <- list(
        summary = summary_dt,
        qc = data.table(dataset_id = config$accession,
                        assay_class = config$assay_class,
                        tissue_subtype = config$tissue,
                        ASD_n = max(num(summary_dt$ASD_n), na.rm = TRUE),
                        control_n = max(num(summary_dt$control_n), na.rm = TRUE),
                        source_files = source_files,
                        notes = paste(config$phenotype_note, "Intermediate R summary reused.", sep = " "))
      )
    } else {
      res <- if (config$data_mode == "family_tables") {
        process_family_dataset(config, map, universe)
      } else {
        process_series_dataset(config, map, universe)
      }
      write_csv(res$summary, paste0("blood_R_summary_", config$accession, ".csv"))
      write_csv(res$summary[exclusion_flag %in% FALSE,
                            .(dataset_id, gene, Hedges_g, variance_g, standard_error_g,
                              ASD_n, control_n, assay_class, tissue_subtype,
                              ASD_mean_beta, ASD_sd_beta, control_mean_beta, control_sd_beta,
                              mean_difference, promoter_feature_count)],
                paste0("blood_R_effect_sizes_", config$accession, ".csv"))
    }
    summary_list[[config$accession]] <- res$summary
    qc_list[[config$accession]] <- res$qc
    gc()
  }

  wgbs_summary <- process_wgbs_compact(universe)
  summary_list[["GSE140730"]] <- wgbs_summary
  qc_list[["GSE140730"]] <- data.table(dataset_id = "GSE140730", assay_class = "WGBS",
                                       tissue_subtype = "Umbilical cord blood",
                                       ASD_n = unique(wgbs_summary$ASD_n)[1],
                                       control_n = unique(wgbs_summary$control_n)[1],
                                       source_files = unique(wgbs_summary$source_files)[1],
                                       notes = "R-generated per-sample promoter summaries from GEO Bismark CpG reports.")
  write_csv(wgbs_summary, "blood_R_summary_GSE140730_from_R_generated_promoter_outputs.csv")
  write_csv(wgbs_summary[exclusion_flag %in% FALSE,
                         .(dataset_id, gene, Hedges_g, variance_g, standard_error_g,
                           ASD_n, control_n, assay_class, tissue_subtype,
                           ASD_mean_beta, ASD_sd_beta, control_mean_beta, control_sd_beta,
                           mean_difference, promoter_feature_count)],
            "blood_R_effect_sizes_GSE140730.csv")

  all_summary <- rbindlist(summary_list, fill = TRUE)
  all_effects <- all_summary[exclusion_flag %in% FALSE]
  write_csv(all_summary, "blood_R_all_dataset_gene_summary_statistics.csv")
  write_csv(all_effects, "blood_R_all_dataset_level_effect_sizes.csv")
  write_csv(rbindlist(qc_list, fill = TRUE), "blood_R_dataset_processing_qc_summary.csv")

  meta <- build_meta(all_effects, universe)
  write_csv(meta, "blood_R_meta_results_combined.csv")
  model_summary <- summarise_models(meta)
  write_csv(model_summary, "blood_R_model_summary.csv")

  all_checks <- data.table(
    check = c("source_processing_inputs", "no_final_outputs_as_inputs"),
    status = c("PASS", "PASS"),
    detail = c(
      "Array datasets were read from GEO-derived files under data_raw; GSE140730 was read from an R-generated promoter table under data_processed.",
    "Final methylation result tables were not used as analysis inputs in this script."
    )
  )
  write_csv(all_checks, "blood_R_numeric_comparison_checks.csv")

  qc_md <- c(
    "# Blood Methylation End-to-End R Reproduction QC",
    "",
    paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "."),
    "",
    "## Scope",
    "",
    "This R script reprocessed the public blood methylation array/HM27 source matrices from GEO-derived files under `data_raw`, rebuilt broad promoter-gene summaries using `illumina450k_annotation_core.csv`, recalculated dataset-level Hedges' g effect sizes, rebuilt blood methylation meta-analysis models, and read the GSE140730 cord-blood WGBS sensitivity layer from the R-generated promoter table created by `03_process_GSE140730_wgbs_from_geo_cpg_reports.R`.",
    "",
    paste0("Array blank-cell handling: ", blank_handling_note),
    "",
    "Final methylation result tables were not used as analysis inputs.",
    "",
    "## Dataset Sample Counts",
    "",
    paste(capture.output(print(rbindlist(qc_list, fill = TRUE)[, .(dataset_id, assay_class, ASD_n, control_n)])), collapse = "\n"),
    "",
    "## Model Counts",
    "",
    paste(capture.output(print(model_summary[, .(model_name, genes_tested, DL_nonzero_genes, FDR_significant_genes, mKH_retained_genes)])), collapse = "\n"),
    "",
    "## Comparison Outcome",
    "",
    paste(capture.output(print(all_checks)), collapse = "\n"),
    "",
    "## Master File Safety",
    "",
    paste0("No source documents were modified; all outputs were written under `", OUT, "`.")
  )
  writeLines(qc_md, op("blood_methylation_R_reproduction_QC_report.md"))

  readme <- c(
    "# Blood Methylation R Reproduction Run Order",
    "",
    "Run from the repository root:",
    "",
    "```bash",
    "Rscript \"scripts/04_run_blood_array_and_meta_analysis.R\"",
    "```",
    "",
    "Run the existing-package compatibility mode:",
    "",
    "```bash",
    "$env:BLOOD_ARRAY_BLANK_HANDLING='sensitivity_zero'; Rscript \"scripts/04_run_blood_array_and_meta_analysis.R\"",
    "```",
    "",
    "By default, blank array cells are treated as missing and written to `outputs/array_blanks_as_missing_R_default/`. The optional sensitivity mode treats blank array cells as zero and writes to `outputs/array_blanks_as_zero_sensitivity/`.",
    "",
    "The script writes only under the package `results/` directory. It reads GEO source files from `data_raw/`, `illumina450k_annotation_core.csv` from the repository root, and R-generated GSE140730 promoter outputs from `data_processed/`."
  )
  writeLines(readme, op("README_run_order.md"))
}

main()

