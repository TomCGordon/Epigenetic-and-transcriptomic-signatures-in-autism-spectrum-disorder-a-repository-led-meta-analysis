#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})
if (identical(Sys.getenv("R_WARNINGS_AS_THEY_OCCUR"), "true")) {
  options(warn = 1)
}

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "lib", "placenta_lcl_methylation_functions.R"))

lcl_raw_dir <- file.path(package_root, "data_raw", "lcl")
processed_dir <- file.path(package_root, "data_processed", "lcl")
annotation_dir <- file.path(package_root, "data_processed", "annotation")
qc_dir <- file.path(package_root, "qc")
dir_create(processed_dir)
dir_create(qc_dir)

summarise_sample_matrix <- function(accession, gene_dt, sample_meta, feature_counts,
                                    assay_class, platform, value_label, notes) {
  long <- melt(gene_dt, id.vars = "gene", variable.name = "sample_id", value.name = "value")
  long <- merge(long, sample_meta[, .(sample_id, group, subject_id)], by = "sample_id", all.x = TRUE)
  long <- long[!is.na(group) & is.finite(value)]
  long <- long[, .(value = mean(value, na.rm = TRUE)), by = .(gene, group, subject_id)]
  asd <- long[group == "ASD", .(ASD_n = uniqueN(subject_id), ASD_mean = mean(value), ASD_SD = sd(value)), by = gene]
  ctl <- long[group == "Control", .(control_n = uniqueN(subject_id), control_mean = mean(value), control_SD = sd(value)), by = gene]
  out <- merge(asd, ctl, by = "gene", all = TRUE)
  out <- merge(out, feature_counts, by = "gene", all.x = TRUE)
  out[, `:=`(
    accession = accession,
    tissue_family = "LCL",
    assay_class = assay_class,
    platform = platform,
    finite_summary = fifelse(is.finite(ASD_mean) & is.finite(control_mean) &
                               is.finite(ASD_SD) & is.finite(control_SD) &
                               ASD_n >= 2 & control_n >= 2, "yes", "no"),
    notes = notes
  )]
  out[is.na(feature_count), feature_count := 0L]
  setcolorder(out, c("accession", "tissue_family", "assay_class", "platform", "gene"))
  effects <- hedges_effects(out, value_label = value_label)
  list(summary = out, effects = effects)
}

process_gse34099 <- function() {
  matrix_path <- file.path(lcl_raw_dir, "GSE34099_series_matrix.txt.gz")
  beta_dt <- read_geo_matrix_table(matrix_path)
  beta_dt <- beta_dt[grepl("^cg", ID_REF)]
  sample_cols <- setdiff(names(beta_dt), "ID_REF")
  beta_dt <- normalise_numeric_matrix(beta_dt, sample_cols)

  meta <- fread(file.path(qc_dir, "GSE34099_GEO_sample_manifest.csv"))
  meta <- meta[include == TRUE & group %in% c("ASD", "Control")]
  meta[, sample_id := sample_id]
  meta[, subject_id := sub("(_[0-9]+)? genomic DNA from LCL$", "", title)]
  meta[, subject_id := sub("_1$", "", subject_id)]
  meta <- meta[sample_id %in% sample_cols]

  hm27_map <- fread(file.path(annotation_dir, "GPL8490_HM27_probe_gene_map.csv"))
  keep <- hm27_map[probe %in% beta_dt$ID_REF]
  keep <- unique(keep[, .(probe, gene)])
  idx <- match(keep$probe, beta_dt$ID_REF)
  keep <- keep[is.finite(idx)]
  idx <- idx[is.finite(idx)]
  split_idx <- split(idx, keep$gene)
  feature_counts <- data.table(gene = names(split_idx), feature_count = lengths(split_idx))
  selected_samples <- meta$sample_id
  mat <- as.matrix(beta_dt[, ..selected_samples])
  storage.mode(mat) <- "numeric"
  gene_values <- lapply(names(split_idx), function(g) {
    vals <- colMeans(mat[split_idx[[g]], , drop = FALSE], na.rm = TRUE)
    vals[is.nan(vals)] <- NA_real_
    data.table(gene = g, t(vals))
  })
  gene_dt <- rbindlist(gene_values, fill = TRUE)
  names(gene_dt)[-1] <- selected_samples
  res <- summarise_sample_matrix(
    "GSE34099", gene_dt, meta, feature_counts,
    assay_class = "HM27",
    platform = "Illumina HumanMethylation27 BeadChip",
    value_label = "HM27 beta-like methylation",
    notes = "R workflow from public GEO series matrix; duplicate subject IDs collapsed before group summary."
  )
  write_csv(meta, file.path(processed_dir, "GSE34099_sample_classification.csv"))
  write_csv(res$summary, file.path(processed_dir, "GSE34099_gene_summary_statistics.csv"))
  write_csv(res$effects, file.path(processed_dir, "GSE34099_effect_sizes.csv"))
  data.table(accession = "GSE34099",
             ASD_n = max(res$summary$ASD_n, na.rm = TRUE),
             control_n = max(res$summary$control_n, na.rm = TRUE),
             finite_gene_summaries = sum(res$summary$finite_summary == "yes", na.rm = TRUE),
             effect_size_rows = nrow(res$effects))
}

process_gse99935 <- function() {
  matrix_path <- file.path(lcl_raw_dir, "GSE99935_Matrix_Transposed_normalized_MeDIP_data_bkgd_subtracted.txt.gz")
  sample_meta <- fread(file.path(qc_dir, "GSE99935_GEO_sample_manifest.csv"))
  sample_meta <- sample_meta[include == TRUE & group %in% c("ASD", "Control")]
  sample_meta[, matrix_col := sprintf("%03d", as.integer(subject_code))]
  sample_meta[, sample_id := matrix_col]
  sample_meta[, subject_id := matrix_col]

  promoters <- fread(file.path(annotation_dir, "hg18_refGene_promoter_coordinates.csv"))
  promoters <- promoters[!is.na(chr) & is.finite(start) & is.finite(end)]
  promoters[, `:=`(start = as.integer(start), end = as.integer(end))]
  setkey(promoters, chr, start, end)
  gene_list <- sort(unique(promoters$gene))
  gene_index <- setNames(seq_along(gene_list), gene_list)

  con <- gzfile(matrix_path, "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  header <- clean_quotes(header)
  header[1] <- "ID"
  selected <- intersect(sample_meta$matrix_col, header[-1])
  sample_meta <- sample_meta[matrix_col %in% selected]
  if (length(selected) != nrow(sample_meta)) stop("GSE99935 matrix/sample metadata mismatch.")
  sum_mat <- matrix(0, nrow = length(gene_list), ncol = length(selected),
                    dimnames = list(gene_list, selected))
  count_vec <- integer(length(gene_list))
  rows_read <- 0L
  rows_mapped <- 0L
  repeat {
    lines <- readLines(con, n = 50000, warn = FALSE)
    if (!length(lines)) break
    rows_read <- rows_read + length(lines)
    dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t",
                col.names = header, showProgress = FALSE)
    setnames(dt, names(dt)[1], "ID")
    dt <- normalise_numeric_matrix(dt, selected)
    dt <- copy(dt)
    dt[, chr := sub("\\.[0-9]+$", "", ID)]
    dt[, pos := suppressWarnings(as.integer(sub("^.*\\.", "", ID)))]
    dt <- dt[grepl("^chr", chr) & is.finite(pos)]
    if (!nrow(dt)) next
    dt[, row_id := .I]
    coords <- dt[, .(row_id, chr, start = pos, end = pos)]
    setkey(coords, chr, start, end)
    hits <- foverlaps(coords, promoters, nomatch = 0L)
    if (!nrow(hits)) next
    vals <- cbind(dt[, .(row_id)], dt[, ..selected])
    hv <- merge(hits[, .(row_id, gene)], vals, by = "row_id", allow.cartesian = TRUE)
    rows_mapped <- rows_mapped + uniqueN(hv$row_id)
    sums <- hv[, lapply(.SD, sum, na.rm = TRUE), by = gene, .SDcols = selected]
    counts <- hv[, .N, by = gene]
    gi <- gene_index[sums$gene]
    sum_mat[gi, ] <- sum_mat[gi, , drop = FALSE] + as.matrix(sums[, ..selected])
    count_vec[gene_index[counts$gene]] <- count_vec[gene_index[counts$gene]] + counts$N
    rm(dt, coords, hits, vals, hv, sums, counts)
    gc()
  }
  gene_dt <- as.data.table(sum_mat)
  gene_dt[, gene := gene_list]
  setcolorder(gene_dt, "gene")
  for (col in selected) {
    gene_dt[[col]] <- ifelse(count_vec > 0, gene_dt[[col]] / count_vec, NA_real_)
  }
  feature_counts <- data.table(gene = gene_list, feature_count = count_vec)
  res <- summarise_sample_matrix(
    "GSE99935", gene_dt, sample_meta, feature_counts,
    assay_class = "MeDIP-chip",
    platform = "Affymetrix Human Promoter 1.0R GeneChip",
    value_label = "MeDIP promoter-enrichment score",
    notes = "R workflow from public GEO normalized MeDIP-minus-input matrix; promoter-level coordinate aggregation uses hg18 refGene promoter windows."
  )
  write_csv(sample_meta, file.path(processed_dir, "GSE99935_sample_classification.csv"))
  write_csv(res$summary, file.path(processed_dir, "GSE99935_gene_summary_statistics.csv"))
  write_csv(res$effects, file.path(processed_dir, "GSE99935_effect_sizes.csv"))
  write_csv(data.table(accession = "GSE99935", rows_read = rows_read, rows_mapped = rows_mapped,
                       selected_samples = length(selected)),
            file.path(qc_dir, "GSE99935_coordinate_processing_QC.csv"))
  data.table(accession = "GSE99935",
             ASD_n = max(res$summary$ASD_n, na.rm = TRUE),
             control_n = max(res$summary$control_n, na.rm = TRUE),
             finite_gene_summaries = sum(res$summary$finite_summary == "yes", na.rm = TRUE),
             effect_size_rows = nrow(res$effects))
}

qc <- rbind(process_gse34099(), process_gse99935(), fill = TRUE)
write_csv(qc, file.path(qc_dir, "lcl_dataset_processing_QC.csv"))
all_summary <- rbindlist(lapply(c("GSE34099", "GSE99935"), function(id) {
  fread(file.path(processed_dir, paste0(id, "_gene_summary_statistics.csv")))
}), fill = TRUE)
all_effects <- rbindlist(lapply(c("GSE34099", "GSE99935"), function(id) {
  fread(file.path(processed_dir, paste0(id, "_effect_sizes.csv")))
}), fill = TRUE)
write_csv(all_summary, file.path(processed_dir, "lcl_all_dataset_gene_summary_statistics.csv"))
write_csv(all_effects, file.path(processed_dir, "lcl_all_dataset_effect_sizes.csv"))
message("LCL HM27 and MeDIP processing completed.")
