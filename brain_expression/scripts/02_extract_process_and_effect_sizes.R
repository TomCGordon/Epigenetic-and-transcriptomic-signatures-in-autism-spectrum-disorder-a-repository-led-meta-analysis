source(file.path("scripts", "lib", "brain_expression_functions.R"))
cfg <- brain_expression_config()

series <- list(
  list(dataset = "GSE113834", matrix = "GSE113834_series_matrix.txt.gz", platform_file = "GPL15207_family.soft.gz", platform = "GPL15207"),
  list(dataset = "GSE28475", matrix = "GSE28475-GPL13388_series_matrix.txt.gz", platform_file = "GPL13388_family.soft.gz", platform = "GPL13388"),
  list(dataset = "GSE28475", matrix = "GSE28475-GPL6883_series_matrix.txt.gz", platform_file = "GPL6883.annot.gz", platform = "GPL6883"),
  list(dataset = "GSE28521", matrix = "GSE28521_series_matrix.txt.gz", platform_file = "GPL6883.annot.gz", platform = "GPL6883"),
  list(dataset = "GSE36315", matrix = "GSE36315_series_matrix.txt.gz", platform_file = "GPL15314_family.soft.gz", platform = "GPL15314"),
  list(dataset = "GSE38322", matrix = "GSE38322_series_matrix.txt.gz", platform_file = "GPL10558.annot.gz", platform = "GPL10558")
)

results <- list()
failed <- list()
for (s in series) {
  message("Processing ", s$matrix)
  one <- try(process_series_dataset(cfg, s$dataset, s$matrix, s$platform_file, platform = s$platform), silent = TRUE)
  if (inherits(one, "try-error")) {
    failed[[length(failed) + 1]] <- data.frame(dataset = s$dataset, source_file = s$matrix, status = "failed", reason = as.character(one))
  } else {
    results[[length(results) + 1]] <- one
  }
}

message("Processing supplemental RNA-seq/source tables")
supp <- try(process_supplemental_datasets(cfg), silent = TRUE)
if (inherits(supp, "try-error")) {
  failed[[length(failed) + 1]] <- data.frame(dataset = "supplemental_brain_routes", source_file = "multiple", status = "failed", reason = as.character(supp))
} else {
  results <- c(results, supp)
}

labels <- rbindlist_fill(lapply(results, `[[`, "labels"))
effects_raw <- rbindlist_fill(lapply(results, `[[`, "effects"))
processed <- rbindlist_fill(lapply(results, `[[`, "processed"))

if (nrow(effects_raw)) {
  effects_collapsed <- collapse_within_study(effects_raw)
} else {
  effects_collapsed <- data.frame()
}

write_csv_safe(labels, file.path(cfg$out_dir, "02_sample_metadata", "brain_expression_R_sample_phenotype_table.csv"))
write_csv_safe(processed, file.path(cfg$out_dir, "03_processed_expression", "brain_expression_R_processed_matrix_summary.csv"))
write_csv_safe(effects_raw, file.path(cfg$out_dir, "05_effect_sizes", "brain_expression_R_dataset_level_effect_sizes_raw_strata.csv"))
write_csv_safe(effects_collapsed, file.path(cfg$out_dir, "05_effect_sizes", "brain_expression_R_dataset_level_effect_sizes_grouped_collapsed.csv"))
failed_out <- rbindlist_fill(failed)
if (nrow(failed_out)) {
  failed_out$status[grepl("reference|fewer than two ASD or control samples", failed_out$reason, ignore.case = TRUE)] <- "excluded_reference_or_not_case_control"
  failed_out$status[grepl("no valid probe-to-gene symbol mapping", failed_out$reason, ignore.case = TRUE)] <- "excluded_no_valid_gene_mapping"
}
write_csv_safe(failed_out, file.path(cfg$out_dir, "00_manifest", "brain_expression_R_processing_failures.csv"))

summary <- data.frame(
  dataset = processed$dataset,
  genes_with_effects = vapply(processed$dataset, function(ds) sum(effects_raw$dataset == ds), integer(1)),
  samples = processed$samples,
  matrix_genes = processed$genes,
  stringsAsFactors = FALSE
)
write_csv_safe(summary, file.path(cfg$out_dir, "04_dataset_gene_summaries", "brain_expression_R_dataset_summary_counts.csv"))
message("Effect-size rows written: raw=", nrow(effects_raw), "; collapsed=", nrow(effects_collapsed))
