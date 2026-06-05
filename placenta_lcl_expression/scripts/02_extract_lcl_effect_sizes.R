source(file.path("scripts", "lib", "placenta_lcl_expression_functions.R"))
cfg <- plcl_expression_config()

datasets <- c("GSE15402", "GSE15451", "GSE29918", "GSE37772", "GSE4187", "GSE7329")
results <- list()
failed <- list()
for (ds in datasets) {
  message("Processing ", ds)
  one <- try(process_lcl_series(cfg, ds), silent = TRUE)
  if (inherits(one, "try-error")) {
    failed[[length(failed) + 1]] <- data.frame(dataset = ds, branch = "LCL", status = "failed", reason = as.character(one))
  } else {
    results[[length(results) + 1]] <- one
  }
}

labels <- rbindlist_fill(lapply(results, `[[`, "labels"))
raw_labels <- rbindlist_fill(lapply(results, `[[`, "raw_labels"))
effects <- rbindlist_fill(lapply(results, `[[`, "effects"))
processed <- rbindlist_fill(lapply(results, `[[`, "processed"))

exclusions <- data.frame(
  dataset = c("GSE43076", "GSE29919", "GSE57802", "GSE285666", "GSE154829", "GSE178205", "GSE178206"),
  branch = c("LCL", "LCL/blood", "LCL", "LCL", "placenta", "placenta", "placenta"),
  status = c("excluded_not_ASD_control", "excluded_parent", "excluded_context", "excluded_context",
             "excluded_no_public_ASD_control_labels", "excluded_wrong_tissue_for_placenta_branch", "excluded_parent_multiomic"),
  reason = c(
    "Public raw route contains ASD and intellectual-disability arrays but no labelled normal/control arrays in the accessible sample set; not an ASD-control comparison.",
    "Superseries/parent route for GSE29918 and GSE29691; GSE29918 LCL child processed directly, parent not double-counted.",
    "LCL CNV/carrier route with controls but no clean ASD-control structure for broad ASD-control expression meta-analysis.",
    "Williams syndrome / parental control LCL route, not ASD-control.",
    "Placenta RNA-seq route lacks public ASD/control phenotype labels for this contrast in accessible GEO metadata.",
    "RNA-seq child matrix is labelled brain plus HEK293T transient overexpression samples, not placenta ASD-control expression.",
    "Superseries/multi-omic parent route; no public broad placenta ASD-control expression matrix was available for safe processing."
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(labels, file.path(cfg$out_dir, "02_sample_metadata", "placenta_lcl_expression_R_biological_sample_table.csv"))
write_csv_safe(raw_labels, file.path(cfg$out_dir, "02_sample_metadata", "placenta_lcl_expression_R_raw_sample_label_table.csv"))
write_csv_safe(processed, file.path(cfg$out_dir, "03_processed_expression", "lcl_expression_R_processed_matrix_summary.csv"))
write_csv_safe(effects, file.path(cfg$out_dir, "05_effect_sizes", "lcl_expression_R_dataset_level_effect_sizes.csv"))
write_csv_safe(data.frame(), file.path(cfg$out_dir, "05_effect_sizes", "placenta_expression_R_dataset_level_effect_sizes.csv"))
write_csv_safe(rbindlist_fill(c(failed, list(exclusions))), file.path(cfg$out_dir, "00_manifest", "placenta_lcl_expression_R_failed_or_excluded_datasets.csv"))

processing <- data.frame(
  dataset = processed$dataset,
  status = "complete",
  biological_ASD_n = vapply(processed$dataset, function(ds) length(unique(labels$sample[labels$dataset == ds & labels$group == "ASD"])), integer(1)),
  biological_control_n = vapply(processed$dataset, function(ds) length(unique(labels$sample[labels$dataset == ds & labels$group == "control"])), integer(1)),
  genes_with_effects = vapply(processed$dataset, function(ds) sum(effects$dataset == ds), integer(1)),
  platform = processed$platform,
  model_role = vapply(processed$dataset, role_for_dataset, character(1)),
  caveat = vapply(processed$dataset, caveat_for_dataset, character(1)),
  stringsAsFactors = FALSE
)
write_csv_safe(processing, file.path(cfg$out_dir, "04_dataset_gene_summaries", "lcl_expression_R_dataset_summary_counts.csv"))
message("LCL effect-size rows written: ", nrow(effects))
