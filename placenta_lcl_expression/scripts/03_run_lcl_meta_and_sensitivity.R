# Build branch-level LCL QC summaries.
#
# The reported LCL expression model is fitted in
# scripts/04_meta_analysis/ using metafor. The outputs from this script are
# retained for branch-level inspection and sensitivity checks.

source(file.path("scripts", "lib", "placenta_lcl_expression_functions.R"))
cfg <- plcl_expression_config()
effects <- data.table::fread(file.path(cfg$out_dir, "05_effect_sizes", "lcl_expression_R_dataset_level_effect_sizes.csv"), data.table = FALSE)

core <- effects[effects$dataset %in% c("GSE15402", "GSE15451", "GSE37772", "GSE4187"), ]
expanded <- effects[effects$dataset %in% c("GSE15402", "GSE15451", "GSE29918", "GSE37772", "GSE4187", "GSE7329"), ]
no_pooled <- effects[effects$dataset %in% c("GSE15402", "GSE15451", "GSE29918", "GSE37772", "GSE4187"), ]

models <- list(
  core = run_meta_model(core, "lcl_expression_core_public_primary_R", "primary", "Individual-level public LCL ASD/control model. Family/twin non-independence remains a caveat."),
  expanded = run_meta_model(expanded, "lcl_expression_expanded_public_sensitivity_R", "sensitivity", "Expanded LCL sensitivity adding ASD-features and pooled/syndromic routes."),
  no_pooled = run_meta_model(no_pooled, "lcl_expression_no_pooled_syndromic_sensitivity_R", "sensitivity", "Excludes pooled syndromic GSE7329 but retains ASD-features route.")
)

summaries <- rbindlist_fill(lapply(models, `[[`, "summary"))
placenta_summary <- data.frame(model = "placenta_expression_public_R", model_role = "not_run_no_eligible_public_dataset",
                               datasets = "", genes_meta_analysed = 0, k1_descriptive = 0,
                               DL_nonzero = 0, FDR_significant = 0, mKH_interval_supported = 0,
                               FDR_mKH_overlap = 0, median_k = NA_real_, median_I2 = NA_real_,
                               caveat = "No public broad placenta ASD-control expression dataset with validated sample labels was available.",
                               stringsAsFactors = FALSE)
write_csv_safe(models$core$meta, file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_core_primary_meta_results.csv"))
write_csv_safe(models$core$k1, file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_core_primary_k1_descriptive_rows.csv"))
write_csv_safe(models$expanded$meta, file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_expanded_sensitivity_meta_results.csv"))
write_csv_safe(models$no_pooled$meta, file.path(cfg$out_dir, "08_sensitivity_analyses", "lcl_expression_R_no_pooled_syndromic_sensitivity_meta_results.csv"))
write_csv_safe(summaries, file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_model_summary.csv"))
write_csv_safe(placenta_summary, file.path(cfg$out_dir, "07_placenta_meta_analysis", "placenta_expression_R_model_summary.csv"))
write_csv_safe(rbindlist_fill(list(summaries, placenta_summary)), file.path(cfg$out_dir, "00_manifest", "placenta_lcl_expression_R_model_summary.csv"))
write_csv_safe(models$core$meta, file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_ranked_gene_results.csv"))
write_csv_safe(models$core$meta[models$core$meta$FDR_significant, ], file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_FDR_significant_genes.csv"))
write_csv_safe(models$core$meta[models$core$meta$mKH_interval_excludes_zero, ], file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_mKH_interval_supported_genes.csv"))
write_csv_safe(data.frame(), file.path(cfg$out_dir, "07_placenta_meta_analysis", "placenta_expression_R_meta_results_combined.csv"))

openxlsx::write.xlsx(
  list(
    lcl_model_summary = summaries,
    placenta_model_summary = placenta_summary,
    lcl_core_primary = models$core$meta,
    lcl_expanded_sensitivity = models$expanded$meta,
    lcl_no_pooled_syndromic = models$no_pooled$meta
  ),
  file.path(cfg$out_dir, "placenta_lcl_expression_R_results_summary.xlsx"),
  overwrite = TRUE
)
message("LCL/placenta model summaries written.")
