# Build branch-level grouped-brain, platform, and region QC summaries.
#
# Reported pooled expression models are fitted in
# scripts/04_meta_analysis/ using metafor. The tables written here are retained
# to inspect route inclusion and sensitivity structure within the branch.

source(file.path("scripts", "lib", "brain_expression_functions.R"))
cfg <- brain_expression_config()

raw <- data.table::fread(file.path(cfg$out_dir, "05_effect_sizes", "brain_expression_R_dataset_level_effect_sizes_raw_strata.csv"), data.table = FALSE)
collapsed <- data.table::fread(file.path(cfg$out_dir, "05_effect_sizes", "brain_expression_R_dataset_level_effect_sizes_grouped_collapsed.csv"), data.table = FALSE)

models <- list()
models[["brain_expression_grouped_primary_public_R"]] <- run_meta_model(
  collapsed,
  "brain_expression_grouped_primary_public_R", "primary",
  "Grouped post-mortem brain public expression model; multi-region and assay heterogeneity interpreted cautiously"
)
models[["brain_expression_microarray_only_sensitivity_R"]] <- run_meta_model(
  collapsed[grepl("microarray", collapsed$assay, ignore.case = TRUE), ],
  "brain_expression_microarray_only_sensitivity_R", "sensitivity",
  "Microarray-only platform sensitivity"
)
models[["brain_expression_RNAseq_only_sensitivity_R"]] <- run_meta_model(
  collapsed[grepl("RNA-seq", collapsed$assay, ignore.case = TRUE), ],
  "brain_expression_RNAseq_only_sensitivity_R", "sensitivity",
  "RNA-seq-only platform sensitivity"
)
models[["brain_expression_cortex_only_sensitivity_R"]] <- run_meta_model(
  raw[grepl("cortex", raw$region_group, ignore.case = TRUE), ],
  "brain_expression_cortex_only_sensitivity_R", "sensitivity",
  "Cortex-only region sensitivity"
)
models[["brain_expression_cerebellum_only_sensitivity_R"]] <- run_meta_model(
  raw[grepl("cerebellum", raw$region_group, ignore.case = TRUE), ],
  "brain_expression_cerebellum_only_sensitivity_R", "sensitivity",
  "Cerebellum/Purkinje sensitivity"
)
models[["brain_expression_prefrontal_cortex_sensitivity_R"]] <- run_meta_model(
  raw[grepl("prefrontal", raw$brain_region, ignore.case = TRUE), ],
  "brain_expression_prefrontal_cortex_sensitivity_R", "sensitivity",
  "Prefrontal cortex sensitivity"
)
models[["brain_expression_BA19_occipital_sensitivity_R"]] <- run_meta_model(
  raw[grepl("BA19|occipital", raw$brain_region, ignore.case = TRUE), ],
  "brain_expression_BA19_occipital_sensitivity_R", "sensitivity",
  "BA19/occipital cortex sensitivity"
)
models[["brain_expression_non_cortical_sensitivity_R"]] <- run_meta_model(
  raw[grepl("non_cortical", raw$region_group, ignore.case = TRUE), ],
  "brain_expression_non_cortical_sensitivity_R", "sensitivity",
  "Non-cortical corpus-callosum sensitivity; expected to be descriptive only if one route contributes"
)

all_meta <- rbindlist_fill(lapply(models, `[[`, "meta"))
all_k1 <- rbindlist_fill(lapply(models, `[[`, "k1"))
all_summary <- rbindlist_fill(lapply(models, `[[`, "summary"))

write_csv_safe(models[["brain_expression_grouped_primary_public_R"]]$meta, file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_grouped_primary_meta_results.csv"))
write_csv_safe(models[["brain_expression_grouped_primary_public_R"]]$k1, file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_grouped_primary_k1_descriptive_rows.csv"))
write_csv_safe(all_summary[grepl("grouped|microarray|RNAseq", all_summary$model), ], file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_model_summary.csv"))
write_csv_safe(rbindlist_fill(lapply(models[c("brain_expression_cortex_only_sensitivity_R", "brain_expression_cerebellum_only_sensitivity_R", "brain_expression_prefrontal_cortex_sensitivity_R", "brain_expression_BA19_occipital_sensitivity_R", "brain_expression_non_cortical_sensitivity_R")], `[[`, "meta")), file.path(cfg$out_dir, "07_region_subtissue_sensitivity", "brain_expression_R_region_sensitivity_results_combined.csv"))
write_csv_safe(all_summary[grepl("cortex|cerebellum|BA19|non_cortical", all_summary$model), ], file.path(cfg$out_dir, "07_region_subtissue_sensitivity", "brain_expression_R_region_sensitivity_model_summary.csv"))
write_csv_safe(rbindlist_fill(lapply(models[c("brain_expression_microarray_only_sensitivity_R", "brain_expression_RNAseq_only_sensitivity_R")], `[[`, "meta")), file.path(cfg$out_dir, "08_platform_sensitivity", "brain_expression_R_platform_sensitivity_results_combined.csv"))
write_csv_safe(all_summary[grepl("microarray|RNAseq", all_summary$model), ], file.path(cfg$out_dir, "08_platform_sensitivity", "brain_expression_R_platform_sensitivity_model_summary.csv"))
write_csv_safe(all_meta, file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_all_model_results_combined.csv"))
write_csv_safe(all_k1, file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_all_k1_descriptive_rows_combined.csv"))

openxlsx::write.xlsx(
  list(
    model_summary = all_summary,
    grouped_primary = models[["brain_expression_grouped_primary_public_R"]]$meta,
    region_summary = all_summary[grepl("cortex|cerebellum|BA19|non_cortical", all_summary$model), ],
    platform_summary = all_summary[grepl("microarray|RNAseq", all_summary$model), ]
  ),
  file.path(cfg$out_dir, "brain_expression_R_results_summary.xlsx"),
  overwrite = TRUE
)
message("Model summaries written: ", nrow(all_summary), " models.")
