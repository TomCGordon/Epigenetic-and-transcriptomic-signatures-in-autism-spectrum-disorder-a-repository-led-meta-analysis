# Calculate dataset-level Hedges' g and branch-level QC model summaries.
#
# The reported pooled models are fitted in
# scripts/04_meta_analysis/ using metafor. The branch-level model tables written
# here are retained as local QC summaries for the blood-expression pipeline.

message("05: calculating effect sizes and blood-expression QC summaries")

sample_metadata <- read.csv(file.path(cfg$out_dir, "02_sample_metadata", "blood_expression_sample_phenotype_table.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)

micro <- read.csv(file.path(cfg$out_dir, "03_processed_expression", "blood_expression_microarray_gene_values_long.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
raw <- read.csv(file.path(cfg$out_dir, "03_processed_expression", "blood_expression_raw_route_gene_values_long.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)

all_values <- rbindlist_fill(list(micro, raw))
all_values <- all_values[is.finite(all_values$expression_value), ]
write_csv_safe(all_values, file.path(cfg$out_dir, "03_processed_expression", "blood_expression_all_gene_values_long.csv"))

raw_effects <- calculate_expression_effect_sizes(all_values)
write_csv_safe(raw_effects, file.path(cfg$out_dir, "05_effect_sizes", "blood_expression_dataset_level_effect_sizes_raw_strata.csv"))

collapsed_effects <- collapse_within_study_platforms(raw_effects)
write_csv_safe(collapsed_effects, file.path(cfg$out_dir, "05_effect_sizes", "blood_expression_dataset_level_effect_sizes.csv"))

primary_effects <- collapsed_effects[collapsed_effects$model_role == "peripheral_blood_primary", ]
cord_effects <- collapsed_effects[collapsed_effects$model_role == "cord_blood_developmental_sensitivity", ]
plus_cord_effects <- rbindlist_fill(list(primary_effects, cord_effects))

primary_meta <- build_meta_results(primary_effects, "blood_expression_peripheral_primary")
plus_cord_meta <- build_meta_results(plus_cord_effects, "blood_expression_plus_cord_blood_sensitivity")

combined <- rbindlist_fill(list(primary_meta$meta, primary_meta$descriptive, plus_cord_meta$meta, plus_cord_meta$descriptive))
write_csv_safe(combined, file.path(cfg$out_dir, "06_meta_analysis", "blood_expression_meta_results_combined.csv"))
write_csv_safe(primary_meta$meta, file.path(cfg$out_dir, "06_meta_analysis", "blood_expression_primary_meta_results.csv"))
write_csv_safe(plus_cord_meta$meta, file.path(cfg$out_dir, "07_sensitivity_analyses", "blood_expression_plus_cord_blood_sensitivity_meta_results.csv"))

model_summary <- rbindlist_fill(list(
  summarise_expression_model(primary_meta$meta, primary_meta$descriptive, primary_effects, "blood_expression_peripheral_primary",
                             "Primary peripheral/postnatal blood-family expression model"),
  summarise_expression_model(plus_cord_meta$meta, plus_cord_meta$descriptive, plus_cord_effects, "blood_expression_plus_cord_blood_sensitivity",
                             "Developmental sensitivity adding public cord-blood expression route where recoverable")
))
write_csv_safe(model_summary, file.path(cfg$out_dir, "06_meta_analysis", "blood_expression_model_summary.csv"))

top_fdr <- primary_meta$meta[primary_meta$meta$FDR_significant == TRUE, ]
top_fdr <- top_fdr[order(top_fdr$FDR, -abs(top_fdr$pooled_g)), ]
top_mkh <- primary_meta$meta[primary_meta$meta$mKH_interval_excludes_zero == TRUE, ]
top_mkh <- top_mkh[order(-abs(top_mkh$pooled_g), top_mkh$p_value), ]
write_csv_safe(top_fdr, file.path(cfg$out_dir, "06_meta_analysis", "blood_expression_primary_FDR_significant_genes.csv"))
write_csv_safe(head(top_mkh, 500), file.path(cfg$out_dir, "06_meta_analysis", "blood_expression_primary_top_mKH_interval_genes.csv"))
