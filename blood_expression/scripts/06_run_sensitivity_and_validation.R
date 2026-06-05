# Branch-level sensitivity/status summaries.
#
# Final reported models are fitted in scripts/04_meta_analysis/ using
# metafor. This step checks that branch-level QC summaries remain consistent
# with metafor estimates and writes a status-transition table for review.

message("06: running branch QC checks and sensitivity summaries")

effects <- read.csv(file.path(cfg$out_dir, "05_effect_sizes", "blood_expression_dataset_level_effect_sizes.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
primary_meta <- read.csv(file.path(cfg$out_dir, "06_meta_analysis", "blood_expression_primary_meta_results.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
plus_cord_meta <- read.csv(file.path(cfg$out_dir, "07_sensitivity_analyses", "blood_expression_plus_cord_blood_sensitivity_meta_results.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)

validation <- validate_meta_with_metafor(effects, primary_meta, "blood_expression_peripheral_primary", n_genes = 100)
write_csv_safe(validation, file.path(cfg$out_dir, "08_quality_control", "blood_expression_metafor_validation_random_100.csv"))

transition <- compare_model_status(primary_meta, plus_cord_meta,
                                   "blood_expression_peripheral_primary",
                                   "blood_expression_plus_cord_blood_sensitivity")
write_csv_safe(transition, file.path(cfg$out_dir, "07_sensitivity_analyses", "blood_expression_primary_vs_plus_cord_status_transition.csv"))
