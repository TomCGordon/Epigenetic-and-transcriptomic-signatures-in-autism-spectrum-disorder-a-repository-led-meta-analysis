#!/usr/bin/env Rscript

# Blood gene-expression public-data processing workflow.
# Run from this pipeline folder.

source(file.path("scripts", "lib", "blood_expression_functions.R"))

cfg <- blood_expression_config()

message("Blood expression output folder: ", cfg$out_dir)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
for (d in cfg$subdirs) dir.create(file.path(cfg$out_dir, d), recursive = TRUE, showWarnings = FALSE)

source(file.path(cfg$scripts_dir, "01_stage_public_sources.R"), local = TRUE)
source(file.path(cfg$scripts_dir, "02_build_sample_metadata.R"), local = TRUE)
source(file.path(cfg$scripts_dir, "03_process_microarray_series_matrices.R"), local = TRUE)
source(file.path(cfg$scripts_dir, "04_process_raw_count_and_signal_routes.R"), local = TRUE)
source(file.path(cfg$scripts_dir, "05_calculate_effect_sizes_and_models.R"), local = TRUE)
source(file.path(cfg$scripts_dir, "06_run_sensitivity_and_validation.R"), local = TRUE)
source(file.path(cfg$scripts_dir, "07_write_qc_and_reproducibility.R"), local = TRUE)

message("Blood expression processing workflow completed.")
