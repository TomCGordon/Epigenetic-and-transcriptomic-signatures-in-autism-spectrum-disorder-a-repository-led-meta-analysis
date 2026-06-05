#!/usr/bin/env Rscript

# Brain gene-expression public-data processing workflow.
# Run from this pipeline folder.

source(file.path("scripts", "lib", "brain_expression_functions.R"))
cfg <- brain_expression_config()
for (d in cfg$subdirs) dir.create(file.path(cfg$out_dir, d), recursive = TRUE, showWarnings = FALSE)

message("01: staging public source files")
source(file.path(cfg$scripts_dir, "01_stage_public_sources.R"), local = TRUE)

message("02: extracting matrices and calculating dataset-level effect sizes")
source(file.path(cfg$scripts_dir, "02_extract_process_and_effect_sizes.R"), local = TRUE)

message("03: writing branch-level grouped-brain, platform, and region/subtissue QC summaries")
source(file.path(cfg$scripts_dir, "03_run_meta_and_sensitivity_models.R"), local = TRUE)

message("04: writing QC and reproducibility reports")
source(file.path(cfg$scripts_dir, "04_write_qc_and_reproducibility.R"), local = TRUE)

message("Brain expression processing workflow completed: ", cfg$out_dir)
