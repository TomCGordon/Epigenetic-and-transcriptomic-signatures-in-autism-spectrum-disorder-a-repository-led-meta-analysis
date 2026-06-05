#!/usr/bin/env Rscript

# LCL expression public-data processing workflow.
# Placenta expression is screened here, but no public broad ASD-control placenta
# expression dataset met criteria for completed analysis.

source(file.path("scripts", "lib", "placenta_lcl_expression_functions.R"))
cfg <- plcl_expression_config()
for (d in cfg$subdirs) dir.create(file.path(cfg$out_dir, d), recursive = TRUE, showWarnings = FALSE)

message("01: staging public source files")
source(file.path(cfg$scripts_dir, "01_stage_public_sources.R"), local = TRUE)

message("02: extracting LCL expression and calculating dataset-level effects")
source(file.path(cfg$scripts_dir, "02_extract_lcl_effect_sizes.R"), local = TRUE)

message("03: writing LCL branch-level QC and sensitivity summaries")
source(file.path(cfg$scripts_dir, "03_run_lcl_meta_and_sensitivity.R"), local = TRUE)

message("04: writing QC and reproducibility reports")
source(file.path(cfg$scripts_dir, "04_write_qc_and_reproducibility.R"), local = TRUE)

message("Placenta/LCL expression processing workflow completed: ", cfg$out_dir)
