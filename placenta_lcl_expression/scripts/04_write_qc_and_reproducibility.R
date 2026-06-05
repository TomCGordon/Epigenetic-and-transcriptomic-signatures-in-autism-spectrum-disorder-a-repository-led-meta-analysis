source(file.path("scripts", "lib", "placenta_lcl_expression_functions.R"))
cfg <- plcl_expression_config()

summary <- data.table::fread(file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_model_summary.csv"), data.table = FALSE)
reference_summary_file <- if (nzchar(cfg$validation_reference_dir)) {
  file.path(cfg$validation_reference_dir, "lcl_expression_model_summary.csv")
} else {
  ""
}
reference_summary <- if (nzchar(reference_summary_file) && file.exists(reference_summary_file)) data.table::fread(reference_summary_file, data.table = FALSE) else data.frame()
comparison <- data.frame()
if (nrow(reference_summary)) {
  comparison <- data.frame(
    metric = c("genes_meta_analysed", "DL_nonzero", "FDR_significant", "mKH_interval_supported"),
    reference_output = c(reference_summary$genes_meta_analysed[1], reference_summary$DL_nonzero[1], reference_summary$FDR_significant[1], reference_summary$mKH_retained[1]),
    R_workflow_output = c(summary$genes_meta_analysed[1], summary$DL_nonzero[1], summary$FDR_significant[1], summary$mKH_interval_supported[1]),
    stringsAsFactors = FALSE
  )
}
write_csv_safe(comparison, file.path(cfg$out_dir, "09_quality_control", "lcl_expression_reference_primary_count_comparison.csv"))

validation <- data.frame()
primary <- data.table::fread(file.path(cfg$out_dir, "06_lcl_meta_analysis", "lcl_expression_R_core_primary_meta_results.csv"), data.table = FALSE)
effects <- data.table::fread(file.path(cfg$out_dir, "05_effect_sizes", "lcl_expression_R_dataset_level_effect_sizes.csv"), data.table = FALSE)
set.seed(20260518)
if (nrow(primary)) {
  sample_genes <- sample(primary$gene, min(100, nrow(primary)))
  validation <- do.call(rbind, lapply(sample_genes, function(g) {
    sub <- effects[effects$gene == g & effects$dataset %in% c("GSE15402", "GSE15451", "GSE37772", "GSE4187"), ]
    fit <- metafor::rma.uni(yi = sub$hedges_g, vi = sub$variance, method = "DL", test = "z")
    reported <- primary[primary$gene == g, ][1, ]
    data.frame(gene = g, metafor_g = as.numeric(fit$b), reported_g = reported$pooled_g,
               abs_diff_g = abs(as.numeric(fit$b) - reported$pooled_g),
               metafor_se = fit$se, reported_se = reported$se,
               abs_diff_se = abs(fit$se - reported$se),
               pass = abs(as.numeric(fit$b) - reported$pooled_g) < 1e-8 & abs(fit$se - reported$se) < 1e-8)
  }))
}
write_csv_safe(validation, file.path(cfg$out_dir, "09_quality_control", "lcl_expression_R_metafor_validation_sample.csv"))

script_manifest <- data.frame(
  script = c("00_run_all_placenta_lcl_expression.R", "01_stage_public_sources.R", "02_extract_lcl_effect_sizes.R", "03_run_lcl_meta_and_sensitivity.R", "04_write_qc_and_reproducibility.R", "lib/placenta_lcl_expression_functions.R"),
  role = c("controller", "public source staging", "LCL extraction and effect sizes", "LCL meta-analysis and placenta no-run summary", "QC and reproducibility reporting", "shared functions"),
  production_script = TRUE,
  notes = c(
    "Runs the full package in order.",
    "Copies from archived public source cache or downloads GEO files if absent.",
    "Processes public LCL GEO series matrices and calculates dataset-level Hedges g.",
    "Runs core LCL primary and expanded sensitivity models; writes placenta no-eligible-public-dataset summary.",
    "Compares headline counts with a user-supplied reference table, when provided, and validates DL pooling against metafor.",
    "Contains source parsers, platform maps, label rules, replicate collapse, Hedges g, DL, FDR and mKH logic."
  )
)
write_csv_safe(script_manifest, file.path(cfg$out_dir, "00_manifest", "placenta_lcl_expression_R_script_role_manifest.csv"))

file_index <- data.frame(file_path = normalizePath(list.files(cfg$out_dir, recursive = TRUE, full.names = TRUE), winslash = "/", mustWork = FALSE),
                         bytes = file.info(list.files(cfg$out_dir, recursive = TRUE, full.names = TRUE))$size,
                         stringsAsFactors = FALSE)
write_csv_safe(file_index, file.path(cfg$out_dir, "00_manifest", "placenta_lcl_expression_R_file_index.csv"))

qc <- c(
  "# Placenta/LCL Expression R Workflow QC Report",
  "",
  paste0("Created: ", Sys.time()),
  "",
  "## Scope",
  "",
  "This package runs the public LCL expression extraction, effect-size calculation, primary meta-analysis, and expanded sensitivity models in R. Placenta expression routes are retained as formal exclusions because no public broad ASD-control placenta expression dataset with validated sample labels was available.",
  "",
  "## Validation",
  "",
  paste0("Metafor validation rows: ", nrow(validation), "; passing rows: ", if (nrow(validation)) sum(validation$pass) else 0, "."),
  "",
  "## Reference-output comparison",
  "",
  if (nrow(comparison)) paste(capture.output(print(comparison)), collapse = "\n") else "Reference-output comparison was not available.",
  "",
  "## Output scope",
  "",
  "Outputs were written only to the configured placenta/LCL expression results directory."
)
writeLines(qc, file.path(cfg$out_dir, "09_quality_control", "placenta_lcl_expression_R_QC_report.md"))

report <- c(
  "# Placenta/LCL Expression R Workflow Completion Report",
  "",
  "## LCL primary model",
  "",
  paste0("Core public LCL model genes with k >= 2: ", summary$genes_meta_analysed[1], "."),
  paste0("DL non-zero genes: ", summary$DL_nonzero[1], "."),
  paste0("FDR-significant genes: ", summary$FDR_significant[1], "."),
  paste0("Genes with modified Knapp-Hartung confidence intervals excluding zero: ", summary$mKH_interval_supported[1], "."),
  "",
  "## Placenta expression",
  "",
  "No public broad ASD-control placenta expression dataset with validated sample labels was available for completed open-data meta-analysis.",
  "",
  "## Output scope",
  "",
  "Outputs were written only to the configured placenta/LCL expression results directory."
)
writeLines(report, file.path(cfg$out_dir, "00_manifest", "placenta_lcl_expression_R_workflow_completion_report.md"))

readme <- c(
  "# Placenta/LCL Expression R Workflow",
  "",
  "Run from the package root with:",
  "",
  "```r",
  "source('scripts/00_run_all_placenta_lcl_expression.R')",
  "```",
  "",
  "Set `PLACENTA_LCL_EXPRESSION_PUBLIC_SOURCE_CACHE` to a folder containing downloaded public source files if you do not want the script to download GEO files."
)
writeLines(readme, file.path(cfg$out_dir, "README_reproducibility.md"))
message("QC and reproducibility files written.")

