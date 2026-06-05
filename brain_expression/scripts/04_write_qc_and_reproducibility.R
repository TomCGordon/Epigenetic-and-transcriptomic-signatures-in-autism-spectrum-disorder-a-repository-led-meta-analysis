source(file.path("scripts", "lib", "brain_expression_functions.R"))
cfg <- brain_expression_config()

summary_file <- file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_model_summary.csv")
new_summary <- if (file.exists(summary_file)) data.table::fread(summary_file, data.table = FALSE) else data.frame()
reference_summary_file <- if (nzchar(cfg$validation_reference_dir)) {
  file.path(cfg$validation_reference_dir, "brain_expression_model_summary.csv")
} else {
  ""
}
reference_summary <- if (nzchar(reference_summary_file) && file.exists(reference_summary_file)) data.table::fread(reference_summary_file, data.table = FALSE) else data.frame()

comparison <- data.frame()
if (nrow(new_summary) && nrow(reference_summary)) {
  primary_reference <- reference_summary[reference_summary$model == "brain_expression_grouped_primary_public", ]
  primary_new <- new_summary[new_summary$model == "brain_expression_grouped_primary_public_R", ]
  comparison <- data.frame(
    metric = c("genes_meta_analysed", "DL_nonzero", "FDR_significant", "mKH_interval_supported"),
    reference_output = c(primary_reference$genes_meta_analysed, primary_reference$DL_nonzero, primary_reference$FDR_significant, primary_reference$mKH_retained),
    R_workflow_output = c(primary_new$genes_meta_analysed, primary_new$DL_nonzero, primary_new$FDR_significant, primary_new$mKH_interval_supported)
  )
}
write_csv_safe(comparison, file.path(cfg$out_dir, "09_quality_control", "brain_expression_reference_primary_count_comparison.csv"))

validation <- data.frame()
primary <- file.path(cfg$out_dir, "06_grouped_brain_meta_analysis", "brain_expression_R_grouped_primary_meta_results.csv")
effects <- file.path(cfg$out_dir, "05_effect_sizes", "brain_expression_R_dataset_level_effect_sizes_grouped_collapsed.csv")
if (file.exists(primary) && file.exists(effects)) {
  meta <- data.table::fread(primary, data.table = FALSE)
  eff <- data.table::fread(effects, data.table = FALSE)
  set.seed(20260518)
  sample_genes <- sample(meta$gene, min(100, nrow(meta)))
  validation <- do.call(rbind, lapply(sample_genes, function(g) {
    sub <- eff[eff$gene == g, ]
    fit <- metafor::rma.uni(yi = sub$hedges_g, vi = sub$variance, method = "DL", test = "z")
    reported <- meta[meta$gene == g, ][1, ]
    data.frame(gene = g,
               metafor_g = as.numeric(fit$b), reported_g = reported$pooled_g,
               abs_diff_g = abs(as.numeric(fit$b) - reported$pooled_g),
               metafor_se = fit$se, reported_se = reported$se,
               abs_diff_se = abs(fit$se - reported$se),
               pass = abs(as.numeric(fit$b) - reported$pooled_g) < 1e-8 & abs(fit$se - reported$se) < 1e-8)
  }))
}
write_csv_safe(validation, file.path(cfg$out_dir, "09_quality_control", "brain_expression_R_metafor_validation_sample.csv"))

file_index <- data.frame(
  file_path = normalizePath(list.files(cfg$out_dir, recursive = TRUE, full.names = TRUE), winslash = "/", mustWork = FALSE),
  bytes = file.info(list.files(cfg$out_dir, recursive = TRUE, full.names = TRUE))$size,
  stringsAsFactors = FALSE
)
write_csv_safe(file_index, file.path(cfg$out_dir, "00_manifest", "brain_expression_R_file_index.csv"))

script_manifest <- data.frame(
  script = c("00_run_all_brain_expression.R", "01_stage_public_sources.R", "02_extract_process_and_effect_sizes.R", "03_run_meta_and_sensitivity_models.R", "04_write_qc_and_reproducibility.R", "lib/brain_expression_functions.R"),
  role = c("controller", "public source staging", "source extraction, processing, effect sizes", "meta-analysis and sensitivity models", "QC and reproducibility reporting", "shared functions"),
  production_script = TRUE,
  notes = c(
    "Runs the full brain expression R workflow in order.",
    "Copies from archived public source cache or downloads public GEO files if absent.",
    "Processes public source matrices/tables and calculates dataset-level Hedges g.",
    "Runs grouped post-mortem brain, platform, and region/subtissue models.",
    "Compares headline counts with a user-supplied reference table, when provided, and validates random-effects calculations against metafor.",
    "Contains parsers, label rules, transformations, Hedges g, DL random effects, FDR, and modified Knapp-Hartung interval logic."
  )
)
write_csv_safe(script_manifest, file.path(cfg$out_dir, "00_manifest", "brain_expression_R_script_role_manifest.csv"))

qc <- c(
  "# Brain Expression R Workflow QC Report",
  "",
  paste0("Created: ", Sys.time()),
  "",
  "## Scope",
  "",
  "This package runs the public post-mortem brain gene-expression extraction, processing, effect-size calculation, grouped-brain meta-analysis, and key region/platform sensitivity analyses in R.",
  "",
  "## Source handling",
  "",
  "Public GEO series matrices, supplemental count/expression tables, platform annotation files, and the GSE62098 public FPKM workbook are staged in `01_source_files/`. Files are copied from the archived public-source cache when present, or downloaded from public repositories if absent.",
  "",
  "## Important caveats",
  "",
  "Some supplemental matrices use public donor/sample IDs rather than GSM IDs. The R scripts encode the public sample-design mappings explicitly for those datasets, and these mappings are documented in the sample metadata output.",
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
  "Outputs were written only to the configured brain-expression results directory."
)
writeLines(qc, file.path(cfg$out_dir, "09_quality_control", "brain_expression_R_QC_report.md"))

readme <- c(
  "# Brain Expression R Workflow",
  "",
  "Run from the package root with:",
  "",
  "```r",
  "source('scripts/00_run_all_brain_expression.R')",
  "```",
  "",
  "The workflow stages public source files, processes microarray and RNA-seq matrices, computes dataset-level Hedges' g, runs DerSimonian-Laird random-effects meta-analyses, applies FDR and modified Knapp-Hartung interval support, and writes grouped-brain plus region/platform sensitivity outputs.",
  "",
  "Set `BRAIN_EXPRESSION_PUBLIC_SOURCE_CACHE` to a folder containing downloaded public source files if you do not want the script to download files from GEO."
)
writeLines(readme, file.path(cfg$out_dir, "README_reproducibility.md"))
message("QC and reproducibility files written.")

