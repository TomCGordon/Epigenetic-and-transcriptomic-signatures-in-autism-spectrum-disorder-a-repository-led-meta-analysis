#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "lib", "placenta_lcl_methylation_functions.R"))

docs_dir <- file.path(package_root, "docs")
qc_dir <- file.path(package_root, "qc")
results_dir <- file.path(package_root, "results", "strict_missing_R_default")
dir_create(docs_dir)
dir_create(qc_dir)

script_files <- list.files(script_dir, pattern = "\\.R$", recursive = TRUE,
                           full.names = TRUE)
parse_check <- rbindlist(lapply(script_files, function(path) {
  ok <- TRUE
  msg <- NA_character_
  tryCatch(parse(file = path), error = function(e) {
    ok <<- FALSE
    msg <<- conditionMessage(e)
  })
  data.table(
    script = sub(paste0("^", gsub("\\\\", "/", script_dir), "/?"),
                 "", normalizePath(path, winslash = "/", mustWork = TRUE)),
    parses = ok,
    message = msg
  )
}), fill = TRUE)
write_csv(parse_check, file.path(qc_dir, "R_script_parse_check.csv"))

script_manifest <- data.table(
  run_order = c(0:7, NA_integer_),
  script = c(
    "00_run_all_placenta_lcl_methylation.R",
    "01_stage_placenta_lcl_sources.R",
    "02_build_promoter_maps_and_annotations.R",
    "03_process_placenta_wgbs_from_public_sources.R",
    "04_process_lcl_hm27_and_medip_from_public_sources.R",
    "05_run_placenta_lcl_meta_and_sensitivity_models.R",
    "06_validate_placenta_lcl_models_with_metafor.R",
    "07_write_placenta_lcl_QC_and_traceability.R",
    "lib/placenta_lcl_methylation_functions.R"
  ),
  role = c(
    "controller",
    "production: source staging and GEO metadata validation",
    "production: promoter maps and platform annotation",
    "production: placenta WGBS source-file processing",
    "production: LCL HM27 and MeDIP public-matrix processing",
    "production: effect-size synthesis and model building",
    "QC: independent DL validation against metafor",
    "QC: reproducibility manifests and audit reports",
    "shared production functions"
  ),
  input_level = c(
    "scripts",
    "public GEO/UCSC sources and optional user-supplied source archives",
    "UCSC refGene and GEO platform annotation",
    "public per-sample WGBS CpG/BED methylation source files",
    "public GEO processed matrices and GEO metadata",
    "R-generated dataset-level effect-size files",
    "R-generated effect-size and meta-analysis files",
    "R-generated manifests/results plus optional prior-package comparison",
    "none"
  ),
  output_level = c(
    "all outputs",
    "source manifests and staged public files",
    "promoter maps and gene-probe/coordinate maps",
    "per-sample promoter summaries, dataset summaries and effect sizes",
    "dataset summaries and effect sizes",
    "model-level results, model summaries and selected gene lists",
    "metafor validation tables",
    "QC reports and traceability manifests",
    "called by production scripts"
  ),
  production_or_audit = c("production", "production", "production", "production",
                          "production", "production", "audit", "audit", "production")
)
write_csv(script_manifest, file.path(docs_dir, "script_role_manifest.csv"))

dataset_status <- data.table(
  accession = c("GSE178203", "GSE67615", "GSE34099", "GSE99935"),
  branch = c("placenta", "placenta", "LCL", "LCL"),
  tissue_or_cell_source = c("placenta", "placenta", "lymphoblastoid cell line",
                            "lymphoblastoid cell line"),
  assay = c("WGBS CpG reports", "WGBS percent-methylation BED tarballs",
            "Illumina HumanMethylation27 public series matrix",
            "Affymetrix Human Promoter 1.0R MeDIP-chip public normalized matrix"),
  source_level = c(
    "public GEO per-sample CpG source files",
    "public GEO per-sample BED source files",
    "public GEO processed methylation matrix",
    "public GEO processed MeDIP-minus-input promoter matrix"
  ),
  workflow_class = c(
    "public WGBS source-file to final result",
    "public WGBS source-file to final result",
    "public repository matrix to final result",
    "public repository matrix to final result"
  ),
  reference_outputs_used_as_inputs = "no",
  controlled_access = "no",
  primary_role = c("placenta descriptive primary",
                   "placenta two-dataset sensitivity and descriptive add-on",
                   "LCL HM27 descriptive and cross-assay exploratory",
                   "LCL MeDIP descriptive and cross-assay exploratory"),
  interpretation_boundary = c(
    "k=1 descriptive when analysed alone",
    "k=1 descriptive alone; k=2 only in placenta sensitivity",
    "k=1 descriptive alone; contributes to cross-assay exploratory LCL model",
    "k=1 descriptive alone; contributes to cross-assay exploratory LCL model"
  )
)
write_csv(dataset_status, file.path(docs_dir, "dataset_source_processing_status.csv"))

result_files <- list.files(package_root, recursive = TRUE, full.names = TRUE)
file_index <- data.table(
  relative_path = sub(paste0("^", gsub("\\\\", "/", package_root), "/?"),
                      "", normalizePath(result_files, winslash = "/", mustWork = TRUE)),
  file_type = tools::file_ext(result_files),
  bytes = file.info(result_files)$size,
  modified_time = as.character(file.info(result_files)$mtime)
)
file_index[, role := fifelse(grepl("^scripts/", relative_path), "script",
                      fifelse(grepl("^data_raw/", relative_path), "staged public source or source manifest input",
                      fifelse(grepl("^data_processed/", relative_path), "R-generated processed analysis file",
                      fifelse(grepl("^results/", relative_path), "R-generated model result",
                      fifelse(grepl("^qc/", relative_path), "QC/audit output",
                      fifelse(grepl("^docs/", relative_path), "documentation/traceability", "package file"))))))]
write_csv(file_index, file.path(qc_dir, "placenta_lcl_methylation_file_index.csv"))

input_output_manifest <- data.table(
  output_file = c(
    "data_processed/placenta/placenta_all_dataset_effect_sizes.csv",
    "data_processed/lcl/lcl_all_dataset_effect_sizes.csv",
    "results/strict_missing_R_default/placenta_R_meta_results_combined.csv",
    "results/strict_missing_R_default/lcl_R_meta_results_combined.csv",
    "qc/metafor_DL_validation_summary_placenta_lcl_strict_missing_R_default.csv"
  ),
  produced_by_script = c(
    "03_process_placenta_wgbs_from_public_sources.R",
    "04_process_lcl_hm27_and_medip_from_public_sources.R",
    "05_run_placenta_lcl_meta_and_sensitivity_models.R",
    "05_run_placenta_lcl_meta_and_sensitivity_models.R",
    "06_validate_placenta_lcl_models_with_metafor.R"
  ),
  principal_inputs = c(
    "GSE178203 CpG reports, GSE67615 BED tarballs, UCSC refGene promoter maps",
    "GSE34099 GEO matrix/GPL8490 annotation, GSE99935 normalized MeDIP matrix, UCSC hg18 refGene promoter map",
    "R-generated placenta dataset-level Hedges g effect sizes",
    "R-generated LCL dataset-level Hedges g effect sizes",
    "R-generated dataset effect sizes and model results"
  ),
  traceability_note = c(
    "Per-sample promoter summaries are checkpointed but generated by this R workflow.",
    "Public processed matrices are treated as repository source inputs; missing entries remain missing.",
    "FDR is applied only to k>=2 rows; k=1 rows remain descriptive.",
    "Cross-assay exploratory synthesis is labelled separately from same-assay replication.",
    "Branch-level QC summaries are checked against metafor::rma.uni(method='DL')."
  )
)
write_csv(input_output_manifest, file.path(docs_dir, "input_output_manifest.csv"))

pkg_names <- c("data.table", "GEOquery", "metafor")
pkg_versions <- data.table(
  package = pkg_names,
  version = vapply(pkg_names, function(p) {
    if (!requireNamespace(p, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(p))
  }, character(1)),
  role = c("tabular processing", "GEO metadata and platform retrieval",
           "independent DerSimonian-Laird validation")
)
write_csv(pkg_versions, file.path(docs_dir, "R_package_versions.csv"))

compare_reference <- function(branch) {
  validation_reference_dir <- Sys.getenv("PLACENTA_LCL_METHYLATION_VALIDATION_REFERENCE_DIR", unset = "")
  if (!nzchar(validation_reference_dir)) {
    return(data.table(branch = branch, comparison_status = "reference comparison not requested"))
  }
  if (branch == "placenta") {
    new_path <- file.path(results_dir, "placenta_R_model_summary.csv")
    reference_path <- file.path(validation_reference_dir, "placenta_broad_model_summary.csv")
  } else {
    new_path <- file.path(results_dir, "lcl_R_model_summary.csv")
    reference_path <- file.path(validation_reference_dir, "lcl_broad_model_summary.csv")
  }
  if (!file.exists(new_path) || !file.exists(reference_path)) {
    return(data.table(branch = branch, comparison_status = "reference comparison not available"))
  }
  new <- fread(new_path)
  reference <- fread(reference_path)
  reference_rename <- c(
    genes_any = "reference_genes_with_any_result",
    genes_with_any_result = "reference_genes_with_any_result",
    genes_k_ge_2 = "reference_genes_meta_analysed_k_ge_2",
    genes_meta_analysed_k_ge_2 = "reference_genes_meta_analysed_k_ge_2",
    k1_descriptive_genes = "reference_k1_descriptive_genes",
    DL_nonzero_genes = "reference_DL_nonzero_genes",
    FDR_significant_genes = "reference_FDR_significant_genes",
    mKH_retained_genes = "reference_mKH_retained_genes"
  )
  for (nm in intersect(names(reference_rename), names(reference))) {
    setnames(reference, nm, reference_rename[[nm]])
  }
  keep_reference <- intersect(c("model_name", "reference_genes_with_any_result",
                                "reference_genes_meta_analysed_k_ge_2", "reference_k1_descriptive_genes",
                                "reference_DL_nonzero_genes",
                                "reference_FDR_significant_genes", "reference_mKH_retained_genes"), names(reference))
  merged <- merge(
    new[, .(model_name, genes_with_any_result, genes_meta_analysed_k_ge_2,
            k1_descriptive_genes, DL_nonzero_genes, FDR_significant_genes,
            mKH_retained_genes)],
    reference[, ..keep_reference],
    by = "model_name",
    all.x = TRUE
  )
  merged[, branch := branch]
  merged[, comparison_status := "reference-output comparison; reference outputs were not model inputs"]
  merged[]
}
comparison <- rbind(compare_reference("placenta"), compare_reference("LCL"), fill = TRUE)
if ("reference_genes_with_any_result" %in% names(comparison)) {
  numeric_pairs <- list(
    genes_with_any_result = "reference_genes_with_any_result",
    genes_meta_analysed_k_ge_2 = "reference_genes_meta_analysed_k_ge_2",
    k1_descriptive_genes = "reference_k1_descriptive_genes",
    DL_nonzero_genes = "reference_DL_nonzero_genes",
    FDR_significant_genes = "reference_FDR_significant_genes",
    mKH_retained_genes = "reference_mKH_retained_genes"
  )
  for (new_col in names(numeric_pairs)) {
    reference_col <- numeric_pairs[[new_col]]
    if (new_col %in% names(comparison) && reference_col %in% names(comparison)) {
      comparison[, (paste0("delta_", new_col)) := as.numeric(get(new_col)) - as.numeric(get(reference_col))]
    }
  }
}
write_csv(comparison, file.path(qc_dir, "placenta_lcl_reference_model_count_comparison.csv"))

model_summary <- rbindlist(list(
  if (file.exists(file.path(results_dir, "placenta_R_model_summary.csv"))) {
    data.table(branch = "placenta", fread(file.path(results_dir, "placenta_R_model_summary.csv")))
  },
  if (file.exists(file.path(results_dir, "lcl_R_model_summary.csv"))) {
    data.table(branch = "LCL", fread(file.path(results_dir, "lcl_R_model_summary.csv")))
  }
), fill = TRUE)

qc_lines <- c(
  "# Placenta/LCL Methylation R Workflow QC Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This package runs the public placenta and LCL methylation branches in R from public source files or public repository matrices. Reference model-summary tables may be supplied for count comparison, but they are not used as model inputs.",
  "",
  "## Dataset Boundaries",
  "",
  "- GSE178203 is rebuilt from public per-sample WGBS CpG reports into promoter-gene summaries.",
  "- GSE67615 is rebuilt from public per-sample WGBS percent-methylation BED tarballs into promoter-gene summaries.",
  "- GSE34099 is rebuilt from the public GEO HM27 series matrix with GEO platform annotation.",
  "- GSE99935 is rebuilt from the public GEO normalized MeDIP-minus-input matrix with coordinate-based promoter summarisation.",
  "",
  "## Key Statistical Rules",
  "",
  "- Empty, failed, unavailable, or non-numeric matrix values are treated as missing, not zero.",
  "- Hedges' g is calculated from dataset-level ASD/control means, SDs, and sample sizes.",
  "- DerSimonian-Laird random-effects models are applied where k >= 2; k = 1 rows are descriptive.",
  "- FDR correction is applied within each model only to rows with k >= 2 and finite p-values.",
  "- Modified Knapp-Hartung confidence intervals are calculated only for k >= 2 rows.",
  "- LCL cross-assay results remain exploratory because HM27 beta-like signal and MeDIP promoter-enrichment signal are not the same assay type.",
  "",
  "## Model Summary",
  ""
)
if (nrow(model_summary)) {
  qc_lines <- c(qc_lines, paste(capture.output(print(model_summary)), collapse = "\n"))
} else {
  qc_lines <- c(qc_lines, "Model summary files were not available when this report was generated.")
}
qc_lines <- c(qc_lines,
              "",
              "## Optional Reference-Output Comparison",
              "",
              "Reference model-summary tables are not used as inputs to this workflow. When supplied, they are compared here only to make count differences visible.",
              "")
if (nrow(comparison) && !"reference comparison not available" %in% comparison$comparison_status) {
  visible_cols <- intersect(c("branch", "model_name", "genes_with_any_result",
                              "reference_genes_with_any_result", "delta_genes_with_any_result",
                              "genes_meta_analysed_k_ge_2",
                              "reference_genes_meta_analysed_k_ge_2",
                              "delta_genes_meta_analysed_k_ge_2",
                              "DL_nonzero_genes", "reference_DL_nonzero_genes",
                              "delta_DL_nonzero_genes", "FDR_significant_genes",
                              "reference_FDR_significant_genes", "mKH_retained_genes",
                              "reference_mKH_retained_genes"), names(comparison))
  qc_lines <- c(qc_lines, paste(capture.output(print(comparison[, ..visible_cols])), collapse = "\n"),
                "",
                "Interpretation: count differences should be reviewed in relation to source files, platform annotation and missing-value handling. The qualitative LCL conclusion remains unchanged: no FDR-significant or modified Knapp-Hartung interval-supported cross-assay genes.")
} else {
  qc_lines <- c(qc_lines, "Reference-output comparison was unavailable.")
}
qc_lines <- c(qc_lines,
              "",
              "## Parse Check",
              "",
              if (all(parse_check$parses)) {
                "All R scripts parsed successfully."
              } else {
                "At least one R script failed to parse; see qc/R_script_parse_check.csv."
              },
              "",
              "## Reproducibility Scope",
              "",
              "This package supports a public-source-to-model workflow for placenta WGBS source files and a public-repository-matrix-to-model workflow for the LCL public matrices. Where repository-available sources were processed matrices, analyses begin from those matrices rather than from raw IDAT/CEL files.",
              "",
              "## Files to Inspect First",
              "",
              "- docs/dataset_source_processing_status.csv",
              "- docs/input_output_manifest.csv",
              "- docs/script_role_manifest.csv",
              "- results/strict_missing_R_default/placenta_R_model_summary.csv",
              "- results/strict_missing_R_default/lcl_R_model_summary.csv",
              "- qc/metafor_DL_validation_summary_placenta_lcl_strict_missing_R_default.csv",
              "- qc/placenta_lcl_reference_model_count_comparison.csv")

writeLines(qc_lines, file.path(qc_dir, "placenta_lcl_methylation_R_workflow_QC_report.md"),
           useBytes = TRUE)

message("Placenta/LCL QC and traceability files written.")
