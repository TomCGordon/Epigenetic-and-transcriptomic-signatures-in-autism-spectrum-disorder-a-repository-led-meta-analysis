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
qc_dir <- file.path(package_root, "qc")
docs_dir <- file.path(package_root, "docs")
results_dir <- file.path(package_root, "results", "strict_missing_R_default")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

scripts <- data.table(
  run_order = 0:7,
  script = c(
    "00_run_all_brain_methylation.R",
    "01_download_brain_sources.R",
    "02_build_brain_gene_universe_and_promoters.R",
    "03_process_brain_array_datasets.R",
    "04_process_brain_wgbs_from_public_cpg_files.R",
    "05_run_brain_meta_and_sensitivity_models.R",
    "06_validate_brain_meta_models_with_metafor.R",
    "07_write_brain_QC_and_traceability.R"
  ),
  role = c(
    "top-level controller",
    "repository retrieval and public source-file staging",
    "harmonised promoter universe and WGBS promoter coordinate generation",
    "array/HM27 matrix processing, gene summaries and effect sizes",
    "WGBS public CpG/BED file processing to promoter-gene summaries and effect sizes",
    "grouped brain primary model and region/platform/WGBS sensitivity models",
    "independent DerSimonian-Laird model validation against metafor",
    "QC, file index and traceability documentation"
  ),
  production_status = "production"
)
fwrite(scripts, file.path(docs_dir, "script_role_manifest.csv"))

dataset_status <- data.table(
  dataset_id = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608", "GSE109875", "GSE81541"),
  tissue_source = c("post-mortem brain mixed route", "post-mortem BA10/BA24 cortex", "post-mortem prefrontal cortex",
                    "post-mortem subventricular zone", "post-mortem dorsal raphe", "post-mortem cerebellum",
                    "post-mortem cerebellum/occipital cortex", "post-mortem BA9 cortex", "post-mortem brain"),
  assay = c(rep("Illumina HumanMethylation450", 6), "Illumina HumanMethylation27", "WGBS", "WGBS"),
  public_repository_source = c(rep("GEO matrix/processed matrix", 7), "GEO public CpG report files", "GEO public BED methylation files"),
  current_R_workflow_start_point = c(rep("repository processed methylation matrix or table", 7),
                                    "public CpG methylation report files", "public BED methylation files"),
  current_R_workflow_status = c(rep("repository-matrix-to-final-results in R", 7),
                               "CpG-report-to-promoter-to-final-results in R",
                               "BED-file-to-promoter-to-final-results in R"),
  notes = c(
    "All regions retained as the broad GSE53162 brain route; region-specific subroutes are generated for sensitivity models only.",
    "Repeated BA10/BA24 measurements are collapsed to independent subject-level values for the grouped parent route; region-specific subroutes are generated for sensitivity models only.",
    "Prefrontal cortex 450K route.",
    "Subventricular zone 450K route.",
    "Dorsal raphe 450K route.",
    "Large cerebellum-only 450K route.",
    "Technical/multi-region rows collapsed to independent subject-level values; region-specific subroutes are generated for sensitivity models only.",
    "Direct BA9 ASD/control WGBS subseries only.",
    "Idiopathic autism/control brain WGBS samples only; duplicate-region, syndromic and cell-culture samples excluded."
  )
)
fwrite(dataset_status, file.path(docs_dir, "dataset_source_processing_status.csv"))

input_manifest <- data.table(
  input = c(
    "GEO series or processed matrices for seven array/HM27 brain routes",
    "GEO sample metadata for WGBS routes",
    "GEO public per-sample CpG/BED WGBS files",
    "illumina450k_annotation_core.csv",
    "UCSC hg38 refGene table"
  ),
  produced_or_retrieved_by = c(
    "01_download_brain_sources.R",
    "01_download_brain_sources.R",
    "01_download_brain_sources.R",
    "01_download_brain_sources.R",
    "02_build_brain_gene_universe_and_promoters.R"
  ),
  used_by = c(
    "03_process_brain_array_datasets.R",
    "04_process_brain_wgbs_from_public_cpg_files.R",
    "04_process_brain_wgbs_from_public_cpg_files.R",
    "02_build_brain_gene_universe_and_promoters.R and 03_process_brain_array_datasets.R",
    "02_build_brain_gene_universe_and_promoters.R"
  ),
  notes = c(
    "Array/HM27 inputs are public repository processed matrices or tables rather than raw IDAT normalisation.",
    "Sample inclusion/exclusion is rebuilt from GEO metadata.",
    "Large source files are staged into data_raw/WGBS before processing.",
    "Defines the harmonised 20,960 promoter-gene universe and promoter-probe map.",
    "Defines hg38 strand-aware promoter windows for WGBS extraction."
  )
)
fwrite(input_manifest, file.path(docs_dir, "input_traceability_manifest.csv"))

packages_to_record <- c("data.table", "GEOquery", "metafor")
package_versions <- rbindlist(lapply(packages_to_record, function(pkg) {
  version <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  data.table(package = pkg, version = version)
}), fill = TRUE)
package_versions <- rbind(data.table(package = "R", version = paste(R.version$major, R.version$minor, sep = ".")), package_versions)
fwrite(package_versions, file.path(docs_dir, "R_package_versions.csv"))

script_files <- list.files(script_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
parse_check <- rbindlist(lapply(script_files, function(path) {
  ok <- tryCatch({ parse(file = path); TRUE }, error = function(e) FALSE)
  data.table(script = sub(paste0("^", gsub("\\\\", "/", script_dir), "/?"), "", gsub("\\\\", "/", path)), parses = ok)
}), fill = TRUE)
fwrite(parse_check, file.path(qc_dir, "R_script_parse_check.csv"))

model_text <- if (file.exists(file.path(results_dir, "brain_R_model_summary.csv"))) {
  ms <- fread(file.path(results_dir, "brain_R_model_summary.csv"))
  paste(capture.output(print(ms[, .(model_name, genes_with_any_result, genes_meta_analysed_k_ge_2,
                                    DL_nonzero_genes, FDR_significant_genes, mKH_retained_genes)])), collapse = "\n")
} else {
  "Model summary not found. Run 05_run_brain_meta_and_sensitivity_models.R."
}

validation_text <- if (file.exists(file.path(qc_dir, "metafor_DL_validation_summary_brain_strict_missing_R_default.csv"))) {
  paste(capture.output(print(fread(file.path(qc_dir, "metafor_DL_validation_summary_brain_strict_missing_R_default.csv")))), collapse = "\n")
} else {
  "metafor validation summary not found."
}

qc_md <- c(
  "# Brain Methylation R Workflow QC",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This package runs the brain methylation branch from public repository inputs using R scripts. Array/HM27 routes begin from public repository processed matrices or tables. WGBS routes begin from public per-sample CpG/BED methylation files.",
  "",
  "## Important Boundaries",
  "",
    "- Final brain methylation result tables are not used as analysis inputs.",
  "- Array/HM27 analyses are not raw-IDAT-to-final re-normalisations; they are repository-matrix-to-final analyses.",
  "- WGBS analyses are rebuilt from public CpG/BED methylation files into promoter-gene summaries.",
  "- Region/subtissue sensitivity models use R-generated region-specific subroutes from the same source data, not anatomically mixed parent accessions.",
  "- Repeated-region measurements are collapsed to subject-level values for GSE53924 and GSE38608 before dataset-level effect-size estimation.",
  "- Blank matrix cells are treated as missing; numeric zeros are retained.",
  "",
  "## Model Summary",
  "",
  "```text",
  model_text,
  "```",
  "",
  "## Metafor Validation",
  "",
  "```text",
  validation_text,
  "```",
  "",
  "## Script Parse Check",
  "",
  "```text",
  paste(capture.output(print(parse_check)), collapse = "\n"),
  "```",
  "",
  "## R Package Versions",
  "",
  "```text",
  paste(capture.output(print(package_versions)), collapse = "\n"),
  "```"
)
writeLines(qc_md, file.path(qc_dir, "brain_methylation_R_workflow_QC_report.md"))

file_index <- data.table(path = list.files(package_root, recursive = TRUE, full.names = TRUE))
file_index[, `:=`(
  relative_path = sub(paste0("^", gsub("\\\\", "/", package_root), "/?"), "", gsub("\\\\", "/", path)),
  bytes = file.info(path)$size
)]
fwrite(file_index[, .(relative_path, bytes)], file.path(qc_dir, "brain_methylation_file_index.csv"))
message("Brain methylation QC report written.")
