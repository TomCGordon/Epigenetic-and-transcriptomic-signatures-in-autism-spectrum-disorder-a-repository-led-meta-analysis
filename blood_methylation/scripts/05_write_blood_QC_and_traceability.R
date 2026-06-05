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
raw_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_RAW_DIR", unset = file.path(package_root, "data_raw")),
                         winslash = "/", mustWork = FALSE)
processed_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_PROCESSED_DIR", unset = file.path(package_root, "data_processed")),
                               winslash = "/", mustWork = FALSE)
results_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_R_OUTPUT_DIR", unset = file.path(package_root, "results")),
                             winslash = "/", mustWork = FALSE)
qc_dir <- normalizePath(file.path(package_root, "qc"), winslash = "/", mustWork = FALSE)
docs_dir <- normalizePath(file.path(package_root, "docs"), winslash = "/", mustWork = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

scripts <- data.table(
  run_order = 0:6,
  script = c(
    "00_run_all_blood_methylation.R",
    "01_download_blood_sources.R",
    "02_build_wgbs_promoter_coordinates.R",
    "03_process_GSE140730_wgbs_from_geo_cpg_reports.R",
    "04_run_blood_array_and_meta_analysis.R",
    "06_validate_meta_models_with_metafor.R",
    "05_write_blood_QC_and_traceability.R"
  ),
  role = c(
    "top-level controller",
    "repository retrieval and source-file staging",
    "harmonised promoter universe and WGBS coordinate generation",
    "GSE140730 WGBS CpG-report to promoter-gene extraction",
    "array processing, effect sizes, random-effects meta-analysis and sensitivity models",
    "independent DerSimonian-Laird model validation against metafor",
    "QC, file index and traceability documentation"
  ),
  production_status = "production"
)
fwrite(scripts, file.path(docs_dir, "script_role_manifest.csv"))

input_manifest <- data.table(
  input = c(
    "GEO series matrices: GSE109905, GSE113967, GSE108785, GSE27044",
    "GEO MINiML family archive: GSE83424",
    "GEO MINiML family archive plus Bismark CpG reports: GSE140730",
    "illumina450k_annotation_core.csv",
    "UCSC hg38 refGene table"
  ),
  produced_or_retrieved_by = c(
    "01_download_blood_sources.R",
    "01_download_blood_sources.R",
    "01_download_blood_sources.R and 03_process_GSE140730_wgbs_from_geo_cpg_reports.R",
    "01_download_blood_sources.R copies the repository snapshot into data_raw/annotation",
    "02_build_wgbs_promoter_coordinates.R downloads from UCSC"
  ),
  used_by = c(
    "04_run_blood_array_and_meta_analysis.R",
    "04_run_blood_array_and_meta_analysis.R",
    "03_process_GSE140730_wgbs_from_geo_cpg_reports.R and 04_run_blood_array_and_meta_analysis.R",
    "02_build_wgbs_promoter_coordinates.R and 04_run_blood_array_and_meta_analysis.R",
    "02_build_wgbs_promoter_coordinates.R"
  ),
  notes = c(
    "Repository processed beta matrices; raw IDAT preprocessing is not possible for all public blood routes because the analysis harmonises repository-derived matrices across platforms.",
    "Sample-level table files are extracted from the MINiML archive.",
    "CpG reports are large. Set BLOOD_WGBS_DOWNLOAD_REPORTS=TRUE to download missing reports; otherwise local reports or direct GEO URLs are used.",
    "Defines the harmonised Illumina promoter-probe universe.",
    "Defines hg38 strand-aware promoter windows for WGBS sensitivity extraction."
  )
)
fwrite(input_manifest, file.path(docs_dir, "input_traceability_manifest.csv"))

packages_to_record <- c("data.table", "GEOquery", "xml2", "metafor")
package_versions <- rbindlist(lapply(packages_to_record, function(pkg) {
  version <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  data.table(package = pkg, version = version)
}), fill = TRUE)
package_versions <- rbind(
  data.table(package = "R", version = paste(R.version$major, R.version$minor, sep = ".")),
  package_versions
)
fwrite(package_versions, file.path(docs_dir, "R_package_versions.csv"))

model_summary_path <- file.path(results_dir, "array_blanks_as_missing_R_default", "blood_R_model_summary.csv")
download_manifest <- file.path(qc_dir, "blood_source_download_manifest.csv")
wgbs_qc <- file.path(qc_dir, "GSE140730_wgbs_processing_QC.csv")
script_files <- list.files(script_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
parse_check <- rbindlist(lapply(script_files, function(path) {
  ok <- tryCatch({ parse(file = path); TRUE }, error = function(e) FALSE)
  data.table(
    script = sub(paste0("^", gsub("\\\\", "/", script_dir), "/?"), "", gsub("\\\\", "/", path)),
    parses = ok
  )
}), fill = TRUE)
fwrite(parse_check, file.path(qc_dir, "R_script_parse_check.csv"))

model_text <- if (file.exists(model_summary_path)) {
  ms <- fread(model_summary_path)
  paste(capture.output(print(ms[, .(model_name, genes_tested, DL_nonzero_genes,
                                    FDR_significant_genes, mKH_retained_genes)])), collapse = "\n")
} else {
  "Model summary was not found. Run 04_run_blood_array_and_meta_analysis.R."
}

download_text <- if (file.exists(download_manifest)) {
  paste(capture.output(print(fread(download_manifest)[, .(status, file, bytes)])), collapse = "\n")
} else {
  "Download manifest was not found. Run 01_download_blood_sources.R."
}

wgbs_text <- if (file.exists(wgbs_qc)) {
  paste(capture.output(print(fread(wgbs_qc))), collapse = "\n")
} else {
  "GSE140730 WGBS QC was not found. Run 03_process_GSE140730_wgbs_from_geo_cpg_reports.R."
}

qc_md <- c(
  "# Blood Methylation R Workflow QC",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This package stages public repository files, builds the harmonised promoter universe, processes blood array/HM27 datasets from GEO source matrices or MINiML sample tables, processes GSE140730 from public Bismark CpG reports, recalculates Hedges' g and random-effects meta-analysis models, and writes sensitivity models.",
  "",
  "## Important Boundaries",
  "",
  "- The array routes use repository-supplied processed beta-value matrices/sample tables rather than re-normalising IDATs, because the completed analysis harmonises public processed matrices across studies and platforms.",
  "- The WGBS route uses public Bismark CpG reports and R-generated promoter summaries.",
  "- Blank array cells are treated as missing by default; a blank-as-zero option exists only as an explicit missing-data sensitivity setting.",
  "- Final methylation result tables are not used as analysis inputs.",
  "",
  "## Source Retrieval",
  "",
  "```text",
  download_text,
  "```",
  "",
  "## GSE140730 WGBS Processing",
  "",
  "```text",
  wgbs_text,
  "```",
  "",
  "## Model Summary",
  "",
  "```text",
  model_text,
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
  "```",
  "",
  "## File Locations",
  "",
  paste0("- Raw/source data staged under: `", raw_dir, "`"),
  paste0("- Processed data written under: `", processed_dir, "`"),
  paste0("- Results written under: `", results_dir, "`"),
  paste0("- QC written under: `", qc_dir, "`")
)
writeLines(qc_md, file.path(qc_dir, "blood_methylation_R_workflow_QC_report.md"))

file_index <- data.table(
  path = list.files(package_root, recursive = TRUE, full.names = TRUE)
)[, `:=`(
  relative_path = sub(paste0("^", gsub("\\\\", "/", package_root), "/?"), "", gsub("\\\\", "/", path)),
  bytes = file.info(path)$size
)]
fwrite(file_index[, .(relative_path, bytes)], file.path(qc_dir, "blood_methylation_file_index.csv"))

message("Blood methylation QC report written.")
