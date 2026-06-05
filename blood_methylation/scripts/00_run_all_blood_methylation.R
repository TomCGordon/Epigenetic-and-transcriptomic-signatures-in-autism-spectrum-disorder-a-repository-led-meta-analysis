#!/usr/bin/env Rscript

# Blood methylation public-data processing workflow.
#
# Run from this pipeline folder:
#   Rscript scripts/00_run_all_blood_methylation.R
#
# For the large GSE140730 WGBS step, set BLOOD_WGBS_DOWNLOAD_REPORTS=TRUE if
# the public GEO CpG reports are not already present locally. This may download
# tens of gigabytes, so the setting is explicit rather than hidden.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(Sys.getenv("ASD_REPO_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
Sys.setenv(ASD_REPO_ROOT = repo_root)

rscript <- Sys.getenv("RSCRIPT_EXE", unset = file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"))

run_step <- function(script_name) {
  path <- file.path(script_dir, script_name)
  message("\n=== Running ", script_name, " ===")
  status <- system2(rscript, shQuote(path))
  if (!identical(status, 0L)) stop("Step failed: ", script_name, call. = FALSE)
}

dir.create(file.path(package_root, "qc"), recursive = TRUE, showWarnings = FALSE)

run_step("01_download_blood_sources.R")
run_step("02_build_wgbs_promoter_coordinates.R")
run_step("03_process_GSE140730_wgbs_from_geo_cpg_reports.R")
run_step("04_run_blood_array_and_meta_analysis.R")
run_step("06_validate_meta_models_with_metafor.R")
run_step("05_write_blood_QC_and_traceability.R")

message("\nBlood methylation processing workflow completed.")
