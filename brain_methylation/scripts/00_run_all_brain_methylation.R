#!/usr/bin/env Rscript

# Brain methylation public-data processing workflow.
#
# Run from this pipeline folder:
#   Rscript scripts/00_run_all_brain_methylation.R
#
# The WGBS source files are large. By default, the download/staging script first
# uses public-source files already present in the repository working cache. If
# those files are unavailable, it downloads them from GEO sample supplementary
# file URLs.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repo_root <- normalizePath(Sys.getenv("ASD_REPO_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
Sys.setenv(ASD_REPO_ROOT = repo_root)

rscript <- Sys.getenv("RSCRIPT_EXE", unset = file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"))

run_step <- function(script_name) {
  path <- file.path(script_dir, script_name)
  message("\n=== Running ", script_name, " ===")
  status <- system2(rscript, shQuote(path))
  if (!identical(status, 0L)) stop("Step failed: ", script_name, call. = FALSE)
}

run_step("01_download_brain_sources.R")
run_step("02_build_brain_gene_universe_and_promoters.R")
run_step("03_process_brain_array_datasets.R")
run_step("04_process_brain_wgbs_from_public_cpg_files.R")
run_step("05_run_brain_meta_and_sensitivity_models.R")
run_step("06_validate_brain_meta_models_with_metafor.R")
run_step("07_write_brain_QC_and_traceability.R")

message("\nBrain methylation processing workflow completed.")
