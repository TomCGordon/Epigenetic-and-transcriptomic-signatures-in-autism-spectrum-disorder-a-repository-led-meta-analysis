#!/usr/bin/env Rscript

# Run methylation and expression random-effects meta-analysis models.
#
# This wrapper calls the two public-facing model-fitting scripts in this folder.
# Run the methylation and expression processing pipelines first so that the
# required dataset-level effect-size tables exist.

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_file)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

run_script <- function(file_name) {
  message("\n=== ", file_name, " ===")
  status <- system2(rscript, file.path(script_dir, file_name))
  if (!identical(status, 0L)) stop("Meta-analysis step failed: ", file_name, call. = FALSE)
}

run_script("01_fit_methylation_models_metafor.R")
run_script("02_fit_expression_models_metafor.R")

message("All meta-analysis models completed.")

