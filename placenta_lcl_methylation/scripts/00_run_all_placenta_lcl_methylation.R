#!/usr/bin/env Rscript

# Placenta and LCL methylation public-data processing workflow.
# Run from this pipeline folder.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

scripts <- c(
  "01_stage_placenta_lcl_sources.R",
  "02_build_promoter_maps_and_annotations.R",
  "03_process_placenta_wgbs_from_public_sources.R",
  "04_process_lcl_hm27_and_medip_from_public_sources.R",
  "05_run_placenta_lcl_meta_and_sensitivity_models.R",
  "06_validate_placenta_lcl_models_with_metafor.R",
  "07_write_placenta_lcl_QC_and_traceability.R"
)

for (script in scripts) {
  message("\n==== Running ", script, " ====")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(rscript, args = shQuote(file.path(script_dir, script)),
                    stdout = "", stderr = "")
  if (!identical(status, 0L)) stop("Script failed: ", script)
}

message("Placenta/LCL methylation processing workflow completed.")
