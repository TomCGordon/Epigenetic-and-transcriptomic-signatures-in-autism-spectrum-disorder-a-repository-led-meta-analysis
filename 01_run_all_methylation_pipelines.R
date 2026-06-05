#!/usr/bin/env Rscript

# Convenience wrapper for the public methylation source-to-effect-size pipelines.
# Individual branch scripts can also be run directly; see docs/RUN_ORDER.md.

run_from <- function(pipeline_dir, script) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(pipeline_dir)
  message("\n=== ", basename(pipeline_dir), " ===")
  status <- system2(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
                    file.path("scripts", script))
  if (!identical(status, 0L)) stop("Pipeline failed: ", pipeline_dir, call. = FALSE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

run_from(file.path(root, "pipelines", "blood_methylation"),
         "00_run_all_blood_methylation.R")
run_from(file.path(root, "pipelines", "brain_methylation"),
         "00_run_all_brain_methylation.R")
run_from(file.path(root, "pipelines", "placenta_lcl_methylation"),
         "00_run_all_placenta_lcl_methylation.R")

message("All methylation pipelines completed.")

