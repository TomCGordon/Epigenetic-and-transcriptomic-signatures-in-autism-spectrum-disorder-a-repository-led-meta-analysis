#!/usr/bin/env Rscript

# Convenience wrapper for the final model fitting, cross-omic convergence,
# pathway-enrichment, and figure-generation scripts.

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

run_script <- function(path) {
  message("\n=== ", path, " ===")
  status <- system2(rscript, path)
  if (!identical(status, 0L)) stop("Script failed: ", path, call. = FALSE)
}

run_script(file.path(root, "scripts", "04_meta_analysis", "00_run_all_meta_analysis_models.R"))
run_script(file.path(root, "scripts", "05_convergence_pathway", "run_cross_omic_pathway_enrichment.R"))
run_script(file.path(root, "scripts", "05_convergence_pathway", "postprocess_pathway_enrichment_outputs.R"))
run_script(file.path(root, "scripts", "05_convergence_pathway", "run_continuous_cross_omic_convergence_enrichment.R"))
run_script(file.path(root, "scripts", "05_convergence_pathway", "run_brain_subtissue_cross_omic_convergence.R"))
run_script(file.path(root, "scripts", "06_figures", "build_dataset_search_flow_diagram_prisma_style.R"))
run_script(file.path(root, "scripts", "06_figures", "build_enhanced_pathway_figures.R"))

message("Downstream analyses completed.")
