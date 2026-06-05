#!/usr/bin/env Rscript

# Placenta/LCL methylation branch-level QC summaries.
#
# This script combines dataset-level Hedges' g inputs and writes branch-level
# descriptive/sensitivity summaries. The reported pooled models are
# fitted separately in scripts/04_meta_analysis/ using metafor.

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

placenta_processed <- file.path(package_root, "data_processed", "placenta")
lcl_processed <- file.path(package_root, "data_processed", "lcl")
results_dir <- file.path(package_root, "results", "strict_missing_R_default")
qc_dir <- file.path(package_root, "qc")
dir_create(results_dir)
dir_create(qc_dir)

placenta_effects <- fread(file.path(placenta_processed, "placenta_all_dataset_effect_sizes.csv"))
lcl_effects <- fread(file.path(lcl_processed, "lcl_all_dataset_effect_sizes.csv"))
write_csv(placenta_effects, file.path(results_dir, "placenta_R_all_dataset_level_effect_sizes.csv"))
write_csv(lcl_effects, file.path(results_dir, "lcl_R_all_dataset_level_effect_sizes.csv"))

placenta_universe <- sort(unique(placenta_effects$gene))
lcl_universe <- sort(unique(lcl_effects$gene))
write_csv(data.table(gene = placenta_universe), file.path(results_dir, "placenta_R_tested_gene_universe.csv"))
write_csv(data.table(gene = lcl_universe), file.path(results_dir, "lcl_R_tested_gene_universe.csv"))

placenta_model_defs <- list(
  placenta_primary_GSE178203_descriptive = "GSE178203",
  placenta_two_dataset_sensitivity = c("GSE178203", "GSE67615"),
  placenta_GSE67615_descriptive = "GSE67615"
)
placenta_model_roles <- data.table(
  model_name = names(placenta_model_defs),
  role = c("primary descriptive", "sensitivity", "descriptive add-on"),
  datasets = vapply(placenta_model_defs, paste, collapse = "; ", FUN.VALUE = character(1)),
  biological_interpretation = c(
    "Primary broad-processable placenta WGBS dataset; k=1 descriptive only.",
    "Two-dataset placenta WGBS sensitivity model; interpret with comparability/independence caveat.",
    "Secondary placenta WGBS dataset retained as descriptive add-on."
  ),
  caveat = c("k=1; not replicated", "k=2 where both datasets cover the gene; independence/comparability caveat",
             "k=1; not primary")
)

lcl_model_defs <- list(
  lcl_cross_assay_exploratory = c("GSE34099", "GSE99935"),
  lcl_HM27_descriptive = "GSE34099",
  lcl_MeDIP_descriptive = "GSE99935"
)
lcl_model_roles <- data.table(
  model_name = names(lcl_model_defs),
  role = c("exploratory primary", "assay-specific descriptive", "assay-specific descriptive"),
  datasets = vapply(lcl_model_defs, paste, collapse = "; ", FUN.VALUE = character(1)),
  biological_interpretation = c(
    "Cross-assay LCL synthesis combining HM27 beta-like methylation and MeDIP promoter-enrichment scores.",
    "HM27 LCL descriptive layer only.",
    "MeDIP-chip LCL descriptive layer only."
  ),
  caveat = c("HM27 beta-like plus MeDIP enrichment; cross-assay caveat", "k=1", "k=1")
)

placenta_meta <- build_meta_results(placenta_effects, placenta_universe, placenta_model_defs, placenta_model_roles)
lcl_meta <- build_meta_results(lcl_effects, lcl_universe, lcl_model_defs, lcl_model_roles)
placenta_summary <- summarise_models(placenta_meta, placenta_model_roles)
lcl_summary <- summarise_models(lcl_meta, lcl_model_roles)

write_csv(placenta_model_roles, file.path(results_dir, "placenta_R_model_definitions.csv"))
write_csv(lcl_model_roles, file.path(results_dir, "lcl_R_model_definitions.csv"))
write_csv(placenta_meta, file.path(results_dir, "placenta_R_meta_results_combined.csv"))
write_csv(lcl_meta, file.path(results_dir, "lcl_R_meta_results_combined.csv"))
write_csv(placenta_summary, file.path(results_dir, "placenta_R_model_summary.csv"))
write_csv(lcl_summary, file.path(results_dir, "lcl_R_model_summary.csv"))

write_csv(placenta_meta[model_name == "placenta_primary_GSE178203_descriptive" & DL_CI_excludes_zero == TRUE],
          file.path(results_dir, "placenta_R_primary_DL_nonzero_genes.csv"))
write_csv(placenta_meta[model_name == "placenta_two_dataset_sensitivity" & DL_CI_excludes_zero == TRUE],
          file.path(results_dir, "placenta_R_two_dataset_DL_nonzero_genes.csv"))
write_csv(placenta_meta[model_name == "placenta_two_dataset_sensitivity" & mKH_CI_excludes_zero == TRUE],
          file.path(results_dir, "placenta_R_two_dataset_mKH_interval_supported_genes.csv"))
write_csv(lcl_meta[model_name == "lcl_cross_assay_exploratory" & DL_CI_excludes_zero == TRUE],
          file.path(results_dir, "lcl_R_cross_assay_DL_nonzero_genes.csv"))
write_csv(lcl_meta[model_name == "lcl_cross_assay_exploratory" & mKH_CI_excludes_zero == TRUE],
          file.path(results_dir, "lcl_R_cross_assay_mKH_interval_supported_genes.csv"))

combined_summary <- rbind(
  data.table(branch = "placenta", placenta_summary),
  data.table(branch = "LCL", lcl_summary),
  fill = TRUE
)
write_csv(combined_summary, file.path(results_dir, "placenta_lcl_R_model_summary_combined.csv"))
message("Placenta/LCL meta-analysis and sensitivity models rebuilt.")
