#!/usr/bin/env Rscript

# Brain methylation branch-level QC summaries.
#
# This script combines dataset-level Hedges' g inputs and writes branch-level
# model/sensitivity summaries for inspection. The reported pooled
# methylation models are fitted separately in scripts/04_meta_analysis/ using
# metafor.

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
source(file.path(script_dir, "lib", "brain_methylation_functions.R"))

array_processed <- file.path(package_root, "data_processed", "arrays")
wgbs_processed <- file.path(package_root, "data_processed", "WGBS")
results_dir <- file.path(package_root, "results", "strict_missing_R_default")
qc_dir <- file.path(package_root, "qc")
dir_create(results_dir)
dir_create(qc_dir)

universe <- fread(file.path(package_root, "data_processed", "annotation", "brain_harmonised_20960_gene_universe.csv"))$gene
array_effects <- fread(file.path(array_processed, "brain_array_all_dataset_effect_sizes.csv"))
wgbs_effects <- fread(file.path(wgbs_processed, "brain_WGBS_all_dataset_effect_sizes.csv"))
region_effect_tables <- list()
array_region_file <- file.path(array_processed, "brain_array_region_subroute_effect_sizes.csv")
wgbs_region_file <- file.path(wgbs_processed, "brain_WGBS_region_subroute_effect_sizes.csv")
if (file.exists(array_region_file)) region_effect_tables$array_region <- fread(array_region_file)
if (file.exists(wgbs_region_file)) region_effect_tables$wgbs_region <- fread(wgbs_region_file)
effects <- rbindlist(c(list(array_effects, wgbs_effects), region_effect_tables), fill = TRUE)
fwrite(effects, file.path(results_dir, "brain_R_all_dataset_level_effect_sizes.csv"))

array_summary <- fread(file.path(array_processed, "brain_array_all_dataset_gene_summary_statistics.csv"))
wgbs_summary <- fread(file.path(wgbs_processed, "brain_WGBS_all_dataset_gene_summary_statistics.csv"))
region_summary_tables <- list()
array_region_summary_file <- file.path(array_processed, "brain_array_region_subroute_gene_summary_statistics.csv")
wgbs_region_summary_file <- file.path(wgbs_processed, "brain_WGBS_region_subroute_gene_summary_statistics.csv")
if (file.exists(array_region_summary_file)) region_summary_tables$array_region <- fread(array_region_summary_file)
if (file.exists(wgbs_region_summary_file)) region_summary_tables$wgbs_region <- fread(wgbs_region_summary_file)
fwrite(rbindlist(c(list(array_summary, wgbs_summary), region_summary_tables), fill = TRUE),
       file.path(results_dir, "brain_R_all_dataset_gene_summary_statistics.csv"))

model_defs <- list(
  brain_array_HM27_primary = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608"),
  brain_450k_only_sensitivity = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285"),
  brain_grouped_primary_with_WGBS = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608", "GSE109875", "GSE81541"),
  brain_array_plus_GSE109875_sensitivity = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608", "GSE109875"),
  brain_array_plus_GSE81541_sensitivity = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608", "GSE81541"),
  cortex_only_sensitivity = c("GSE53162_prefrontal_cortex", "GSE53162_temporal_cortex", "GSE53924", "GSE80017", "GSE109875", "GSE81541_BA9", "GSE38608_occipital_cortex"),
  prefrontal_BA9_cortex_sensitivity = c("GSE53162_prefrontal_cortex", "GSE53924_frontal_cortex_ba10", "GSE80017", "GSE109875", "GSE81541_BA9"),
  cerebellum_only_sensitivity = c("GSE53162_cerebellum", "GSE278285", "GSE38608_cerebellum"),
  non_cortical_sensitivity = c("GSE53162_cerebellum", "GSE278285", "GSE38608_cerebellum", "GSE131706", "GSE242427"),
  without_large_cerebellum_sensitivity = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE38608", "GSE109875", "GSE81541"),
  WGBS_only_sensitivity = c("GSE109875", "GSE81541"),
  BA9_WGBS_only_sensitivity = c("GSE109875", "GSE81541_BA9"),
  WGBS_excluded_sensitivity = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE38608"),
  HM27_excluded_with_WGBS_sensitivity = c("GSE53162", "GSE53924", "GSE80017", "GSE131706", "GSE242427", "GSE278285", "GSE109875", "GSE81541")
)

model_roles <- data.table(
  model_name = names(model_defs),
  role = c("array/HM27 benchmark", "platform sensitivity", "primary", "WGBS sensitivity",
           "WGBS sensitivity", "region sensitivity", "region sensitivity",
           "region sensitivity", "region sensitivity", "influence sensitivity",
           "platform sensitivity", "region/WGBS sensitivity", "platform sensitivity",
           "platform sensitivity"),
  datasets = vapply(model_defs, paste, collapse = "; ", FUN.VALUE = character(1)),
  biological_interpretation = c(
    "Broad array/HM27 grouped-brain benchmark.",
    "450K-only grouped-brain sensitivity.",
    "Grouped post-mortem brain synthesis including broad-processable array/HM27 and WGBS datasets.",
    "Array/HM27 plus BA9 cortical WGBS only.",
    "Array/HM27 plus GSE81541 brain WGBS only.",
    "Cortex-only model using R-generated region-specific subroutes where parent accessions include multiple brain regions.",
    "Prefrontal/BA9-compatible cortex model using R-generated region-specific subroutes.",
    "Cerebellum-only sensitivity model using R-generated region-specific subroutes plus GSE278285.",
    "Non-cortical sensitivity model including cerebellar, subventricular-zone and dorsal-raphe evidence.",
    "WGBS-integrated model excluding large cerebellum-only GSE278285.",
    "WGBS-only coordinate-based model.",
    "BA9-compatible WGBS-only sensitivity model.",
    "WGBS-excluded array/HM27 benchmark.",
    "Excludes HM27 but includes WGBS."
  )
)
fwrite(model_roles, file.path(results_dir, "brain_R_model_definitions.csv"))

meta_results <- build_meta_results(effects, universe, model_defs, model_roles)
model_summary <- summarise_models(meta_results, model_roles)

fwrite(meta_results, file.path(results_dir, "brain_R_meta_results_combined.csv"))
fwrite(model_summary, file.path(results_dir, "brain_R_model_summary.csv"))

fwrite(meta_results[model_name == "brain_grouped_primary_with_WGBS" & DL_CI_excludes_zero == TRUE],
       file.path(results_dir, "brain_R_primary_DL_nonzero_genes.csv"))
fwrite(meta_results[model_name == "brain_grouped_primary_with_WGBS" & FDR_significant == TRUE],
       file.path(results_dir, "brain_R_primary_FDR_significant_genes.csv"))
fwrite(meta_results[model_name == "brain_grouped_primary_with_WGBS" & mKH_CI_excludes_zero == TRUE],
       file.path(results_dir, "brain_R_primary_mKH_interval_supported_genes.csv"))

changed <- merge(
  meta_results[model_name == "brain_array_HM27_primary", .(gene, array_DL = DL_CI_excludes_zero, array_mKH = mKH_CI_excludes_zero, array_FDR = FDR_significant)],
  meta_results[model_name == "brain_grouped_primary_with_WGBS", .(gene, wgbs_DL = DL_CI_excludes_zero, wgbs_mKH = mKH_CI_excludes_zero, wgbs_FDR = FDR_significant)],
  by = "gene", all = TRUE
)
changed[, `:=`(
  DL_status_changed = array_DL != wgbs_DL,
  mKH_status_changed = array_mKH != wgbs_mKH,
  FDR_status_changed = array_FDR != wgbs_FDR
)]
fwrite(changed, file.path(results_dir, "brain_R_WGBS_status_change_table.csv"))

message("Brain meta-analysis and sensitivity models rebuilt.")
