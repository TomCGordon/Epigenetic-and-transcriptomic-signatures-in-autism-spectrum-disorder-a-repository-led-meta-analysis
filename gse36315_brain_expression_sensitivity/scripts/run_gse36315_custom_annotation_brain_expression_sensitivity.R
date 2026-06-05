#!/usr/bin/env Rscript

# GSE36315 custom-annotated brain expression sensitivity analyses.
#
# GSE36315 is not included in the primary post-mortem brain expression model
# because the public platform annotation is incomplete for broad gene-level
# synthesis. This script adds it only as a sensitivity analysis after remapping
# public GPL15314 probe sequences to GENCODE v19 transcripts.

suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_file)
pipeline_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(pipeline_root, "..", ".."), winslash = "/", mustWork = TRUE)

base_pkg <- file.path(package_root, "pipelines", "brain_expression")
remap_pkg <- file.path(package_root, "pipelines", "gse36315_custom_probe_remapping")
out_dir <- pipeline_root

for (d in c("00_manifest", "01_gse36315_inputs", "02_effect_size_inputs",
            "03_grouped_brain_sensitivity", "04_region_subtissue_sensitivity",
            "05_platform_sensitivity", "06_status_transitions", "07_quality_control",
            "08_reports", "scripts")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

source(file.path(base_pkg, "scripts/lib/brain_expression_functions.R"))

raw_file <- file.path(base_pkg, "05_effect_sizes/brain_expression_R_dataset_level_effect_sizes_raw_strata.csv")
collapsed_file <- file.path(base_pkg, "05_effect_sizes/brain_expression_R_dataset_level_effect_sizes_grouped_collapsed.csv")
gse_file <- file.path(remap_pkg, "04_recoverability_assessment/gse36315_custom_annotated_region_effect_sizes.csv")
remap_summary_file <- file.path(remap_pkg, "04_recoverability_assessment/GSE36315_custom_bowtie_remapping_summary.csv")
gene_recovery_file <- file.path(remap_pkg, "04_recoverability_assessment/GSE36315_gene_level_recoverability_summary.csv")

stopifnot(file.exists(raw_file), file.exists(collapsed_file), file.exists(gse_file))

raw_base <- fread(raw_file)
collapsed_base <- fread(collapsed_file)
gse_raw <- fread(gse_file)

gse_effects <- copy(gse_raw)
gse_effects[, dataset := fifelse(region == "prefrontal cortex",
                                 "GSE36315_prefrontal_cortex_custom_annotation",
                                 "GSE36315_cerebellum_custom_annotation")]
gse_effects[, study_id := "GSE36315"]
gse_effects[, source_accession := "GSE36315"]
gse_effects[, assay := "microarray"]
gse_effects[, platform := "GPL15314 custom GENCODEv19 remap"]
gse_effects[, brain_region := region]
gse_effects[, region_group := fifelse(region == "prefrontal cortex", "cortex", "cerebellum")]
gse_effects[, model_role := "custom_annotation_sensitivity_only"]
gse_effects[, ASD_n := n_asd]
gse_effects[, control_n := n_control]
gse_effects[, ASD_mean := NA_real_]
gse_effects[, control_mean := NA_real_]
gse_effects[, ASD_sd := NA_real_]
gse_effects[, control_sd := NA_real_]
gse_effects[, direction := fifelse(hedges_g > 0, "higher_in_ASD", fifelse(hedges_g < 0, "lower_in_ASD", "zero"))]
gse_effects[, source_file := "GSE36315_series_matrix.txt.gz + GPL15314 public sequence remap"]
gse_effects[, transform_note := "GEO series matrix values collapsed to genes after custom GPL15314 probe remapping to GENCODE v19 GRCh37 with Rbowtie; custom-annotated sensitivity only"]
gse_effects <- gse_effects[, names(raw_base), with = FALSE]

write_csv_safe(gse_effects, file.path(out_dir, "01_gse36315_inputs/gse36315_effect_sizes.csv"))
if (file.exists(remap_summary_file)) file.copy(remap_summary_file, file.path(out_dir, "01_gse36315_inputs/GSE36315_custom_bowtie_remapping_summary.csv"), overwrite = TRUE)
if (file.exists(gene_recovery_file)) file.copy(gene_recovery_file, file.path(out_dir, "01_gse36315_inputs/GSE36315_gene_level_recoverability_summary.csv"), overwrite = TRUE)

gse_pfc <- gse_effects[brain_region == "prefrontal cortex"]
gse_cere <- gse_effects[brain_region == "cerebellum"]

raw_plus <- rbindlist(list(raw_base, gse_effects), fill = TRUE)
collapsed_plus_both <- as.data.table(collapse_within_study(raw_plus))
collapsed_plus_pfc <- rbindlist(list(collapsed_base, gse_pfc), fill = TRUE)
collapsed_plus_cere <- rbindlist(list(collapsed_base, gse_cere), fill = TRUE)

write_csv_safe(raw_plus, file.path(out_dir, "02_effect_size_inputs/raw_plus_GSE36315.csv"))
write_csv_safe(collapsed_plus_both, file.path(out_dir, "02_effect_size_inputs/collapsed_plus_GSE36315_both.csv"))
write_csv_safe(collapsed_plus_pfc, file.path(out_dir, "02_effect_size_inputs/collapsed_plus_GSE36315_pfc.csv"))
write_csv_safe(collapsed_plus_cere, file.path(out_dir, "02_effect_size_inputs/collapsed_plus_GSE36315_cere.csv"))

models <- list()
models[["baseline_grouped_recalculated"]] <- run_meta_model(
  collapsed_base,
  "brain_expression_grouped_primary_public_R_recalculated_for_GSE36315_sensitivity",
  "primary_reference",
  "Reference grouped post-mortem brain expression model recalculated from existing collapsed effect sizes"
)
models[["plus_gse36315_prefrontal_only"]] <- run_meta_model(
  collapsed_plus_pfc,
  "brain_expression_grouped_plus_GSE36315_prefrontal_custom_annotation_sensitivity_R",
  "custom_annotation_sensitivity",
  "Grouped post-mortem brain expression sensitivity adding GSE36315 prefrontal cortex only; custom probe remapping caveat"
)
models[["plus_gse36315_cerebellum_only"]] <- run_meta_model(
  collapsed_plus_cere,
  "brain_expression_grouped_plus_GSE36315_cerebellum_custom_annotation_sensitivity_R",
  "custom_annotation_sensitivity",
  "Grouped post-mortem brain expression sensitivity adding GSE36315 cerebellum only; custom probe remapping caveat"
)
models[["plus_gse36315_two_region_collapsed"]] <- run_meta_model(
  collapsed_plus_both,
  "brain_expression_grouped_plus_GSE36315_two_region_collapsed_custom_annotation_sensitivity_R",
  "custom_annotation_sensitivity",
  "Grouped post-mortem brain expression sensitivity adding GSE36315 after within-study prefrontal/cerebellum collapse; same-donor region caveat"
)
models[["microarray_plus_gse36315_two_region_collapsed"]] <- run_meta_model(
  collapsed_plus_both[grepl("microarray", assay, ignore.case = TRUE)],
  "brain_expression_microarray_plus_gse36315_custom_annotated_sensitivity_R",
  "custom_annotation_platform_sensitivity",
  "Microarray-only sensitivity including custom-annotated GSE36315 after within-study region collapse"
)
models[["cortex_plus_gse36315_prefrontal"]] <- run_meta_model(
  rbindlist(list(raw_base[grepl("cortex", region_group, ignore.case = TRUE)], gse_pfc), fill = TRUE),
  "brain_expression_cortex_plus_GSE36315_prefrontal_custom_annotation_sensitivity_R",
  "custom_annotation_region_sensitivity",
  "Cortex-only sensitivity adding custom-annotated GSE36315 prefrontal cortex"
)
models[["prefrontal_plus_gse36315_prefrontal"]] <- run_meta_model(
  rbindlist(list(raw_base[grepl("prefrontal", brain_region, ignore.case = TRUE)], gse_pfc), fill = TRUE),
  "brain_expression_prefrontal_plus_gse36315_custom_annotated_sensitivity_R",
  "custom_annotation_region_sensitivity",
  "Prefrontal cortex sensitivity adding custom-annotated GSE36315 prefrontal cortex"
)
models[["cerebellum_plus_gse36315_cerebellum"]] <- run_meta_model(
  rbindlist(list(raw_base[grepl("cerebellum", region_group, ignore.case = TRUE)], gse_cere), fill = TRUE),
  "brain_expression_cerebellum_plus_gse36315_custom_annotated_sensitivity_R",
  "custom_annotation_region_sensitivity",
  "Cerebellum/Purkinje sensitivity adding custom-annotated GSE36315 cerebellum"
)
models[["gse36315_prefrontal_descriptive"]] <- run_meta_model(
  gse_pfc,
  "GSE36315_prefrontal_custom_annotation_descriptive_R",
  "descriptive",
  "Single custom-annotated GSE36315 prefrontal cortex route; no pooled meta-analysis"
)
models[["gse36315_cerebellum_descriptive"]] <- run_meta_model(
  gse_cere,
  "GSE36315_cerebellum_custom_annotation_descriptive_R",
  "descriptive",
  "Single custom-annotated GSE36315 cerebellum route; no pooled meta-analysis"
)

all_summary <- rbindlist_fill(lapply(models, `[[`, "summary"))
all_meta <- rbindlist_fill(lapply(models, `[[`, "meta"))
all_k1 <- rbindlist_fill(lapply(models, `[[`, "k1"))

write_csv_safe(all_summary, file.path(out_dir, "03_grouped_brain_sensitivity/model_summary.csv"))
write_csv_safe(all_meta, file.path(out_dir, "03_grouped_brain_sensitivity/all_meta_results.csv"))
write_csv_safe(all_k1, file.path(out_dir, "03_grouped_brain_sensitivity/all_k1_rows.csv"))

region_models <- c("cortex_plus_gse36315_prefrontal", "prefrontal_plus_gse36315_prefrontal", "cerebellum_plus_gse36315_cerebellum",
                   "gse36315_prefrontal_descriptive", "gse36315_cerebellum_descriptive")
write_csv_safe(rbindlist_fill(lapply(models[region_models], `[[`, "meta")), file.path(out_dir, "04_region_subtissue_sensitivity/region_results.csv"))
write_csv_safe(rbindlist_fill(lapply(models[region_models], `[[`, "summary")), file.path(out_dir, "04_region_subtissue_sensitivity/region_summary.csv"))
write_csv_safe(models[["microarray_plus_gse36315_two_region_collapsed"]]$meta, file.path(out_dir, "05_platform_sensitivity/microarray_results.csv"))
write_csv_safe(models[["microarray_plus_gse36315_two_region_collapsed"]]$summary, file.path(out_dir, "05_platform_sensitivity/microarray_summary.csv"))

status_cols <- c("gene", "pooled_g", "p_value", "FDR", "DL_nonzero", "FDR_significant", "mKH_interval_excludes_zero", "direction", "k", "I2")
transition_table <- function(reference, comparison, comparison_name) {
  ref <- as.data.table(reference)[, ..status_cols]
  cmp <- as.data.table(comparison)[, ..status_cols]
  setnames(ref, setdiff(names(ref), "gene"), paste0("reference_", setdiff(names(ref), "gene")))
  setnames(cmp, setdiff(names(cmp), "gene"), paste0("comparison_", setdiff(names(cmp), "gene")))
  out <- merge(ref, cmp, by = "gene", all = TRUE)
  out[, comparison_model := comparison_name]
  for (field in c("DL_nonzero", "FDR_significant", "mKH_interval_excludes_zero")) {
    r <- paste0("reference_", field)
    c <- paste0("comparison_", field)
    out[is.na(get(r)), (r) := FALSE]
    out[is.na(get(c)), (c) := FALSE]
    out[, paste0(field, "_transition") := fifelse(!get(r) & get(c), "gained",
                                                  fifelse(get(r) & !get(c), "lost",
                                                          fifelse(get(r) & get(c), "retained", "neither")))]
  }
  out[]
}

transitions <- rbindlist(list(
  transition_table(models[["baseline_grouped_recalculated"]]$meta, models[["plus_gse36315_prefrontal_only"]]$meta, "plus_GSE36315_prefrontal_only"),
  transition_table(models[["baseline_grouped_recalculated"]]$meta, models[["plus_gse36315_cerebellum_only"]]$meta, "plus_GSE36315_cerebellum_only"),
  transition_table(models[["baseline_grouped_recalculated"]]$meta, models[["plus_gse36315_two_region_collapsed"]]$meta, "plus_GSE36315_two_region_collapsed")
), fill = TRUE)
write_csv_safe(transitions, file.path(out_dir, "06_status_transitions/status_transitions.csv"))

transition_summary <- transitions[, .(
  genes_in_comparison = .N,
  DL_gained = sum(DL_nonzero_transition == "gained", na.rm = TRUE),
  DL_lost = sum(DL_nonzero_transition == "lost", na.rm = TRUE),
  DL_retained = sum(DL_nonzero_transition == "retained", na.rm = TRUE),
  FDR_gained = sum(FDR_significant_transition == "gained", na.rm = TRUE),
  FDR_lost = sum(FDR_significant_transition == "lost", na.rm = TRUE),
  FDR_retained = sum(FDR_significant_transition == "retained", na.rm = TRUE),
  mKH_gained = sum(mKH_interval_excludes_zero_transition == "gained", na.rm = TRUE),
  mKH_lost = sum(mKH_interval_excludes_zero_transition == "lost", na.rm = TRUE),
  mKH_retained = sum(mKH_interval_excludes_zero_transition == "retained", na.rm = TRUE)
), by = comparison_model]
write_csv_safe(transition_summary, file.path(out_dir, "06_status_transitions/status_transition_summary.csv"))

set.seed(36315)
validation_rows <- list()
validation_models <- c("plus_gse36315_prefrontal_only", "plus_gse36315_cerebellum_only", "plus_gse36315_two_region_collapsed")
for (mname in validation_models) {
  effects <- switch(mname,
                    plus_gse36315_prefrontal_only = collapsed_plus_pfc,
                    plus_gse36315_cerebellum_only = collapsed_plus_cere,
                    plus_gse36315_two_region_collapsed = collapsed_plus_both)
  meta <- as.data.table(models[[mname]]$meta)
  if (!nrow(meta)) next
  sample_genes <- sample(meta$gene, min(50, nrow(meta)))
  validation_rows[[mname]] <- rbindlist(lapply(sample_genes, function(gene_id) {
    e <- effects[effects$gene == gene_id & is.finite(hedges_g) & is.finite(variance) & variance > 0]
    fit <- tryCatch(metafor::rma.uni(yi = e$hedges_g, vi = e$variance, method = "DL", test = "z"), error = function(err) NULL)
    branch_row <- meta[meta$gene == gene_id][1]
    if (is.null(fit)) {
      return(data.table(model_key = mname, gene = gene_id, metafor_ok = FALSE))
    }
    data.table(model_key = mname, gene = gene_id, metafor_ok = TRUE,
               branch_pooled_g = branch_row$pooled_g,
               metafor_pooled_g = as.numeric(fit$b),
               abs_diff_pooled_g = abs(branch_row$pooled_g - as.numeric(fit$b)),
               branch_tau2 = branch_row$tau2,
               metafor_tau2 = fit$tau2,
               abs_diff_tau2 = abs(branch_row$tau2 - fit$tau2))
  }), fill = TRUE)
}
validation <- rbindlist(validation_rows, fill = TRUE)
validation[, pass := metafor_ok & abs_diff_pooled_g < 1e-8 & abs_diff_tau2 < 1e-8]
write_csv_safe(validation, file.path(out_dir, "07_quality_control/metafor_validation_sample.csv"))

baseline_source_summary <- fread(file.path(base_pkg, "06_grouped_brain_meta_analysis/brain_expression_R_model_summary.csv"))
baseline_recalc_summary <- models[["baseline_grouped_recalculated"]]$summary
source_row <- baseline_source_summary[model == "brain_expression_grouped_primary_public_R"]
recalc_row <- as.data.table(baseline_recalc_summary)
compare_metrics <- c("genes_meta_analysed", "k1_descriptive", "DL_nonzero", "FDR_significant",
                     "mKH_interval_supported", "FDR_mKH_overlap", "median_k", "median_I2",
                     "high_I2_gt50", "high_I2_gt75")
baseline_comparison <- data.table(
  metric = compare_metrics,
  source_value = vapply(compare_metrics, function(m) as.character(source_row[[m]][1]), character(1)),
  recalculated_value = vapply(compare_metrics, function(m) as.character(recalc_row[[m]][1]), character(1))
)
baseline_comparison[, matches_source := source_value == recalculated_value]
write_csv_safe(baseline_comparison, file.path(out_dir, "07_quality_control/baseline_recalc_check.csv"))

pkg_versions <- data.table(
  package = c("data.table", "metafor", "openxlsx"),
  version = vapply(c("data.table", "metafor", "openxlsx"), function(p) as.character(utils::packageVersion(p)), character(1))
)
write_csv_safe(pkg_versions, file.path(out_dir, "00_manifest/R_package_versions.csv"))

file_index <- data.table(
  file_path = list.files(out_dir, recursive = TRUE, full.names = TRUE),
  bytes = file.info(list.files(out_dir, recursive = TRUE, full.names = TRUE))$size
)
file_index[, file_path := normalizePath(file_path, winslash = "/", mustWork = FALSE)]
write_csv_safe(file_index, file.path(out_dir, "00_manifest/file_index.csv"))

openxlsx::write.xlsx(
  list(
    model_summary = as.data.frame(all_summary),
    status_transition_summary = as.data.frame(transition_summary),
    grouped_plus_two_region = as.data.frame(models[["plus_gse36315_two_region_collapsed"]]$meta),
    grouped_plus_prefrontal = as.data.frame(models[["plus_gse36315_prefrontal_only"]]$meta),
    grouped_plus_cerebellum = as.data.frame(models[["plus_gse36315_cerebellum_only"]]$meta),
    region_summary = as.data.frame(rbindlist_fill(lapply(models[region_models], `[[`, "summary"))),
    metafor_validation = as.data.frame(validation)
  ),
  file.path(out_dir, "GSE36315_sensitivity_results.xlsx"),
  overwrite = TRUE
)

qc_lines <- c(
  "# GSE36315 Custom-Annotated Brain Expression Sensitivity QC Report",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Base brain expression package: `", normalizePath(base_pkg, winslash = "/", mustWork = FALSE), "`"),
  paste0("- Custom GSE36315 remapping package: `", normalizePath(remap_pkg, winslash = "/", mustWork = FALSE), "`"),
  "- GSE36315 was not added to the primary model; all new models are labelled custom-annotated sensitivity or descriptive.",
  "",
  "## Key Model Summary",
  "",
  paste(capture.output(print(all_summary)), collapse = "\n"),
  "",
  "## Status Transition Summary",
  "",
  paste(capture.output(print(transition_summary)), collapse = "\n"),
  "",
  "## QC",
  "",
  paste0("- Metafor validation rows: ", nrow(validation), ". Passed: ", sum(validation$pass, na.rm = TRUE), ". Failed: ", sum(!validation$pass, na.rm = TRUE), "."),
  "- Baseline grouped model was recalculated from existing collapsed effect sizes to provide a reference comparison.",
  "- GSE36315 same-donor prefrontal/cerebellum structure was handled by reporting prefrontal-only, cerebellum-only, and within-study-collapsed sensitivity variants rather than treating both regions as independent primary datasets.",
    "- Outputs were written only to the configured GSE36315 sensitivity results directory.",
  "",
  "## Recommendation",
  "",
  "Keep the original grouped post-mortem brain expression model as the primary brain expression result. Report GSE36315, if used, as a custom-annotated sensitivity analysis only."
)
writeLines(qc_lines, file.path(out_dir, "07_quality_control/QC_report.md"))

readme_lines <- c(
  "# Brain Expression GSE36315 Custom-Annotated Sensitivity Package",
  "",
  "This package adds GSE36315 only as a sensitivity analysis after custom public-sequence remapping of GPL15314 probes to GENCODE v19 transcripts.",
  "",
  "Run from the repository root:",
  "",
  "```r",
  "source(\"pipelines/gse36315_brain_expression_sensitivity/scripts/run_gse36315_custom_annotation_brain_expression_sensitivity.R\")",
  "```",
  "",
  "or:",
  "",
  "```bash",
  "Rscript \"pipelines/gse36315_brain_expression_sensitivity/scripts/run_gse36315_custom_annotation_brain_expression_sensitivity.R\" .",
  "```",
  "",
  "Required prior inputs:",
  "",
  "- `brain_expression`",
  "- `gse36315_custom_probe_remapping`",
  "",
  "Interpretation: sensitivity/descriptive only; do not replace the primary grouped-brain expression model without official Arraystar/author annotation validation."
)
writeLines(readme_lines, file.path(out_dir, "README_reproducibility.md"))

message("Finished GSE36315 custom-annotated sensitivity package: ", out_dir)
