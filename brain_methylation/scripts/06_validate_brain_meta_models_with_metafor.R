#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

results_dir <- file.path(package_root, "results", "strict_missing_R_default")
qc_dir <- file.path(package_root, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

effects <- fread(file.path(results_dir, "brain_R_all_dataset_level_effect_sizes.csv"))
meta_results <- fread(file.path(results_dir, "brain_R_meta_results_combined.csv"))
model_defs <- fread(file.path(results_dir, "brain_R_model_definitions.csv"))

set.seed(20260515)
eligible <- meta_results[k >= 2 & is.finite(pooled_Hedges_g) & is.finite(SE)]
if (nrow(eligible) > 100) eligible <- eligible[sample(.N, 100)]

checks <- rbindlist(lapply(seq_len(nrow(eligible)), function(i) {
  row <- eligible[i]
  datasets <- strsplit(model_defs[model_name == row$model_name]$datasets, ";\\s*")[[1]]
  e <- effects[gene == row$gene & accession %in% datasets & is.finite(Hedges_g) & is.finite(variance_g) & variance_g > 0]
  fit <- metafor::rma.uni(yi = e$Hedges_g, vi = e$variance_g, method = "DL", test = "z")
  data.table(
    model_name = row$model_name,
    gene = row$gene,
    k = row$k,
    branch_effect = row$pooled_Hedges_g,
    metafor_effect = as.numeric(fit$b[1, 1]),
    abs_effect_diff = abs(row$pooled_Hedges_g - as.numeric(fit$b[1, 1])),
    branch_SE = row$SE,
    metafor_SE = fit$se,
    abs_SE_diff = abs(row$SE - fit$se),
    branch_tau2 = row$tau2,
    metafor_tau2 = fit$tau2,
    abs_tau2_diff = abs(row$tau2 - fit$tau2),
    branch_I2 = row$I2,
    metafor_I2 = fit$I2,
    abs_I2_diff = abs(row$I2 - fit$I2)
  )
}), fill = TRUE)

summary <- checks[, .(
  sampled_models = .N,
  max_abs_effect_diff = max(abs_effect_diff, na.rm = TRUE),
  max_abs_SE_diff = max(abs_SE_diff, na.rm = TRUE),
  max_abs_tau2_diff = max(abs_tau2_diff, na.rm = TRUE),
  max_abs_I2_diff = max(abs_I2_diff, na.rm = TRUE)
)]

fwrite(checks, file.path(qc_dir, "metafor_DL_validation_sample_brain_strict_missing_R_default.csv"))
fwrite(summary, file.path(qc_dir, "metafor_DL_validation_summary_brain_strict_missing_R_default.csv"))
message("Brain DL validation against metafor completed.")
