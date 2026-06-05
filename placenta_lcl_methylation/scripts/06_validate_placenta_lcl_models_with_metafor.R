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
source(file.path(script_dir, "lib", "placenta_lcl_methylation_functions.R"))

results_dir <- file.path(package_root, "results", "strict_missing_R_default")
qc_dir <- file.path(package_root, "qc")
dir_create(qc_dir)

validate_branch <- function(branch, effects_file, meta_file) {
  effects <- fread(effects_file)
  meta <- fread(meta_file)
  candidates <- meta[k >= 2 & is.finite(pooled_Hedges_g)]
  if (!nrow(candidates)) return(list(sample = data.table(), summary = data.table(branch = branch, sampled_models = 0L)))
  set.seed(20260515)
  sample_rows <- candidates[sample(.N, min(.N, 100L))]
  rows <- rbindlist(lapply(seq_len(nrow(sample_rows)), function(i) {
    r <- sample_rows[i]
    dat <- effects[gene == r$gene & accession %in% strsplit(r$contributing_datasets, ";", fixed = TRUE)[[1]]]
    fit <- metafor::rma.uni(yi = dat$Hedges_g, vi = dat$variance_g, method = "DL", test = "z")
    data.table(
      branch = branch,
      model_name = r$model_name,
      gene = r$gene,
      branch_pooled_g = r$pooled_Hedges_g,
      metafor_pooled_g = as.numeric(fit$b[1]),
      pooled_abs_diff = abs(r$pooled_Hedges_g - as.numeric(fit$b[1])),
      branch_SE = r$SE,
      metafor_SE = fit$se,
      SE_abs_diff = abs(r$SE - fit$se),
      branch_tau2 = r$tau2,
      metafor_tau2 = fit$tau2,
      tau2_abs_diff = abs(r$tau2 - fit$tau2),
      branch_I2 = r$I2,
      metafor_I2 = fit$I2,
      I2_abs_diff = abs(r$I2 - fit$I2)
    )
  }), fill = TRUE)
  summary <- rows[, .(
    sampled_models = .N,
    max_abs_effect_diff = max(pooled_abs_diff, na.rm = TRUE),
    max_abs_SE_diff = max(SE_abs_diff, na.rm = TRUE),
    max_abs_tau2_diff = max(tau2_abs_diff, na.rm = TRUE),
    max_abs_I2_diff = max(I2_abs_diff, na.rm = TRUE)
  )]
  summary[, branch := branch]
  list(sample = rows, summary = summary)
}

placenta_val <- validate_branch("placenta",
                                file.path(results_dir, "placenta_R_all_dataset_level_effect_sizes.csv"),
                                file.path(results_dir, "placenta_R_meta_results_combined.csv"))
lcl_val <- validate_branch("LCL",
                           file.path(results_dir, "lcl_R_all_dataset_level_effect_sizes.csv"),
                           file.path(results_dir, "lcl_R_meta_results_combined.csv"))

write_csv(rbind(placenta_val$sample, lcl_val$sample, fill = TRUE),
          file.path(qc_dir, "metafor_DL_validation_sample_placenta_lcl_strict_missing_R_default.csv"))
write_csv(rbind(placenta_val$summary, lcl_val$summary, fill = TRUE),
          file.path(qc_dir, "metafor_DL_validation_summary_placenta_lcl_strict_missing_R_default.csv"))
message("Placenta/LCL DL validation against metafor completed.")
