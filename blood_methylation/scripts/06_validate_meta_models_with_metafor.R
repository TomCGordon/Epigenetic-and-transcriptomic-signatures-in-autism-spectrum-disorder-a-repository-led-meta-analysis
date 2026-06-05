#!/usr/bin/env Rscript

# Sampled validation that the explicit DerSimonian-Laird calculations in the
# blood workflow agree with metafor::rma.uni(method = "DL").

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
results_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_R_OUTPUT_DIR", unset = file.path(package_root, "results")),
                             winslash = "/", mustWork = FALSE)
qc_dir <- normalizePath(file.path(package_root, "qc"), winslash = "/", mustWork = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

mode <- Sys.getenv("BLOOD_VALIDATION_MODE", unset = "array_blanks_as_missing_R_default")
result_dir <- file.path(results_dir, mode)
effects <- fread(file.path(result_dir, "blood_R_all_dataset_level_effect_sizes.csv"), showProgress = FALSE)
meta <- fread(file.path(result_dir, "blood_R_meta_results_combined.csv"), showProgress = FALSE)

set.seed(20260514)
eligible <- meta[k >= 2 & is.finite(pooled_Hedges_g) & is.finite(SE)]
n_sample <- min(100L, nrow(eligible))
sampled <- eligible[sample(.N, n_sample)]

validate_one <- function(row) {
  datasets <- strsplit(row$contributing_datasets, ";", fixed = TRUE)[[1]]
  dat <- effects[gene == row$gene & dataset_id %in% datasets & is.finite(Hedges_g) & is.finite(variance_g) & variance_g > 0]
  fit <- metafor::rma.uni(yi = dat$Hedges_g, vi = dat$variance_g, method = "DL", test = "z")
  data.table(
    model_name = row$model_name,
    gene = row$gene,
    k = nrow(dat),
    branch_pooled_Hedges_g = row$pooled_Hedges_g,
    metafor_pooled_Hedges_g = as.numeric(fit$b[1, 1]),
    diff_pooled_Hedges_g = row$pooled_Hedges_g - as.numeric(fit$b[1, 1]),
    branch_SE = row$SE,
    metafor_SE = fit$se,
    diff_SE = row$SE - fit$se,
    branch_tau2 = row$tau2,
    metafor_tau2 = fit$tau2,
    diff_tau2 = row$tau2 - fit$tau2,
    branch_I2 = row$I2,
    metafor_I2 = fit$I2,
    diff_I2 = row$I2 - fit$I2
  )
}

checks <- rbindlist(lapply(seq_len(nrow(sampled)), function(i) validate_one(sampled[i])), fill = TRUE)
fwrite(checks, file.path(qc_dir, paste0("metafor_DL_validation_sample_", mode, ".csv")))

summary <- data.table(
  validation_mode = mode,
  sampled_models = nrow(checks),
  max_abs_effect_diff = max(abs(checks$diff_pooled_Hedges_g), na.rm = TRUE),
  max_abs_SE_diff = max(abs(checks$diff_SE), na.rm = TRUE),
  max_abs_tau2_diff = max(abs(checks$diff_tau2), na.rm = TRUE),
  max_abs_I2_diff = max(abs(checks$diff_I2), na.rm = TRUE)
)
fwrite(summary, file.path(qc_dir, paste0("metafor_DL_validation_summary_", mode, ".csv")))
print(summary)
