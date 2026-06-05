#!/usr/bin/env Rscript

# Shared functions for tissue-stratified random-effects meta-analysis.
#
# These functions are used by the methylation and expression model-fitting
# scripts. They take dataset-level Hedges' g estimates and sampling variances as
# input and fit DerSimonian-Laird random-effects models using metafor.

suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(openxlsx)
})

options(stringsAsFactors = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

safe_name <- function(x, n = 55) {
  sx <- gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", x))
  if (nchar(sx) > n) sx <- substr(sx, 1, n)
  sx
}

split_dataset_string <- function(x) {
  trimws(unlist(strsplit(x, ";", fixed = TRUE)))
}

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

q_safe <- function(x, p) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}

read_required <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Required input not found: ", file_path, call. = FALSE)
  }
  fread(file_path)
}

standardise_effects <- function(dt, source_label, dataset_col, yi_col, vi_col,
                                n_asd_col = NULL, n_ctrl_col = NULL,
                                extra_cols = character()) {
  required <- c(dataset_col, yi_col, vi_col, "gene")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop("Input table is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  out <- data.table(
    source_label = source_label,
    dataset = as.character(dt[[dataset_col]]),
    gene = as.character(dt[["gene"]]),
    yi = suppressWarnings(as.numeric(dt[[yi_col]])),
    vi = suppressWarnings(as.numeric(dt[[vi_col]]))
  )
  out[, se := sqrt(vi)]

  if (!is.null(n_asd_col) && n_asd_col %in% names(dt)) {
    out[, ASD_n := suppressWarnings(as.numeric(dt[[n_asd_col]]))]
  } else {
    out[, ASD_n := NA_real_]
  }

  if (!is.null(n_ctrl_col) && n_ctrl_col %in% names(dt)) {
    out[, control_n := suppressWarnings(as.numeric(dt[[n_ctrl_col]]))]
  } else {
    out[, control_n := NA_real_]
  }

  for (cc in extra_cols) {
    if (cc %in% names(dt)) out[, (cc) := dt[[cc]]]
  }

  out[is.finite(yi) & is.finite(vi) & vi > 0 & !is.na(gene) & gene != ""]
}

collapse_within_study_platforms <- function(dt) {
  dt <- as.data.table(dt)
  required <- c("dataset", "study_id", "model_role", "platform_id", "gene",
                "feature_count", "n_asd", "n_control", "hedges_g", "variance_g")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop("Blood-expression platform-stratified table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dt <- dt[is.finite(hedges_g) & is.finite(variance_g) & variance_g > 0]
  multi <- dt[, .N, by = .(study_id, gene)][N > 1]
  if (!nrow(multi)) return(dt)

  multi_key <- paste(multi$study_id, multi$gene, sep = "||")
  dt[, collapse_key := paste(study_id, gene, sep = "||")]
  keep <- dt[!collapse_key %in% multi_key]

  collapsed <- dt[collapse_key %in% multi_key, {
    w <- 1 / variance_g
    g <- sum(w * hedges_g) / sum(w)
    v <- 1 / sum(w)
    .(
      dataset = study_id[1],
      study_id = study_id[1],
      model_role = model_role[1],
      platform_id = paste(sort(unique(platform_id)), collapse = ";"),
      gene = gene[1],
      feature_count = sum(feature_count, na.rm = TRUE),
      n_asd = sum(unique(n_asd)),
      n_control = sum(unique(n_control)),
      ASD_mean = NA_real_,
      control_mean = NA_real_,
      ASD_sd = NA_real_,
      control_sd = NA_real_,
      hedges_g = g,
      variance_g = v,
      se_g = sqrt(v),
      ci_lower = g - 1.96 * sqrt(v),
      ci_upper = g + 1.96 * sqrt(v),
      direction = ifelse(g > 0, "ASD_higher", "ASD_lower"),
      notes = paste("within-study platform strata collapsed by inverse-variance fixed effect:",
                    paste(sort(unique(dataset)), collapse = ";"))
    )
  }, by = collapse_key]

  keep[, collapse_key := NULL]
  collapsed[, collapse_key := NULL]
  rbindlist(list(keep, collapsed), fill = TRUE)
}

fit_one_gene <- function(gdat) {
  if (nrow(gdat) < 2) return(NULL)

  dl <- tryCatch(
    metafor::rma.uni(yi = gdat$yi, vi = gdat$vi, method = "DL", test = "z"),
    error = function(e) NULL
  )
  if (is.null(dl)) return(NULL)

  mkh <- tryCatch(
    metafor::rma.uni(yi = gdat$yi, vi = gdat$vi, method = "DL", test = "adhoc"),
    error = function(e) dl
  )

  contributor <- if ("study_id" %in% names(gdat) && any(!is.na(gdat$study_id))) {
    gdat$study_id
  } else {
    gdat$dataset
  }

  pooled <- as.numeric(dl$b[1])
  data.table(
    k = nrow(gdat),
    contributing_datasets = paste(unique(contributor), collapse = ";"),
    total_ASD_n = sum(gdat$ASD_n, na.rm = TRUE),
    total_control_n = sum(gdat$control_n, na.rm = TRUE),
    pooled_Hedges_g = pooled,
    DL_SE = as.numeric(dl$se),
    DL_CI_lower = as.numeric(dl$ci.lb),
    DL_CI_upper = as.numeric(dl$ci.ub),
    DL_z_or_stat = as.numeric(dl$zval),
    p_value = as.numeric(dl$pval),
    tau2 = as.numeric(dl$tau2),
    I2 = as.numeric(dl$I2),
    H2 = as.numeric(dl$H2),
    Q = as.numeric(dl$QE),
    Q_p_value = as.numeric(dl$QEp),
    mKH_SE = as.numeric(mkh$se),
    mKH_CI_lower = as.numeric(mkh$ci.lb),
    mKH_CI_upper = as.numeric(mkh$ci.ub),
    mKH_p_value = as.numeric(mkh$pval),
    direction = ifelse(pooled > 0, "ASD_higher", ifelse(pooled < 0, "ASD_lower", "zero"))
  )
}

run_meta_model <- function(effects, model_name, omic_layer, tissue_branch, model_role,
                           datasets, notes = "", sensitivity_family = "") {
  selected <- effects[dataset %in% datasets]
  if ("study_id" %in% names(effects)) {
    selected <- unique(rbindlist(list(selected, effects[study_id %in% datasets]), fill = TRUE))
  }
  selected <- selected[is.finite(yi) & is.finite(vi) & vi > 0]

  if (!nrow(selected)) {
    return(list(
      results = data.table(),
      k1 = data.table(),
      summary = data.table(
        omic_layer, tissue_branch, model_name, model_role, sensitivity_family,
        datasets_included = paste(datasets, collapse = ";"),
        genes_with_any_result = 0L,
        genes_meta_analysed = 0L,
        k1_descriptive = 0L,
        DL_nonzero = 0L,
        FDR_significant = 0L,
        mKH_interval_supported = 0L,
        FDR_mKH_overlap = 0L,
        median_k = NA_real_,
        median_I2 = NA_real_,
        IQR_I2_lower = NA_real_,
        IQR_I2_upper = NA_real_,
        I2_gt50 = 0L,
        I2_gt75 = 0L,
        median_tau2 = NA_real_,
        max_k = 0L,
        max_total_ASD_n = NA_real_,
        max_total_control_n = NA_real_,
        notes = notes
      )
    ))
  }

  gene_counts <- selected[, .(k = .N), by = gene]
  k1_genes <- gene_counts[k == 1, gene]
  meta_genes <- gene_counts[k >= 2, gene]

  k1 <- selected[gene %in% k1_genes]
  if (nrow(k1)) {
    k1 <- k1[, .(
      omic_layer = omic_layer,
      tissue_branch = tissue_branch,
      model_name = model_name,
      model_role = model_role,
      gene,
      dataset,
      Hedges_g = yi,
      variance = vi,
      ASD_n,
      control_n
    )]
  }

  res <- data.table()
  if (length(meta_genes)) {
    setkey(selected, gene)
    pieces <- vector("list", length(meta_genes))
    for (i in seq_along(meta_genes)) {
      gene_i <- meta_genes[i]
      fit <- fit_one_gene(selected[J(gene_i)])
      if (!is.null(fit)) {
        fit[, gene := gene_i]
        pieces[[i]] <- fit
      }
      if (i %% 5000L == 0L) {
        message(format(Sys.time(), "%H:%M:%S"), " | ", model_name, " | fitted ", i, " / ", length(meta_genes), " genes")
      }
    }
    res <- rbindlist(pieces, fill = TRUE)
  }

  if (nrow(res)) {
    res[, FDR := p.adjust(p_value, method = "BH")]
    res[, DL_nonzero := (DL_CI_lower > 0 | DL_CI_upper < 0)]
    res[, mKH_interval_supported := (mKH_CI_lower > 0 | mKH_CI_upper < 0)]
    res[, FDR_significant := FDR < 0.05]
    res[, FDR_mKH_overlap := FDR_significant & mKH_interval_supported]
    res[, `:=`(
      omic_layer = omic_layer,
      tissue_branch = tissue_branch,
      model_name = model_name,
      model_role = model_role,
      sensitivity_family = sensitivity_family,
      datasets_included = paste(datasets, collapse = ";"),
      notes = notes
    )]
    setcolorder(res, c("omic_layer", "tissue_branch", "model_name", "model_role", "sensitivity_family", "gene"))
  }

  summary <- data.table(
    omic_layer = omic_layer,
    tissue_branch = tissue_branch,
    model_name = model_name,
    model_role = model_role,
    sensitivity_family = sensitivity_family,
    datasets_included = paste(datasets, collapse = ";"),
    genes_with_any_result = uniqueN(selected$gene),
    genes_meta_analysed = nrow(res),
    k1_descriptive = length(k1_genes),
    DL_nonzero = if (nrow(res)) sum(res$DL_nonzero, na.rm = TRUE) else 0L,
    FDR_significant = if (nrow(res)) sum(res$FDR_significant, na.rm = TRUE) else 0L,
    mKH_interval_supported = if (nrow(res)) sum(res$mKH_interval_supported, na.rm = TRUE) else 0L,
    FDR_mKH_overlap = if (nrow(res)) sum(res$FDR_mKH_overlap, na.rm = TRUE) else 0L,
    median_k = if (nrow(res)) median(res$k, na.rm = TRUE) else NA_real_,
    median_I2 = if (nrow(res)) median(res$I2, na.rm = TRUE) else NA_real_,
    IQR_I2_lower = if (nrow(res)) q_safe(res$I2, 0.25) else NA_real_,
    IQR_I2_upper = if (nrow(res)) q_safe(res$I2, 0.75) else NA_real_,
    I2_gt50 = if (nrow(res)) sum(res$I2 > 50, na.rm = TRUE) else 0L,
    I2_gt75 = if (nrow(res)) sum(res$I2 > 75, na.rm = TRUE) else 0L,
    median_tau2 = if (nrow(res)) median(res$tau2, na.rm = TRUE) else NA_real_,
    max_k = if (nrow(res)) max(res$k, na.rm = TRUE) else 0L,
    max_total_ASD_n = if (nrow(res)) max(res$total_ASD_n, na.rm = TRUE) else NA_real_,
    max_total_control_n = if (nrow(res)) max(res$total_control_n, na.rm = TRUE) else NA_real_,
    notes = notes
  )

  list(results = res, k1 = k1, summary = summary)
}

write_model_outputs <- function(model_specs, effects_by_input, output_dir, prefix) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  by_model_dir <- file.path(output_dir, "by_model")
  by_model_k1_dir <- file.path(output_dir, "by_model_k1")
  by_model_summary_dir <- file.path(output_dir, "by_model_summary")
  dir.create(by_model_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(by_model_k1_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(by_model_summary_dir, recursive = TRUE, showWarnings = FALSE)

  all_results <- list()
  all_k1 <- list()
  all_summary <- list()

  for (i in seq_along(model_specs)) {
    spec <- model_specs[[i]]
    message("\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", spec$model_name)
    ans <- run_meta_model(
      effects = effects_by_input[[spec$input_id]],
      model_name = spec$model_name,
      omic_layer = spec$omic_layer,
      tissue_branch = spec$tissue_branch,
      model_role = spec$model_role,
      datasets = spec$datasets,
      notes = spec$notes,
      sensitivity_family = spec$sensitivity_family
    )

    file_stub <- paste0(sprintf("%02d_", i), safe_name(spec$model_name))
    fwrite(ans$results, file.path(by_model_dir, paste0(file_stub, "_meta_results.csv")))
    fwrite(ans$k1, file.path(by_model_k1_dir, paste0(file_stub, "_k1_descriptive_rows.csv")))
    fwrite(ans$summary, file.path(by_model_summary_dir, paste0(file_stub, "_summary.csv")))

    all_results[[i]] <- ans$results
    all_k1[[i]] <- ans$k1
    all_summary[[i]] <- ans$summary
  }

  combined_results <- rbindlist(all_results, fill = TRUE)
  combined_k1 <- rbindlist(all_k1, fill = TRUE)
  combined_summary <- rbindlist(all_summary, fill = TRUE)

  fwrite(combined_results, file.path(output_dir, paste0(prefix, "_gene_level_meta_results.csv")))
  fwrite(combined_k1, file.path(output_dir, paste0(prefix, "_k1_descriptive_rows.csv")))
  fwrite(combined_summary, file.path(output_dir, paste0(prefix, "_model_summary.csv")))

  wb <- createWorkbook()
  addWorksheet(wb, "model_summary")
  writeData(wb, "model_summary", combined_summary)
  saveWorkbook(wb, file.path(output_dir, paste0(prefix, "_model_summary.xlsx")), overwrite = TRUE)

  invisible(list(results = combined_results, k1 = combined_k1, summary = combined_summary))
}

