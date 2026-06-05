suppressPackageStartupMessages({
  library(metafor)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

clean_quotes <- function(x) {
  x <- gsub('^"', "", x)
  x <- gsub('"$', "", x)
  x
}

dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_copy <- function(from, to, overwrite = FALSE) {
  dir_create(dirname(to))
  if (file.exists(to) && !overwrite) return(TRUE)
  if (!file.exists(from)) return(FALSE)
  file.copy(from, to, overwrite = overwrite)
}

download_or_copy <- function(url = NA_character_, local_fallback = NA_character_, dest, overwrite = FALSE) {
  dir_create(dirname(dest))
  status <- "missing"
  source_used <- NA_character_
  if (file.exists(dest) && !overwrite) {
    return(data.table::data.table(dest = dest, status = "already_present", source_used = dest, bytes = file.info(dest)$size))
  }
  if (!is.na(local_fallback) && file.exists(local_fallback)) {
    ok <- file.copy(local_fallback, dest, overwrite = TRUE)
    status <- if (ok) "copied_from_public_source_cache" else "copy_failed"
    source_used <- local_fallback
  } else if (!is.na(url) && nzchar(url)) {
    ok <- tryCatch({
      utils::download.file(url, dest, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    status <- if (ok) "downloaded" else "download_failed"
    source_used <- url
  }
  data.table::data.table(dest = dest, status = status, source_used = source_used,
                         bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_)
}

read_geo_sample_metadata <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  lines <- lines[seq_len(which(lines == "!series_matrix_table_begin")[1] - 1L)]
  sample_lines <- lines[grepl("^!Sample_", lines)]
  rows <- lapply(sample_lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    data.table::data.table(field = sub("^!Sample_", "", parts[1]), sample_index = seq_along(parts[-1]),
                           value = clean_quotes(parts[-1]))
  })
  long <- data.table::rbindlist(rows, fill = TRUE)
  wide <- data.table::dcast(long, sample_index ~ field, value.var = "value", fun.aggregate = function(z) paste(z, collapse = " | "))
  data.table::setnames(wide, old = intersect(names(wide), c("geo_accession", "title", "source_name_ch1", "description")),
                       new = intersect(names(wide), c("geo_accession", "title", "source_name_ch1", "description")))
  wide
}

read_geo_matrix_table <- function(path) {
  # GEO series matrix files are not completely consistent about whether fread()
  # treats the table-boundary marker or the following ID_REF row as the header.
  # Extracting the table block explicitly keeps the reader deterministic across
  # the brain array/HM27 routes.
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  start <- which(lines == "!series_matrix_table_begin")[1]
  end <- which(lines == "!series_matrix_table_end")[1]
  if (!is.finite(start) || !is.finite(end) || end <= start + 1L) {
    stop("Could not identify GEO series matrix table block in: ", path)
  }
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines[(start + 1L):(end - 1L)], tmp, useBytes = TRUE)
  dt <- data.table::fread(tmp, showProgress = FALSE)
  if (!"ID_REF" %in% names(dt)) data.table::setnames(dt, names(dt)[1], "ID_REF")
  dt <- dt[grepl("^cg", ID_REF)]
  dt
}

normalise_numeric_matrix <- function(dt, sample_cols, blanks_as_zero = FALSE) {
  for (col in sample_cols) {
    raw <- dt[[col]]
    if (is.character(raw)) {
      if (!blanks_as_zero) raw[trimws(raw) == ""] <- NA_character_
      suppressWarnings(dt[[col]] <- as.numeric(raw))
    } else {
      dt[[col]] <- as.numeric(raw)
    }
  }
  dt
}

extract_after_colon <- function(x, pattern) {
  hit <- grep(pattern, x, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub("^[^:]+:\\s*", "", hit[1]))
}

infer_group <- function(text) {
  text <- tolower(paste(text, collapse = " "))
  if (grepl("autism spectrum disorder|asd/autism|\\basd\\b|autism|idioaut", text)) return("ASD")
  if (grepl("control|disease state:\\s*normal|\\bnormal\\b|unaffected|neurotypical|typically developing|healthy|no diagnosed neurological disorders", text)) return("Control")
  NA_character_
}

infer_region_group <- function(region) {
  r <- tolower(region %||% "")
  if (grepl("cortex|ba9|ba10|ba19|prefrontal|temporal|brodmann", r)) return("cortex")
  if (grepl("cerebell", r)) return("cerebellum")
  if (grepl("subventricular", r)) return("subventricular zone")
  if (grepl("raphe", r)) return("dorsal raphe")
  if (grepl("mixed|brain", r)) return("mixed brain")
  "other"
}

row_sds <- function(mat) {
  n <- rowSums(is.finite(mat))
  s <- rowSums(mat, na.rm = TRUE)
  ss <- rowSums(mat^2, na.rm = TRUE)
  out <- sqrt(pmax(0, (ss - (s^2 / n)) / pmax(1, n - 1)))
  out[n < 2] <- NA_real_
  out
}

build_gene_matrix <- function(beta_dt, probe_map, sample_cols) {
  probes <- beta_dt$ID_REF
  keep <- probe_map[probe %in% probes]
  keep <- keep[!duplicated(paste(probe, gene))]
  idx <- match(keep$probe, probes)
  split_idx <- split(idx, keep$gene)
  feature_counts <- data.table::data.table(gene = names(split_idx), feature_probe_count = lengths(split_idx))
  mat <- as.matrix(beta_dt[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  gene_values <- lapply(names(split_idx), function(gene) {
    rows <- split_idx[[gene]]
    vals <- colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
    vals[is.nan(vals)] <- NA_real_
    data.table::data.table(gene = gene, t(vals))
  })
  gene_dt <- data.table::rbindlist(gene_values, fill = TRUE)
  names(gene_dt)[-1] <- sample_cols
  list(gene_matrix = gene_dt, feature_counts = feature_counts)
}

summarise_gene_matrix <- function(dataset_id, gene_dt, sample_meta, feature_counts,
                                  brain_region, broader_region_group, platform, assay_class,
                                  collapse_subject = FALSE) {
  long <- data.table::melt(gene_dt, id.vars = "gene", variable.name = "sample_id", value.name = "value")
  sample_meta <- data.table::as.data.table(sample_meta)
  long <- merge(long, sample_meta[, .(sample_id, group, subject_id)], by = "sample_id", all.x = TRUE)
  if (collapse_subject) {
    long <- long[!is.na(group) & is.finite(value), .(value = mean(value, na.rm = TRUE)),
                 by = .(gene, group, subject_id)]
    group_id <- "subject_id"
  } else {
    long <- long[!is.na(group) & is.finite(value)]
    group_id <- "sample_id"
  }
  asd <- long[group == "ASD", .(
    ASD_n = data.table::uniqueN(get(group_id)),
    ASD_mean_methylation = mean(value),
    ASD_SD = stats::sd(value)
  ), by = gene]
  ctl <- long[group == "Control", .(
    control_n = data.table::uniqueN(get(group_id)),
    control_mean_methylation = mean(value),
    control_SD = stats::sd(value)
  ), by = gene]
  out <- merge(asd, ctl, by = "gene", all = TRUE)
  out <- merge(out, feature_counts, by = "gene", all.x = TRUE)
  out[, `:=`(
    accession = dataset_id,
    brain_region = brain_region,
    broader_region_group = broader_region_group,
    platform = platform,
    assay_class = assay_class,
    mean_difference = ASD_mean_methylation - control_mean_methylation,
    finite_summary = fifelse(is.finite(ASD_mean_methylation) & is.finite(control_mean_methylation) &
                               is.finite(ASD_SD) & is.finite(control_SD) &
                               ASD_n >= 2 & control_n >= 2, "yes", "no"),
    missingness_reason = fifelse(is.na(ASD_n) | is.na(control_n), "missing group", ""),
    phenotype_labels_used = "ASD vs Control",
    promoter_definition_used = "Illumina annotation promoter groups: TSS200, TSS1500, 5UTR, 1stExon",
    aggregation_method = if (collapse_subject) {
      "Mean promoter feature value per sample; collapsed technical/multi-region measurements to independent subject-level values before group summary."
    } else {
      "Mean promoter feature value per sample; no subject collapse required."
    },
    notes = if (collapse_subject) {
      "Autism/control labels and brain region parsed from title/source; technical and multi-region rows collapsed to independent subject values."
    } else {
      "Autism/control labels parsed from repository metadata."
    }
  )]
  data.table::setcolorder(out, c("accession", "brain_region", "broader_region_group", "platform", "assay_class", "gene"))
  out[]
}

hedges_effects <- function(summary_dt) {
  dt <- data.table::copy(summary_dt)
  dt <- dt[finite_summary == "yes"]
  dt[, `:=`(
    ASD_mean = ASD_mean_methylation,
    control_mean = control_mean_methylation
  )]
  dt[, pooled_sd := sqrt(((ASD_n - 1) * ASD_SD^2 + (control_n - 1) * control_SD^2) / (ASD_n + control_n - 2))]
  dt <- dt[is.finite(pooled_sd) & pooled_sd > 0]
  dt[, cohen_d := (ASD_mean - control_mean) / pooled_sd]
  dt[, J := 1 - (3 / (4 * (ASD_n + control_n) - 9))]
  dt[, Hedges_g := J * cohen_d]
  dt[, variance_g := ((ASD_n + control_n) / (ASD_n * control_n)) + (Hedges_g^2 / (2 * (ASD_n + control_n - 2)))]
  dt[, SE_g := sqrt(variance_g)]
  dt[, CI_lower := Hedges_g - 1.96 * SE_g]
  dt[, CI_upper := Hedges_g + 1.96 * SE_g]
  dt[, feature_count := feature_probe_count]
  dt[, .(accession, brain_region, broader_region_group, platform, assay_class, gene, ASD_n, control_n,
         ASD_mean, ASD_SD, control_mean, control_SD, Hedges_g, variance_g, SE_g, CI_lower, CI_upper,
         feature_count, notes)]
}

ci_excludes_zero <- function(lo, hi) is.finite(lo) & is.finite(hi) & ((lo > 0 & hi > 0) | (lo < 0 & hi < 0))

random_effects_dl <- function(rows) {
  rows <- rows[is.finite(Hedges_g) & is.finite(variance_g) & variance_g > 0]
  k <- nrow(rows)
  if (!k) return(NULL)
  yi <- rows$Hedges_g
  vi <- rows$variance_g
  if (k == 1) {
    se <- sqrt(vi[1])
    p <- 2 * stats::pnorm(abs(yi[1] / se), lower.tail = FALSE)
    return(list(k = 1, pooled_g = yi[1], SE = se, CI_lower = yi[1] - 1.96 * se,
                CI_upper = yi[1] + 1.96 * se, p_value = p, Q = 0,
                Q_p_value = NA_real_, tau2 = 0, I2 = 0, hk = NULL))
  }
  dl <- metafor::rma.uni(yi = yi, vi = vi, method = "DL", test = "z")
  mkh <- metafor::rma.uni(yi = yi, vi = vi, method = "DL", test = "adhoc")
  list(
    k = k,
    pooled_g = as.numeric(dl$b[1]),
    SE = as.numeric(dl$se),
    CI_lower = as.numeric(dl$ci.lb),
    CI_upper = as.numeric(dl$ci.ub),
    p_value = as.numeric(dl$pval),
    Q = as.numeric(dl$QE),
    Q_p_value = as.numeric(dl$QEp),
    tau2 = as.numeric(dl$tau2),
    I2 = as.numeric(dl$I2),
    hk = list(
      mKH_CI_lower = as.numeric(mkh$ci.lb),
      mKH_CI_upper = as.numeric(mkh$ci.ub)
    )
  )
}

build_meta_results <- function(effects, universe, model_defs, model_roles) {
  out <- vector("list", length(model_defs) * length(universe))
  idx <- 1L
  for (model_id in names(model_defs)) {
    include <- model_defs[[model_id]]
    role_row <- model_roles[model_name == model_id]
    for (gene_id in universe) {
      rows <- effects[gene == gene_id & accession %in% include]
      re <- random_effects_dl(rows)
      if (is.null(re)) {
        out[[idx]] <- data.table::data.table(model_name = model_id, gene = gene_id, k = 0L)
      } else {
        hk <- re$hk
        datasets <- paste(rows$accession, collapse = ";")
        out[[idx]] <- data.table::data.table(
          model_name = model_id,
          role = role_row$role %||% NA_character_,
          gene = gene_id,
          k = re$k,
          contributing_datasets = datasets,
          pooled_Hedges_g = re$pooled_g,
          SE = re$SE,
          DL_CI_lower = re$CI_lower,
          DL_CI_upper = re$CI_upper,
          p_value = re$p_value,
          Q = re$Q,
          Q_p_value = re$Q_p_value,
          tau2 = re$tau2,
          I2 = re$I2,
          mKH_CI_lower = if (is.null(hk)) NA_real_ else hk$mKH_CI_lower,
          mKH_CI_upper = if (is.null(hk)) NA_real_ else hk$mKH_CI_upper,
          DL_CI_excludes_zero = ci_excludes_zero(re$CI_lower, re$CI_upper),
          mKH_CI_excludes_zero = if (is.null(hk)) FALSE else ci_excludes_zero(hk$mKH_CI_lower, hk$mKH_CI_upper),
          total_ASD_n = sum(rows$ASD_n, na.rm = TRUE),
          total_control_n = sum(rows$control_n, na.rm = TRUE)
        )
      }
      idx <- idx + 1L
    }
  }
  res <- data.table::rbindlist(out, fill = TRUE)
  res[, FDR := NA_real_]
  res[k >= 2 & is.finite(p_value), FDR := stats::p.adjust(p_value, method = "BH"), by = model_name]
  res[, FDR_significant := is.finite(FDR) & FDR < 0.05]
  res[]
}

summarise_models <- function(meta_results, model_roles) {
  ms <- meta_results[, .(
    genes_with_any_result = sum(k >= 1, na.rm = TRUE),
    genes_meta_analysed_k_ge_2 = sum(k >= 2, na.rm = TRUE),
    k1_descriptive_genes = sum(k == 1, na.rm = TRUE),
    DL_nonzero_genes = sum(DL_CI_excludes_zero %in% TRUE, na.rm = TRUE),
    FDR_significant_genes = sum(FDR_significant %in% TRUE, na.rm = TRUE),
    mKH_retained_genes = sum(mKH_CI_excludes_zero %in% TRUE, na.rm = TRUE),
    FDR_and_mKH_overlap = sum(FDR_significant %in% TRUE & mKH_CI_excludes_zero %in% TRUE, na.rm = TRUE),
    max_k = max(k, na.rm = TRUE),
    max_total_ASD_n = max(total_ASD_n, na.rm = TRUE),
    max_total_control_n = max(total_control_n, na.rm = TRUE)
  ), by = model_name]
  merge(model_roles, ms, by = "model_name", all.x = TRUE)
}
