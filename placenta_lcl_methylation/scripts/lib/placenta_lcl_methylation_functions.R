suppressPackageStartupMessages({
  library(metafor)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

clean_quotes <- function(x) {
  x <- gsub('^"', "", x)
  x <- gsub('"$', "", x)
  x
}

normalise_public_url <- function(url) {
  ifelse(is.na(url) | !nzchar(url), url, sub("^ftp://", "https://", url))
}

download_or_copy <- function(url = NA_character_, local_fallback = NA_character_, dest,
                             overwrite = FALSE) {
  dir_create(dirname(dest))
  if (file.exists(dest) && !overwrite) {
    return(data.table::data.table(dest = dest, status = "already_present",
                                  source_used = dest, bytes = file.info(dest)$size))
  }
  if (!is.na(local_fallback) && file.exists(local_fallback)) {
    ok <- file.copy(local_fallback, dest, overwrite = TRUE)
    status <- if (ok) "copied_from_user_supplied_source_archive" else "copy_failed"
    source_used <- local_fallback
  } else if (!is.na(url) && nzchar(url)) {
    tmp <- paste0(dest, ".download_tmp")
    if (file.exists(tmp)) unlink(tmp)
    ok <- tryCatch({
      utils::download.file(normalise_public_url(url), tmp, mode = "wb", quiet = TRUE)
      file.rename(tmp, dest)
    }, error = function(e) FALSE)
    if (!ok && file.exists(tmp)) unlink(tmp)
    status <- if (ok) "downloaded" else "download_failed"
    source_used <- normalise_public_url(url)
  } else {
    status <- "missing"
    source_used <- NA_character_
  }
  data.table::data.table(dest = dest, status = status, source_used = source_used,
                         bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_)
}

read_geo_sample_metadata <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  table_start <- which(lines == "!series_matrix_table_begin")[1]
  if (is.finite(table_start)) lines <- lines[seq_len(table_start - 1L)]
  sample_lines <- lines[grepl("^!Sample_", lines)]
  rows <- lapply(sample_lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    data.table::data.table(field = sub("^!Sample_", "", parts[1]),
                           sample_index = seq_along(parts[-1]),
                           value = clean_quotes(parts[-1]))
  })
  long <- data.table::rbindlist(rows, fill = TRUE)
  data.table::dcast(long, sample_index ~ field, value.var = "value",
                    fun.aggregate = function(z) paste(z, collapse = " | "))
}

read_geo_matrix_table <- function(path) {
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
  dt
}

normalise_numeric_matrix <- function(dt, sample_cols) {
  for (col in sample_cols) {
    raw <- dt[[col]]
    if (is.character(raw)) {
      raw[trimws(raw) == ""] <- NA_character_
      suppressWarnings(dt[[col]] <- as.numeric(raw))
    } else {
      dt[[col]] <- as.numeric(raw)
    }
  }
  dt
}

split_gene_symbols <- function(x) {
  out <- unlist(strsplit(x, "\\s*(;|,|///|//|\\|)\\s*"))
  out <- trimws(out)
  out[nzchar(out)]
}

load_refgene_promoters <- function(path, build) {
  dt <- data.table::fread(path, header = FALSE, showProgress = FALSE)
  if (ncol(dt) >= 16) {
    data.table::setnames(dt, seq_len(16),
                         c("bin", "name", "chr", "strand", "txStart", "txEnd",
                           "cdsStart", "cdsEnd", "exonCount", "exonStarts",
                           "exonEnds", "score", "name2", "cdsStartStat",
                           "cdsEndStat", "exonFrames"))
    dt[, gene := name2]
  } else {
    data.table::setnames(dt, seq_len(13),
                         c("bin", "name", "chr", "strand", "txStart", "txEnd",
                           "cdsStart", "cdsEnd", "exonCount", "exonStarts",
                           "exonEnds", "score", "gene"))
  }
  dt <- dt[grepl("^chr([0-9]+|X|Y)$", chr) & nzchar(gene)]
  dt[, tss := ifelse(strand == "+", txStart, txEnd)]
  dt[, start := ifelse(strand == "+", pmax(0L, tss - 1500L), pmax(0L, tss - 200L))]
  dt[, end := ifelse(strand == "+", tss + 200L, tss + 1500L)]
  dt[, `:=`(build = build, start = as.integer(start), end = as.integer(end))]
  unique(dt[, .(gene, build, transcript = name, chr, strand, tss, start, end,
                rule = "strand-aware promoter: upstream 1500 bp/downstream 200 bp")])
}

write_csv <- function(x, path) {
  dir_create(dirname(path))
  data.table::fwrite(x, path)
  invisible(path)
}

ci_excludes_zero <- function(lo, hi) {
  is.finite(lo) & is.finite(hi) & ((lo > 0 & hi > 0) | (lo < 0 & hi < 0))
}

hedges_effects <- function(summary_dt, value_label = "methylation") {
  dt <- data.table::copy(summary_dt)
  dt <- dt[finite_summary == "yes"]
  dt[, pooled_sd := sqrt(((ASD_n - 1) * ASD_SD^2 + (control_n - 1) * control_SD^2) /
                           (ASD_n + control_n - 2))]
  dt <- dt[is.finite(pooled_sd) & pooled_sd > 0]
  dt[, cohen_d := (ASD_mean - control_mean) / pooled_sd]
  dt[, J := 1 - (3 / (4 * (ASD_n + control_n) - 9))]
  dt[, Hedges_g := J * cohen_d]
  dt[, variance_g := ((ASD_n + control_n) / (ASD_n * control_n)) +
       (Hedges_g^2 / (2 * (ASD_n + control_n - 2)))]
  dt[, SE_g := sqrt(variance_g)]
  dt[, CI_lower := Hedges_g - 1.96 * SE_g]
  dt[, CI_upper := Hedges_g + 1.96 * SE_g]
  dt[, signal_type := value_label]
  dt[, .(accession, tissue_family, assay_class, platform, gene, ASD_n, control_n,
         ASD_mean, ASD_SD, control_mean, control_SD, Hedges_g, variance_g, SE_g,
         CI_lower, CI_upper, feature_count, signal_type, notes)]
}

random_effects_dl <- function(rows) {
  rows <- rows[is.finite(Hedges_g) & is.finite(variance_g) & variance_g > 0]
  k <- nrow(rows)
  if (!k) return(NULL)
  yi <- rows$Hedges_g
  vi <- rows$variance_g
  if (k == 1) {
    se <- sqrt(vi[1])
    p <- 2 * stats::pnorm(abs(yi[1] / se), lower.tail = FALSE)
    return(list(k = 1, pooled_g = yi[1], SE = se,
                CI_lower = yi[1] - 1.96 * se, CI_upper = yi[1] + 1.96 * se,
                p_value = p, Q = 0, Q_p_value = NA_real_, tau2 = 0, I2 = 0,
                hk = NULL))
  }
  dl <- metafor::rma.uni(yi = yi, vi = vi, method = "DL", test = "z")
  mkh <- metafor::rma.uni(yi = yi, vi = vi, method = "DL", test = "adhoc")
  list(k = k,
       pooled_g = as.numeric(dl$b[1]),
       SE = as.numeric(dl$se),
       CI_lower = as.numeric(dl$ci.lb),
       CI_upper = as.numeric(dl$ci.ub),
       p_value = as.numeric(dl$pval),
       Q = as.numeric(dl$QE),
       Q_p_value = as.numeric(dl$QEp),
       tau2 = as.numeric(dl$tau2),
       I2 = as.numeric(dl$I2),
       hk = list(mKH_CI_lower = as.numeric(mkh$ci.lb),
                 mKH_CI_upper = as.numeric(mkh$ci.ub)))
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
        out[[idx]] <- data.table::data.table(model_name = model_id,
                                             role = role_row$role %||% NA_character_,
                                             gene = gene_id, k = 0L)
      } else {
        hk <- re$hk
        out[[idx]] <- data.table::data.table(
          model_name = model_id,
          role = role_row$role %||% NA_character_,
          gene = gene_id,
          k = re$k,
          contributing_datasets = paste(rows$accession, collapse = ";"),
          contributing_assays = paste(unique(rows$assay_class), collapse = ";"),
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
