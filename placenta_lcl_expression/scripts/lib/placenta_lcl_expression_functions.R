suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(huex10sttranscriptcluster.db)
  library(metafor)
  library(openxlsx)
})

plcl_expression_config <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(root) == "scripts") root <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = TRUE)
  if (basename(root) == "placenta_lcl_expression") {
    out_dir <- root
    repo_root <- normalizePath(file.path(out_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  } else {
    repo_root <- root
    out_dir <- file.path(repo_root, "pipelines", "placenta_lcl_expression")
  }
  source_cache <- Sys.getenv("PLACENTA_LCL_EXPRESSION_PUBLIC_SOURCE_CACHE", unset = "")
  list(
    repo_root = repo_root,
    out_dir = out_dir,
    scripts_dir = file.path(out_dir, "scripts"),
    source_cache = source_cache,
    validation_reference_dir = Sys.getenv("PLACENTA_LCL_EXPRESSION_VALIDATION_REFERENCE_DIR", unset = ""),
    subdirs = c("00_manifest", "01_source_files", "02_sample_metadata", "03_processed_expression",
                "04_dataset_gene_summaries", "05_effect_sizes", "06_lcl_meta_analysis",
                "07_placenta_meta_analysis", "08_sensitivity_analyses", "09_quality_control")
  )
}

write_csv_safe <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.frame(x), file)
}

rbindlist_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  x <- x[vapply(x, nrow, integer(1)) > 0]
  if (!length(x)) return(data.frame())
  data.table::rbindlist(x, fill = TRUE)
}

geo_series_prefix <- function(gse) paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
geo_platform_prefix <- function(gpl) paste0(substr(gpl, 1, nchar(gpl) - 3), "nnn")
geo_matrix_url <- function(gse) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz", geo_series_prefix(gse), gse, gse)
geo_platform_url <- function(gpl) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz", geo_platform_prefix(gpl), gpl, gpl)

stage_one_file <- function(source_cache, dest_dir, filename, url) {
  dest <- file.path(dest_dir, filename)
  cache_file <- file.path(source_cache, filename)
  status <- "present"
  source_used <- "destination already present"
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
      file.copy(cache_file, dest, overwrite = TRUE)
      status <- "copied_from_public_source_cache"
      source_used <- "local public source cache"
    } else {
      old_timeout <- getOption("timeout")
      options(timeout = max(1200, old_timeout))
      on.exit(options(timeout = old_timeout), add = TRUE)
      ok <- try(download.file(url, dest, mode = "wb", quiet = TRUE), silent = TRUE)
      if (inherits(ok, "try-error")) {
        status <- "not_available_at_standard_geo_url"
        source_used <- "not staged"
        if (file.exists(dest) && file.info(dest)$size == 0) unlink(dest)
      } else {
        status <- "downloaded"
        source_used <- "downloaded from GEO"
      }
    }
  }
  data.frame(file_name = filename, local_path = normalizePath(dest, winslash = "/", mustWork = FALSE),
             url = url, source_used = source_used, status = status,
             bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_, stringsAsFactors = FALSE)
}

stage_plcl_sources <- function(cfg) {
  dest_dir <- file.path(cfg$out_dir, "01_source_files")
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  datasets <- c("GSE15402", "GSE15451", "GSE29918", "GSE37772", "GSE4187", "GSE7329",
                "GSE154829", "GSE178205", "GSE285666", "GSE57802")
  platforms <- c("GPL3427", "GPL6883", "GPL1708", "GPL13158")
  files <- lapply(datasets, function(ds) c(ds, paste0(ds, "_series_matrix.txt.gz"), geo_matrix_url(ds), "series_matrix"))
  files <- c(files, lapply(platforms, function(gpl) c(gpl, paste0(gpl, ".annot.gz"), geo_platform_url(gpl), "platform_annotation")))
  out <- lapply(files, function(z) {
    s <- stage_one_file(cfg$source_cache, dest_dir, z[[2]], z[[3]])
    s$dataset_id <- z[[1]]
    s$route_type <- z[[4]]
    s
  })
  rbindlist_fill(out)
}

read_series_matrix <- function(file) {
  lines <- readLines(gzfile(file, "rt"), warn = FALSE)
  begin <- match("!series_matrix_table_begin", lines)
  end <- match("!series_matrix_table_end", lines)
  if (is.na(begin) || is.na(end)) stop("No series matrix table found in ", file)
  sample_meta <- list()
  for (ln in lines[seq_len(begin - 1)][grepl("^!Sample_", lines[seq_len(begin - 1)])]) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    sample_meta[[sub("^!", "", parts[1])]] <- gsub('^"|"$', "", parts[-1])
  }
  dt <- data.table::fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                          data.table = FALSE, check.names = FALSE, na.strings = c("NA", "null", "NULL", ""))
  names(dt)[1] <- "probe_id"
  list(values = dt, sample_meta = sample_meta)
}

sample_meta_frame <- function(sm) {
  n <- length(sm$Sample_geo_accession)
  field <- function(k) if (!is.null(sm[[k]])) sm[[k]] else rep("", n)
  ch <- field("Sample_characteristics_ch1")
  char <- if (is.matrix(ch)) apply(ch, 2, paste, collapse = "; ") else ch
  if (length(char) != n) char <- rep(paste(ch, collapse = "; "), n)
  data.frame(sample = field("Sample_geo_accession"), title = field("Sample_title"),
             source = field("Sample_source_name_ch1"), characteristics = char,
             platform = field("Sample_platform_id"), stringsAsFactors = FALSE)
}

split_gene_field <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  parts <- unlist(strsplit(x, "///|//|;|,|\\|", perl = TRUE), use.names = FALSE)
  parts <- trimws(gsub("\\s*\\(.*?\\)", "", parts))
  parts <- parts[nzchar(parts)]
  parts <- parts[!toupper(parts) %in% c("NA", "N/A", "---", "--", "NULL", "NONE")]
  unique(parts[grepl("^[A-Za-z0-9_.-]+$", parts)])
}

read_annot_platform_map <- function(file) {
  cache_dir <- file.path(dirname(dirname(file)), "03_processed_expression", "platform_maps")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(cache_dir, paste0(basename(file), ".platform_map.csv.gz"))
  if (file.exists(cache) && file.info(cache)$size > 0) return(data.table::fread(cache, data.table = FALSE))
  lines <- readLines(gzfile(file, "rt"), warn = FALSE)
  header <- grep("^(ID|\"ID\")\t", lines)[1]
  fields <- strsplit(lines[header], "\t", fixed = TRUE)[[1]]
  id_i <- which(tolower(fields) %in% c("id", "id_ref"))[1]
  gene_i <- which(tolower(fields) %in% c("gene symbol", "gene_symbol", "symbol"))[1]
  rows <- lapply(lines[(header + 1):length(lines)], function(ln) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < max(id_i, gene_i, na.rm = TRUE)) return(NULL)
    genes <- split_gene_field(parts[gene_i])
    if (!length(genes)) return(NULL)
    data.frame(probe_id = parts[id_i], gene = genes, stringsAsFactors = FALSE)
  })
  out <- unique(rbindlist_fill(rows))
  write_csv_safe(out, cache)
  out
}

read_platform_map <- function(cfg, platform) {
  if (platform == "GPL5175") {
    keys <- AnnotationDbi::keys(huex10sttranscriptcluster.db, keytype = "PROBEID")
    syms <- AnnotationDbi::mapIds(huex10sttranscriptcluster.db, keys = keys, keytype = "PROBEID",
                                  column = "SYMBOL", multiVals = "first")
    out <- data.frame(probe_id = keys, gene = as.character(syms), stringsAsFactors = FALSE)
    return(out[!is.na(out$gene) & nzchar(out$gene), ])
  }
  read_annot_platform_map(file.path(cfg$out_dir, "01_source_files", paste0(platform, ".annot.gz")))
}

label_lcl_sample <- function(dataset, title, source, characteristics) {
  txt <- tolower(paste(title, source, characteristics))
  if (dataset == "GSE15451") {
    if (grepl("^control", title, ignore.case = TRUE)) return(c("control", "Control title label"))
    if (grepl("^autistic", title, ignore.case = TRUE)) return(c("ASD", "Autistic title label"))
    return(c("exclude", "unclear GSE15451 title label"))
  }
  if (dataset == "GSE4187") {
    s <- tolower(source)
    if (grepl("nonautistic", s)) return(c("control", "nonautistic source label"))
    if (grepl("autistic|mildly autistic", s)) return(c("ASD", "autistic source label"))
    return(c("exclude", "unclear GSE4187 source label"))
  }
  if (dataset == "GSE29918") {
    if (grepl("^Control", title, ignore.case = TRUE) || grepl("healthy", txt)) return(c("control", "Control/healthy label"))
    if (grepl("^Index", title, ignore.case = TRUE) || grepl("autistic features|affected", txt)) return(c("ASD", "Index/ASD-features label"))
  }
  if (dataset == "GSE7329") {
    if (grepl("autism", txt)) return(c("ASD", "pooled autism/syndromic autism label"))
    if (grepl("control|normal", txt)) return(c("control", "pooled control label"))
  }
  if (grepl("nonautistic|control|normal", txt)) return(c("control", "public control label"))
  if (grepl("autistic|autism|asd", txt)) return(c("ASD", "public ASD/autism label"))
  c("exclude", "unclear ASD/control label")
}

bio_id_for_sample <- function(dataset, sample, title, source) {
  if (dataset == "GSE29918") {
    m <- regmatches(title, regexpr("^(Control|Index)[0-9]+", title, ignore.case = TRUE))
    if (length(m) && nzchar(m)) return(m)
    return(sample)
  }
  if (dataset == "GSE4187") {
    m <- regmatches(source, regexpr("HI[0-9]+", source, ignore.case = TRUE))
    if (length(m) && nzchar(m)) return(toupper(m))
    m <- regmatches(title, regexpr("HI[0-9]+", title, ignore.case = TRUE))
    if (length(m) && nzchar(m)) return(toupper(m))
  }
  sample
}

role_for_dataset <- function(dataset) {
  if (dataset %in% c("GSE15402", "GSE15451", "GSE37772", "GSE4187")) return("lcl_core_primary_candidate")
  if (dataset == "GSE29918") return("lcl_asd_features_sensitivity_candidate")
  if (dataset == "GSE7329") return("lcl_pooled_syndromic_sensitivity_candidate")
  "excluded_or_context"
}

caveat_for_dataset <- function(dataset) {
  if (dataset == "GSE29918") return("Global developmental delay with autistic features; expanded sensitivity only.")
  if (dataset == "GSE7329") return("Pooled syndromic LCL route; expanded sensitivity only.")
  if (dataset == "GSE15451") return("Sibling/family LCL comparison; family non-independence caveat.")
  if (dataset == "GSE4187") return("Twin/sibling LCL design; technical replicate columns collapsed by HI identifier.")
  ""
}

collapse_gene_matrix <- function(mat) {
  rows <- split(seq_len(nrow(mat)), rownames(mat))
  out <- vapply(rows, function(idx) {
    if (length(idx) == 1) as.numeric(mat[idx, ]) else colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(mat)))
  out <- t(out)
  colnames(out) <- colnames(mat)
  out
}

collapse_replicates <- function(mat, labels) {
  labels <- labels[labels$group %in% c("ASD", "control") & labels$sample %in% colnames(mat), ]
  labels$key <- paste(labels$bio_id, labels$group, sep = "||")
  groups <- split(labels$sample, labels$key)
  out <- vapply(groups, function(cols) {
    if (length(cols) == 1) as.numeric(mat[, cols]) else rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(mat)))
  out <- as.matrix(out)
  rownames(out) <- rownames(mat)
  new_labels <- labels[match(names(groups), labels$key), ]
  new_labels$sample <- names(groups)
  new_labels$source_samples <- vapply(groups, paste, character(1), collapse = ";")
  list(mat = out, labels = new_labels)
}

calculate_effect_table <- function(mat, labels, dataset, platform, source_file, transform_note) {
  collapsed <- collapse_replicates(mat, labels)
  mat <- collapsed$mat
  labels <- collapsed$labels
  asd_cols <- labels$sample[labels$group == "ASD"]
  ctl_cols <- labels$sample[labels$group == "control"]
  if (length(asd_cols) < 2 || length(ctl_cols) < 2) stop(dataset, " has fewer than two ASD/control biological samples after replicate collapse")
  res <- lapply(seq_len(nrow(mat)), function(i) {
    x1 <- as.numeric(mat[i, asd_cols]); x0 <- as.numeric(mat[i, ctl_cols])
    x1 <- x1[is.finite(x1)]; x0 <- x0[is.finite(x0)]
    n1 <- length(x1); n0 <- length(x0)
    if (n1 < 2 || n0 < 2) return(NULL)
    sd1 <- sd(x1); sd0 <- sd(x0)
    if (!is.finite(sd1) || !is.finite(sd0) || (sd1 == 0 && sd0 == 0)) return(NULL)
    sp <- sqrt(((n1 - 1) * sd1^2 + (n0 - 1) * sd0^2) / (n1 + n0 - 2))
    if (!is.finite(sp) || sp == 0) return(NULL)
    d <- (mean(x1) - mean(x0)) / sp
    J <- 1 - 3 / (4 * (n1 + n0) - 9)
    g <- J * d
    var_g <- ((n1 + n0) / (n1 * n0)) + (g^2 / (2 * (n1 + n0 - 2)))
    data.frame(gene = rownames(mat)[i], dataset = dataset, platform = platform,
               model_role = role_for_dataset(dataset), caveat = caveat_for_dataset(dataset),
               ASD_n = n1, control_n = n0, ASD_mean = mean(x1), control_mean = mean(x0),
               ASD_sd = sd1, control_sd = sd0, hedges_g = g, variance = var_g, se = sqrt(var_g),
               direction = ifelse(g > 0, "higher_in_ASD", ifelse(g < 0, "lower_in_ASD", "zero")),
               source_file = source_file, transform_note = transform_note, stringsAsFactors = FALSE)
  })
  list(effects = rbindlist_fill(res), labels = labels)
}

process_lcl_series <- function(cfg, dataset) {
  sm <- read_series_matrix(file.path(cfg$out_dir, "01_source_files", paste0(dataset, "_series_matrix.txt.gz")))
  md <- sample_meta_frame(sm$sample_meta)
  platform <- unique(md$platform)[1]
  pmap <- read_platform_map(cfg, platform)
  labels <- rbindlist_fill(lapply(seq_len(nrow(md)), function(i) {
    lab <- label_lcl_sample(dataset, md$title[i], md$source[i], md$characteristics[i])
    data.frame(dataset = dataset, sample = md$sample[i], title = md$title[i], source = md$source[i],
               characteristics = md$characteristics[i], group = lab[1], label_reason = lab[2],
               bio_id = bio_id_for_sample(dataset, md$sample[i], md$title[i], md$source[i]),
               model_role = role_for_dataset(dataset), caveat = caveat_for_dataset(dataset),
               stringsAsFactors = FALSE)
  }))
  vals <- sm$values
  value_cols <- intersect(labels$sample[labels$group %in% c("ASD", "control")], names(vals))
  long <- merge(data.frame(probe_id = as.character(vals$probe_id), vals[, value_cols, drop = FALSE], check.names = FALSE),
                pmap, by = "probe_id")
  mat <- as.matrix(long[, value_cols, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- long$gene
  mat <- collapse_gene_matrix(mat)
  eff <- calculate_effect_table(mat, labels, dataset, platform, paste0(dataset, "_series_matrix.txt.gz"), "GEO series matrix values used as supplied; probe rows collapsed to gene means")
  eff$all_labels <- labels
  list(labels = eff$labels, raw_labels = labels, effects = eff$effects, processed = data.frame(dataset = dataset, platform = platform, genes = nrow(mat), samples = ncol(mat), stringsAsFactors = FALSE))
}

dl_meta_one <- function(g, v) {
  ok <- is.finite(g) & is.finite(v) & v > 0
  g <- g[ok]
  v <- v[ok]
  k <- length(g)
  if (k < 2) return(NULL)
  dl <- metafor::rma.uni(yi = g, vi = v, method = "DL", test = "z")
  mkh <- metafor::rma.uni(yi = g, vi = v, method = "DL", test = "adhoc")
  data.frame(
    k = k,
    pooled_g = as.numeric(dl$b[1]),
    se = as.numeric(dl$se),
    tau2 = as.numeric(dl$tau2),
    Q = as.numeric(dl$QE),
    I2 = as.numeric(dl$I2),
    p_value = as.numeric(dl$pval),
    DL_CI_lower = as.numeric(dl$ci.lb),
    DL_CI_upper = as.numeric(dl$ci.ub),
    mKH_CI_lower = as.numeric(mkh$ci.lb),
    mKH_CI_upper = as.numeric(mkh$ci.ub)
  )
}

run_meta_model <- function(effects, model, role, caveat) {
  genes <- split(effects, effects$gene)
  meta <- lapply(names(genes), function(gene) {
    m <- dl_meta_one(genes[[gene]]$hedges_g, genes[[gene]]$variance)
    if (is.null(m)) return(NULL)
    cbind(data.frame(gene = gene, model = model, model_role = role,
                     contributing_datasets = paste(unique(genes[[gene]]$dataset), collapse = ";"),
                     stringsAsFactors = FALSE), m)
  })
  out <- rbindlist_fill(meta)
  gene_counts <- table(effects$gene)
  k1 <- data.frame(gene = names(gene_counts)[gene_counts == 1], model = model, model_role = "descriptive_k1", stringsAsFactors = FALSE)
  if (!nrow(out)) {
    summary <- data.frame(model = model, model_role = role, datasets = paste(unique(effects$dataset), collapse = ";"),
                          genes_meta_analysed = 0, k1_descriptive = nrow(k1), DL_nonzero = 0,
                          FDR_significant = 0, mKH_interval_supported = 0, FDR_mKH_overlap = 0,
                          median_k = NA_real_, median_I2 = NA_real_, caveat = caveat)
    return(list(meta = out, k1 = k1, summary = summary))
  }
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out$DL_nonzero <- out$DL_CI_lower > 0 | out$DL_CI_upper < 0
  out$FDR_significant <- out$FDR < 0.05
  out$mKH_interval_excludes_zero <- out$mKH_CI_lower > 0 | out$mKH_CI_upper < 0
  out$FDR_mKH_overlap <- out$FDR_significant & out$mKH_interval_excludes_zero
  out$direction <- ifelse(out$pooled_g > 0, "higher_in_ASD", "lower_in_ASD")
  out <- out[order(out$FDR, -abs(out$pooled_g)), ]
  summary <- data.frame(model = model, model_role = role, datasets = paste(unique(effects$dataset), collapse = ";"),
                        genes_meta_analysed = nrow(out), k1_descriptive = nrow(k1),
                        DL_nonzero = sum(out$DL_nonzero), FDR_significant = sum(out$FDR_significant),
                        mKH_interval_supported = sum(out$mKH_interval_excludes_zero),
                        FDR_mKH_overlap = sum(out$FDR_mKH_overlap),
                        median_k = median(out$k), median_I2 = median(out$I2, na.rm = TRUE),
                        caveat = caveat, stringsAsFactors = FALSE)
  list(meta = out, k1 = k1, summary = summary)
}

package_version_table <- function() {
  pkgs <- c("data.table", "AnnotationDbi", "org.Hs.eg.db", "huex10sttranscriptcluster.db", "metafor", "openxlsx")
  data.frame(package = pkgs,
             installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
             version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
             stringsAsFactors = FALSE)
}
