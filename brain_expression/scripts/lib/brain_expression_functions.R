suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(metafor)
  library(openxlsx)
})

brain_expression_config <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(root) == "scripts") root <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = TRUE)
  if (basename(root) == "brain_expression") {
    out_dir <- root
    repo_root <- normalizePath(file.path(out_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  } else {
    repo_root <- root
    out_dir <- file.path(repo_root, "pipelines", "brain_expression")
  }
  source_cache <- Sys.getenv("BRAIN_EXPRESSION_PUBLIC_SOURCE_CACHE", unset = "")
  list(
    repo_root = repo_root,
    out_dir = out_dir,
    scripts_dir = file.path(out_dir, "scripts"),
    source_cache = source_cache,
    validation_reference_dir = Sys.getenv("BRAIN_EXPRESSION_VALIDATION_REFERENCE_DIR", unset = ""),
    subdirs = c("00_manifest", "01_source_files", "02_sample_metadata", "03_processed_expression",
                "04_dataset_gene_summaries", "05_effect_sizes", "06_grouped_brain_meta_analysis",
                "07_region_subtissue_sensitivity", "08_platform_sensitivity", "09_quality_control")
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
geo_matrix_url <- function(gse, file = paste0(gse, "_series_matrix.txt.gz")) {
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", geo_series_prefix(gse), gse, file)
}
geo_platform_url <- function(gpl) {
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz", geo_platform_prefix(gpl), gpl, gpl)
}
geo_platform_soft_url <- function(gpl) {
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s_family.soft.gz", geo_platform_prefix(gpl), gpl, gpl)
}
geo_supp_url <- function(gse, file) {
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", geo_series_prefix(gse), gse, file)
}

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
      options(timeout = max(1800, old_timeout))
      on.exit(options(timeout = old_timeout), add = TRUE)
      download.file(url, dest, mode = "wb", quiet = TRUE)
      status <- "downloaded"
      source_used <- "downloaded from public repository"
    }
  }
  data.frame(file_name = filename,
             local_path = normalizePath(dest, winslash = "/", mustWork = FALSE),
             url = url,
             source_used = source_used,
             status = status,
             bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_,
             stringsAsFactors = FALSE)
}

stage_brain_expression_sources <- function(cfg) {
  dest_dir <- file.path(cfg$out_dir, "01_source_files")
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  routes <- list(
    c("GSE113834", "GSE113834_series_matrix.txt.gz", geo_matrix_url("GSE113834"), "series_matrix"),
    c("GSE28475", "GSE28475-GPL13388_series_matrix.txt.gz", geo_matrix_url("GSE28475", "GSE28475-GPL13388_series_matrix.txt.gz"), "series_matrix"),
    c("GSE28475", "GSE28475-GPL6883_series_matrix.txt.gz", geo_matrix_url("GSE28475", "GSE28475-GPL6883_series_matrix.txt.gz"), "series_matrix"),
    c("GSE28521", "GSE28521_series_matrix.txt.gz", geo_matrix_url("GSE28521"), "series_matrix"),
    c("GSE36315", "GSE36315_series_matrix.txt.gz", geo_matrix_url("GSE36315"), "series_matrix"),
    c("GSE38322", "GSE38322_series_matrix.txt.gz", geo_matrix_url("GSE38322"), "series_matrix"),
    c("GSE102741", "GSE102741_log2RPKMcounts.xlsx", geo_supp_url("GSE102741", "GSE102741_log2RPKMcounts.xlsx"), "supplemental_matrix"),
    c("GSE211154", "GSE211154_counts.txt.gz", geo_supp_url("GSE211154", "GSE211154_counts.txt.gz"), "supplemental_counts"),
    c("GSE236761", "GSE236761_raw_counts.txt.gz", geo_supp_url("GSE236761", "GSE236761_raw_counts.txt.gz"), "supplemental_counts"),
    c("GSE269105", "GSE269105_Processed_data_gene_expression_read_and_fpkm.txt.gz", geo_supp_url("GSE269105", "GSE269105_Processed_data_gene_expression_read_and_fpkm.txt.gz"), "supplemental_counts"),
    c("GSE59288", "GSE59288_exp_mRNA.txt.gz", geo_supp_url("GSE59288", "GSE59288_exp_mRNA.txt.gz"), "supplemental_matrix"),
    c("GSE62098", "GSE62098_processedFPKMS.xls.gz", geo_supp_url("GSE62098", "GSE62098_processedFPKMS.xls.gz"), "supplemental_workbook"),
    c("GSE64018", "GSE64018_countlevel_12asd_12ctl.txt.gz", geo_supp_url("GSE64018", "GSE64018_countlevel_12asd_12ctl.txt.gz"), "supplemental_counts"),
    c("HGNC", "hgnc_complete_set.txt", "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt", "annotation")
  )
  platforms <- list(
    c("GPL10558", "GPL10558.annot.gz", geo_platform_url("GPL10558"), "platform_annotation"),
    c("GPL6883", "GPL6883.annot.gz", geo_platform_url("GPL6883"), "platform_annotation"),
    c("GPL13388", "GPL13388_family.soft.gz", geo_platform_soft_url("GPL13388"), "platform_annotation"),
    c("GPL15314", "GPL15314_family.soft.gz", geo_platform_soft_url("GPL15314"), "platform_annotation"),
    c("GPL15207", "GPL15207_family.soft.gz", geo_platform_soft_url("GPL15207"), "platform_annotation")
  )
  all <- c(routes, platforms)
  out <- lapply(all, function(z) {
    s <- stage_one_file(cfg$source_cache, dest_dir, z[[2]], z[[3]])
    s$dataset_id <- z[[1]]
    s$route_type <- z[[4]]
    s
  })
  rbindlist_fill(out)
}

safe_split_gene <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("\\s*\\(.*?\\)", "", x)
  parts <- unlist(strsplit(x, "///|//|;|,|\\|", perl = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  parts <- parts[!toupper(parts) %in% c("NA", "N/A", "---", "--", "NULL", "NONE")]
  unique(parts[grepl("^[A-Za-z0-9_.-]+$", parts)])
}

map_ensembl_to_symbol <- function(ids) {
  base <- sub("\\..*$", "", ids)
  base_unique <- unique(base[!is.na(base) & grepl("^ENSG", base)])
  if (!length(base_unique)) return(rep(NA_character_, length(ids)))
  out <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = base_unique, keytype = "ENSEMBL",
                               column = "SYMBOL", multiVals = "first")
  as.character(out[base])
}

map_entrez_to_symbol <- function(ids) {
  out <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(as.character(ids)), keytype = "ENTREZID",
                               column = "SYMBOL", multiVals = "first")
  as.character(out[as.character(ids)])
}

read_series_matrix <- function(file) {
  lines <- readLines(gzfile(file, open = "rt"), warn = FALSE)
  begin <- match("!series_matrix_table_begin", lines)
  end <- match("!series_matrix_table_end", lines)
  if (is.na(begin) || is.na(end)) stop("No series matrix table found in ", file)
  sample_meta <- list()
  for (ln in lines[seq_len(begin - 1)][grepl("^!Sample_", lines[seq_len(begin - 1)])]) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    key <- sub("^!", "", parts[1])
    sample_meta[[key]] <- gsub('^"|"$', "", parts[-1])
  }
  dt <- data.table::fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                          data.table = FALSE, check.names = FALSE)
  names(dt)[1] <- "probe_id"
  list(values = dt, sample_meta = sample_meta)
}

sample_meta_frame <- function(sm) {
  n <- length(sm$Sample_geo_accession)
  field <- function(k) if (!is.null(sm[[k]])) sm[[k]] else rep("", n)
  ch <- field("Sample_characteristics_ch1")
  char <- if (is.matrix(ch)) apply(ch, 2, paste, collapse = "; ") else ch
  if (length(char) != n) char <- rep(paste(ch, collapse = "; "), n)
  data.frame(
    sample = field("Sample_geo_accession"),
    title = field("Sample_title"),
    source_name = field("Sample_source_name_ch1"),
    characteristics = char,
    stringsAsFactors = FALSE
  )
}

stream_soft_platform_map <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  in_table <- FALSE
  header <- NULL
  id_i <- gene_i <- NA_integer_
  gene_kind <- "symbol"
  rows <- list()
  n <- 0
  repeat {
    chunk <- readLines(con, n = 20000, warn = FALSE)
    if (!length(chunk)) break
    for (ln in chunk) {
      if (!in_table) {
        if (identical(ln, "!platform_table_begin")) in_table <- TRUE
        next
      }
      if (identical(ln, "!platform_table_end")) {
        in_table <- FALSE
        break
      }
      if (is.null(header)) {
        header <- strsplit(ln, "\t", fixed = TRUE)[[1]]
        lower <- tolower(header)
        id_i <- which(lower %in% c("id", "id_ref"))[1]
        gene_i <- which(lower %in% c("gene symbol", "gene_symbol", "symbol"))[1]
        if (is.na(gene_i)) {
          gene_i <- which(lower %in% c("ensembl_id", "ensembl"))[1]
          gene_kind <- "ensembl"
        }
        if (is.na(gene_i)) {
          gene_i <- which(lower %in% c("entrez gene", "entrez_gene", "entrezid"))[1]
          gene_kind <- "entrez"
        }
        if (is.na(gene_i)) {
          gene_i <- grep("gene.?symbol|refgene|gene_assignment|transcript", header, ignore.case = TRUE)[1]
          gene_kind <- "symbol"
        }
        next
      }
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(parts) < max(id_i, gene_i, na.rm = TRUE)) next
      genes <- safe_split_gene(parts[gene_i])
      if (!length(genes)) next
      n <- n + 1
      rows[[n]] <- data.frame(probe_id = parts[id_i], gene = genes, gene_kind = gene_kind, stringsAsFactors = FALSE)
    }
  }
  out <- unique(rbindlist_fill(rows))
  if (!nrow(out)) return(out[, c("probe_id", "gene"), drop = FALSE])
  if (any(out$gene_kind == "ensembl")) {
    idx <- out$gene_kind == "ensembl"
    out$gene[idx] <- map_ensembl_to_symbol(out$gene[idx])
  }
  if (any(out$gene_kind == "entrez")) {
    idx <- out$gene_kind == "entrez"
    out$gene[idx] <- map_entrez_to_symbol(out$gene[idx])
  }
  out <- out[!is.na(out$gene) & nzchar(out$gene), c("probe_id", "gene"), drop = FALSE]
  unique(out)
}

read_platform_map <- function(file) {
  cache_dir <- file.path(dirname(dirname(file)), "03_processed_expression", "platform_maps")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0(basename(file), ".platform_map.csv.gz"))
  if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
    return(data.table::fread(cache_file, data.table = FALSE))
  }
  if (basename(file) == "GPL15207_family.soft.gz") {
    dt <- data.table::fread(file, skip = "ID\tGeneChip", nrows = 49359,
                            select = c("ID", "Gene Symbol"), data.table = FALSE,
                            showProgress = FALSE)
    out <- unique(data.frame(probe_id = as.character(dt[["ID"]]),
                             gene = as.character(dt[["Gene Symbol"]]),
                             stringsAsFactors = FALSE))
    out <- out[!is.na(out$gene) & nzchar(out$gene), ]
    write_csv_safe(out, cache_file)
    return(out)
  }
  if (basename(file) == "GPL15314_family.soft.gz") {
    dt <- data.table::fread(file, skip = "ID\tSPOT_ID", nrows = 60756,
                            select = c("ID", "ENSEMBL_ID"), data.table = FALSE,
                            showProgress = FALSE, fill = TRUE)
    ens <- as.character(dt[["ENSEMBL_ID"]])
    ens[!grepl("^ENSG", ens)] <- NA_character_
    genes <- map_ensembl_to_symbol(stats::na.omit(ens))
    gene_map <- setNames(genes, stats::na.omit(ens))
    genes <- unname(gene_map[ens])
    out <- unique(data.frame(probe_id = as.character(dt[["ID"]]),
                             gene = genes,
                             stringsAsFactors = FALSE))
    out <- out[!is.na(out$gene) & nzchar(out$gene), ]
    write_csv_safe(out, cache_file)
    return(out)
  }
  if (grepl("_family\\.soft\\.gz$", file, ignore.case = TRUE)) {
    out <- stream_soft_platform_map(file)
    write_csv_safe(out, cache_file)
    return(out)
  }
  lines <- readLines(gzfile(file, "rt"), warn = FALSE)
  header <- grep("^(ID|\"ID\")\t", lines)[1]
  if (is.na(header)) return(data.frame(probe_id = character(), gene = character()))
  header_fields <- strsplit(lines[header], "\t", fixed = TRUE)[[1]]
  id_i <- which(tolower(header_fields) %in% c("id", "id_ref"))[1]
  lower <- tolower(header_fields)
  gene_i <- which(lower %in% c("gene symbol", "gene_symbol", "symbol"))[1]
  if (is.na(gene_i)) gene_i <- grep("gene.?symbol|refgene|gene_assignment|transcript|orf", header_fields, ignore.case = TRUE)[1]
  rows <- lapply(lines[(header + 1):length(lines)], function(ln) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < max(id_i, gene_i, na.rm = TRUE)) return(NULL)
    genes <- safe_split_gene(parts[gene_i])
    if (!length(genes)) return(NULL)
    data.frame(probe_id = parts[id_i], gene = genes, stringsAsFactors = FALSE)
  })
  out <- unique(rbindlist_fill(rows))
  write_csv_safe(out, cache_file)
  out
}

guess_region <- function(dataset, title = "", source = "", characteristics = "") {
  txt <- tolower(paste(dataset, title, source, characteristics))
  if (grepl("corpus callosum", txt)) return(c("corpus callosum", "non_cortical"))
  if (grepl("purkinje|cerebell", txt)) return(c("cerebellum/Purkinje", "cerebellum"))
  if (grepl("ba19|brodmann area 19|occipital", txt)) return(c("BA19/occipital cortex", "cortex"))
  if (grepl("ba41|ba42|ba22|temporal", txt)) return(c("temporal cortex", "cortex"))
  if (grepl("prefrontal|dlpfc|ba8|ba9", txt)) return(c("prefrontal cortex", "cortex"))
  c("post-mortem brain", "mixed_brain")
}

label_series_sample <- function(dataset, sample, title, source = "", characteristics = "") {
  txt <- tolower(paste(title, source, characteristics))
  if (dataset == "GSE113834") {
    if (!grepl("input", txt)) return(c("exclude", "WASH/ELUTED fraction excluded; INPUT retained"))
    if (grepl("ctrl|control", txt)) return(c("control", "CTRL INPUT sample"))
    if (grepl("asd|autism", txt)) return(c("ASD", "ASD INPUT sample"))
  }
  if (dataset == "GSE28521") {
    if (grepl("^a_", title, ignore.case = TRUE)) return(c("ASD", "A-coded ASD sample"))
    if (grepl("^c_", title, ignore.case = TRUE)) return(c("control", "C-coded control sample"))
  }
  if (grepl("reference", txt)) return(c("exclude", "reference/non-case-control sample excluded"))
  if (grepl("control|ctrl|healthy", txt)) return(c("control", "public metadata control label"))
  if (grepl("autism|asd", txt)) return(c("ASD", "public metadata ASD label"))
  c("exclude", "unclear ASD/control label")
}

long_to_gene_matrix <- function(dt, id_col, value_cols, id_type = c("symbol", "ensembl", "entrez"), log_transform = FALSE, cpm_transform = FALSE) {
  id_type <- match.arg(id_type)
  genes <- as.character(dt[[id_col]])
  if (id_type == "ensembl") genes <- map_ensembl_to_symbol(genes)
  if (id_type == "entrez") genes <- map_entrez_to_symbol(genes)
  keep <- !is.na(genes) & nzchar(genes)
  mat <- as.matrix(dt[keep, value_cols, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- genes[keep]
  if (cpm_transform) mat <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  if (log_transform) mat <- log2(mat + 1)
  collapse_gene_matrix(mat)
}

collapse_gene_matrix <- function(mat) {
  if (!nrow(mat)) return(mat)
  rows <- split(seq_len(nrow(mat)), rownames(mat))
  out <- vapply(rows, function(idx) {
    if (length(idx) == 1) as.numeric(mat[idx, ]) else colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(mat)))
  out <- t(out)
  colnames(out) <- colnames(mat)
  out
}

calculate_effect_table <- function(mat, labels, dataset, study_id, assay, platform, region, region_group, model_role, source_file, transform_note) {
  labels <- labels[labels$sample %in% colnames(mat), ]
  labels <- labels[labels$group %in% c("ASD", "control"), ]
  mat <- mat[, labels$sample, drop = FALSE]
  asd_cols <- labels$sample[labels$group == "ASD"]
  ctl_cols <- labels$sample[labels$group == "control"]
  if (length(asd_cols) < 2 || length(ctl_cols) < 2) stop(dataset, " has fewer than two ASD or control samples")
  res <- lapply(seq_len(nrow(mat)), function(i) {
    x1 <- as.numeric(mat[i, asd_cols])
    x0 <- as.numeric(mat[i, ctl_cols])
    x1 <- x1[is.finite(x1)]
    x0 <- x0[is.finite(x0)]
    n1 <- length(x1); n0 <- length(x0)
    if (n1 < 2 || n0 < 2) return(NULL)
    sd1 <- stats::sd(x1); sd0 <- stats::sd(x0)
    if (!is.finite(sd1) || !is.finite(sd0) || (sd1 == 0 && sd0 == 0)) return(NULL)
    sp <- sqrt(((n1 - 1) * sd1^2 + (n0 - 1) * sd0^2) / (n1 + n0 - 2))
    if (!is.finite(sp) || sp == 0) return(NULL)
    d <- (mean(x1) - mean(x0)) / sp
    J <- 1 - 3 / (4 * (n1 + n0) - 9)
    g <- J * d
    var_g <- ((n1 + n0) / (n1 * n0)) + (g^2 / (2 * (n1 + n0 - 2)))
    data.frame(
      gene = rownames(mat)[i], dataset = dataset, study_id = study_id, source_accession = dataset,
      assay = assay, platform = platform, brain_region = region, region_group = region_group,
      model_role = model_role, ASD_n = n1, control_n = n0,
      ASD_mean = mean(x1), control_mean = mean(x0), ASD_sd = sd1, control_sd = sd0,
      hedges_g = g, variance = var_g, se = sqrt(var_g),
      direction = ifelse(g > 0, "higher_in_ASD", ifelse(g < 0, "lower_in_ASD", "zero")),
      source_file = source_file, transform_note = transform_note, stringsAsFactors = FALSE)
  })
  rbindlist_fill(res)
}

process_series_dataset <- function(cfg, dataset, matrix_file, platform_file, assay = "microarray", platform = NA_character_) {
  sm <- read_series_matrix(file.path(cfg$out_dir, "01_source_files", matrix_file))
  meta <- sample_meta_frame(sm$sample_meta)
  labels <- lapply(seq_len(nrow(meta)), function(i) {
    lab <- label_series_sample(dataset, meta$sample[i], meta$title[i], meta$source_name[i], meta$characteristics[i])
    reg <- guess_region(dataset, meta$title[i], meta$source_name[i], meta$characteristics[i])
    data.frame(dataset = dataset, study_id = dataset, sample = meta$sample[i], group = lab[1],
               label_reason = lab[2], title = meta$title[i], tissue = meta$source_name[i],
               brain_region = reg[1], region_group = reg[2], assay = assay, platform = platform,
               model_role = "brain_grouped_primary_candidate", stringsAsFactors = FALSE)
  })
  labels <- rbindlist_fill(labels)
  dataset_id <- if (dataset == "GSE28475" && grepl("GPL6883", matrix_file)) "GSE28475-GPL6883" else if (dataset == "GSE28475") "GSE28475-GPL13388" else dataset
  if (dataset_id %in% c("GSE28475-GPL6883", "GSE28521", "GSE38322")) {
    # These public series matrices do not consistently encode all diagnosis
    # and/or brain-region fields needed for branch-specific sensitivity models.
    # Use the curated public phenotype map from the repository-validation
    # package, and keep it as an explicit input dependency.
    pheno_file <- file.path(cfg$out_dir, "02_sample_metadata", "brain_expression_sample_phenotype_table.csv")
    if (!file.exists(pheno_file) && nzchar(cfg$validation_reference_dir)) {
      pheno_file <- file.path(cfg$validation_reference_dir, "brain_expression_sample_phenotype_table.csv")
    }
    if (!file.exists(pheno_file)) {
      stop(
        "Curated public phenotype table not found for ", dataset_id, ": ", pheno_file,
        "\nProvide it in 02_sample_metadata/ or set BRAIN_EXPRESSION_VALIDATION_REFERENCE_DIR.",
        call. = FALSE
      )
    }
    pheno <- data.table::fread(pheno_file, data.table = FALSE)
    pheno <- pheno[pheno$dataset == dataset_id, ]
    labels_base <- as.data.frame(labels)[, setdiff(names(labels), c("group", "label_reason", "brain_region", "region_group")), drop = FALSE]
    labels <- merge(labels_base,
                    pheno[, c("sample", "group", "label_reason", "brain_region", "region_group")],
                    by = "sample", all.x = TRUE)
    labels$group[is.na(labels$group)] <- "exclude"
    labels$label_reason[is.na(labels$label_reason)] <- "not present in curated public phenotype map"
  }
  pmap <- read_platform_map(file.path(cfg$out_dir, "01_source_files", platform_file))
  if (!nrow(pmap)) stop(dataset_id, " has no valid probe-to-gene symbol mapping in the public platform annotation")
  vals <- sm$values
  value_cols <- intersect(labels$sample[labels$group %in% c("ASD", "control")], names(vals))
  long <- merge(data.frame(probe_id = as.character(vals$probe_id), vals[, value_cols, drop = FALSE], check.names = FALSE),
                pmap, by = "probe_id")
  mat <- as.matrix(long[, value_cols, drop = FALSE]); mode(mat) <- "numeric"; rownames(mat) <- long$gene
  mat <- collapse_gene_matrix(mat)
  common_region <- unique(labels$brain_region[labels$group %in% c("ASD", "control")])
  common_group <- unique(labels$region_group[labels$group %in% c("ASD", "control")])
  labels$dataset <- dataset_id
  eff <- calculate_effect_table(mat, labels, dataset_id, dataset, assay, platform,
                                paste(common_region, collapse = ";"), paste(common_group, collapse = ";"),
                                "brain_grouped_primary_candidate", matrix_file, "GEO series matrix values used as supplied")
  list(labels = labels, effects = eff, processed = data.frame(dataset = dataset_id, genes = nrow(mat), samples = ncol(mat)))
}

curated_label_table <- function(dataset, samples) {
  if (dataset == "GSE102741") {
    groups <- c(rep("control", 39), rep("ASD", 13))
    samples <- paste0("sample", seq_along(groups))
    return(data.frame(sample = samples, group = groups, title = c(paste("Control", 1:39), paste("ASD", 1:13)),
                      brain_region = "prefrontal cortex", region_group = "cortex", label_reason = "public sample-order labels from GEO/design metadata"))
  }
  if (dataset == "GSE211154") {
    ids <- samples
    asd <- c("4334","4849","4999","5027","5115","5144","5176","5278","5294","5297","5308","5403","5419","5565","5841","5878","5940","5945","5978","6033")
    ctrl <- c("1158","4599","5079","5113","5163","5170","5334","5387","5391","5566","5646","5669","5705","5813","5889","5893","5922","5926","914")
    return(data.frame(sample = ids, group = ifelse(ids %in% asd, "ASD", ifelse(ids %in% ctrl, "control", "exclude")),
                      title = ids, brain_region = "cerebellum/Purkinje", region_group = "cerebellum",
                      label_reason = "public donor-ID diagnosis mapping from GEO sample titles"))
  }
  if (dataset == "GSE269105") {
    ids <- sub("^read_count_", "", samples)
    asd <- c("CBell_B4925","CBell_B5000","CBell_B5144","CBell_B5173","CBell_B6337","CBell_B6401","CBell_B6677")
    ctrl <- c("CBell_B5251","CBell_B5718","CBell_B5873","CBell_B7333","CBell_B7369")
    return(data.frame(sample = samples, group = ifelse(ids %in% asd, "ASD", ifelse(ids %in% ctrl, "control", "exclude")),
                      title = ids, brain_region = "cerebellum", region_group = "cerebellum",
                      label_reason = "public sample-ID diagnosis mapping; CBell_B6994 excluded"))
  }
  if (dataset == "GSE59288") {
    return(data.frame(sample = samples, group = ifelse(grepl("^aut", samples, ignore.case = TRUE), "ASD", "control"),
                      title = samples, brain_region = "prefrontal cortex", region_group = "cortex",
                      label_reason = "public matrix column label: aut versus ptr"))
  }
  if (dataset == "GSE64018") {
    return(data.frame(sample = samples, group = ifelse(seq_along(samples) <= 12, "ASD", "control"),
                      title = samples, brain_region = "temporal cortex", region_group = "cortex",
                      label_reason = "public matrix contains first 12 ASD and next 12 control temporal-cortex columns"))
  }
  if (dataset == "GSE62098") {
    ids <- samples
    asd <- c("5308", "5144", "4899", "4999", "5302", "5403")
    ctrl <- c("4727", "5163", "5391", "5242", "4670", "5407")
    return(data.frame(sample = ids, group = ifelse(ids %in% asd, "ASD", ifelse(ids %in% ctrl, "control", "exclude")),
                      title = ids, brain_region = "corpus callosum", region_group = "non_cortical",
                      label_reason = "public GSE62098 design: six ASD and six matched controls"))
  }
  if (dataset == "GSE236761") {
    return(data.frame(sample = samples, group = ifelse(grepl("^A|^ASD", samples), "ASD", "control"),
                      title = samples, brain_region = "BA19/occipital cortex", region_group = "cortex",
                      label_reason = "public Total-RNA column labels"))
  }
  stop("No curated label rule for ", dataset)
}

process_supplemental_datasets <- function(cfg) {
  src <- file.path(cfg$out_dir, "01_source_files")
  out <- list()
  # GSE102741 log2 RPKM workbook
  x <- as.data.frame(readxl::read_excel(file.path(src, "GSE102741_log2RPKMcounts.xlsx"), sheet = 1))
  names(x)[1] <- "gene_id"
  mat <- long_to_gene_matrix(x, "gene_id", names(x)[-1], "ensembl", log_transform = FALSE)
  labs <- curated_label_table("GSE102741", colnames(mat))
  out[["GSE102741"]] <- list(labels = transform(labs, dataset="GSE102741", study_id="GSE102741", tissue="DLPFC", assay="bulk RNA-seq", platform="log2 RPKM", model_role="brain_grouped_primary_candidate"),
                              effects = calculate_effect_table(mat, transform(labs, dataset="GSE102741"), "GSE102741", "GSE102741", "bulk RNA-seq", "log2 RPKM", "prefrontal cortex", "cortex", "brain_grouped_primary_candidate", "GSE102741_log2RPKMcounts.xlsx", "public log2 RPKM values used as supplied"),
                              processed = data.frame(dataset="GSE102741", genes=nrow(mat), samples=ncol(mat)))
  # GSE211154 counts
  d <- fread(file.path(src, "GSE211154_counts.txt.gz"), header = FALSE, data.table = FALSE)
  colnames(d) <- as.character(d[1, ]); d <- d[-1, ]
  mat <- long_to_gene_matrix(d, "Symbol", names(d)[-1], "ensembl", cpm_transform = TRUE)
  labs <- curated_label_table("GSE211154", colnames(mat))
  out[["GSE211154"]] <- list(labels = transform(labs, dataset="GSE211154", study_id="GSE211154", tissue="Purkinje cells", assay="bulk RNA-seq", platform="gene counts", model_role="brain_grouped_primary_candidate"),
                              effects = calculate_effect_table(mat, transform(labs, dataset="GSE211154"), "GSE211154", "GSE211154", "bulk RNA-seq", "gene counts", "cerebellum/Purkinje", "cerebellum", "brain_grouped_primary_candidate", "GSE211154_counts.txt.gz", "public counts transformed to log2(CPM+1) with edgeR"),
                              processed = data.frame(dataset="GSE211154", genes=nrow(mat), samples=ncol(mat)))
  # GSE236761 total RNA counts
  d <- fread(file.path(src, "GSE236761_raw_counts.txt.gz"), data.table = FALSE)
  value_cols <- grep("_Total$", names(d), value = TRUE)
  mat <- long_to_gene_matrix(d, "geneID", value_cols, "symbol", cpm_transform = TRUE)
  labs <- curated_label_table("GSE236761", colnames(mat))
  out[["GSE236761"]] <- list(labels = transform(labs, dataset="GSE236761", study_id="GSE236761", tissue="BA19/occipital cortex", assay="bulk RNA-seq", platform="total RNA raw counts", model_role="brain_grouped_primary_candidate"),
                              effects = calculate_effect_table(mat, transform(labs, dataset="GSE236761"), "GSE236761", "GSE236761", "bulk RNA-seq", "total RNA raw counts", "BA19/occipital cortex", "cortex", "brain_grouped_primary_candidate", "GSE236761_raw_counts.txt.gz", "Total RNA count columns retained; Poly/translatome excluded; transformed to log2(CPM+1)"),
                              processed = data.frame(dataset="GSE236761", genes=nrow(mat), samples=ncol(mat)))
  # GSE269105 read-count columns
  d <- fread(file.path(src, "GSE269105_Processed_data_gene_expression_read_and_fpkm.txt.gz"), data.table = FALSE)
  value_cols <- grep("^read_count_", names(d), value = TRUE)
  mat <- long_to_gene_matrix(d, "Gene id", value_cols, "entrez", cpm_transform = TRUE)
  labs <- curated_label_table("GSE269105", colnames(mat))
  out[["GSE269105"]] <- list(labels = transform(labs, dataset="GSE269105", study_id="GSE269105", tissue="cerebellum", assay="bulk RNA-seq", platform="processed read counts", model_role="brain_grouped_primary_candidate"),
                              effects = calculate_effect_table(mat, transform(labs, dataset="GSE269105"), "GSE269105", "GSE269105", "bulk RNA-seq", "processed read counts", "cerebellum", "cerebellum", "brain_grouped_primary_candidate", "GSE269105_Processed_data_gene_expression_read_and_fpkm.txt.gz", "read-count columns transformed to log2(CPM+1); one unresolved sample excluded"),
                              processed = data.frame(dataset="GSE269105", genes=nrow(mat), samples=ncol(mat)))
  # GSE59288 mRNA expression
  f59288 <- file.path(src, "GSE59288_exp_mRNA.txt.gz")
  hdr <- scan(gzfile(f59288, "rt"), what = character(), nlines = 1, quiet = TRUE)
  d <- fread(f59288, header = FALSE, skip = 1, data.table = FALSE)
  names(d) <- c("gene_id", hdr)
  mat <- long_to_gene_matrix(d, "gene_id", names(d)[-1], "ensembl", log_transform = TRUE)
  labs <- curated_label_table("GSE59288", colnames(mat))
  out[["GSE59288"]] <- list(labels = transform(labs, dataset="GSE59288", study_id="GSE59288", tissue="prefrontal cortex", assay="bulk RNA-seq", platform="mRNA expression matrix", model_role="brain_grouped_primary_candidate"),
                             effects = calculate_effect_table(mat, transform(labs, dataset="GSE59288"), "GSE59288", "GSE59288", "bulk RNA-seq", "mRNA expression matrix", "prefrontal cortex", "cortex", "brain_grouped_primary_candidate", "GSE59288_exp_mRNA.txt.gz", "public expression values transformed to log2(value+1)"),
                             processed = data.frame(dataset="GSE59288", genes=nrow(mat), samples=ncol(mat)))
  # GSE64018 count matrix
  f64018 <- file.path(src, "GSE64018_countlevel_12asd_12ctl.txt.gz")
  d <- fread(f64018, header = FALSE, data.table = FALSE)
  sample64018 <- c(
    "AN02987_ba41-42-22_8.6","AN04682_ba41-42-22_8.2","UMB5278_ba41-42-22_7.8",
    "AN01570_ba41-42-22_7.1","AN00493_ba41-42-22_7.3","AN12457_ba41-42-22_6",
    "AN08166_ba41-42-22_6.4","AN08792_ba41-42-22_5.1","AN01971_ba41-42-22_2.9",
    "AN03632_ba41-42-22_8.1","AN08043_ba41-42-22_7.8","AN09714_ba41-42-22_7.2",
    "AN17425_ba41-42-22_7.9","AN07444_ba41-42-22_8.2","UMB4590_ba41-42-22_8.3",
    "AN10833_ba41-42-22_1.8","AN14757_ba41-42-22_8","AN19760_ba41-42-22_6.1",
    "AN12137_ba41-42-22_6.4","UMB5079_ba41-42-22_7.9","AN08161_ba41-42-22_8.1",
    "AN08677_ba41-42-22_8","UMB4842_ba41-42-22_8.1","AN13295_ba41-42-22_6.3")
  names(d) <- c("gene_id", sample64018)
  mat <- long_to_gene_matrix(d, "gene_id", names(d)[-1], "ensembl", cpm_transform = TRUE)
  labs <- curated_label_table("GSE64018", colnames(mat))
  out[["GSE64018"]] <- list(labels = transform(labs, dataset="GSE64018", study_id="GSE64018", tissue="temporal cortex", assay="bulk RNA-seq", platform="RNA-seq read counts", model_role="brain_grouped_primary_candidate"),
                             effects = calculate_effect_table(mat, transform(labs, dataset="GSE64018"), "GSE64018", "GSE64018", "bulk RNA-seq", "RNA-seq read counts", "temporal cortex", "cortex", "brain_grouped_primary_candidate", "GSE64018_countlevel_12asd_12ctl.txt.gz", "public count matrix transformed to log2(CPM+1)"),
                             processed = data.frame(dataset="GSE64018", genes=nrow(mat), samples=ncol(mat)))
  # GSE62098 gzipped public Excel workbook
  tmp <- tempfile(fileext = ".xls")
  R.utils::gunzip(file.path(src, "GSE62098_processedFPKMS.xls.gz"), tmp, remove = FALSE, overwrite = TRUE)
  wb <- as.data.frame(readxl::read_excel(tmp, sheet = 1, col_names = FALSE))
  unlink(tmp)
  sample_ids <- as.character(unlist(wb[1, seq(1, ncol(wb), by = 3)]))
  values <- list()
  for (i in seq_along(sample_ids)) {
    c0 <- 1 + (i - 1) * 3
    genes <- as.character(wb[[c0]][-(1:2)])
    vals <- suppressWarnings(as.numeric(wb[[c0 + 1]][-(1:2)]))
    values[[i]] <- data.frame(gene = genes, sample = sample_ids[i], value = vals)
  }
  long <- rbindlist_fill(values)
  wide <- dcast(long[is.finite(long$value) & nzchar(long$gene), ], gene ~ sample, value.var = "value", fun.aggregate = mean)
  mat <- as.matrix(wide[, -1, drop = FALSE]); mode(mat) <- "numeric"; rownames(mat) <- wide$gene
  mat <- collapse_gene_matrix(log2(mat + 1))
  labs <- curated_label_table("GSE62098", colnames(mat))
  out[["GSE62098"]] <- list(labels = transform(labs, dataset="GSE62098", study_id="GSE62098", tissue="corpus callosum", assay="bulk RNA-seq", platform="processed FPKM workbook", model_role="brain_grouped_primary_candidate"),
                             effects = calculate_effect_table(mat, transform(labs, dataset="GSE62098"), "GSE62098", "GSE62098", "bulk RNA-seq", "processed FPKM workbook", "corpus callosum", "non_cortical", "brain_grouped_primary_candidate", "GSE62098_processedFPKMS.xls.gz", "public processed FPKM workbook transformed to log2(FPKM+1)"),
                             processed = data.frame(dataset="GSE62098", genes=nrow(mat), samples=ncol(mat)))
  out
}

collapse_within_study <- function(effects) {
  dt <- data.table::as.data.table(effects)
  out <- dt[, {
    w <- 1 / variance
    g <- sum(w * hedges_g) / sum(w)
    v <- 1 / sum(w)
    data.table::data.table(
      dataset = if (.N == 1) dataset[1] else paste0(study_id[1], "_within_study_collapsed"),
      source_accession = paste(unique(dataset), collapse = ";"),
      assay = paste(unique(assay), collapse = ";"),
      platform = paste(unique(platform), collapse = ";"),
      brain_region = paste(unique(brain_region), collapse = ";"),
      region_group = paste(unique(region_group), collapse = ";"),
      model_role = "brain_grouped_primary_candidate",
      ASD_n = max(ASD_n, na.rm = TRUE),
      control_n = max(control_n, na.rm = TRUE),
      ASD_mean = NA_real_,
      control_mean = NA_real_,
      ASD_sd = NA_real_,
      control_sd = NA_real_,
      hedges_g = g,
      variance = v,
      se = sqrt(v),
      direction = ifelse(g > 0, "higher_in_ASD", ifelse(g < 0, "lower_in_ASD", "zero")),
      source_file = paste(unique(source_file), collapse = ";"),
      transform_note = if (.N == 1) transform_note[1] else "within-study multi-region/platform effects collapsed by inverse-variance fixed effect before grouped-brain model"
    )
  }, by = .(gene, study_id)]
  as.data.frame(out)
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
  if (!nrow(effects)) {
    summary <- data.frame(model = model, model_role = role, datasets = "",
                          genes_meta_analysed = 0, k1_descriptive = 0,
                          DL_nonzero = 0, FDR_significant = 0,
                          mKH_interval_supported = 0, FDR_mKH_overlap = 0,
                          median_k = NA_real_, median_I2 = NA_real_,
                          high_I2_gt50 = 0, high_I2_gt75 = 0,
                          caveat = caveat, stringsAsFactors = FALSE)
    return(list(meta = data.frame(), k1 = data.frame(), summary = summary))
  }
  genes <- split(effects, effects$gene)
  meta <- lapply(names(genes), function(gene) {
    m <- dl_meta_one(genes[[gene]]$hedges_g, genes[[gene]]$variance)
    if (is.null(m)) return(NULL)
    cbind(data.frame(gene = gene, model = model, model_role = role,
                     contributing_datasets = paste(unique(genes[[gene]]$dataset), collapse = ";"),
                     contributing_regions = paste(unique(genes[[gene]]$brain_region), collapse = ";"),
                     contributing_assays = paste(unique(genes[[gene]]$assay), collapse = ";"),
                     stringsAsFactors = FALSE), m)
  })
  out <- rbindlist_fill(meta)
  gene_counts <- table(effects$gene)
  k1 <- data.frame(gene = names(gene_counts)[gene_counts == 1], model = model, model_role = "descriptive_k1", stringsAsFactors = FALSE)
  if (!nrow(out)) {
    summary <- data.frame(model = model, model_role = role,
                          datasets = paste(unique(effects$dataset), collapse = ";"),
                          genes_meta_analysed = 0, k1_descriptive = nrow(k1),
                          DL_nonzero = 0, FDR_significant = 0,
                          mKH_interval_supported = 0, FDR_mKH_overlap = 0,
                          median_k = NA_real_, median_I2 = NA_real_,
                          high_I2_gt50 = 0, high_I2_gt75 = 0,
                          caveat = caveat, stringsAsFactors = FALSE)
    return(list(meta = data.frame(), k1 = k1, summary = summary))
  }
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out$DL_nonzero <- out$DL_CI_lower > 0 | out$DL_CI_upper < 0
  out$FDR_significant <- out$FDR < 0.05
  out$mKH_interval_excludes_zero <- out$mKH_CI_lower > 0 | out$mKH_CI_upper < 0
  out$direction <- ifelse(out$pooled_g > 0, "higher_in_ASD", "lower_in_ASD")
  out <- out[order(out$FDR, -abs(out$pooled_g)), ]
  summary <- data.frame(model = model, model_role = role,
                        datasets = paste(unique(effects$dataset), collapse = ";"),
                        genes_meta_analysed = nrow(out), k1_descriptive = nrow(k1),
                        DL_nonzero = sum(out$DL_nonzero), FDR_significant = sum(out$FDR_significant),
                        mKH_interval_supported = sum(out$mKH_interval_excludes_zero),
                        FDR_mKH_overlap = sum(out$FDR_significant & out$mKH_interval_excludes_zero),
                        median_k = stats::median(out$k), median_I2 = stats::median(out$I2, na.rm = TRUE),
                        high_I2_gt50 = sum(out$I2 > 50, na.rm = TRUE), high_I2_gt75 = sum(out$I2 > 75, na.rm = TRUE),
                        caveat = caveat, stringsAsFactors = FALSE)
  list(meta = out, k1 = k1, summary = summary)
}

package_version_table <- function() {
  pkgs <- c("data.table", "readxl", "edgeR", "AnnotationDbi", "org.Hs.eg.db", "metafor", "openxlsx", "R.utils")
  data.frame(package = pkgs,
             installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
             version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
             stringsAsFactors = FALSE)
}
