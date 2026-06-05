suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(limma)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(metafor)
})

blood_expression_config <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(root) == "scripts") root <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = TRUE)
  if (basename(root) == "blood_expression") {
    out_dir <- root
    repo_root <- normalizePath(file.path(out_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  } else {
    repo_root <- root
    out_dir <- file.path(repo_root, "pipelines", "blood_expression")
  }
  scripts_dir <- file.path(out_dir, "scripts")
  source_cache <- Sys.getenv("BLOOD_EXPRESSION_SOURCE_CACHE", unset = "")
  metadata_table <- Sys.getenv(
    "BLOOD_EXPRESSION_METADATA_TABLE",
    unset = file.path(out_dir, "02_sample_metadata", "geo_sample_metadata_compact_with_soft_fallback.csv")
  )
  list(
    repo_root = repo_root,
    out_dir = out_dir,
    scripts_dir = scripts_dir,
    source_cache = source_cache,
    compact_metadata = metadata_table,
    subdirs = c("00_manifest", "01_source_files", "02_sample_metadata", "03_processed_expression",
                "04_dataset_gene_summaries", "05_effect_sizes", "06_meta_analysis",
                "07_sensitivity_analyses", "08_quality_control")
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

package_version_table <- function(pkgs) {
  data.frame(
    package = pkgs,
    installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(pkgs, function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
    }, character(1)),
    stringsAsFactors = FALSE
  )
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
      source_used <- "local public-source cache"
      status <- "copied_from_cache"
    } else {
      old_timeout <- getOption("timeout")
      options(timeout = max(1200, old_timeout))
      on.exit(options(timeout = old_timeout), add = TRUE)
      if (file.exists(dest) && file.info(dest)$size == 0) unlink(dest)
      download.file(url, dest, mode = "wb", quiet = TRUE)
      source_used <- "downloaded_from_public_repository"
      status <- "downloaded"
    }
  }
  data.frame(
    file_name = filename,
    local_path = normalizePath(dest, winslash = "/", mustWork = FALSE),
    url = url,
    source_used = source_used,
    status = status,
    bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_,
    stringsAsFactors = FALSE
  )
}

stage_blood_expression_sources <- function(cfg) {
  dest_dir <- file.path(cfg$out_dir, "01_source_files")
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list(
    list("GSE111175", "GSE111175_series_matrix.txt.gz", geo_matrix_url("GSE111175"), "GPL10558", "series_matrix", "GSE111175"),
    list("GSE111176", "GSE111176_series_matrix.txt.gz", geo_matrix_url("GSE111176"), "GPL10558", "series_matrix", "GSE111176"),
    list("GSE18123", "GSE18123-GPL570_series_matrix.txt.gz", geo_matrix_url("GSE18123", "GSE18123-GPL570_series_matrix.txt.gz"), "GPL570", "series_matrix", "GSE18123-GPL570"),
    list("GSE18123", "GSE18123-GPL6244_series_matrix.txt.gz", geo_matrix_url("GSE18123", "GSE18123-GPL6244_series_matrix.txt.gz"), "GPL6244", "series_matrix", "GSE18123-GPL6244"),
    list("GSE25507", "GSE25507_series_matrix.txt.gz", geo_matrix_url("GSE25507"), "GPL570", "series_matrix", "GSE25507"),
    list("GSE26415", "GSE26415_series_matrix.txt.gz", geo_matrix_url("GSE26415"), "GPL6480", "series_matrix", "GSE26415"),
    list("GSE123302", "GSE123302_series_matrix.txt.gz", geo_matrix_url("GSE123302"), "GPL16686", "series_matrix", "GSE123302"),
    list("GSE212645", "GSE212645_series_matrix.txt.gz", geo_matrix_url("GSE212645"), "GPL24676", "metadata_only", "GSE212645"),
    list("GSE140702", "GSE140702_RAW.tar", geo_supp_url("GSE140702", "GSE140702_RAW.tar"), "GPL20301", "featurecounts_tar", "GSE140702"),
    list("GSE77103", "GSE77103_RAW.tar", geo_supp_url("GSE77103", "GSE77103_RAW.tar"), "GPL17077", "agilent_tar", "GSE77103"),
    list("HGNC", "hgnc_complete_set.txt", "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt", NA_character_, "annotation", "HGNC")
  )
  annot <- unique(c("GPL10558", "GPL570", "GPL6244", "GPL6480"))
  for (gpl in annot) files[[length(files) + 1]] <- list(gpl, paste0(gpl, ".annot.gz"), geo_platform_url(gpl), gpl, "platform_annotation", gpl)
  files[[length(files) + 1]] <- list("GPL16686", "GPL16686_family.soft.gz", geo_platform_soft_url("GPL16686"), "GPL16686", "platform_annotation", "GPL16686")
  out <- lapply(files, function(z) {
    s <- stage_one_file(cfg$source_cache, dest_dir, z[[2]], z[[3]])
    s$dataset_id <- z[[1]]
    s$platform_id <- z[[4]]
    s$route_type <- z[[5]]
    s$route_id <- z[[6]]
    s
  })
  rbindlist_fill(out)
}

read_series_matrix <- function(file) {
  lines <- readLines(gzfile(file, open = "rt"), warn = FALSE)
  begin <- match("!series_matrix_table_begin", lines)
  end <- match("!series_matrix_table_end", lines)
  if (is.na(begin) || is.na(end)) stop("No series matrix table found in ", file)
  meta_lines <- lines[seq_len(begin - 1)]
  sample_meta <- list()
  for (ln in meta_lines[grepl("^!Sample_", meta_lines)]) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    key <- sub("^!", "", parts[1])
    sample_meta[[key]] <- gsub('^"|"$', "", parts[-1])
  }
  table_lines <- lines[(begin + 1):(end - 1)]
  dt <- data.table::fread(text = paste(table_lines, collapse = "\n"), data.table = FALSE, check.names = FALSE)
  names(dt)[1] <- "probe_id"
  list(values = dt, sample_meta = sample_meta)
}

find_header_line <- function(file, patterns = c("^ID\t", "^\"ID\"\t", "^ID,")) {
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  i <- 0
  repeat {
    x <- readLines(con, n = 1, warn = FALSE)
    if (!length(x)) break
    i <- i + 1
    if (any(grepl(paste(patterns, collapse = "|"), x))) return(i)
  }
  NA_integer_
}

read_soft_platform_table <- function(file) {
  lines <- readLines(gzfile(file, open = "rt"), warn = FALSE)
  begin <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  if (is.na(begin) || is.na(end)) return(data.frame())
  data.table::fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                    data.table = FALSE, fill = TRUE, quote = "\"")
}

split_gene_field <- function(x) {
  x <- ifelse(is.na(x), "", x)
  genes <- unlist(strsplit(x, "///|//|;|,|\\|", perl = TRUE), use.names = FALSE)
  genes <- trimws(gsub("\\s*\\(.*?\\)\\s*", "", genes))
  genes <- genes[grepl("^[A-Za-z0-9_.-]+$", genes)]
  genes <- genes[!toupper(genes) %in% c("", "NA", "N/A", "---", "--", "NULL", "NONE")]
  unique(genes)
}

read_platform_map <- function(file) {
  if (grepl("_family\\.soft\\.gz$", file, ignore.case = TRUE)) {
    dt <- read_soft_platform_table(file)
    if (!nrow(dt)) return(data.frame(probe_id = character(), gene = character()))
  } else {
    header <- find_header_line(file)
    if (is.na(header)) return(data.frame(probe_id = character(), gene = character()))
    lines <- readLines(gzfile(file, open = "rt"), warn = FALSE)
    header_fields <- strsplit(lines[header], "\t", fixed = TRUE)[[1]]
    id_i <- which(tolower(header_fields) %in% c("id", "id_ref"))[1]
    gene_i <- which(tolower(header_fields) %in% tolower(c("Gene symbol", "GENE_SYMBOL", "Symbol", "Gene symbol", "UCSC_REFGENE_NAME", "ILMN_Gene", "ORF")))[1]
    if (is.na(gene_i)) gene_i <- grep("gene symbol|symbol|refgene|ilmn_gene|orf", header_fields, ignore.case = TRUE)[1]
    if (is.na(id_i) || is.na(gene_i)) return(data.frame(probe_id = character(), gene = character()))
    rows <- lapply(lines[(header + 1):length(lines)], function(ln) {
      if (!nzchar(ln) || grepl("^#", ln)) return(NULL)
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(parts) < max(id_i, gene_i)) return(NULL)
      genes <- split_gene_field(parts[gene_i])
      if (!length(genes)) return(NULL)
      data.frame(probe_id = parts[id_i], gene = genes, stringsAsFactors = FALSE)
    })
    out <- unique(rbindlist_fill(rows))
    out <- out[out$probe_id != "" & out$gene != "", ]
    return(harmonise_platform_gene_tokens(out))
  }
  id_col <- names(dt)[tolower(names(dt)) %in% c("id", "id_ref")]
  if (!length(id_col)) id_col <- names(dt)[1]
  gene_cols <- names(dt)[tolower(names(dt)) %in% tolower(c(
    "Gene Symbol", "GENE_SYMBOL", "Symbol", "Gene symbol", "gene_assignment",
    "Gene assignment", "UCSC_REFGENE_NAME", "ILMN_Gene", "ORF", "GB_ACC",
    "gene", "genesymbol", "transcript_assignment"
  ))]
  if (!length(gene_cols)) {
    gene_cols <- names(dt)[grepl("symbol|gene|assignment|refgene|orf|gb_acc|transcript", names(dt), ignore.case = TRUE)]
  }
  rows <- list()
  for (gc in gene_cols) {
    vals <- dt[[gc]]
    tmp <- lapply(seq_along(vals), function(i) {
      g <- split_gene_field(as.character(vals[i]))
      if (!length(g)) return(NULL)
      data.frame(probe_id = as.character(dt[[id_col]][i]), gene = g, stringsAsFactors = FALSE)
    })
    rows <- c(rows, tmp)
  }
  out <- unique(rbindlist_fill(rows))
  out <- out[out$probe_id != "" & out$gene != "", ]
  out <- harmonise_platform_gene_tokens(out)
  out
}

harmonise_platform_gene_tokens <- function(map) {
  if (!nrow(map)) return(map)
  map$gene <- as.character(map$gene)
  refseq_like <- grepl("^(NM|NR|XM|XR)_[0-9]+", map$gene)
  entrez_like <- grepl("^[0-9]+$", map$gene)
  converted <- map[!refseq_like & !entrez_like, c("probe_id", "gene"), drop = FALSE]
  if (any(refseq_like)) {
    keys <- unique(map$gene[refseq_like])
    suppressMessages({
      sel <- AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = "REFSEQ", columns = "SYMBOL")
    })
    sel <- sel[!is.na(sel$SYMBOL) & sel$SYMBOL != "", ]
    if (nrow(sel)) {
      tmp <- merge(map[refseq_like, c("probe_id", "gene"), drop = FALSE],
                   sel[, c("REFSEQ", "SYMBOL")], by.x = "gene", by.y = "REFSEQ")
      converted <- rbind(converted, data.frame(probe_id = tmp$probe_id, gene = tmp$SYMBOL, stringsAsFactors = FALSE))
    }
  }
  if (any(entrez_like)) {
    keys <- unique(map$gene[entrez_like])
    suppressMessages({
      sel <- AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = "ENTREZID", columns = "SYMBOL")
    })
    sel <- sel[!is.na(sel$SYMBOL) & sel$SYMBOL != "", ]
    if (nrow(sel)) {
      tmp <- merge(map[entrez_like, c("probe_id", "gene"), drop = FALSE],
                   sel[, c("ENTREZID", "SYMBOL")], by.x = "gene", by.y = "ENTREZID")
      converted <- rbind(converted, data.frame(probe_id = tmp$probe_id, gene = tmp$SYMBOL, stringsAsFactors = FALSE))
    }
  }
  unique(converted)
}

hgnc_maps <- function(cfg) {
  h <- data.table::fread(file.path(cfg$out_dir, "01_source_files", "hgnc_complete_set.txt"), data.table = FALSE)
  symbol <- h$symbol
  names(symbol) <- h$symbol
  ensembl <- data.frame(
    ensembl_id = sub("\\..*$", "", h$ensembl_gene_id),
    gene = h$symbol,
    stringsAsFactors = FALSE
  )
  ensembl <- ensembl[!is.na(ensembl$ensembl_id) & ensembl$ensembl_id != "" & ensembl$gene != "", ]
  refseq_cols <- intersect(c("refseq_accession", "rna_central_ids", "vega_id"), names(h))
  refseq <- list()
  for (cc in refseq_cols) {
    refseq <- c(refseq, lapply(seq_len(nrow(h)), function(i) {
      vals <- split_gene_field(as.character(h[[cc]][i]))
      if (!length(vals)) return(NULL)
      data.frame(refseq = vals, gene = h$symbol[i], stringsAsFactors = FALSE)
    }))
  }
  list(
    ensembl = unique(ensembl),
    refseq = unique(rbindlist_fill(refseq))
  )
}

build_blood_expression_sample_metadata <- function(cfg, source_inventory) {
  if (!file.exists(cfg$compact_metadata)) {
    stop(
      "Blood-expression sample metadata table not found: ", cfg$compact_metadata,
      "\nProvide the accession-validated metadata table with BLOOD_EXPRESSION_METADATA_TABLE, ",
      "or place it in 02_sample_metadata/ before running this step.",
      call. = FALSE
    )
  }
  compact <- read.csv(cfg$compact_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  compact <- compact[compact$series_accession %in% c("GSE111175", "GSE111176", "GSE18123", "GSE25507", "GSE26415", "GSE77103", "GSE140702", "GSE123302", "GSE212645"), ]
  rows <- lapply(seq_len(nrow(compact)), function(i) {
    r <- compact[i, ]
    lab <- label_blood_expression_sample(r)
    platform <- extract_platform_id(r$all_metadata_compact)
    data.frame(
      dataset_id = r$series_accession,
      sample_id = r$sample_accession,
      title = r$title,
      source_name = r$source_name_ch1,
      disease_or_status = r$disease_or_status,
      inferred_group = r$inferred_group,
      platform_id = platform,
      group = lab$group,
      include = lab$include,
      model_role = lab$model_role,
      exclusion_reason = lab$reason,
      condition = lab$condition,
      stringsAsFactors = FALSE
    )
  })
  out <- rbindlist_fill(rows)
  out[order(out$dataset_id, out$sample_id), ]
}

extract_platform_id <- function(x) {
  m <- regmatches(x, regexpr("platform_id=[^ |]+", x))
  if (!length(m) || m == "") return(NA_character_)
  sub("platform_id=", "", m)
}

label_blood_expression_sample <- function(r) {
  dataset <- r$series_accession
  txt <- tolower(paste(r$title, r$source_name_ch1, r$disease_or_status, r$inferred_group, r$all_metadata_compact, sep = " "))
  title <- tolower(r$title)
  if (dataset == "GSE140702") {
    if (!grepl("nt", title, fixed = TRUE)) {
      return(list(group = "exclude", include = FALSE, model_role = NA_character_, reason = "stimulated sample excluded; NT only retained", condition = "stimulated"))
    }
    if (grepl("full autistic|pddnos|asperger", txt)) return(list(group = "ASD", include = TRUE, model_role = "peripheral_blood_primary", reason = "ASD-spectrum NT monocyte", condition = "NT"))
    if (grepl("typical", txt)) return(list(group = "control", include = TRUE, model_role = "peripheral_blood_primary", reason = "Typical NT monocyte", condition = "NT"))
    return(list(group = "exclude", include = FALSE, model_role = NA_character_, reason = "unclear GSE140702 phenotype", condition = "NT"))
  }
  if (dataset == "GSE123302") {
    if (grepl("non-td|nontd", txt)) return(list(group = "exclude", include = FALSE, model_role = NA_character_, reason = "Non-TD comparison group excluded", condition = "cord_blood"))
    if (grepl("cord_blood_td|\\btd\\b|typical", txt)) return(list(group = "control", include = TRUE, model_role = "cord_blood_developmental_sensitivity", reason = "TD cord-blood control", condition = "cord_blood"))
    if (grepl("cord_blood_asd|\\basd\\b|autism", txt)) return(list(group = "ASD", include = TRUE, model_role = "cord_blood_developmental_sensitivity", reason = "ASD cord blood", condition = "cord_blood"))
    return(list(group = "exclude", include = FALSE, model_role = NA_character_, reason = "unclear cord-blood phenotype", condition = "cord_blood"))
  }
  if (dataset == "GSE212645") {
    return(list(group = "exclude", include = FALSE, model_role = NA_character_, reason = "sample-level ASD/control labels not exposed in public compact metadata", condition = NA_character_))
  }
  if (dataset == "GSE26415") {
    if (grepl("control|nonautistic", txt)) return(list(group = "control", include = TRUE, model_role = "peripheral_blood_primary", reason = "nonautistic control", condition = NA_character_))
    if (grepl("autism spectrum", txt)) return(list(group = "ASD", include = TRUE, model_role = "peripheral_blood_primary", reason = "ASD", condition = NA_character_))
  }
  if (grepl("control|typical|healthy|\\btd\\b", txt)) return(list(group = "control", include = TRUE, model_role = "peripheral_blood_primary", reason = "metadata control label", condition = NA_character_))
  if (grepl("autism|\\basd\\b|pdd|asperger", txt)) return(list(group = "ASD", include = TRUE, model_role = "peripheral_blood_primary", reason = "metadata ASD-spectrum label", condition = NA_character_))
  list(group = "exclude", include = FALSE, model_role = NA_character_, reason = "unclear phenotype label", condition = NA_character_)
}

collapse_probe_to_gene <- function(values_dt, probe_map) {
  values_dt$probe_id <- as.character(values_dt$probe_id)
  probe_map$probe_id <- as.character(probe_map$probe_id)
  long <- data.table::melt(data.table::as.data.table(values_dt), id.vars = "probe_id",
                           variable.name = "sample_id", value.name = "expression_value")
  long[, expression_value := suppressWarnings(as.numeric(expression_value))]
  long <- long[is.finite(expression_value)]
  pm <- data.table::as.data.table(probe_map)
  merged <- merge(long, pm, by = "probe_id", allow.cartesian = TRUE)
  merged <- merged[gene != "" & !is.na(gene)]
  merged[, .(
    expression_value = mean(expression_value, na.rm = TRUE),
    feature_count = .N
  ), by = .(gene, sample_id)]
}

process_series_matrix_route <- function(cfg, route, sample_metadata) {
  file <- route$local_path
  sm <- read_series_matrix(file)
  gpl <- route$platform_id
  annot_file <- file.path(cfg$out_dir, "01_source_files", paste0(gpl, ".annot.gz"))
  if (!file.exists(annot_file)) {
    annot_file <- file.path(cfg$out_dir, "01_source_files", paste0(gpl, "_family.soft.gz"))
  }
  probe_map <- read_platform_map(annot_file)
  if (!nrow(probe_map)) {
    return(list(gene_expression = data.frame(), log = data.frame(
      dataset_id = route$dataset_id, route_id = route$route_id, status = "failed",
      reason = paste("no usable platform gene-symbol map for", gpl), samples_total = 0,
      ASD_n = 0, control_n = 0, excluded_samples = 0, genes_with_values = 0,
      stringsAsFactors = FALSE
    )))
  }
  values <- sm$values
  samples <- intersect(names(values), sample_metadata$sample_id)
  values <- values[, c("probe_id", samples), drop = FALSE]
  sample_rows <- sample_metadata[sample_metadata$sample_id %in% samples & sample_metadata$include == TRUE, ]
  sample_rows <- sample_rows[sample_rows$dataset_id == route$dataset_id, ]
  if (!nrow(sample_rows)) {
    return(list(gene_expression = data.frame(), log = data.frame(
      dataset_id = route$dataset_id, route_id = route$route_id, status = "excluded",
      reason = "no included samples after phenotype rules", samples_total = length(samples),
      ASD_n = 0, control_n = 0, excluded_samples = length(samples), genes_with_values = 0,
      stringsAsFactors = FALSE
    )))
  }
  keep <- c("probe_id", sample_rows$sample_id)
  keep <- keep[keep %in% names(values)]
  gene_values <- collapse_probe_to_gene(values[, keep, drop = FALSE], probe_map)
  gene_values <- merge(gene_values, sample_rows[, c("sample_id", "dataset_id", "group", "model_role")],
                       by = "sample_id")
  gene_values$route_id <- route$route_id
  gene_values$platform_id <- gpl
  gene_values$processing_method <- "GEO series matrix collapsed to gene by mean across mapped probes"
  list(
    gene_expression = as.data.frame(gene_values),
    log = data.frame(
      dataset_id = route$dataset_id,
      route_id = route$route_id,
      status = "complete",
      reason = "",
      samples_total = length(samples),
      ASD_n = sum(sample_rows$group == "ASD"),
      control_n = sum(sample_rows$group == "control"),
      excluded_samples = length(samples) - nrow(sample_rows),
      genes_with_values = length(unique(gene_values$gene)),
      platform_id = gpl,
      source_file = file,
      stringsAsFactors = FALSE
    )
  )
}

extract_tar_if_needed <- function(tar_file, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(list.files(out_dir, full.names = TRUE))) utils::untar(tar_file, exdir = out_dir)
  list.files(out_dir, full.names = TRUE, recursive = TRUE)
}

process_gse140702_featurecounts <- function(cfg, sample_metadata) {
  tar_file <- file.path(cfg$out_dir, "01_source_files", "GSE140702_RAW.tar")
  raw_dir <- file.path(cfg$out_dir, "01_source_files", "GSE140702_RAW")
  files <- extract_tar_if_needed(tar_file, raw_dir)
  files <- files[grepl("\\.counts\\.txt\\.gz$", files)]
  sample_ids <- sub("_.*$", "", basename(files))
  included <- sample_metadata[sample_metadata$dataset_id == "GSE140702" & sample_metadata$include == TRUE, ]
  files <- files[sample_ids %in% included$sample_id]
  sample_ids <- sub("_.*$", "", basename(files))
  maps <- hgnc_maps(cfg)
  ens_map <- data.table::as.data.table(maps$ensembl)
  count_list <- lapply(seq_along(files), function(i) {
    dt <- data.table::fread(files[i], data.table = FALSE)
    gene_col <- if ("Geneid" %in% names(dt)) "Geneid" else names(dt)[1]
    count_col <- names(dt)[ncol(dt)]
    data.frame(ensembl_id = sub("\\..*$", "", dt[[gene_col]]), sample_id = sample_ids[i],
               count = suppressWarnings(as.numeric(dt[[count_col]])), stringsAsFactors = FALSE)
  })
  counts <- rbindlist_fill(count_list)
  counts <- merge(data.table::as.data.table(counts), ens_map, by = "ensembl_id")
  counts <- counts[is.finite(count)]
  counts <- counts[, .(count = sum(count)), by = .(gene, sample_id)]
  mat <- data.table::dcast(counts, gene ~ sample_id, value.var = "count", fill = 0)
  rn <- mat$gene
  count_mat <- as.matrix(mat[, -1, drop = FALSE])
  rownames(count_mat) <- rn
  y <- edgeR::DGEList(counts = count_mat)
  keep <- edgeR::filterByExpr(y, group = included$group[match(colnames(count_mat), included$sample_id)])
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  long <- data.table::melt(data.table::as.data.table(logcpm, keep.rownames = "gene"),
                           id.vars = "gene", variable.name = "sample_id",
                           value.name = "expression_value")
  long <- merge(long, included[, c("sample_id", "dataset_id", "group", "model_role")], by = "sample_id")
  long$route_id <- "GSE140702"
  long$platform_id <- "GPL20301"
  long$feature_count <- 1
  long$processing_method <- "edgeR TMM-normalised logCPM from public featureCounts files; NT unstimulated samples only"
  list(gene_expression = as.data.frame(long), log = data.frame(
    dataset_id = "GSE140702", route_id = "GSE140702", status = "complete",
    reason = "", samples_total = length(unique(sample_ids)),
    ASD_n = sum(included$group == "ASD"), control_n = sum(included$group == "control"),
    excluded_samples = nrow(sample_metadata[sample_metadata$dataset_id == "GSE140702", ]) - nrow(included),
    genes_with_values = length(unique(long$gene)), platform_id = "GPL20301",
    source_file = tar_file, stringsAsFactors = FALSE
  ))
}

process_gse77103_agilent <- function(cfg, sample_metadata) {
  tar_file <- file.path(cfg$out_dir, "01_source_files", "GSE77103_RAW.tar")
  raw_dir <- file.path(cfg$out_dir, "01_source_files", "GSE77103_RAW")
  files <- extract_tar_if_needed(tar_file, raw_dir)
  files <- files[grepl("\\.txt\\.gz$", files)]
  sample_ids <- sub("_.*$", "", basename(files))
  included <- sample_metadata[sample_metadata$dataset_id == "GSE77103" & sample_metadata$include == TRUE, ]
  files <- files[sample_ids %in% included$sample_id]
  sample_ids <- sub("_.*$", "", basename(files))
  maps <- hgnc_maps(cfg)
  ref_map <- data.table::as.data.table(maps$refseq)
  parse_one <- function(file, sample_id) {
    lines <- readLines(gzfile(file, open = "rt"), warn = FALSE)
    start <- grep("^FEATURES\t", lines)[1]
    if (is.na(start)) return(data.frame())
    feature_lines <- lines[start:length(lines)]
    feature_lines[1] <- sub("^FEATURES\t", "", feature_lines[1])
    feature_lines[-1] <- sub("^DATA\t", "", feature_lines[-1])
    dt <- data.table::fread(text = paste(feature_lines, collapse = "\n"), data.table = FALSE, fill = TRUE)
    nm <- tolower(names(dt))
    sig_col <- names(dt)[which(nm %in% tolower(c("gProcessedSignal", "gMeanSignal", "rProcessedSignal")))[1]]
    ref_col <- names(dt)[grep("refseq|accession|gb_acc|gene|systematicname", nm)[1]]
    if (is.na(sig_col) || is.na(ref_col)) {
      return(data.frame())
    }
    vals <- data.frame(refseq = as.character(dt[[ref_col]]),
                       sample_id = sample_id,
                       signal = suppressWarnings(as.numeric(dt[[sig_col]])),
                       stringsAsFactors = FALSE)
    vals <- vals[is.finite(vals$signal), ]
    vals <- merge(data.table::as.data.table(vals), ref_map, by = "refseq")
    vals[, .(expression_value = mean(log2(signal + 1), na.rm = TRUE), feature_count = .N), by = .(gene, sample_id)]
  }
  long <- rbindlist_fill(Map(parse_one, files, sample_ids))
  long <- merge(long, included[, c("sample_id", "dataset_id", "group", "model_role")], by = "sample_id")
  long$route_id <- "GSE77103"
  long$platform_id <- "GPL17077"
  long$processing_method <- "log2(gProcessedSignal + 1) from public Agilent one-color signal files mapped through HGNC RefSeq accessions"
  list(gene_expression = as.data.frame(long), log = data.frame(
    dataset_id = "GSE77103", route_id = "GSE77103", status = "complete",
    reason = "", samples_total = length(unique(sample_ids)),
    ASD_n = sum(included$group == "ASD"), control_n = sum(included$group == "control"),
    excluded_samples = nrow(sample_metadata[sample_metadata$dataset_id == "GSE77103", ]) - nrow(included),
    genes_with_values = length(unique(long$gene)), platform_id = "GPL17077",
    source_file = tar_file, stringsAsFactors = FALSE
  ))
}

hedges_g <- function(asd, ctrl) {
  asd <- as.numeric(asd[is.finite(asd)])
  ctrl <- as.numeric(ctrl[is.finite(ctrl)])
  n1 <- length(asd); n0 <- length(ctrl)
  if (n1 < 2 || n0 < 2) return(NULL)
  m1 <- mean(asd); m0 <- mean(ctrl)
  s1 <- stats::sd(asd); s0 <- stats::sd(ctrl)
  if (!is.finite(s1) || !is.finite(s0)) return(NULL)
  sp2 <- ((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2)
  if (!is.finite(sp2) || sp2 <= 0) return(NULL)
  d <- (m1 - m0) / sqrt(sp2)
  J <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- J * d
  v <- (n1 + n0) / (n1 * n0) + (g^2) / (2 * (n1 + n0 - 2))
  data.frame(
    n_asd = n1, n_control = n0,
    ASD_mean = m1, control_mean = m0,
    ASD_sd = s1, control_sd = s0,
    hedges_g = g, variance_g = v, se_g = sqrt(v),
    ci_lower = g - 1.96 * sqrt(v), ci_upper = g + 1.96 * sqrt(v),
    direction = ifelse(g > 0, "ASD_higher", "ASD_lower"),
    stringsAsFactors = FALSE
  )
}

calculate_expression_effect_sizes <- function(values) {
  dt <- data.table::as.data.table(values)
  split_dt <- split(dt, paste(dt$route_id, dt$gene, sep = "||"))
  rows <- lapply(split_dt, function(x) {
    eff <- hedges_g(x$expression_value[x$group == "ASD"], x$expression_value[x$group == "control"])
    if (is.null(eff)) return(NULL)
    data.frame(
      dataset = x$route_id[1],
      study_id = x$dataset_id[1],
      model_role = x$model_role[1],
      platform_id = x$platform_id[1],
      gene = x$gene[1],
      feature_count = sum(x$feature_count, na.rm = TRUE),
      eff,
      notes = x$processing_method[1],
      stringsAsFactors = FALSE
    )
  })
  rbindlist_fill(rows)
}

collapse_within_study_platforms <- function(effects) {
  dt <- data.table::as.data.table(effects)
  multi <- dt[, .N, by = .(study_id, gene)][N > 1]
  if (!nrow(multi)) return(as.data.frame(dt))
  key <- paste(multi$study_id, multi$gene, sep = "||")
  dt[, collapse_key := paste(study_id, gene, sep = "||")]
  keep <- dt[!collapse_key %in% key]
  collapsed <- dt[collapse_key %in% key, {
    w <- 1 / variance_g
    g <- sum(w * hedges_g) / sum(w)
    v <- 1 / sum(w)
    .(
      dataset = study_id[1],
      study_id = study_id[1],
      model_role = model_role[1],
      platform_id = paste(unique(platform_id), collapse = ";"),
      gene = gene[1],
      feature_count = sum(feature_count, na.rm = TRUE),
      n_asd = sum(unique(n_asd)),
      n_control = sum(unique(n_control)),
      ASD_mean = NA_real_, control_mean = NA_real_, ASD_sd = NA_real_, control_sd = NA_real_,
      hedges_g = g,
      variance_g = v,
      se_g = sqrt(v),
      ci_lower = g - 1.96 * sqrt(v),
      ci_upper = g + 1.96 * sqrt(v),
      direction = ifelse(g > 0, "ASD_higher", "ASD_lower"),
      notes = paste("within-study platform strata collapsed by inverse-variance fixed effect:", paste(unique(dataset), collapse = ";"))
    )
  }, by = collapse_key]
  collapsed$collapse_key <- NULL
  keep$collapse_key <- NULL
  as.data.frame(rbindlist_fill(list(keep, collapsed)))
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
    ci_lower = as.numeric(dl$ci.lb),
    ci_upper = as.numeric(dl$ci.ub),
    p_value = as.numeric(dl$pval),
    tau2 = as.numeric(dl$tau2),
    Q = as.numeric(dl$QE),
    I2 = as.numeric(dl$I2),
    mKH_se = as.numeric(mkh$se),
    mKH_ci_lower = as.numeric(mkh$ci.lb),
    mKH_ci_upper = as.numeric(mkh$ci.ub),
    stringsAsFactors = FALSE
  )
}

build_meta_results <- function(effects, model_name) {
  dt <- data.table::as.data.table(effects)
  rows <- dt[, {
    if (.N >= 2) dl_meta_one(hedges_g, variance_g) else data.frame()
  }, by = gene]
  if (nrow(rows)) {
    rows$model <- model_name
    rows$datasets_contributing <- dt[, .(datasets_contributing = paste(unique(study_id), collapse = ";")), by = gene]$datasets_contributing[match(rows$gene, dt[, unique(gene)])]
    rows$DL_interval_excludes_zero <- with(rows, ci_lower > 0 | ci_upper < 0)
    rows$mKH_interval_excludes_zero <- with(rows, mKH_ci_lower > 0 | mKH_ci_upper < 0)
    rows$direction <- ifelse(rows$pooled_g > 0, "ASD_higher", "ASD_lower")
    rows$FDR <- p.adjust(rows$p_value, method = "BH")
    rows$FDR_significant <- rows$FDR < 0.05
    rows$status <- ifelse(rows$FDR_significant & rows$mKH_interval_excludes_zero, "FDR_and_mKH",
                          ifelse(rows$FDR_significant, "FDR_only",
                                 ifelse(rows$mKH_interval_excludes_zero, "mKH_interval_only",
                                        ifelse(rows$DL_interval_excludes_zero, "DL_only", "null"))))
  }
  desc <- dt[, .N, by = gene][N == 1, .(gene)]
  if (nrow(desc)) {
    desc <- merge(desc, dt, by = "gene")
    desc$model <- model_name
    desc$status <- "k1_descriptive"
  }
  list(meta = as.data.frame(rows), descriptive = as.data.frame(desc))
}

summarise_expression_model <- function(meta, desc, effects, model, role) {
  data.frame(
    model = model,
    role = role,
    datasets = paste(sort(unique(effects$study_id)), collapse = ";"),
    genes_meta_analysed = nrow(meta),
    k1_descriptive = nrow(desc),
    DL_nonzero = sum(meta$DL_interval_excludes_zero, na.rm = TRUE),
    FDR_significant = sum(meta$FDR_significant, na.rm = TRUE),
    mKH_interval_supported = sum(meta$mKH_interval_excludes_zero, na.rm = TRUE),
    FDR_mKH_overlap = sum(meta$FDR_significant & meta$mKH_interval_excludes_zero, na.rm = TRUE),
    median_I2 = median(meta$I2, na.rm = TRUE),
    median_tau2 = median(meta$tau2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

validate_meta_with_metafor <- function(effects, meta, model, n_genes = 100) {
  dt <- data.table::as.data.table(effects)
  if (model == "blood_expression_peripheral_primary") {
    dt <- dt[model_role == "peripheral_blood_primary"]
  } else if (model == "blood_expression_plus_cord_blood_sensitivity") {
    dt <- dt[model_role %in% c("peripheral_blood_primary", "cord_blood_developmental_sensitivity")]
  }
  genes <- intersect(meta$gene, dt[, .N, by = gene][N >= 2, gene])
  set.seed(20260518)
  genes <- sample(genes, min(length(genes), n_genes))
  out <- lapply(genes, function(gene_i) {
    x <- dt[gene == gene_i]
    m <- metafor::rma.uni(yi = x$hedges_g, vi = x$variance_g, method = "DL", test = "knha")
    old <- meta[meta$gene == gene_i, ]
    data.frame(
      model = model,
      gene = gene_i,
      branch_pooled_g = old$pooled_g,
      metafor_pooled_g = as.numeric(m$b[1]),
      branch_tau2 = old$tau2,
      metafor_tau2 = m$tau2,
      diff_pooled_g = abs(old$pooled_g - as.numeric(m$b[1])),
      diff_tau2 = abs(old$tau2 - m$tau2),
      within_tolerance = abs(old$pooled_g - as.numeric(m$b[1])) < 1e-8 & abs(old$tau2 - m$tau2) < 1e-8,
      stringsAsFactors = FALSE
    )
  })
  rbindlist_fill(out)
}

compare_model_status <- function(primary, sensitivity, primary_name, sensitivity_name) {
  p <- primary[, c("gene", "status", "FDR_significant", "mKH_interval_excludes_zero", "DL_interval_excludes_zero")]
  s <- sensitivity[, c("gene", "status", "FDR_significant", "mKH_interval_excludes_zero", "DL_interval_excludes_zero")]
  names(p)[-1] <- paste0("primary_", names(p)[-1])
  names(s)[-1] <- paste0("sensitivity_", names(s)[-1])
  m <- merge(p, s, by = "gene", all = TRUE)
  m$primary_model <- primary_name
  m$sensitivity_model <- sensitivity_name
  m$status_transition <- paste(ifelse(is.na(m$primary_status), "absent", m$primary_status),
                               "to",
                               ifelse(is.na(m$sensitivity_status), "absent", m$sensitivity_status))
  m
}
