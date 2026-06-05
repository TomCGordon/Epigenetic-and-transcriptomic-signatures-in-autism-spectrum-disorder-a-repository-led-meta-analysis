#!/usr/bin/env Rscript

# Custom GPL15314 probe remapping for GSE36315.
#
# Purpose:
#   Construct a public, sequence-based probe-to-gene annotation for the
#   Arraystar Human LncRNA Microarray V2.0 (GPL15314) using public probe
#   sequences and GENCODE GRCh37 transcript sequences. This is intended as a
#   sensitivity check only; it does not modify the primary brain
#   expression analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_file)
work_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(work_dir, "..", ".."), winslash = "/", mustWork = TRUE)

brain_source_dir <- file.path(package_root, "pipelines", "brain_expression", "01_source_files")
optional_annotation_dir <- Sys.getenv("GSE36315_OPTIONAL_ANNOTATION_SOURCE", unset = "")

dir.create(file.path(work_dir, "01_source"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "02_probe_sequences"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "03_gencode_mapping"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "04_recoverability_assessment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "05_reports"), recursive = TRUE, showWarnings = FALSE)

message("Working directory: ", work_dir)

copy_if_available <- function(from, to) {
  if (file.exists(from) && !file.exists(to)) file.copy(from, to, overwrite = FALSE)
}

copy_if_available(file.path(brain_source_dir, "GPL15314_family.soft.gz"),
                  file.path(work_dir, "01_source/GPL15314_family.soft.gz"))
copy_if_available(file.path(brain_source_dir, "GSE36315_series_matrix.txt.gz"),
                  file.path(work_dir, "01_source/GSE36315_series_matrix.txt.gz"))
if (nzchar(optional_annotation_dir)) {
  copy_if_available(file.path(optional_annotation_dir, "GPL15314_AnnoProbe_pipe.csv"),
                    file.path(work_dir, "01_source/GPL15314_AnnoProbe_pipe.csv"))
}

gencode_urls <- c(
  pc = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz",
  lncRNA = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz"
)
gencode_fastas <- file.path(work_dir, "01_source", basename(gencode_urls))
for (i in seq_along(gencode_urls)) {
  if (!file.exists(gencode_fastas[[i]])) {
    message("Downloading ", names(gencode_urls)[[i]], " GENCODE v19 GRCh37 transcript FASTA...")
    download.file(gencode_urls[[i]], gencode_fastas[[i]], mode = "wb", quiet = FALSE)
  }
}

read_gpl_probe_table <- function(soft_gz) {
  message("Reading GPL15314 SOFT platform file...")
  lines <- readLines(gzfile(soft_gz), warn = FALSE)
  header_i <- grep("^!platform_table_begin", lines)
  end_i <- grep("^!platform_table_end", lines)
  if (length(header_i) != 1 || length(end_i) != 1) {
    stop("Could not locate unique platform table in ", soft_gz)
  }
  tab_txt <- lines[(header_i + 1):(end_i - 1)]
  dt <- fread(text = paste(tab_txt, collapse = "\n"), sep = "\t", quote = "")
  setnames(dt, make.names(trimws(names(dt)), unique = TRUE))
  if (!all(c("ID", "SEQUENCE") %in% names(dt))) {
    stop("GPL table lacks ID/SEQUENCE columns. Columns: ", paste(names(dt), collapse = ", "))
  }
  dt[, probe_id := as.character(ID)]
  dt[, sequence := toupper(gsub("[^ACGT]", "", as.character(SEQUENCE)))]
  dt <- dt[nchar(sequence) >= 40]
  unique(dt[, .(probe_id, sequence)])
}

read_series_probe_ids <- function(series_gz) {
  message("Reading GSE36315 series matrix probe IDs...")
  lines <- readLines(gzfile(series_gz), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) stop("Could not locate series matrix table.")
  mat <- fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"), sep = "\t", quote = "")
  names(mat)[1] <- "probe_id"
  unique(as.character(mat$probe_id))
}

parse_gencode_headers <- function(x) {
  hdr <- names(x)
  pipe_fields <- strsplit(hdr, "\\|", fixed = FALSE)
  get_pos <- function(i) vapply(pipe_fields, function(z) if (length(z) >= i) z[[i]] else NA_character_, character(1))
  get_label <- function(label, default = NA_character_) {
    pattern <- paste0(".*\\|", label, ":([^|]+).*")
    out <- sub(pattern, "\\1", hdr, perl = TRUE)
    out[out == hdr] <- default
    out
  }
  data.table(
    transcript_id = sub("\\|.*$", "", hdr),
    gene_id = fifelse(!is.na(get_label("gene")), get_label("gene"), get_pos(2)),
    gene_symbol = fifelse(!is.na(get_label("gene_symbol")), get_label("gene_symbol"), get_pos(5)),
    gene_type = fifelse(!is.na(get_label("gene_type")), get_label("gene_type"), get_pos(6)),
    transcript_type = fifelse(!is.na(get_label("transcript_type")), get_label("transcript_type"), get_pos(6))
  )
}

gpl <- read_gpl_probe_table(file.path(work_dir, "01_source/GPL15314_family.soft.gz"))
fwrite(gpl, file.path(work_dir, "02_probe_sequences/GPL15314_public_probe_sequences.csv.gz"))

series_probes <- read_series_probe_ids(file.path(work_dir, "01_source/GSE36315_series_matrix.txt.gz"))
series_probe_dt <- data.table(probe_id = series_probes)
fwrite(series_probe_dt, file.path(work_dir, "02_probe_sequences/GSE36315_series_matrix_probe_ids.csv"))

message("GPL probes with sequence: ", nrow(gpl))
message("GSE36315 series probes: ", length(series_probes))

message("Reading GENCODE transcript FASTA files...")
tx <- do.call(c, lapply(gencode_fastas, readDNAStringSet))
tx_meta <- parse_gencode_headers(tx)
names(tx) <- tx_meta$transcript_id
message("GENCODE transcripts: ", length(tx))
message("GENCODE unique gene symbols: ", uniqueN(tx_meta$gene_symbol[!is.na(tx_meta$gene_symbol) & tx_meta$gene_symbol != ""]))

map_probes <- function(probe_dt, transcript_set, transcript_meta, max_mismatch = 0, chunk_size = 2000) {
  out <- vector("list", ceiling(nrow(probe_dt) / chunk_size))
  chunk_id <- 0L
  for (start in seq(1, nrow(probe_dt), by = chunk_size)) {
    chunk_id <- chunk_id + 1L
    end <- min(start + chunk_size - 1, nrow(probe_dt))
    message("Mapping probes ", start, "-", end, " with max.mismatch=", max_mismatch)
    chunk <- probe_dt[start:end]
    patterns <- DNAStringSet(chunk$sequence)
    names(patterns) <- chunk$probe_id
    rc_patterns <- reverseComplement(patterns)
    names(rc_patterns) <- paste0(names(patterns), "__RC")

    hits_fwd <- vmatchPattern(patterns, transcript_set, max.mismatch = max_mismatch, fixed = FALSE)
    hits_rev <- vmatchPattern(rc_patterns, transcript_set, max.mismatch = max_mismatch, fixed = FALSE)

    to_dt <- function(hit_list, orientation) {
      nonempty <- which(elementNROWS(hit_list) > 0)
      if (!length(nonempty)) return(data.table())
      rbindlist(lapply(nonempty, function(i) {
        probe_name <- names(hit_list)[i]
        if (orientation == "reverse_complement") {
          probe_name <- sub("__RC$", "", probe_name)
        }
        data.table(
          probe_id = probe_name,
          transcript_id = names(hit_list)[i],
          orientation = orientation,
          hit_count_within_transcript = length(hit_list[[i]])
        )
      }))
    }

    dt <- rbindlist(list(
      to_dt(hits_fwd, "as_reported"),
      to_dt(hits_rev, "reverse_complement")
    ), fill = TRUE)
    if (nrow(dt)) {
      dt <- merge(dt, transcript_meta, by = "transcript_id", all.x = TRUE)
      dt[, max_mismatch := max_mismatch]
    }
    out[[chunk_id]] <- dt
    gc()
  }
  rbindlist(out, fill = TRUE)
}

# First attempt: exact transcript matches. If this is sparse, try one mismatch
# only for probes not mapped exactly. This follows the published precedent of
# aligning Arraystar V2.0 60-mer probes to transcript databases while avoiding
# unnecessary multi-match inflation.
exact_hits <- map_probes(gpl, tx, tx_meta, max_mismatch = 0, chunk_size = 2000)
fwrite(exact_hits, file.path(work_dir, "03_gencode_mapping/GPL15314_GENCODEv19_exact_transcript_hits.csv.gz"))

exact_mapped <- unique(exact_hits$probe_id)
unmapped_after_exact <- gpl[!probe_id %in% exact_mapped]
message("Exact mapped probes: ", length(exact_mapped), " / ", nrow(gpl))

mismatch_hits <- data.table()
if (nrow(unmapped_after_exact)) {
  # One-mismatch matching is substantially slower; for a first custom annotation pass,
  # apply it only to GSE36315 probes not mapped exactly.
  gse_unmapped <- unmapped_after_exact[probe_id %in% series_probes]
  message("One-mismatch mapping will be attempted for GSE36315 unmapped probes only: ", nrow(gse_unmapped))
  if (nrow(gse_unmapped)) {
    mismatch_hits <- map_probes(gse_unmapped, tx, tx_meta, max_mismatch = 1, chunk_size = 500)
    fwrite(mismatch_hits, file.path(work_dir, "03_gencode_mapping/GSE36315_GENCODEv19_one_mismatch_transcript_hits.csv.gz"))
  }
}

all_hits <- rbindlist(list(exact_hits, mismatch_hits), fill = TRUE)
all_hits <- all_hits[!is.na(gene_symbol) & gene_symbol != ""]
setorder(all_hits, probe_id, max_mismatch, gene_symbol)

summarise_probe_map <- function(hits) {
  if (!nrow(hits)) {
    return(data.table(
      probe_id = character(), mapping_tier = character(), n_transcripts = integer(),
      n_genes = integer(), gene_symbols = character(), gene_types = character(),
      orientations = character()
    ))
  }
  hits[, .(
    mapping_tier = ifelse(min(max_mismatch, na.rm = TRUE) == 0, "exact", "one_mismatch"),
    n_transcripts = uniqueN(transcript_id),
    n_genes = uniqueN(gene_symbol),
    gene_symbols = paste(sort(unique(gene_symbol)), collapse = "|"),
    gene_types = paste(sort(unique(gene_type)), collapse = "|"),
    orientations = paste(sort(unique(orientation)), collapse = "|")
  ), by = probe_id]
}

probe_map <- summarise_probe_map(all_hits)
probe_map[, mapping_class := fifelse(n_genes == 1, "unique_gene",
                                     fifelse(n_genes > 1, "multi_gene", "unmapped"))]
fwrite(probe_map, file.path(work_dir, "03_gencode_mapping/GPL15314_GENCODEv19_probe_to_gene_map.csv.gz"))

gse_probe_map <- merge(series_probe_dt, probe_map, by = "probe_id", all.x = TRUE)
gse_probe_map[is.na(mapping_class), mapping_class := "unmapped"]
gse_probe_map[is.na(n_genes), `:=`(n_genes = 0L, n_transcripts = 0L)]
fwrite(gse_probe_map, file.path(work_dir, "04_recoverability_assessment/GSE36315_GENCODEv19_probe_mapping_assessment.csv"))

annoprobe_path <- file.path(work_dir, "01_source/GPL15314_AnnoProbe_pipe.csv")
annoprobe_overlap <- data.table()
if (file.exists(annoprobe_path)) {
  ap <- fread(annoprobe_path)
  setnames(ap, names(ap), tolower(names(ap)))
  if (all(c("probe_id", "symbol") %in% names(ap))) {
    ap[, symbol := toupper(trimws(symbol))]
    annoprobe_overlap <- merge(gse_probe_map, ap, by = "probe_id", all.x = TRUE, suffixes = c("_gencode", "_annoprobe"))
    annoprobe_overlap[, gencode_symbol_set := toupper(gene_symbols)]
    annoprobe_overlap[, annoprobe_in_gencode_set := !is.na(symbol) & grepl(paste0("(^|\\|)", symbol, "(\\||$)"), gencode_symbol_set)]
    fwrite(annoprobe_overlap, file.path(work_dir, "04_recoverability_assessment/GSE36315_GENCODEv19_vs_AnnoProbe_overlap.csv"))
  }
}

summary_rows <- rbindlist(list(
  data.table(scope = "GPL15314_all_public_probes", metric = "probes_with_public_sequence", value = nrow(gpl)),
  data.table(scope = "GPL15314_all_public_probes", metric = "exact_mapped_probe_ids", value = length(exact_mapped)),
  data.table(scope = "GPL15314_all_public_probes", metric = "unique_gene_probe_ids_any_tier", value = nrow(probe_map[mapping_class == "unique_gene"])),
  data.table(scope = "GPL15314_all_public_probes", metric = "multi_gene_probe_ids_any_tier", value = nrow(probe_map[mapping_class == "multi_gene"])),
  data.table(scope = "GPL15314_all_public_probes", metric = "unique_gene_symbols_any_tier", value = uniqueN(unlist(strsplit(probe_map[mapping_class != "unmapped"]$gene_symbols, "\\|", fixed = FALSE)))),
  data.table(scope = "GSE36315_series_matrix", metric = "series_probe_ids", value = length(series_probes)),
  data.table(scope = "GSE36315_series_matrix", metric = "unique_gene_probe_ids_any_tier", value = nrow(gse_probe_map[mapping_class == "unique_gene"])),
  data.table(scope = "GSE36315_series_matrix", metric = "multi_gene_probe_ids_any_tier", value = nrow(gse_probe_map[mapping_class == "multi_gene"])),
  data.table(scope = "GSE36315_series_matrix", metric = "unmapped_probe_ids", value = nrow(gse_probe_map[mapping_class == "unmapped"])),
  data.table(scope = "GSE36315_series_matrix", metric = "unique_gene_symbols_any_tier", value = uniqueN(unlist(strsplit(gse_probe_map[mapping_class != "unmapped"]$gene_symbols, "\\|", fixed = FALSE))))
))
fwrite(summary_rows, file.path(work_dir, "04_recoverability_assessment/GSE36315_custom_remapping_summary.csv"))

report <- c(
  "# GSE36315 / GPL15314 Custom Probe Remapping Report",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  "- GPL15314 public GEO SOFT platform table, including probe IDs and 60-mer sequences.",
  "- GSE36315 public GEO series matrix probe IDs.",
  "- GENCODE v19 GRCh37 transcript FASTA.",
  "- Optional comparison against the local AnnoProbe GPL15314 probe-symbol map, if present.",
  "",
  "## Mapping Strategy",
  "",
  "- Probe sequences were matched to GENCODE transcript sequences in R using Biostrings.",
  "- Both the reported probe sequence and reverse complement were tested.",
  "- Exact transcript matches were prioritised.",
  "- One-mismatch matching was attempted only for GSE36315 probes not recovered by exact matching.",
  "- Probe-to-gene assignments were labelled unique-gene or multi-gene, rather than forcing ambiguous probes into one symbol.",
  "",
  "## Key Counts",
  "",
  paste(capture.output(print(summary_rows)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "This remapping provides a public sequence-derived annotation route for GPL15314. It is more defensible than using an undocumented proprietary Arraystar map, but it remains a custom annotation and should be labelled as a sensitivity analysis unless independently validated against an official manufacturer or author-provided full platform annotation.",
  "",
  "Recommended use: custom-annotated GSE36315 brain-expression sensitivity analysis, not primary grouped-brain expression model inclusion."
)
writeLines(report, file.path(work_dir, "05_reports/GSE36315_custom_probe_remapping_report.md"))

message("Done. Summary written to: ", file.path(work_dir, "04_recoverability_assessment/GSE36315_custom_remapping_summary.csv"))
