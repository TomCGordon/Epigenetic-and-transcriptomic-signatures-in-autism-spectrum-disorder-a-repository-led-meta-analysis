#!/usr/bin/env Rscript

# GPL15314 public sequence remapping with Rbowtie.
#
# This script attempts to recover a broad probe-to-gene annotation for
# GSE36315/GPL15314 by aligning public 60-mer Arraystar probe sequences against
# public GENCODE v19 GRCh37 transcript sequences. It is intended for an
# explicitly labelled sensitivity analysis, not for primary-model
# inclusion.

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(Rbowtie)
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

for (d in c("01_source", "02_probe_sequences", "03_gencode_mapping", "04_recoverability_assessment", "05_reports")) {
  dir.create(file.path(work_dir, d), recursive = TRUE, showWarnings = FALSE)
}

copy_if_available <- function(from, to) {
  if (file.exists(from) && !file.exists(to)) invisible(file.copy(from, to, overwrite = FALSE))
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
    message("Downloading ", names(gencode_urls)[[i]], " GENCODE v19 transcript FASTA...")
    download.file(gencode_urls[[i]], gencode_fastas[[i]], mode = "wb")
  }
}

read_gpl_probe_table <- function(soft_gz) {
  lines <- readLines(gzfile(soft_gz), warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) stop("Could not locate GPL platform table.")
  dt <- fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"), sep = "\t", quote = "")
  setnames(dt, make.names(trimws(names(dt)), unique = TRUE))
  if (!all(c("ID", "SEQUENCE") %in% names(dt))) stop("GPL table lacks ID/SEQUENCE columns.")
  dt[, probe_id := as.character(ID)]
  dt[, sequence := toupper(gsub("[^ACGT]", "", as.character(SEQUENCE)))]
  unique(dt[nchar(sequence) >= 40, .(probe_id, sequence)])
}

read_series_probe_ids <- function(series_gz) {
  lines <- readLines(gzfile(series_gz), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) stop("Could not locate series matrix table.")
  dt <- fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"), sep = "\t", quote = "")
  names(dt)[1] <- "probe_id"
  unique(gsub('^"+|"+$', "", as.character(dt$probe_id)))
}

parse_gencode_headers <- function(x) {
  hdr <- names(x)
  fields <- strsplit(hdr, "\\|")
  get_pos <- function(i) vapply(fields, function(z) if (length(z) >= i) z[[i]] else NA_character_, character(1))
  data.table(
    transcript_id = sub("\\|.*$", "", hdr),
    gene_id = get_pos(2),
    transcript_name = get_pos(5),
    gene_symbol = get_pos(6),
    gene_type = NA_character_
  )
}

message("Extracting GPL15314 public probe sequences...")
gpl <- read_gpl_probe_table(file.path(work_dir, "01_source/GPL15314_family.soft.gz"))
fwrite(gpl, file.path(work_dir, "02_probe_sequences/GPL15314_public_probe_sequences.csv.gz"))
probe_fasta <- file.path(work_dir, "02_probe_sequences/GPL15314_public_probe_sequences.fa")
probe_set <- DNAStringSet(gpl$sequence)
names(probe_set) <- gpl$probe_id
writeXStringSet(probe_set, probe_fasta)

series_probes <- read_series_probe_ids(file.path(work_dir, "01_source/GSE36315_series_matrix.txt.gz"))
fwrite(data.table(probe_id = series_probes), file.path(work_dir, "02_probe_sequences/GSE36315_series_matrix_probe_ids.csv"))

message("Reading GENCODE v19 transcript FASTA files...")
tx <- do.call(c, lapply(gencode_fastas, readDNAStringSet))
tx_meta <- parse_gencode_headers(tx)
tx_meta <- tx_meta[!is.na(transcript_id) & transcript_id != "" & !is.na(gene_symbol) & gene_symbol != ""]
tx <- tx[match(tx_meta$transcript_id, sub("\\|.*$", "", names(tx)))]
names(tx) <- tx_meta$transcript_id
tx <- tx[!duplicated(names(tx))]
tx_meta <- tx_meta[!duplicated(transcript_id)]

tx_fasta <- file.path(work_dir, "03_gencode_mapping/gencode_v19_pc_lncRNA_transcripts_simplified.fa")
if (!file.exists(tx_fasta)) writeXStringSet(tx, tx_fasta)
fwrite(tx_meta, file.path(work_dir, "03_gencode_mapping/gencode_v19_transcript_metadata.csv.gz"))

index_dir <- file.path(work_dir, "03_gencode_mapping/bowtie_index")
dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
index_prefix <- file.path(index_dir, "gencode_v19_pc_lncRNA")
if (!file.exists(paste0(index_prefix, ".1.ebwt"))) {
  message("Building Rbowtie transcriptome index...")
  bowtie_build(tx_fasta, outdir = index_dir, prefix = "gencode_v19_pc_lncRNA", force = TRUE)
}

alignment_file <- file.path(work_dir, "03_gencode_mapping/GPL15314_to_GENCODEv19_bowtie_v1_beststrata.txt")
message("Aligning GPL15314 probes to GENCODE v19 transcripts with <=1 mismatch...")
if (!file.exists(alignment_file) || file.info(alignment_file)$size == 0) {
  bowtie(
    sequences = probe_fasta,
    index = index_prefix,
    f = TRUE, v = 1, a = TRUE, best = TRUE, strata = TRUE,
    p = max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L)),
    type = "single",
    outfile = alignment_file,
    force = TRUE
  )
} else {
  message("Existing bowtie alignment found; reusing ", alignment_file)
}

message("Parsing bowtie output...")
if (file.info(alignment_file)$size == 0) {
  hits <- data.table()
} else {
  hits <- fread(alignment_file, header = FALSE, sep = "\t", fill = TRUE)
  setnames(hits, c("probe_id", "strand", "transcript_id", "offset0", "aligned_sequence", "quality", "other_hits", "mismatch_string")[seq_along(names(hits))])
  if (!"mismatch_string" %in% names(hits)) hits[, mismatch_string := ""]
  hits[, mismatch_count := fifelse(is.na(mismatch_string) | mismatch_string == "", 0L, lengths(strsplit(mismatch_string, ",")))]
  hits <- merge(hits, tx_meta, by = "transcript_id", all.x = TRUE)
}
fwrite(hits, file.path(work_dir, "03_gencode_mapping/GPL15314_GENCODEv19_bowtie_transcript_hits.csv.gz"))

probe_map <- hits[!is.na(gene_symbol) & gene_symbol != "", .(
  n_transcripts = uniqueN(transcript_id),
  n_genes = uniqueN(gene_symbol),
  min_mismatch = min(mismatch_count, na.rm = TRUE),
  gene_symbols = paste(sort(unique(gene_symbol)), collapse = "|"),
  gene_types = paste(sort(unique(gene_type)), collapse = "|"),
  strands = paste(sort(unique(strand)), collapse = "|")
), by = probe_id]
probe_map[, mapping_class := fifelse(n_genes == 1, "unique_gene", "multi_gene")]

all_probe_map <- merge(gpl[, .(probe_id)], probe_map, by = "probe_id", all.x = TRUE)
all_probe_map[is.na(mapping_class), mapping_class := "unmapped"]
all_probe_map[is.na(n_genes), `:=`(n_genes = 0L, n_transcripts = 0L)]
fwrite(all_probe_map, file.path(work_dir, "03_gencode_mapping/GPL15314_GENCODEv19_bowtie_probe_to_gene_map.csv.gz"))

gse_probe_map <- merge(data.table(probe_id = series_probes), all_probe_map, by = "probe_id", all.x = TRUE)
gse_probe_map[is.na(mapping_class), mapping_class := "unmapped"]
gse_probe_map[is.na(n_genes), `:=`(n_genes = 0L, n_transcripts = 0L)]
fwrite(gse_probe_map, file.path(work_dir, "04_recoverability_assessment/GSE36315_GENCODEv19_bowtie_probe_mapping_assessment.csv"))

annoprobe_path <- file.path(work_dir, "01_source/GPL15314_AnnoProbe_pipe.csv")
if (file.exists(annoprobe_path)) {
  ap <- fread(annoprobe_path)
  setnames(ap, names(ap), tolower(names(ap)))
  if (all(c("probe_id", "symbol") %in% names(ap))) {
    ap[, symbol := toupper(trimws(symbol))]
    ap <- unique(ap[!is.na(symbol) & symbol != "", .(probe_id, annoprobe_symbol = symbol)])
    overlap <- merge(gse_probe_map, ap, by = "probe_id", all.x = TRUE)
    overlap[, gencode_symbol_set := toupper(fifelse(is.na(gene_symbols), "", gene_symbols))]
    overlap[, annoprobe_in_gencode_set := mapply(
      function(sym, set) !is.na(sym) && nzchar(sym) && grepl(paste0("(^|\\|)", sym, "(\\||$)"), set),
      annoprobe_symbol, gencode_symbol_set
    )]
    fwrite(overlap, file.path(work_dir, "04_recoverability_assessment/GSE36315_GENCODEv19_bowtie_vs_AnnoProbe_overlap.csv"))
  }
}

split_symbols <- function(x) unique(unlist(strsplit(x[!is.na(x) & x != ""], "\\|")))

summary_rows <- rbindlist(list(
  data.table(scope = "GPL15314_all_public_probes", metric = "probes_with_public_sequence", value = nrow(gpl)),
  data.table(scope = "GPL15314_all_public_probes", metric = "mapped_probe_ids_any", value = nrow(all_probe_map[mapping_class != "unmapped"])),
  data.table(scope = "GPL15314_all_public_probes", metric = "unique_gene_probe_ids", value = nrow(all_probe_map[mapping_class == "unique_gene"])),
  data.table(scope = "GPL15314_all_public_probes", metric = "multi_gene_probe_ids", value = nrow(all_probe_map[mapping_class == "multi_gene"])),
  data.table(scope = "GPL15314_all_public_probes", metric = "unique_gene_symbols", value = length(split_symbols(all_probe_map[mapping_class != "unmapped"]$gene_symbols))),
  data.table(scope = "GSE36315_series_matrix", metric = "series_probe_ids", value = length(series_probes)),
  data.table(scope = "GSE36315_series_matrix", metric = "mapped_probe_ids_any", value = nrow(gse_probe_map[mapping_class != "unmapped"])),
  data.table(scope = "GSE36315_series_matrix", metric = "unique_gene_probe_ids", value = nrow(gse_probe_map[mapping_class == "unique_gene"])),
  data.table(scope = "GSE36315_series_matrix", metric = "multi_gene_probe_ids", value = nrow(gse_probe_map[mapping_class == "multi_gene"])),
  data.table(scope = "GSE36315_series_matrix", metric = "unmapped_probe_ids", value = nrow(gse_probe_map[mapping_class == "unmapped"])),
  data.table(scope = "GSE36315_series_matrix", metric = "unique_gene_symbols", value = length(split_symbols(gse_probe_map[mapping_class != "unmapped"]$gene_symbols)))
))
fwrite(summary_rows, file.path(work_dir, "04_recoverability_assessment/GSE36315_custom_bowtie_remapping_summary.csv"))

report <- c(
  "# GSE36315 / GPL15314 Custom Probe Remapping Report",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## What Was Done",
  "",
  "GPL15314 public 60-mer probe sequences were aligned to public GENCODE v19 GRCh37 protein-coding and lncRNA transcript sequences using Rbowtie. Alignments allowed up to one mismatch and reported best-stratum hits. Both forward and reverse-complement matches are handled by Bowtie's default strand-aware search.",
  "",
  "## Key Counts",
  "",
  paste(capture.output(print(summary_rows)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "This is a public sequence-derived annotation route. It is substantially more transparent than an undocumented proprietary platform map and does not depend on AnnoProbe for the mapping itself. However, because it is a custom remapping to GENCODE v19 rather than an official Arraystar full annotation, GSE36315 should still be treated as a custom-annotated sensitivity dataset unless an official manufacturer or author annotation is obtained.",
  "",
  "Recommended next step: if the mapped coverage is adequate, run GSE36315 as a clearly labelled brain-expression sensitivity route, ideally region-specific because the dataset contains paired prefrontal cortex and cerebellum samples from the same donors."
)
writeLines(report, file.path(work_dir, "05_reports/GSE36315_custom_bowtie_probe_remapping_report.md"))

message("Finished custom remapping. Summary: ", file.path(work_dir, "04_recoverability_assessment/GSE36315_custom_bowtie_remapping_summary.csv"))
