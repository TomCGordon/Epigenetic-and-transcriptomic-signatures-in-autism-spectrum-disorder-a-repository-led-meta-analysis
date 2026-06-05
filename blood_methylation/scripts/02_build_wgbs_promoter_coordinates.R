#!/usr/bin/env Rscript

# Build hg38 strand-aware promoter coordinates for the harmonised 20,960-gene
# blood/brain methylation universe. These coordinates are used only for the
# GSE140730 WGBS developmental sensitivity analysis.

suppressPackageStartupMessages({
  library(data.table)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
raw_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_RAW_DIR", unset = file.path(package_root, "data_raw")),
                         winslash = "/", mustWork = FALSE)
processed_dir <- normalizePath(Sys.getenv("BLOOD_METHYLATION_PROCESSED_DIR", unset = file.path(package_root, "data_processed")),
                               winslash = "/", mustWork = FALSE)
qc_dir <- normalizePath(file.path(package_root, "qc"), winslash = "/", mustWork = FALSE)
dir.create(file.path(raw_dir, "annotation"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(processed_dir, "annotation"), recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

annotation_path <- file.path(raw_dir, "annotation", "illumina450k_annotation_core.csv")
if (!file.exists(annotation_path)) {
  stop("Missing annotation file: ", annotation_path,
       "\nRun 01_download_blood_sources.R first.", call. = FALSE)
}

promoter_terms <- c("TSS200", "TSS1500", "5'UTR", "1stExon")
anno <- fread(annotation_path, showProgress = FALSE)
needed <- c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")
if (!all(needed %in% names(anno))) stop("Annotation file lacks required columns: ", paste(needed, collapse = ", "))

map_rows <- vector("list", nrow(anno))
for (i in seq_len(nrow(anno))) {
  genes <- strsplit(as.character(anno$UCSC_RefGene_Name[[i]]), ";", fixed = TRUE)[[1]]
  groups <- strsplit(as.character(anno$UCSC_RefGene_Group[[i]]), ";", fixed = TRUE)[[1]]
  genes <- trimws(genes)
  groups <- trimws(groups)
  n <- max(length(genes), length(groups))
  rows <- vector("list", n)
  for (j in seq_len(n)) {
    gene <- if (j <= length(genes)) genes[[j]] else genes[[1]]
    group <- if (j <= length(groups)) groups[[j]] else groups[[1]]
    if (nzchar(gene) && group %in% promoter_terms) {
      rows[[j]] <- data.table(probe = anno$Name[[i]], gene_symbol = gene, promoter_group = group)
    }
  }
  map_rows[[i]] <- rbindlist(rows, fill = TRUE)
}
probe_gene_map <- unique(rbindlist(map_rows, fill = TRUE))
universe <- sort(unique(probe_gene_map$gene_symbol))
fwrite(data.table(gene_symbol = universe), file.path(processed_dir, "annotation", "blood_harmonised_20960_gene_universe.csv"))
fwrite(probe_gene_map, file.path(processed_dir, "annotation", "blood_450k_promoter_probe_gene_map.csv"))

refgene_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
refgene_path <- file.path(raw_dir, "annotation", "UCSC_hg38_refGene.txt.gz")
if (!file.exists(refgene_path) || toupper(Sys.getenv("FORCE_DOWNLOAD", "FALSE")) == "TRUE") {
  message("Downloading UCSC hg38 refGene table")
  utils::download.file(refgene_url, refgene_path, mode = "wb", quiet = FALSE)
}

ref_cols <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
              "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
              "cdsEndStat", "exonFrames")
ref <- fread(refgene_path, header = FALSE, col.names = ref_cols, showProgress = FALSE)
ref <- ref[name2 %in% universe & chrom %chin% paste0("chr", c(1:22, "X", "Y"))]
ref[, tss_coordinate := ifelse(strand == "+", txStart, txEnd)]
ref[, promoter_start := ifelse(strand == "+", pmax(1L, tss_coordinate - 1500L), pmax(1L, tss_coordinate - 200L))]
ref[, promoter_end := ifelse(strand == "+", tss_coordinate + 200L, tss_coordinate + 1500L)]
coord <- unique(ref[, .(
  gene_symbol = name2,
  transcript_id = name,
  chromosome = chrom,
  strand,
  tss_coordinate,
  promoter_start,
  promoter_end,
  genome_build = "hg38 / GRCh38",
  source_annotation = "UCSC hg38 refGene",
  promoter_definition = "strand-aware TSS1500/TSS200 interval: + strand TSS-1500 to TSS+200; - strand TSS-200 to TSS+1500"
)])
setorder(coord, gene_symbol, chromosome, promoter_start, promoter_end, transcript_id)
coord[, promoter_interval_id := paste0(gene_symbol, "_", seq_len(.N)), by = gene_symbol]

out_path <- file.path(processed_dir, "annotation", "GSE140730_hg38_promoter_coordinates_from_UCSC_refGene.csv")
fwrite(coord, out_path)

qc <- data.table(
  metric = c("harmonised_universe_genes", "genes_with_refGene_promoter_coordinate", "promoter_intervals", "source_url"),
  value = c(length(universe), uniqueN(coord$gene_symbol), nrow(coord), refgene_url)
)
fwrite(qc, file.path(qc_dir, "wgbs_promoter_coordinate_QC.csv"))

message("WGBS promoter coordinate table written: ", out_path)
