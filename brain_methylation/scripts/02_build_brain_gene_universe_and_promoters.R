#!/usr/bin/env Rscript

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
source(file.path(script_dir, "lib", "brain_methylation_functions.R"))

raw_annotation <- file.path(package_root, "data_raw", "annotation", "illumina450k_annotation_core.csv")
processed_annotation <- file.path(package_root, "data_processed", "annotation")
qc_dir <- file.path(package_root, "qc")
dir_create(processed_annotation)
dir_create(qc_dir)

promoter_terms <- c("TSS200", "TSS1500", "5'UTR", "5UTR", "1stExon")
ann <- fread(raw_annotation)
setnames(ann, old = intersect(names(ann), c("Name", "IlmnID")), new = "probe")
if (!"probe" %in% names(ann)) stop("Annotation file must contain Name or IlmnID probe column.")
if (!all(c("UCSC_RefGene_Name", "UCSC_RefGene_Group") %in% names(ann))) {
  stop("Annotation file must contain UCSC_RefGene_Name and UCSC_RefGene_Group.")
}

probe_map <- ann[!is.na(UCSC_RefGene_Name) & nzchar(UCSC_RefGene_Name) &
                   !is.na(UCSC_RefGene_Group) & nzchar(UCSC_RefGene_Group),
                 .(probe, gene_names = UCSC_RefGene_Name, groups = UCSC_RefGene_Group)]
probe_map <- probe_map[, {
  genes <- unlist(strsplit(gene_names, ";", fixed = TRUE))
  groups <- unlist(strsplit(groups, ";", fixed = TRUE))
  len <- max(length(genes), length(groups))
  data.table(gene = rep(genes, length.out = len), group = rep(groups, length.out = len))
}, by = probe]
probe_map[, `:=`(gene = trimws(gene), group = trimws(group))]
probe_map <- probe_map[nzchar(gene) & group %in% promoter_terms]
probe_map <- unique(probe_map[, .(probe, gene)])
gene_universe <- sort(unique(probe_map$gene))

fwrite(data.table(gene = gene_universe), file.path(processed_annotation, "brain_harmonised_20960_gene_universe.csv"))
fwrite(probe_map, file.path(processed_annotation, "brain_450k_promoter_probe_gene_map.csv"))

refgene_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
refgene_gz <- file.path(package_root, "data_raw", "annotation", "UCSC_hg38_refGene.txt.gz")
if (!file.exists(refgene_gz)) {
  download.file(refgene_url, refgene_gz, mode = "wb", quiet = TRUE)
}
refgene <- fread(refgene_gz, header = FALSE)
setnames(refgene, c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                    "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
                    "cdsEndStat", "exonFrames"))
refgene <- refgene[name2 %in% gene_universe & grepl("^chr[0-9XYM]+$", chrom)]
refgene[, tss := fifelse(strand == "+", txStart + 1L, txEnd)]
refgene[, start := fifelse(strand == "+", pmax(1L, tss - 1500L), pmax(1L, tss - 200L))]
refgene[, end := fifelse(strand == "+", tss + 200L, tss + 1500L)]
promoters <- unique(refgene[, .(chr = chrom, start, end, gene = name2, strand)])
setorder(promoters, chr, start, end, gene)
fwrite(promoters, file.path(processed_annotation, "brain_hg38_refGene_promoter_coordinates.csv"))

qc <- data.table(
  metric = c("harmonised_array_promoter_gene_universe", "promoter_probe_gene_pairs",
             "genes_with_hg38_refGene_promoter_interval", "promoter_intervals_total", "genome_build"),
  value = c(length(gene_universe), nrow(probe_map), uniqueN(promoters$gene), nrow(promoters), "hg38 / GRCh38")
)
fwrite(qc, file.path(qc_dir, "brain_gene_universe_and_promoter_QC.csv"))
message("Brain gene universe and promoter coordinate files written.")
