#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "lib", "placenta_lcl_methylation_functions.R"))

raw_annotation_dir <- file.path(package_root, "data_raw", "annotation")
processed_annotation_dir <- file.path(package_root, "data_processed", "annotation")
qc_dir <- file.path(package_root, "qc")
dir_create(processed_annotation_dir)
dir_create(qc_dir)

promoters_hg18 <- load_refgene_promoters(file.path(raw_annotation_dir, "hg18_refGene.txt.gz"), "hg18")
promoters_hg19 <- load_refgene_promoters(file.path(raw_annotation_dir, "hg19_refGene.txt.gz"), "hg19")
promoters_hg38 <- load_refgene_promoters(file.path(raw_annotation_dir, "hg38_refGene.txt.gz"), "hg38")
write_csv(promoters_hg18, file.path(processed_annotation_dir, "hg18_refGene_promoter_coordinates.csv"))
write_csv(promoters_hg19, file.path(processed_annotation_dir, "hg19_refGene_promoter_coordinates.csv"))
write_csv(promoters_hg38, file.path(processed_annotation_dir, "hg38_refGene_promoter_coordinates.csv"))

gpl8490 <- getGEO("GPL8490", AnnotGPL = FALSE)
hm27_tab <- as.data.table(Table(gpl8490))
hm27_tab[, Symbol := trimws(Symbol)]
hm27_map <- hm27_tab[nzchar(Symbol), .(gene = split_gene_symbols(Symbol)), by = .(probe = ID)]
hm27_map <- unique(hm27_map[nzchar(gene)])
write_csv(hm27_map, file.path(processed_annotation_dir, "GPL8490_HM27_probe_gene_map.csv"))

qc <- data.table(
  resource = c("hg18_refGene_promoters", "hg19_refGene_promoters", "hg38_refGene_promoters",
               "GPL8490_HM27_probe_gene_map"),
  rows = c(nrow(promoters_hg18), nrow(promoters_hg19), nrow(promoters_hg38), nrow(hm27_map)),
  unique_genes = c(uniqueN(promoters_hg18$gene), uniqueN(promoters_hg19$gene),
                   uniqueN(promoters_hg38$gene), uniqueN(hm27_map$gene)),
  role = c("LCL MeDIP hg18 coordinate overlap", "GSE67615 hg19 WGBS coordinate overlap",
           "GSE178203 hg38 WGBS coordinate overlap", "GSE34099 HM27 probe-to-gene aggregation")
)
write_csv(qc, file.path(qc_dir, "placenta_lcl_annotation_QC.csv"))
message("Placenta/LCL promoter maps and HM27 annotation written.")
