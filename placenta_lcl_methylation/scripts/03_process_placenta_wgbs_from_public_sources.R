#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})
options(timeout = max(3600, getOption("timeout")))

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "lib", "placenta_lcl_methylation_functions.R"))

raw_dir <- file.path(package_root, "data_raw", "placenta")
processed_dir <- file.path(package_root, "data_processed", "placenta")
annotation_dir <- file.path(package_root, "data_processed", "annotation")
qc_dir <- file.path(package_root, "qc")
dir_create(raw_dir)
dir_create(processed_dir)
dir_create(qc_dir)

manifest <- fread(file.path(qc_dir, "placenta_GEO_sample_manifest.csv"))
manifest <- manifest[include == TRUE & group %in% c("ASD", "Control")]

load_promoters <- function(build) {
  p <- fread(file.path(annotation_dir, paste0(build, "_refGene_promoter_coordinates.csv")))
  p <- p[!is.na(chr) & is.finite(start) & is.finite(end)]
  p[, `:=`(start = as.integer(start), end = as.integer(end))]
  setkey(p, chr, start, end)
  p
}

promoters_by_build <- list(
  hg38 = load_promoters("hg38"),
  hg19 = load_promoters("hg19")
)
genes_by_build <- lapply(promoters_by_build, function(p) sort(unique(p$gene)))

process_cpg_report <- function(source_file, promoters) {
  cpg <- fread(source_file, header = FALSE, select = c(1, 2, 4, 5),
               col.names = c("chr", "pos", "methylated", "unmethylated"),
               showProgress = FALSE)
  cpg[, coverage := methylated + unmethylated]
  cpg <- cpg[coverage > 0 & grepl("^chr", chr)]
  cpg[, methylation := methylated / coverage]
  cpg[, `:=`(start = as.integer(pos), end = as.integer(pos))]
  cpg <- cpg[is.finite(start) & is.finite(methylation), .(chr, start, end, methylation)]
  setkey(cpg, chr, start, end)
  hits <- foverlaps(cpg, promoters, nomatch = 0L)
  out <- hits[, .(mean_methylation = mean(methylation, na.rm = TRUE),
                  feature_count = .N), by = gene]
  rm(cpg, hits)
  gc()
  out
}

process_bed_tarball <- function(source_file, promoters) {
  tmp <- tempfile("gse67615_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(source_file, exdir = tmp)
  bed_files <- list.files(tmp, pattern = "\\.BED$", full.names = TRUE, recursive = TRUE)
  parts <- lapply(bed_files, function(bed) {
    dt <- fread(bed, header = FALSE, skip = 1, select = c(1, 2, 3, 4),
                col.names = c("chr", "start0", "end0", "methylation"),
                showProgress = FALSE)
    dt <- dt[grepl("^chr", chr) & is.finite(methylation)]
    dt[, `:=`(start = as.integer(start0 + 1L), end = as.integer(start0 + 1L))]
    dt[, .(chr, start, end, methylation)]
  })
  cpg <- rbindlist(parts, fill = TRUE)
  setkey(cpg, chr, start, end)
  hits <- foverlaps(cpg, promoters, nomatch = 0L)
  out <- hits[, .(mean_methylation = mean(methylation, na.rm = TRUE),
                  feature_count = .N), by = gene]
  rm(cpg, hits)
  gc()
  out
}

process_sample <- function(row) {
  per_dir <- file.path(processed_dir, row$accession, "per_sample_promoter")
  dir_create(per_dir)
  out_file <- file.path(per_dir, paste0(row$sample_id, "_promoter_values.csv"))
  if (file.exists(out_file)) {
    return(data.table(accession = row$accession, sample_id = row$sample_id,
                      group = row$group, status = "already_present",
                      output_file = out_file, rows = NA_integer_))
  }
  dir_create(dirname(row$staged_file))
  source_file <- row$staged_file
  if (!file.exists(source_file)) {
    ans <- download_or_copy(row$source_file_url, NA_character_, source_file)
    if (!file.exists(source_file)) stop("Could not stage WGBS source file for ", row$sample_id, ": ", ans$status)
  }
  message("Processing placenta WGBS ", row$accession, " ", row$sample_id)
  promoters <- promoters_by_build[[row$build]]
  gene_list <- genes_by_build[[row$build]]
  vals <- if (row$accession == "GSE178203") {
    process_cpg_report(source_file, promoters)
  } else if (row$accession == "GSE67615") {
    process_bed_tarball(source_file, promoters)
  } else {
    stop("Unsupported placenta WGBS accession: ", row$accession)
  }
  vals <- merge(data.table(gene = gene_list), vals, by = "gene", all.x = TRUE)
  vals[, `:=`(
    accession = row$accession,
    sample_id = row$sample_id,
    group = row$group,
    build = row$build,
    finite = is.finite(mean_methylation)
  )]
  vals[is.na(feature_count), feature_count := 0L]
  setcolorder(vals, c("accession", "sample_id", "group", "build", "gene",
                      "mean_methylation", "feature_count", "finite"))
  write_csv(vals, out_file)
  data.table(accession = row$accession, sample_id = row$sample_id,
             group = row$group, status = "processed",
             output_file = out_file, rows = nrow(vals))
}

log <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) process_sample(manifest[i])), fill = TRUE)
write_csv(log, file.path(qc_dir, "placenta_WGBS_R_processing_log.csv"))

build_dataset_summary <- function(accession, build, platform_note) {
  acc_id <- accession
  sm <- manifest[accession == acc_id]
  per <- rbindlist(lapply(sm$sample_id, function(id) {
    f <- file.path(processed_dir, accession, "per_sample_promoter", paste0(id, "_promoter_values.csv"))
    fread(f)
  }), fill = TRUE)
  long <- per[finite == TRUE & is.finite(mean_methylation)]
  asd <- long[group == "ASD", .(
    ASD_n = uniqueN(sample_id),
    ASD_mean = mean(mean_methylation),
    ASD_SD = sd(mean_methylation)
  ), by = gene]
  ctl <- long[group == "Control", .(
    control_n = uniqueN(sample_id),
    control_mean = mean(mean_methylation),
    control_SD = sd(mean_methylation)
  ), by = gene]
  feat <- long[, .(feature_count = mean(feature_count, na.rm = TRUE)), by = gene]
  out <- merge(data.table(gene = genes_by_build[[build]]), merge(asd, ctl, by = "gene", all = TRUE),
               by = "gene", all.x = TRUE)
  out <- merge(out, feat, by = "gene", all.x = TRUE)
  out[, `:=`(
    accession = accession,
    tissue_family = "placenta",
    assay_class = "WGBS",
    platform = platform_note,
    finite_summary = fifelse(is.finite(ASD_mean) & is.finite(control_mean) &
                               is.finite(ASD_SD) & is.finite(control_SD) &
                               ASD_n >= 2 & control_n >= 2, "yes", "no"),
    notes = paste0("R workflow from public ", accession, " WGBS source files.")
  )]
  setcolorder(out, c("accession", "tissue_family", "assay_class", "platform", "gene"))
  effects <- hedges_effects(out, value_label = "WGBS promoter methylation")
  write_csv(out, file.path(processed_dir, paste0(accession, "_gene_summary_statistics.csv")))
  write_csv(effects, file.path(processed_dir, paste0(accession, "_effect_sizes.csv")))
  data.table(accession = accession,
             ASD_n = max(out$ASD_n, na.rm = TRUE),
             control_n = max(out$control_n, na.rm = TRUE),
             finite_gene_summaries = sum(out$finite_summary == "yes", na.rm = TRUE),
             effect_size_rows = nrow(effects))
}

summary_qc <- rbind(
  build_dataset_summary("GSE178203", "hg38", "WGBS Bismark CpG reports"),
  build_dataset_summary("GSE67615", "hg19", "WGBS BED percent methylation")
)
write_csv(summary_qc, file.path(qc_dir, "placenta_WGBS_dataset_summary_QC.csv"))

all_summary <- rbindlist(lapply(c("GSE178203", "GSE67615"), function(id) {
  fread(file.path(processed_dir, paste0(id, "_gene_summary_statistics.csv")))
}), fill = TRUE)
all_effects <- rbindlist(lapply(c("GSE178203", "GSE67615"), function(id) {
  fread(file.path(processed_dir, paste0(id, "_effect_sizes.csv")))
}), fill = TRUE)
write_csv(all_summary, file.path(processed_dir, "placenta_all_dataset_gene_summary_statistics.csv"))
write_csv(all_effects, file.path(processed_dir, "placenta_all_dataset_effect_sizes.csv"))
message("Placenta WGBS processing completed.")
