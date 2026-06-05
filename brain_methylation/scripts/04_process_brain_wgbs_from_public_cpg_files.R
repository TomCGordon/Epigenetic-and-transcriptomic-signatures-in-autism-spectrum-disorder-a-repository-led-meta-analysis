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

wgbs_raw_dir <- file.path(package_root, "data_raw", "WGBS")
processed_dir <- file.path(package_root, "data_processed", "WGBS")
qc_dir <- file.path(package_root, "qc")
dir_create(processed_dir)
dir_create(qc_dir)

promoters <- fread(file.path(package_root, "data_processed", "annotation", "brain_hg38_refGene_promoter_coordinates.csv"))
setnames(promoters, c("chr", "start", "end"), c("chr", "start", "end"))
promoters <- promoters[!is.na(chr) & is.finite(start) & is.finite(end)]
promoters[, `:=`(start = as.integer(start), end = as.integer(end))]
setkey(promoters, chr, start, end)

universe <- fread(file.path(package_root, "data_processed", "annotation", "brain_harmonised_20960_gene_universe.csv"))$gene
sample_manifest <- fread(file.path(qc_dir, "brain_WGBS_GEO_sample_manifest.csv"))
sample_manifest <- sample_manifest[include == TRUE & group %in% c("ASD", "Control")]
staged_manifest <- fread(file.path(qc_dir, "brain_WGBS_staged_source_manifest.csv"))
sample_manifest <- merge(sample_manifest, staged_manifest[, .(accession, sample_id, source_file = dest, stage_status = status)],
                         by = c("accession", "sample_id"), all.x = TRUE)
if (any(!file.exists(sample_manifest$source_file))) {
  stop("Missing staged WGBS source files. Re-run 01_download_brain_sources.R.")
}

process_wgbs_sample <- function(row) {
  out_dir <- file.path(processed_dir, row$accession, "per_sample_promoter")
  dir_create(out_dir)
  out_file <- file.path(out_dir, paste0(row$sample_id, "_promoter_values.csv"))
  if (file.exists(out_file)) {
    return(data.table(accession = row$accession, sample_id = row$sample_id, group = row$group,
                      output_file = out_file, status = "already_present"))
  }
  message("Processing WGBS ", row$accession, " ", row$sample_id)
  if (row$accession == "GSE109875") {
    cpg <- fread(row$source_file, header = FALSE, select = c(1, 2, 4, 5),
                 col.names = c("chr", "pos", "methylated", "unmethylated"), showProgress = FALSE)
    cpg[, coverage := methylated + unmethylated]
    cpg <- cpg[coverage > 0]
    cpg[, methylation := methylated / coverage]
    cpg[, `:=`(start = as.integer(pos), end = as.integer(pos))]
  } else if (row$accession == "GSE81541") {
    cpg <- fread(row$source_file, header = FALSE, select = c(1, 2, 3, 4),
                 col.names = c("chr", "start0", "end0", "name"), showProgress = FALSE)
    cpg[, methylation := suppressWarnings(as.numeric(sub("-.*$", "", name)))]
    cpg <- cpg[is.finite(methylation)]
    cpg[, `:=`(start = as.integer(start0 + 1L), end = as.integer(start0 + 1L))]
  } else {
    stop("Unsupported WGBS accession: ", row$accession)
  }
  cpg <- cpg[grepl("^chr", chr) & is.finite(start) & is.finite(end) & is.finite(methylation),
             .(chr, start, end, methylation)]
  setkey(cpg, chr, start, end)
  hits <- foverlaps(cpg, promoters, nomatch = 0L)
  gene_values <- hits[, .(
    promoterMean = mean(methylation, na.rm = TRUE),
    promoterCpGsUsed = .N
  ), by = gene]
  gene_values <- merge(data.table(gene = universe), gene_values, by = "gene", all.x = TRUE)
  fwrite(gene_values, out_file)
  rm(cpg, hits, gene_values)
  gc()
  data.table(accession = row$accession, sample_id = row$sample_id, group = row$group,
             output_file = out_file, status = "processed")
}

log <- rbindlist(lapply(seq_len(nrow(sample_manifest)), function(i) process_wgbs_sample(sample_manifest[i])), fill = TRUE)
fwrite(log, file.path(qc_dir, "brain_WGBS_R_processing_log.csv"))

build_wgbs_dataset_summary <- function(accession, brain_region, broader_region_group,
                                       output_accession = accession,
                                       sample_ids = NULL) {
  acc_id <- accession
  sm <- sample_manifest[accession == acc_id]
  if (!is.null(sample_ids)) sm <- sm[sample_id %in% sample_ids]
  if (sum(sm$group == "ASD") < 2 || sum(sm$group == "Control") < 2) {
    return(NULL)
  }
  per_sample <- rbindlist(lapply(seq_len(nrow(sm)), function(i) {
    f <- file.path(processed_dir, sm$accession[i], "per_sample_promoter", paste0(sm$sample_id[i], "_promoter_values.csv"))
    dt <- fread(f)
    dt[, `:=`(sample_id = sm$sample_id[i], group = sm$group[i])]
    dt
  }), fill = TRUE)
  long <- per_sample[is.finite(promoterMean)]
  asd <- long[group == "ASD", .(
    ASD_n = uniqueN(sample_id),
    ASD_mean_methylation = mean(promoterMean),
    ASD_SD = sd(promoterMean)
  ), by = gene]
  ctl <- long[group == "Control", .(
    control_n = uniqueN(sample_id),
    control_mean_methylation = mean(promoterMean),
    control_SD = sd(promoterMean)
  ), by = gene]
  feat <- long[, .(feature_probe_count = mean(promoterCpGsUsed, na.rm = TRUE)), by = gene]
  out <- merge(data.table(gene = universe), merge(asd, ctl, by = "gene", all = TRUE), by = "gene", all.x = TRUE)
  out <- merge(out, feat, by = "gene", all.x = TRUE)
  out[, `:=`(
    accession = output_accession,
    brain_region = brain_region,
    broader_region_group = broader_region_group,
    platform = "Illumina HiSeq 2000",
    assay_class = "WGBS",
    mean_difference = ASD_mean_methylation - control_mean_methylation,
    finite_summary = fifelse(is.finite(ASD_mean_methylation) & is.finite(control_mean_methylation) &
                               is.finite(ASD_SD) & is.finite(control_SD) &
                               ASD_n >= 2 & control_n >= 2, "yes", "no"),
    missingness_reason = fifelse(is.na(ASD_n) | is.na(control_n), "gene not represented by finite WGBS promoter summaries", ""),
    phenotype_labels_used = "ASD vs Control",
    promoter_definition_used = "Coordinate-based hg38 refGene promoter windows generated in R",
    aggregation_method = "Per-sample promoter methylation values summarised by ASD/control group.",
    notes = if (accession == "GSE109875") {
      "Retain direct BA9 cortical ASD-control WGBS subseries only; umbrella routes GSE119981 and PRJNA490887 are annotations only."
    } else if (output_accession != accession) {
      "Region-specific sensitivity subset generated from public per-sample promoter values; not combined with the parent GSE81541 route in the same model."
    } else {
      "Retain Brain_IdioAut and Brain_Control samples; exclude BraindfBA duplicate-region, syndromic groups and cell-culture samples; PRJNA321909 is umbrella annotation only."
    }
  )]
  setcolorder(out, c("accession", "brain_region", "broader_region_group", "platform", "assay_class", "gene"))
  effects <- hedges_effects(out)
  fwrite(out, file.path(processed_dir, paste0(output_accession, "_WGBS_gene_summary_statistics.csv")))
  fwrite(effects, file.path(processed_dir, paste0(output_accession, "_WGBS_effect_sizes.csv")))
  data.table(
    accession = output_accession,
    parent_accession = accession,
    ASD_n = max(out$ASD_n, na.rm = TRUE),
    control_n = max(out$control_n, na.rm = TRUE),
    finite_gene_summaries = sum(out$finite_summary == "yes", na.rm = TRUE),
    effect_size_rows = nrow(effects)
  )
}

wgbs_qc <- rbindlist(list(
  build_wgbs_dataset_summary("GSE109875", "BA9 dorsal lateral prefrontal cortex", "cortex"),
  build_wgbs_dataset_summary("GSE81541", "Grouped post-mortem brain idiopathic autism/control", "mixed brain")
), fill = TRUE)
fwrite(wgbs_qc, file.path(qc_dir, "brain_WGBS_dataset_summary_QC.csv"))

gse81541_ba9 <- sample_manifest[accession == "GSE81541" & grepl("BA9", source_name, ignore.case = TRUE)]$sample_id
wgbs_region_qc <- rbindlist(list(
  build_wgbs_dataset_summary("GSE81541", "Region-specific subroute: BA9 cortex", "cortex",
                             output_accession = "GSE81541_BA9", sample_ids = gse81541_ba9)
), fill = TRUE)
if (nrow(wgbs_region_qc)) {
  fwrite(wgbs_region_qc, file.path(qc_dir, "brain_WGBS_region_subroute_QC.csv"))
}

all_summary <- rbindlist(lapply(c("GSE109875", "GSE81541"), function(id) {
  fread(file.path(processed_dir, paste0(id, "_WGBS_gene_summary_statistics.csv")))
}), fill = TRUE)
all_effects <- rbindlist(lapply(c("GSE109875", "GSE81541"), function(id) {
  fread(file.path(processed_dir, paste0(id, "_WGBS_effect_sizes.csv")))
}), fill = TRUE)
fwrite(all_summary, file.path(processed_dir, "brain_WGBS_all_dataset_gene_summary_statistics.csv"))
fwrite(all_effects, file.path(processed_dir, "brain_WGBS_all_dataset_effect_sizes.csv"))
region_summary_files <- file.path(processed_dir, paste0(wgbs_region_qc$accession, "_WGBS_gene_summary_statistics.csv"))
region_effect_files <- file.path(processed_dir, paste0(wgbs_region_qc$accession, "_WGBS_effect_sizes.csv"))
if (length(region_summary_files) && all(file.exists(region_summary_files))) {
  fwrite(rbindlist(lapply(region_summary_files, fread), fill = TRUE),
         file.path(processed_dir, "brain_WGBS_region_subroute_gene_summary_statistics.csv"))
}
if (length(region_effect_files) && all(file.exists(region_effect_files))) {
  fwrite(rbindlist(lapply(region_effect_files, fread), fill = TRUE),
         file.path(processed_dir, "brain_WGBS_region_subroute_effect_sizes.csv"))
}
message("Brain WGBS datasets processed from public CpG/BED files.")
