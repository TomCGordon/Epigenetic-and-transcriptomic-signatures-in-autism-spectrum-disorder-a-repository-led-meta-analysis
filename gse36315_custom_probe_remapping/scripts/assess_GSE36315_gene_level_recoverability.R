#!/usr/bin/env Rscript

# Assess whether the custom GPL15314 sequence remapping yields usable
# gene-level expression summaries for GSE36315.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_file)
work_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

series_gz <- file.path(work_dir, "01_source/GSE36315_series_matrix.txt.gz")
map_file <- file.path(work_dir, "04_recoverability_assessment/GSE36315_GENCODEv19_bowtie_probe_mapping_assessment.csv")
out_dir <- file.path(work_dir, "04_recoverability_assessment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_series_matrix <- function(series_gz) {
  lines <- readLines(gzfile(series_gz), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  dt <- fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"), sep = "\t", quote = "")
  names(dt)[1] <- "probe_id"
  names(dt) <- gsub('^"+|"+$', "", names(dt))
  dt[, probe_id := gsub('^"+|"+$', "", probe_id)]
  dt
}

hedges_g <- function(x_case, x_control) {
  x_case <- x_case[is.finite(x_case)]
  x_control <- x_control[is.finite(x_control)]
  n1 <- length(x_case); n0 <- length(x_control)
  if (n1 < 2 || n0 < 2) return(list(g = NA_real_, se = NA_real_, var = NA_real_))
  m1 <- mean(x_case); m0 <- mean(x_control)
  s1 <- sd(x_case); s0 <- sd(x_control)
  df <- n1 + n0 - 2
  sp <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / df)
  if (!is.finite(sp) || sp == 0) return(list(g = NA_real_, se = NA_real_, var = NA_real_))
  d <- (m1 - m0) / sp
  J <- 1 - 3 / (4 * df - 1)
  g <- J * d
  v <- (n1 + n0) / (n1 * n0) + (g^2 / (2 * (n1 + n0 - 2)))
  list(g = g, se = sqrt(v), var = v)
}

message("Reading GSE36315 series matrix and custom probe map...")
mat <- read_series_matrix(series_gz)
probe_map <- fread(map_file)
unique_map <- probe_map[mapping_class == "unique_gene", .(probe_id, gene = gene_symbols)]

mapped <- merge(unique_map, mat, by = "probe_id", all.x = FALSE, all.y = FALSE)
sample_cols <- setdiff(names(mat), "probe_id")
for (cc in sample_cols) mapped[, (cc) := as.numeric(get(cc))]

gene_matrix <- mapped[, lapply(.SD, mean, na.rm = TRUE), by = gene, .SDcols = sample_cols]
for (cc in sample_cols) gene_matrix[is.nan(get(cc)), (cc) := NA_real_]
fwrite(gene_matrix, file.path(out_dir, "GSE36315_GENCODEv19_unique_gene_expression_matrix.csv.gz"))

sample_meta <- data.table(
  sample_id = c("GSM887867", "GSM887868", "GSM887869", "GSM887870", "GSM887871", "GSM887872", "GSM887873", "GSM887874"),
  diagnosis = c("ASD", "ASD", "ASD", "ASD", "Control", "Control", "Control", "Control"),
  region = c("prefrontal cortex", "cerebellum", "prefrontal cortex", "cerebellum", "prefrontal cortex", "cerebellum", "prefrontal cortex", "cerebellum"),
  donor_pair = c("ASD_1", "ASD_1", "ASD_2", "ASD_2", "Control_1", "Control_1", "Control_2", "Control_2")
)
fwrite(sample_meta, file.path(out_dir, "GSE36315_sample_metadata_from_GEO.csv"))

effect_for_region <- function(region_name) {
  sm <- sample_meta[region == region_name]
  asd <- sm[diagnosis == "ASD"]$sample_id
  ctl <- sm[diagnosis == "Control"]$sample_id
  out <- gene_matrix[, {
    eff <- hedges_g(unlist(.SD[, ..asd]), unlist(.SD[, ..ctl]))
    .(k = 1L, n_asd = length(asd), n_control = length(ctl), hedges_g = eff$g, se = eff$se, variance = eff$var)
  }, by = gene]
  out[, region := region_name]
  out[, dataset := "gse36315_custom_annotated"]
  out[]
}

effects <- rbindlist(list(effect_for_region("prefrontal cortex"), effect_for_region("cerebellum")), fill = TRUE)
setcolorder(effects, c("dataset", "region", "gene", "k", "n_asd", "n_control", "hedges_g", "se", "variance"))
fwrite(effects, file.path(out_dir, "gse36315_custom_annotated_region_effect_sizes.csv"))

summary <- rbindlist(list(
  data.table(metric = "series_matrix_probes", value = nrow(mat)),
  data.table(metric = "unique_gene_mapped_probes_used", value = nrow(unique_map)),
  data.table(metric = "mapped_probe_rows_in_series", value = nrow(mapped)),
  data.table(metric = "gene_level_rows_after_probe_collapse", value = nrow(gene_matrix)),
  data.table(metric = "genes_with_all_8_samples_finite", value = nrow(gene_matrix[complete.cases(gene_matrix[, ..sample_cols])])),
  data.table(metric = "prefrontal_effect_rows_finite", value = nrow(effects[region == "prefrontal cortex" & is.finite(hedges_g) & is.finite(se)])),
  data.table(metric = "cerebellum_effect_rows_finite", value = nrow(effects[region == "cerebellum" & is.finite(hedges_g) & is.finite(se)]))
))
fwrite(summary, file.path(out_dir, "GSE36315_gene_level_recoverability_summary.csv"))

report <- c(
  "# GSE36315 Gene-Level Recoverability Assessment",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "The custom GENCODE-remapped annotation is sufficient to construct a broad gene-level expression matrix for many genes in GSE36315. Because the dataset contains two ASD donors and two control donors, each represented in prefrontal cortex and cerebellum, region-specific effect sizes should be treated as very small-sample sensitivity/descriptive evidence. A combined all-sample model would not be independent without explicit donor-level handling."
)
writeLines(report, file.path(work_dir, "05_reports/GSE36315_gene_level_recoverability_report.md"))

message("Finished recoverability assessment.")
