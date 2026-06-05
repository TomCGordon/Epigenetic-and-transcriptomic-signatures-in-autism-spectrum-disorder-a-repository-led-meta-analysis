# Process public GEO series matrices into gene-level expression values.

message("03: processing GEO microarray series matrices")

source_inventory <- read.csv(file.path(cfg$out_dir, "00_manifest", "blood_expression_public_source_inventory.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
sample_metadata <- read.csv(file.path(cfg$out_dir, "02_sample_metadata", "blood_expression_sample_phenotype_table.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)

microarray_routes <- source_inventory[
  source_inventory$route_type == "series_matrix" &
    source_inventory$dataset_id %in% c("GSE111175", "GSE111176", "GSE18123", "GSE25507", "GSE26415", "GSE123302"),
]

processed <- list()
logs <- list()
for (i in seq_len(nrow(microarray_routes))) {
  route <- microarray_routes[i, ]
  message("  processing ", route$route_id)
  res <- process_series_matrix_route(cfg, route, sample_metadata)
  processed[[length(processed) + 1]] <- res$gene_expression
  logs[[length(logs) + 1]] <- res$log
}

gene_expression <- rbindlist_fill(processed)
processing_log <- rbindlist_fill(logs)

write_csv_safe(gene_expression, file.path(cfg$out_dir, "03_processed_expression", "blood_expression_microarray_gene_values_long.csv"))
write_csv_safe(processing_log, file.path(cfg$out_dir, "03_processed_expression", "blood_expression_microarray_processing_log.csv"))
