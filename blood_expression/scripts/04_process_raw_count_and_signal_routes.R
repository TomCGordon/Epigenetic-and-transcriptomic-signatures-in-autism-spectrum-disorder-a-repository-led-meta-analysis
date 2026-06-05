# Process public RAW-derived expression routes:
# - GSE140702: featureCounts files, NT unstimulated monocytes retained.
# - GSE77103: Agilent one-color signal files.

message("04: processing raw count/signal routes")

source_inventory <- read.csv(file.path(cfg$out_dir, "00_manifest", "blood_expression_public_source_inventory.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
sample_metadata <- read.csv(file.path(cfg$out_dir, "02_sample_metadata", "blood_expression_sample_phenotype_table.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)

res_140702 <- process_gse140702_featurecounts(cfg, sample_metadata)
res_77103 <- process_gse77103_agilent(cfg, sample_metadata)

gene_expression <- rbindlist_fill(list(res_140702$gene_expression, res_77103$gene_expression))
processing_log <- rbindlist_fill(list(res_140702$log, res_77103$log))

write_csv_safe(gene_expression, file.path(cfg$out_dir, "03_processed_expression", "blood_expression_raw_route_gene_values_long.csv"))
write_csv_safe(processing_log, file.path(cfg$out_dir, "03_processed_expression", "blood_expression_raw_route_processing_log.csv"))
