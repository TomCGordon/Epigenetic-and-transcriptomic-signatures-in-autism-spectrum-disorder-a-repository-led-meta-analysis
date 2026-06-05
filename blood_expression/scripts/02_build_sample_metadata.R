# Build sample-level phenotype and inclusion metadata.

message("02: building sample metadata")

source_inventory <- read.csv(file.path(cfg$out_dir, "00_manifest", "blood_expression_public_source_inventory.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
sample_metadata <- build_blood_expression_sample_metadata(cfg, source_inventory)

write_csv_safe(sample_metadata, file.path(cfg$out_dir, "02_sample_metadata", "blood_expression_sample_phenotype_table.csv"))

sample_summary <- aggregate(
  sample_id ~ dataset_id + platform_id + model_role + group + include,
  data = sample_metadata,
  FUN = length
)
names(sample_summary)[names(sample_summary) == "sample_id"] <- "n_samples"
write_csv_safe(sample_summary, file.path(cfg$out_dir, "02_sample_metadata", "blood_expression_sample_summary.csv"))
