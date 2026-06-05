# Stage public GEO source files and annotation files.

message("01: staging public source files")

source_inventory <- stage_blood_expression_sources(cfg)
write_csv_safe(source_inventory, file.path(cfg$out_dir, "00_manifest", "blood_expression_public_source_inventory.csv"))

pkg_versions <- package_version_table(c(
  "GEOquery", "Biobase", "limma", "edgeR", "AnnotationDbi", "org.Hs.eg.db",
  "data.table", "readr", "metafor", "meta"
))
write_csv_safe(pkg_versions, file.path(cfg$out_dir, "00_manifest", "blood_expression_R_package_versions.csv"))
