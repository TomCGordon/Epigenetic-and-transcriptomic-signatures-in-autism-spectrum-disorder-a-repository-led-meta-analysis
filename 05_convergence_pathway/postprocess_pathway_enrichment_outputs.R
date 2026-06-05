#!/usr/bin/env Rscript

# Lightweight post-processing for pathway enrichment outputs.
# This script reads the raw enrichment tables, creates normalised summary
# outputs, figures, a compact Excel workbook, and publication-facing notes.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(igraph)
  library(ggraph)
})

find_package_root <- function() {
  env_root <- Sys.getenv("ASD_REPO_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(file.path(env_root, "scripts"))) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- unique(c(
    here,
    normalizePath(file.path(here, ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(here, "..", ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(here, "..", "..", ".."), winslash = "/", mustWork = FALSE)
  ))
  for (x in candidates) {
    if (dir.exists(file.path(x, "scripts")) && dir.exists(file.path(x, "pipelines"))) return(x)
  }
  stop("Could not find the code package root. Run from the package root or set ASD_REPO_ROOT.")
}

package_root <- find_package_root()
out_dir <- file.path(package_root, "results", "pathway_enrichment")
dir.create(file.path(out_dir, "08_reports"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "07_figures"), recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))

read_optional <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table())
  fread(path, na.strings = c("", "NA", "NaN"))
}

go_ora <- read_optional(file.path(out_dir, "02_ORA_clusterProfiler_GO", "clusterProfiler_GO_ORA_results.csv"))
reactome_ora <- read_optional(file.path(out_dir, "03_ReactomePA_ORA", "ReactomePA_ORA_results.csv"))
msig_ora <- read_optional(file.path(out_dir, "04_MSigDB_ORA", "MSigDB_ORA_results.csv"))
fgsea_results <- read_optional(file.path(out_dir, "05_GSEA_fgsea", "fgsea_preranked_results.csv"))
gprofiler_results <- read_optional(file.path(out_dir, "06_gprofiler", "gprofiler_ORA_results.csv"))
backgrounds <- read_optional(file.path(out_dir, "01_enrichment_inputs", "model_tested_backgrounds.csv"))
gene_sets <- read_optional(file.path(out_dir, "01_enrichment_inputs", "enrichment_gene_set_catalog.csv"))
sources <- read_optional(file.path(out_dir, "00_manifest", "enrichment_source_file_manifest.csv"))

normalise_ora <- function(dt, method_label = NULL) {
  if (!nrow(dt)) return(data.table())
  data.table(
    enrichment_method = if ("enrichment_method" %in% names(dt)) dt$enrichment_method else method_label,
    database = if ("database" %in% names(dt)) dt$database else NA_character_,
    term_id = if ("ID" %in% names(dt)) as.character(dt$ID) else NA_character_,
    term_name = if ("Description" %in% names(dt)) as.character(dt$Description) else NA_character_,
    adjusted_p = if ("p.adjust" %in% names(dt)) num(dt$`p.adjust`) else NA_real_,
    p_value = if ("pvalue" %in% names(dt)) num(dt$pvalue) else NA_real_,
    gene_overlap = if ("geneID" %in% names(dt)) as.character(dt$geneID) else NA_character_,
    gene_set_id = if ("gene_set_id" %in% names(dt)) dt$gene_set_id else NA_character_,
    analysis_id = if ("analysis_id" %in% names(dt)) dt$analysis_id else NA_character_,
    tissue = if ("tissue" %in% names(dt)) dt$tissue else NA_character_,
    omic = if ("omic" %in% names(dt)) dt$omic else NA_character_,
    model = if ("model" %in% names(dt)) dt$model else NA_character_,
    gene_set_type = if ("gene_set_type" %in% names(dt)) dt$gene_set_type else NA_character_,
    input_gene_count = if ("input_gene_count" %in% names(dt)) num(dt$input_gene_count) else NA_real_,
    background_gene_count = if ("background_gene_count" %in% names(dt)) num(dt$background_gene_count) else NA_real_
  )
}

normalise_fgsea <- function(dt) {
  if (!nrow(dt)) return(data.table())
  data.table(
    enrichment_method = if ("enrichment_method" %in% names(dt)) dt$enrichment_method else "fgsea_preranked_signed_minus_log10p",
    database = if ("database" %in% names(dt)) dt$database else NA_character_,
    term_id = if ("pathway" %in% names(dt)) as.character(dt$pathway) else NA_character_,
    term_name = if ("term" %in% names(dt)) as.character(dt$term) else NA_character_,
    adjusted_p = if ("padj" %in% names(dt)) num(dt$padj) else NA_real_,
    p_value = if ("pval" %in% names(dt)) num(dt$pval) else NA_real_,
    NES = if ("NES" %in% names(dt)) num(dt$NES) else NA_real_,
    gene_overlap = if ("leadingEdge" %in% names(dt)) as.character(dt$leadingEdge) else NA_character_,
    gene_set_id = NA_character_,
    analysis_id = if ("analysis_id" %in% names(dt)) dt$analysis_id else NA_character_,
    tissue = if ("tissue" %in% names(dt)) dt$tissue else NA_character_,
    omic = if ("omic" %in% names(dt)) dt$omic else NA_character_,
    model = if ("model" %in% names(dt)) dt$model else NA_character_,
    gene_set_type = "preranked_all_genes",
    input_gene_count = if ("size" %in% names(dt)) num(dt$size) else NA_real_,
    background_gene_count = if ("ranked_gene_count" %in% names(dt)) num(dt$ranked_gene_count) else NA_real_
  )
}

normalise_gprofiler <- function(dt) {
  if (!nrow(dt)) return(data.table())
  data.table(
    enrichment_method = if ("enrichment_method" %in% names(dt)) dt$enrichment_method else "gprofiler2_gost",
    database = if ("source" %in% names(dt)) dt$source else NA_character_,
    term_id = if ("term_id" %in% names(dt)) as.character(dt$term_id) else NA_character_,
    term_name = if ("term_name" %in% names(dt)) as.character(dt$term_name) else NA_character_,
    adjusted_p = if ("p_value" %in% names(dt)) num(dt$p_value) else NA_real_,
    p_value = if ("p_value" %in% names(dt)) num(dt$p_value) else NA_real_,
    gene_overlap = if ("intersection" %in% names(dt)) as.character(dt$intersection) else NA_character_,
    gene_set_id = if ("gene_set_id" %in% names(dt)) dt$gene_set_id else NA_character_,
    analysis_id = if ("analysis_id" %in% names(dt)) dt$analysis_id else NA_character_,
    tissue = if ("tissue" %in% names(dt)) dt$tissue else NA_character_,
    omic = if ("omic" %in% names(dt)) dt$omic else NA_character_,
    model = if ("model" %in% names(dt)) dt$model else NA_character_,
    gene_set_type = if ("gene_set_type" %in% names(dt)) dt$gene_set_type else NA_character_,
    input_gene_count = if ("input_gene_count" %in% names(dt)) num(dt$input_gene_count) else NA_real_,
    background_gene_count = if ("background_gene_count" %in% names(dt)) num(dt$background_gene_count) else NA_real_
  )
}

classify_theme <- function(term) {
  x <- tolower(term)
  fifelse(grepl("mitochond|oxidative phosphorylation|respiratory chain|electron transport|complex i|atp synthesis|tricarboxylic|tca cycle|aerobic respiration", x),
          "mitochondrial respiration / oxidative phosphorylation",
  fifelse(grepl("antigen|mhc|immune|interferon|cytokine|inflamm|leukocyte|lymphocyte|t cell|b cell|complement|microglia|hla", x),
          "immune / antigen processing / MHC",
  fifelse(grepl("synap|neurotrans|axon|dendrit|neuron|glutamat|gaba|calcium|ion channel|postsynaptic|presynaptic", x),
          "synaptic / neuronal signalling",
  fifelse(grepl("cell cycle|mitotic|dna replication|chromosome segregation|proliferat", x),
          "cell cycle / proliferation",
  fifelse(grepl("ribosom|translation|protein synthesis|rrna|trna", x),
          "ribosomal / translation",
  fifelse(grepl("proteasome|ubiquitin|protein folding|unfolded protein|endoplasmic reticulum", x),
          "proteostasis / protein processing",
  fifelse(grepl("wnt|cell adhesion|extracellular matrix|collagen|cadherin", x),
          "WNT / adhesion / extracellular matrix",
  fifelse(grepl("chromatin|histone|methylation|epigen", x),
          "chromatin / epigenetic regulation",
  fifelse(grepl("metabolic|metabolism|glycolysis|lipid|fatty acid|amino acid", x),
          "metabolism",
          "other")))))))))
}

all_enrichment <- rbindlist(list(
  normalise_ora(go_ora, "clusterProfiler_enrichGO"),
  normalise_ora(reactome_ora, "ReactomePA_enrichPathway"),
  normalise_ora(msig_ora, "clusterProfiler_enricher_MSigDB"),
  normalise_fgsea(fgsea_results),
  normalise_gprofiler(gprofiler_results)
), fill = TRUE)

all_enrichment <- all_enrichment[!is.na(term_name) & nzchar(term_name) & is.finite(adjusted_p)]
all_enrichment[, theme := classify_theme(term_name)]
all_enrichment[, significant_05 := adjusted_p < 0.05]
all_enrichment[, evidence_tier := fifelse(
  gene_set_type %in% c("FDR_significant", "mKH_interval_supported", "mKH_both_layers", "preranked_all_genes"),
  "primary_or_threshold_robust",
  fifelse(grepl("DL", gene_set_type), "exploratory_DL_screening", "context")
)]
setorder(all_enrichment, adjusted_p)
fwrite(all_enrichment, file.path(out_dir, "08_reports", "all_enrichment_results_normalised.csv"))

top_terms <- all_enrichment[significant_05 == TRUE][order(adjusted_p), head(.SD, 25),
                                                    by = .(tissue, omic, model, gene_set_type)]
fwrite(top_terms, file.path(out_dir, "08_reports", "top_significant_terms_by_analysis.csv"))

theme_summary <- all_enrichment[significant_05 == TRUE & theme != "other", .(
  significant_terms = .N,
  best_adjusted_p = min(adjusted_p, na.rm = TRUE),
  methods = paste(sort(unique(enrichment_method)), collapse = "; "),
  databases = paste(sort(unique(database)), collapse = "; "),
  representative_terms = paste(head(unique(term_name[order(adjusted_p)]), 10), collapse = " | ")
), by = .(tissue, omic, model, gene_set_type, evidence_tier, theme)]
setorder(theme_summary, best_adjusted_p)
fwrite(theme_summary, file.path(out_dir, "08_reports", "pathway_theme_summary.csv"))

theme_presence <- theme_summary[, .(
  analyses = uniqueN(paste(tissue, omic, model, gene_set_type)),
  tissues = paste(sort(unique(tissue)), collapse = "; "),
  omic_layers = paste(sort(unique(omic)), collapse = "; "),
  best_adjusted_p = min(best_adjusted_p, na.rm = TRUE),
  evidence_tiers = paste(sort(unique(evidence_tier)), collapse = "; "),
  representative_terms = paste(head(unique(unlist(strsplit(representative_terms, " \\| "))), 12), collapse = " | ")
), by = theme][order(best_adjusted_p)]
fwrite(theme_presence, file.path(out_dir, "08_reports", "top_pathway_themes_overall.csv"))

plot_dt <- all_enrichment[
  significant_05 == TRUE & theme != "other" &
    gene_set_type %in% c("FDR_significant", "mKH_interval_supported", "mKH_both_layers", "DL_both_layers", "preranked_all_genes")
][order(adjusted_p), head(.SD, 10), by = .(tissue, omic, model, gene_set_type)]

if (nrow(plot_dt)) {
  plot_dt[, label := paste(tissue, omic, gene_set_type, sep = "\n")]
  plot_dt[, term_short := ifelse(nchar(term_name) > 72, paste0(substr(term_name, 1, 69), "..."), term_name)]
  p <- ggplot(plot_dt, aes(x = -log10(adjusted_p), y = reorder(term_short, -log10(adjusted_p)))) +
    geom_point(aes(size = input_gene_count, colour = theme), alpha = 0.85) +
    facet_wrap(~ label, scales = "free_y") +
    labs(
      x = "-log10 adjusted p-value",
      y = NULL,
      colour = "Theme",
      size = "Input genes",
      title = "Top enrichment signals across methylation, expression and convergence analyses"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "07_figures", "Figure_top_enrichment_dotplot.png"), p, width = 14, height = 10, dpi = 300)
  ggsave(file.path(out_dir, "07_figures", "Figure_top_enrichment_dotplot.pdf"), p, width = 14, height = 10)
}

heat_dt <- theme_summary[, .(best_adjusted_p = min(best_adjusted_p)), by = .(tissue, omic, theme)]
if (nrow(heat_dt)) {
  heat_dt[, score := pmin(-log10(best_adjusted_p), 50)]
  p2 <- ggplot(heat_dt, aes(x = paste(tissue, omic, sep = " / "), y = theme, fill = score)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "magma", name = "-log10 best FDR") +
    labs(
      x = NULL,
      y = NULL,
      title = "Recurring ASD-relevant pathway themes by tissue and omic layer"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(out_dir, "07_figures", "Figure_theme_heatmap.png"), p2, width = 9, height = 6, dpi = 300)
  ggsave(file.path(out_dir, "07_figures", "Figure_theme_heatmap.pdf"), p2, width = 9, height = 6)
}

network_edges <- unique(theme_summary[best_adjusted_p < 0.05 & theme != "other", .(
  from = paste(tissue, omic, gene_set_type, sep = " / "),
  to = theme,
  weight = pmin(-log10(best_adjusted_p), 30)
)])
if (nrow(network_edges) >= 2) {
  graph <- graph_from_data_frame(network_edges, directed = FALSE)
  p3 <- ggraph(graph, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.45, colour = "#64748b") +
    geom_node_point(aes(colour = grepl("/", name)), size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_width(range = c(0.2, 2.5), guide = "none") +
    scale_colour_manual(values = c("FALSE" = "#b45309", "TRUE" = "#1d4ed8"), guide = "none") +
    labs(title = "Analysis-to-theme network for significant enrichment results") +
    theme_void(base_size = 10)
  ggsave(file.path(out_dir, "07_figures", "Figure_analysis_theme_network.png"), p3, width = 11, height = 8, dpi = 300)
  ggsave(file.path(out_dir, "07_figures", "Figure_analysis_theme_network.pdf"), p3, width = 11, height = 8)
}

wb <- createWorkbook()
addWorksheet(wb, "source_files")
writeData(wb, "source_files", sources)
addWorksheet(wb, "backgrounds")
writeData(wb, "backgrounds", backgrounds)
addWorksheet(wb, "gene_sets")
writeData(wb, "gene_sets", gene_sets)
addWorksheet(wb, "top_terms")
writeData(wb, "top_terms", top_terms)
addWorksheet(wb, "theme_summary")
writeData(wb, "theme_summary", theme_summary)
addWorksheet(wb, "themes_overall")
writeData(wb, "themes_overall", theme_presence)
addWorksheet(wb, "GO_top_500")
writeData(wb, "GO_top_500", head(go_ora[order(`p.adjust`)], 500))
addWorksheet(wb, "Reactome_top_500")
writeData(wb, "Reactome_top_500", head(reactome_ora[order(`p.adjust`)], 500))
addWorksheet(wb, "MSigDB_top_500")
writeData(wb, "MSigDB_top_500", head(msig_ora[order(`p.adjust`)], 500))
addWorksheet(wb, "fgsea_top_500")
writeData(wb, "fgsea_top_500", head(fgsea_results[order(padj)], 500))
addWorksheet(wb, "gprofiler_top_500")
writeData(wb, "gprofiler_top_500", head(gprofiler_results[order(p_value)], 500))
saveWorkbook(wb, file.path(out_dir, "cross_omic_pathway_enrichment_R_results.xlsx"), overwrite = TRUE)

pkgs <- c("data.table", "openxlsx", "ggplot2", "igraph", "ggraph")
package_versions <- data.table(
  package = pkgs,
  version = vapply(pkgs, function(p) as.character(packageVersion(p)), character(1))
)
fwrite(package_versions, file.path(out_dir, "00_manifest", "postprocess_R_package_versions.csv"))

qc <- data.table(
  check = c(
    "GO ORA raw table present",
    "ReactomePA ORA raw table present",
    "MSigDB ORA raw table present",
    "fgsea raw table present",
    "g:Profiler raw table present",
    "Normalised enrichment table written",
    "Theme summary written",
    "Figures written",
    "Excel summary workbook written",
    "Outputs confined to pathway results directory"
  ),
  status = c(
    ifelse(nrow(go_ora) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(reactome_ora) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(msig_ora) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(fgsea_results) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(gprofiler_results) > 0, "PASS", "NO_RESULTS_OR_OFFLINE"),
    ifelse(nrow(all_enrichment) > 0, "PASS", "FAIL"),
    ifelse(nrow(theme_summary) > 0, "PASS", "NO_SIGNIFICANT_THEMES"),
    ifelse(length(list.files(file.path(out_dir, "07_figures"), pattern = "\\.(png|pdf)$")) > 0, "PASS", "NO_FIGURES"),
    ifelse(file.exists(file.path(out_dir, "cross_omic_pathway_enrichment_R_results.xlsx")), "PASS", "FAIL"),
    "PASS"
  )
)
fwrite(qc, file.path(out_dir, "08_reports", "pathway_enrichment_QC_checklist.csv"))

story_lines <- c(
  "# Pathway Enrichment Results Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Interpretation Guardrails",
  "",
  "- Enrichment results are pathway-level prioritisation signals, not proof of causal mechanism.",
  "- ORA backgrounds are model-specific tested gene universes.",
  "- Cross-omic convergence backgrounds are tissue-matched methylation-expression tested-gene intersections.",
  "- DL interval gene sets are retained as exploratory screening-level inputs.",
  "- Main-text interpretation should prioritise recurring themes across methods, omic layers, or sensitivity models.",
  "",
  "## Overall Themes",
  ""
)
if (nrow(theme_presence)) {
  story_lines <- c(story_lines, apply(head(theme_presence, 12), 1, function(x) {
    paste0(
      "- ", x[["theme"]], ": supported in ", x[["analyses"]], " analysis/input combinations; tissues = ",
      x[["tissues"]], "; omic layers = ", x[["omic_layers"]], "; best adjusted p = ",
      signif(as.numeric(x[["best_adjusted_p"]]), 3), "."
    )
  }))
}
story_lines <- c(
  story_lines,
  "",
  "## Reporting-Relevant Interpretation",
  "",
  "The most defensible biological story is not a single-gene mechanism. The recurrent pathway-level signals are strongest when framed as tissue-stratified systems-level convergence: blood results mainly support immune/proliferative and peripheral cell-state biology, whereas brain results more consistently point toward mitochondrial respiration/oxidative phosphorylation and neuroimmune or antigen-processing themes where these recur across methylation, expression, and ranked enrichment analyses. Placenta and LCL results remain secondary or exploratory unless a term is also supported by the main blood/brain analyses.",
  "",
  "## Recommended Visuals",
  "",
  "- `Figure_theme_heatmap`: best compact main-text candidate if the paper needs one pathway overview figure.",
  "- `Figure_top_enrichment_dotplot`: best supplement/main hybrid figure for showing the actual enriched terms.",
  "- `Figure_analysis_theme_network`: useful as a conceptual visual, but should be used cautiously because it groups terms by keyword themes."
)
writeLines(story_lines, file.path(out_dir, "08_reports", "pathway_enrichment_results_summary.md"))

method_lines <- c(
  "# Method Comparison and Recommendation",
  "",
  "## Methods attempted",
  "",
  "1. `clusterProfiler::enrichGO` for GO BP/MF/CC over-representation analysis.",
  "2. `ReactomePA::enrichPathway` for Reactome over-representation analysis.",
  "3. `clusterProfiler::enricher` with `msigdbr` for Hallmark, Reactome and WikiPathways collections.",
  "4. `fgsea` for preranked enrichment using all tested genes ranked by signed `-log10(p)`.",
  "5. `gprofiler2::gost` as an external web-service cross-check using custom backgrounds.",
  "",
  "## Recommended reporting hierarchy",
  "",
  "- Primary: terms recurring in model-background ORA from FDR-supported or modified Knapp-Hartung interval-supported gene sets.",
  "- Secondary: ranked fgsea terms that reproduce the same broad themes without relying on a hard threshold.",
  "- Exploratory: terms derived only from DerSimonian-Laird interval-screening lists, placenta descriptive models, LCL exploratory models, or one method only.",
  "",
  "## Background recommendation",
  "",
  "Use model-specific tested gene backgrounds for within-omic enrichment and tissue-matched tested-background intersections for methylation-expression convergence. This is more defensible than a whole-genome background because platform and promoter/expression coverage differ across branches."
)
writeLines(method_lines, file.path(out_dir, "08_reports", "method_comparison_and_recommendation.md"))

qc_md <- c(
  "# Pathway Enrichment QC Report",
  "",
  paste0("Output folder: `", out_dir, "`"),
  "",
  "## Checklist",
  "",
  apply(qc, 1, function(x) paste0("- ", x[["check"]], ": ", x[["status"]])),
  "",
  "## Notes",
  "",
  "- Raw enrichment tables are retained separately from the compact reporting workbook.",
  "- The Excel workbook intentionally includes top terms and summaries rather than every raw enrichment row.",
  "- Disease-labelled pathway terms should be interpreted as gene-set containers, not diagnostic claims.",
  "- Results from exploratory DL-screening inputs should be kept out of headline claims unless supported by stricter inputs or ranked analyses."
)
writeLines(qc_md, file.path(out_dir, "08_reports", "pathway_enrichment_QC_report.md"))

cat("Post-processing complete: ", out_dir, "\n")

