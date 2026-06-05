#!/usr/bin/env Rscript

# Continuous/ranked methylation-expression convergence analysis.
#
# This is a sensitivity layer built from final methylation and expression
# meta-analysis outputs. It avoids relying only on hard overlap categories by
# combining per-gene methylation and expression evidence within tissue-matched
# tested backgrounds, then testing pathway enrichment of the ranked combined
# evidence.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(scales)
  library(viridis)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = FALSE)
env_root <- Sys.getenv("ASD_REPO_ROOT", unset = "")
package_root <- if (nzchar(env_root)) normalizePath(env_root, winslash = "/", mustWork = TRUE) else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source_pkg <- file.path(package_root, "results", "pathway_enrichment")
source_input_dir <- file.path(source_pkg, "01_enrichment_inputs")
out_dir <- file.path(package_root, "results", "continuous_convergence_enrichment")
for (d in c("01_inputs", "02_gene_level_convergence", "03_pathway_enrichment", "04_figures", "05_reports", "06_quality_control")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

copy_if <- function(src, dst) if (file.exists(src) && !file.exists(dst)) file.copy(src, dst, overwrite = FALSE)
copy_if(file.path(source_input_dir, "standardised_model_gene_results.csv"),
        file.path(out_dir, "01_inputs/standardised_model_gene_results.csv"))
copy_if(file.path(source_input_dir, "msigdb_term2gene_used.csv"),
        file.path(out_dir, "01_inputs/msigdb_term2gene_used.csv"))

res <- fread(file.path(out_dir, "01_inputs/standardised_model_gene_results.csv"))
term2gene <- fread(file.path(out_dir, "01_inputs/msigdb_term2gene_used.csv"))
term2gene <- unique(term2gene[!is.na(symbol) & symbol != "" & !is.na(term) & term != "", .(database, term, symbol)])

theme_for_term <- function(term) {
  t <- tolower(term)
  fifelse(grepl("mitochond|oxidative|respiratory|electron transport|complex i|complex iii|complex iv|tca|atp synth", t),
          "mitochondrial respiration / oxidative phosphorylation",
  fifelse(grepl("immune|interferon|mhc|antigen|cytokine|complement|microglia|inflamm|t cell|lymphocyte|hla", t),
          "immune / antigen processing / MHC",
  fifelse(grepl("synap|neuro|glutamat|gaba|cholin|calcium|axon|dendrit|postsynaptic|presynaptic", t),
          "synaptic / neuronal signalling",
  fifelse(grepl("ribosom|translation|trna|rrna|protein synthesis", t),
          "ribosomal / translation",
  fifelse(grepl("cell cycle|mitotic|prolifer|chromosome|spindle|dna replication", t),
          "cell cycle / proliferation",
  fifelse(grepl("proteas|ubiquitin|unfolded|protein folding|er stress|endoplasmic", t),
          "proteostasis / protein processing",
  fifelse(grepl("wnt|adhesion|cadherin|extracellular matrix|ecm|integrin", t),
          "WNT / adhesion / extracellular matrix",
  fifelse(grepl("chromatin|histone|methyl|acetyl|epigen", t),
          "chromatin / epigenetic regulation",
  fifelse(grepl("metabolism|metabolic|amino acid|fatty acid|lipid|glycol|nucleotide", t),
          "metabolism", "other")))))))))
}

select_primary_pair <- function(tissue_name) {
  if (tissue_name == "blood") {
    list(
      methylation = res[tissue == "blood" & omic == "methylation" & model == "blood_array_peripheral_primary"],
      expression = res[tissue == "blood" & omic == "expression" & model == "blood_expression_peripheral_primary"]
    )
  } else if (tissue_name == "brain") {
    list(
      methylation = res[tissue == "brain" & omic == "methylation" & model == "brain_grouped_primary_with_WGBS"],
      expression = res[tissue == "brain" & omic == "expression" & model == "brain_expression_grouped_primary_public_R"]
    )
  } else {
    stop("Unsupported tissue: ", tissue_name)
  }
}

signed_z <- function(effect, p) {
  p <- pmin(pmax(p, 1e-300), 1)
  sign(effect) * qnorm(p / 2, lower.tail = FALSE)
}

combine_pair <- function(tissue_name) {
  pair <- select_primary_pair(tissue_name)
  m <- pair$methylation
  e <- pair$expression
  m <- m[!duplicated(gene)]
  e <- e[!duplicated(gene)]
  joined <- merge(
    m[, .(gene, methylation_model = model, methylation_effect = effect, methylation_p = p_value,
          methylation_FDR = FDR, methylation_DL = DL_nonzero, methylation_mKH = mKH_supported)],
    e[, .(gene, expression_model = model, expression_effect = effect, expression_p = p_value,
          expression_FDR = FDR, expression_DL = DL_nonzero, expression_mKH = mKH_supported)],
    by = "gene"
  )
  joined[, tissue := tissue_name]
  joined[, methylation_z := signed_z(methylation_effect, methylation_p)]
  joined[, expression_z := signed_z(expression_effect, expression_p)]
  joined[, fisher_chisq := -2 * (log(pmax(methylation_p, 1e-300)) + log(pmax(expression_p, 1e-300)))]
  joined[, fisher_p := pchisq(fisher_chisq, df = 4, lower.tail = FALSE)]
  joined[, fisher_FDR := p.adjust(fisher_p, method = "BH")]
  joined[, combined_score := -log10(pmax(fisher_p, 1e-300))]
  joined[, combined_FDR_score := -log10(pmax(fisher_FDR, 1e-300))]
  joined[, direction_product := methylation_z * expression_z]
  joined[, abs_direction_product := abs(direction_product)]
  joined[, direction_class := fifelse(direction_product > 0, "same_direction",
                                      fifelse(direction_product < 0, "opposite_direction", "neutral"))]
  joined[, same_direction_score := fifelse(direction_class == "same_direction", combined_score, -combined_score)]
  joined[, opposite_direction_score := fifelse(direction_class == "opposite_direction", combined_score, -combined_score)]
  joined[, stringent_both_FDR := methylation_FDR < 0.05 & expression_FDR < 0.05]
  joined[, both_mKH := methylation_mKH & expression_mKH]
  joined[, both_DL := methylation_DL & expression_DL]
  setorder(joined, fisher_p, -abs_direction_product)
  joined
}

run_ora <- function(genes, background, universe_term2gene, label, direction = "all") {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  background <- unique(background[!is.na(background) & background != ""])
  t2g <- unique(universe_term2gene[symbol %in% background])
  terms <- split(t2g$symbol, paste(t2g$database, t2g$term, sep = "::"))
  out <- rbindlist(lapply(names(terms), function(term_id) {
    gs <- unique(terms[[term_id]])
    a <- length(intersect(genes, gs))
    if (a < 3) return(NULL)
    b <- length(genes) - a
    c <- length(setdiff(gs, genes))
    d <- length(background) - a - b - c
    if (d < 0) return(NULL)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    parts <- strsplit(term_id, "::", fixed = TRUE)[[1]]
    data.table(
      analysis_id = label,
      direction_set = direction,
      database = parts[1],
      term = parts[2],
      overlap = a,
      input_genes = length(genes),
      background_genes = length(background),
      p_value = ft$p.value,
      genes = paste(sort(intersect(genes, gs)), collapse = ";")
    )
  }), fill = TRUE)
  if (!nrow(out)) return(out)
  out[, adjusted_p := p.adjust(p_value, method = "BH")]
  out[, theme := theme_for_term(term)]
  setorder(out, adjusted_p, p_value)
  out
}

run_fgsea_one <- function(stats, label, direction = "combined") {
  stats <- stats[is.finite(stats)]
  stats <- stats[!duplicated(names(stats))]
  pathways <- split(term2gene$symbol, paste(term2gene$database, term2gene$term, sep = "::"))
  pathways <- lapply(pathways, unique)
  fg <- suppressWarnings(fgsea::fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500, eps = 0))
  fg <- as.data.table(fg)
  if (!nrow(fg)) return(fg)
  fg[, c("database", "term") := tstrsplit(pathway, "::", fixed = TRUE)]
  fg[, analysis_id := label]
  fg[, direction_set := direction]
  fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  fg[, theme := theme_for_term(term)]
  setorder(fg, padj, pval)
  fg
}

all_convergence <- rbindlist(lapply(c("blood", "brain"), combine_pair), fill = TRUE)
fwrite(all_convergence, file.path(out_dir, "02_gene_level_convergence/continuous_gene_level_convergence_scores.csv"))

gene_summary <- all_convergence[, .(
  shared_background_genes = .N,
  both_FDR = sum(stringent_both_FDR, na.rm = TRUE),
  both_mKH = sum(both_mKH, na.rm = TRUE),
  both_DL = sum(both_DL, na.rm = TRUE),
  fisher_FDR_05 = sum(fisher_FDR < 0.05, na.rm = TRUE),
  fisher_p_001 = sum(fisher_p < 0.001, na.rm = TRUE),
  same_direction = sum(direction_class == "same_direction", na.rm = TRUE),
  opposite_direction = sum(direction_class == "opposite_direction", na.rm = TRUE),
  median_combined_score = median(combined_score, na.rm = TRUE),
  top_1pct_n = ceiling(.N * 0.01),
  top_5pct_n = ceiling(.N * 0.05)
), by = tissue]
fwrite(gene_summary, file.path(out_dir, "02_gene_level_convergence/continuous_convergence_summary_by_tissue.csv"))

top_genes <- all_convergence[, head(.SD, 100), by = tissue]
fwrite(top_genes, file.path(out_dir, "02_gene_level_convergence/top_100_continuous_convergence_genes_by_tissue.csv"))

ora_all <- list()
fgsea_all <- list()
for (tiss in c("blood", "brain")) {
  dt <- all_convergence[tissue == tiss]
  bg <- dt$gene
  n5 <- ceiling(nrow(dt) * 0.05)
  n1 <- ceiling(nrow(dt) * 0.01)
  dt_order <- dt[order(fisher_p, -abs_direction_product)]
  same_order <- dt[direction_class == "same_direction"][order(fisher_p)]
  opp_order <- dt[direction_class == "opposite_direction"][order(fisher_p)]
  ora_all[[paste0(tiss, "_top5")]] <- run_ora(head(dt_order$gene, n5), bg, term2gene, paste0(tiss, "__continuous_convergence_top5pct"), "top_5pct_any_direction")
  ora_all[[paste0(tiss, "_top1")]] <- run_ora(head(dt_order$gene, n1), bg, term2gene, paste0(tiss, "__continuous_convergence_top1pct"), "top_1pct_any_direction")
  ora_all[[paste0(tiss, "_same5")]] <- run_ora(head(same_order$gene, min(n5, nrow(same_order))), bg, term2gene, paste0(tiss, "__continuous_convergence_same_direction_top5pct"), "same_direction_top_5pct")
  ora_all[[paste0(tiss, "_opp5")]] <- run_ora(head(opp_order$gene, min(n5, nrow(opp_order))), bg, term2gene, paste0(tiss, "__continuous_convergence_opposite_direction_top5pct"), "opposite_direction_top_5pct")
  comb_stats <- dt$combined_score
  names(comb_stats) <- dt$gene
  same_stats <- dt$same_direction_score
  names(same_stats) <- dt$gene
  opp_stats <- dt$opposite_direction_score
  names(opp_stats) <- dt$gene
  fgsea_all[[paste0(tiss, "_combined")]] <- run_fgsea_one(comb_stats, paste0(tiss, "__continuous_convergence_fgsea_combined"), "combined_evidence")
  fgsea_all[[paste0(tiss, "_same")]] <- run_fgsea_one(same_stats, paste0(tiss, "__continuous_convergence_fgsea_same_direction"), "same_direction_signed")
  fgsea_all[[paste0(tiss, "_opposite")]] <- run_fgsea_one(opp_stats, paste0(tiss, "__continuous_convergence_fgsea_opposite_direction"), "opposite_direction_signed")
}

ora <- rbindlist(ora_all, fill = TRUE)
fgsea_res <- rbindlist(fgsea_all, fill = TRUE)
fwrite(ora, file.path(out_dir, "03_pathway_enrichment/continuous_convergence_ORA_results.csv"))
fwrite(fgsea_res, file.path(out_dir, "03_pathway_enrichment/continuous_convergence_fgsea_results.csv"))

theme_summary_ora <- if (nrow(ora)) {
  ora[adjusted_p < 0.05, .(
    significant_terms = .N,
    best_adjusted_p = min(adjusted_p),
    representative_terms = paste(head(term[order(adjusted_p)], 8), collapse = " | ")
  ), by = .(analysis_id, direction_set, theme)]
} else data.table()

theme_summary_fgsea <- if (nrow(fgsea_res)) {
  fgsea_res[padj < 0.05, .(
    significant_terms = .N,
    best_adjusted_p = min(padj),
    max_abs_NES = max(abs(NES), na.rm = TRUE),
    representative_terms = paste(head(term[order(padj)], 8), collapse = " | ")
  ), by = .(analysis_id, direction_set, theme)]
} else data.table()

theme_summary_ora[, method := "top-percentile ORA"]
theme_summary_fgsea[, method := "preranked fgsea"]
theme_summary <- rbindlist(list(theme_summary_ora, theme_summary_fgsea), fill = TRUE)
theme_summary[, tissue := fifelse(grepl("^blood", analysis_id), "blood", "brain")]
setorder(theme_summary, tissue, best_adjusted_p)
fwrite(theme_summary, file.path(out_dir, "03_pathway_enrichment/continuous_convergence_theme_summary.csv"))

# Compact plots
plot_theme <- theme_summary[theme != "other"]
plot_theme[, score := pmin(-log10(pmax(best_adjusted_p, 1e-300)), 80)]
plot_theme[, tissue := factor(tissue, levels = c("blood", "brain"), labels = c("Blood", "Post-mortem brain"))]
plot_theme[, method := factor(method, levels = c("top-percentile ORA", "preranked fgsea"))]
plot_theme[, theme := fct_reorder(theme, score, .fun = max)]

p_theme <- ggplot(plot_theme, aes(x = score, y = theme, fill = direction_set)) +
  geom_col(position = position_dodge2(width = 0.75, preserve = "single"), width = 0.68, colour = "white", linewidth = 0.2) +
  facet_grid(tissue ~ method, scales = "free_y", space = "free_y") +
  scale_fill_viridis_d(option = "turbo", end = 0.86, name = "Convergence input") +
  labs(title = "Continuous cross-omic convergence highlights pathway-level signal beyond strict gene overlap",
       subtitle = "Pathways are tested using Fisher-combined methylation-expression gene evidence and ranked/top-percentile convergence sets.",
       x = expression(-log[10]("adjusted p")), y = NULL,
       caption = "This sensitivity analysis tests summary-level tissue-matched convergence; it is not participant-level regulatory coupling.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = "#4b5563"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "04_figures/Figure_continuous_convergence_theme_summary.png"), p_theme, width = 12.5, height = 8.5, dpi = 320)
ggsave(file.path(out_dir, "04_figures/Figure_continuous_convergence_theme_summary.pdf"), p_theme, width = 12.5, height = 8.5, device = cairo_pdf)

top_term <- rbindlist(list(
  ora[adjusted_p < 0.05, .(method = "top-percentile ORA", analysis_id, direction_set, database, term, adjusted_p, theme, overlap, NES = NA_real_)],
  fgsea_res[padj < 0.05, .(method = "preranked fgsea", analysis_id, direction_set, database, term, adjusted_p = padj, theme, overlap = size, NES)]
), fill = TRUE)
top_term[, tissue := fifelse(grepl("^blood", analysis_id), "Blood", "Post-mortem brain")]
top_term <- top_term[theme != "other"]
top_term[, term_short := str_to_sentence(gsub("_", " ", gsub("REACTOME_|WP_|HALLMARK_", "", term)))]
top_term[, term_short := str_trunc(term_short, 58)]
top_plot <- top_term[order(tissue, adjusted_p), head(.SD, 20), by = tissue]
top_plot[, score := pmin(-log10(pmax(adjusted_p, 1e-300)), 80)]
top_plot[, term_short := fct_reorder(term_short, score)]
fwrite(top_term, file.path(out_dir, "03_pathway_enrichment/continuous_convergence_significant_terms_combined.csv"))

p_terms <- ggplot(top_plot, aes(x = score, y = term_short)) +
  geom_segment(aes(x = 0, xend = score, yend = term_short, colour = theme), alpha = 0.55, linewidth = 0.6) +
  geom_point(aes(size = overlap, fill = theme), shape = 21, colour = "white", stroke = 0.3) +
  facet_wrap(~ tissue, scales = "free_y", ncol = 1) +
  scale_fill_viridis_d(option = "turbo", end = 0.86, name = "Theme") +
  scale_colour_viridis_d(option = "turbo", end = 0.86, guide = "none") +
  scale_size_continuous(range = c(2, 7), name = "Genes / pathway size") +
  labs(title = "Top continuous methylation-expression convergence pathway terms",
       subtitle = "Terms are ranked by adjusted p-value across top-percentile ORA and preranked fgsea sensitivity analyses.",
       x = expression(-log[10]("adjusted p")), y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = "#4b5563"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "04_figures/Figure_top_continuous_convergence_terms.png"), p_terms, width = 11.5, height = 10.5, dpi = 320)
ggsave(file.path(out_dir, "04_figures/Figure_top_continuous_convergence_terms.pdf"), p_terms, width = 11.5, height = 10.5, device = cairo_pdf)

openxlsx::write.xlsx(
  list(
    gene_summary = as.data.frame(gene_summary),
    top_genes = as.data.frame(top_genes),
    ORA = as.data.frame(ora),
    fgsea = as.data.frame(fgsea_res),
    theme_summary = as.data.frame(theme_summary),
    significant_terms = as.data.frame(top_term)
  ),
  file.path(out_dir, "continuous_cross_omic_convergence_enrichment_results.xlsx"),
  overwrite = TRUE
)

report_lines <- c(
  "# Continuous Cross-Omic Convergence Enrichment Report",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Why This Was Run",
  "",
  "The thresholded convergence pathway analysis relies mainly on genes passing support criteria in both methylation and expression layers. That is conservative but can miss pathway-level convergence distributed across many genes with modest same-tissue evidence. This sensitivity analysis uses continuous gene-level evidence by Fisher-combining methylation and expression p-values within tissue-matched tested backgrounds, then tests ranked and top-percentile pathway enrichment.",
  "",
  "## Gene-Level Continuous Convergence Summary",
  "",
  paste(capture.output(print(gene_summary)), collapse = "\n"),
  "",
  "## Significant Pathway Themes",
  "",
  paste(capture.output(print(theme_summary[order(tissue, best_adjusted_p)])), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "This analysis should be treated as a sensitivity/prioritisation layer rather than replacement for strict gene-overlap convergence. Stronger results here mean that pathway-level methylation-expression convergence is distributed across ranked evidence rather than concentrated in many genes passing both omic thresholds.",
  "",
  "## Recommended Reporting",
  "",
  "- Report strict overlap as limited.",
  "- Add that continuous/ranked convergence sensitivity supports pathway-level convergence, especially where themes match the primary expression/methylation enrichment story.",
  "- Avoid claiming direct regulatory coupling because datasets are not paired within individuals."
)
writeLines(report_lines, file.path(out_dir, "05_reports/continuous_cross_omic_convergence_enrichment_report.md"))

qc <- c(
  "# Continuous Cross-Omic Convergence QC",
  "",
  "- Used final standardised model-gene results from results/pathway_enrichment.",
  "- Used tissue-matched primary blood and brain methylation/expression model pairs.",
  "- Used the shared tested genes in each tissue as the convergence background.",
  "- Fisher-combined gene evidence and ranked fgsea were used as sensitivity methods, not replacement primary analyses.",
  "- Top-percentile ORA used the tissue-specific shared background.",
  "- Outputs are written only to the configured convergence-enrichment results folder."
)
writeLines(qc, file.path(out_dir, "06_quality_control/continuous_convergence_QC_report.md"))

pkg_versions <- data.table(
  package = c("data.table", "fgsea", "ggplot2", "forcats", "stringr", "scales", "viridis", "openxlsx"),
  version = vapply(c("data.table", "fgsea", "ggplot2", "forcats", "stringr", "scales", "viridis", "openxlsx"),
                   function(p) as.character(packageVersion(p)), character(1))
)
fwrite(pkg_versions, file.path(out_dir, "06_quality_control/R_package_versions.csv"))

file_index <- data.table(file_path = normalizePath(list.files(out_dir, recursive = TRUE, full.names = TRUE), winslash = "/", mustWork = FALSE))
file_index[, bytes := file.info(file_path)$size]
fwrite(file_index, file.path(out_dir, "06_quality_control/file_index.csv"))

readme <- c(
  "# Continuous Cross-Omic Convergence Enrichment",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript scripts/05_convergence_pathway/run_continuous_cross_omic_convergence_enrichment.R",
  "```",
  "",
  "This is a sensitivity analysis designed to complement, not replace, strict methylation-expression overlap tests."
)
writeLines(readme, file.path(out_dir, "README_reproducibility.md"))

message("Finished continuous convergence enrichment package: ", out_dir)

