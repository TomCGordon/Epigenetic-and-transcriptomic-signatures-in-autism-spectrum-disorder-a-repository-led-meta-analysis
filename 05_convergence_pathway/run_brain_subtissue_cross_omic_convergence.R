#!/usr/bin/env Rscript

# Brain subtissue methylation-expression convergence sensitivity analysis.
#
# This script uses final audited methylation and expression model-gene outputs
# to test whether pathway-level cross-omic convergence is driven only by the
# grouped post-mortem brain model or is also visible within region/subtissue
# sensitivity models. The grouped model remains an umbrella synthesis; the
# subtissue pairs below are secondary/sensitivity comparisons.

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
package_root <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source_pkg <- file.path(package_root, "results", "pathway_enrichment")
source_input_dir <- file.path(source_pkg, "01_enrichment_inputs")
gse36315_pkg <- file.path(package_root, "pipelines", "gse36315_brain_expression_sensitivity")
out_dir <- file.path(package_root, "results", "brain_subtissue_convergence")

for (d in c("01_inputs", "02_gene_level_convergence", "03_pathway_enrichment",
            "04_figures", "05_reports", "06_quality_control")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

copy_if <- function(src, dst) {
  if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}

copy_if(file.path(source_input_dir, "standardised_model_gene_results.csv"),
        file.path(out_dir, "01_inputs/standardised_model_gene_results.csv"))
copy_if(file.path(source_input_dir, "msigdb_term2gene_used.csv"),
        file.path(out_dir, "01_inputs/msigdb_term2gene_used.csv"))

res <- fread(file.path(out_dir, "01_inputs/standardised_model_gene_results.csv"))
term2gene <- fread(file.path(out_dir, "01_inputs/msigdb_term2gene_used.csv"))
term2gene <- unique(term2gene[!is.na(symbol) & symbol != "" & !is.na(term) & term != "",
                              .(database, term, symbol)])

gse36315_region_file <- file.path(gse36315_pkg, "04_region_subtissue_sensitivity/region_results.csv")
if (file.exists(gse36315_region_file)) {
  gse <- fread(gse36315_region_file)
  gse_std <- gse[, .(
    source_file = normalizePath(gse36315_region_file, winslash = "/", mustWork = FALSE),
    source_label = "gse36315_custom_annotated_expression_sensitivity",
    omic = "expression",
    tissue = "brain",
    gene,
    model,
    role = model_role,
    effect = pooled_g,
    p_value,
    FDR,
    k,
    I2,
    DL_nonzero = as.logical(DL_nonzero),
    mKH_supported = as.logical(mKH_interval_excludes_zero),
    direction,
    analysis_id = paste("brain", "expression", model, sep = "__"),
    signed_rank = sign(pooled_g) * -log10(pmax(p_value, 1e-300))
  )]
  res <- rbindlist(list(res, gse_std), fill = TRUE)
}

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

signed_z <- function(effect, p) {
  p <- pmin(pmax(p, 1e-300), 1)
  sign(effect) * qnorm(p / 2, lower.tail = FALSE)
}

pair_defs <- data.table(
  pair_id = c(
    "grouped_postmortem_brain_reference",
    "cortex_sensitivity",
    "prefrontal_BA9_sensitivity",
    "cerebellum_sensitivity",
    "cortex_plus_gse36315_custom_annotated",
    "prefrontal_plus_gse36315_custom_annotated",
    "cerebellum_plus_gse36315_custom_annotated"
  ),
  evidence_layer = c(
    "primary grouped reference",
    "official subtissue sensitivity",
    "official subtissue sensitivity",
    "official subtissue sensitivity",
    "custom-annotated expression sensitivity",
    "custom-annotated expression sensitivity",
    "custom-annotated expression sensitivity"
  ),
  methylation_model = c(
    "brain_grouped_primary_with_WGBS",
    "cortex_only_sensitivity",
    "prefrontal_BA9_cortex_sensitivity",
    "cerebellum_only_sensitivity",
    "cortex_only_sensitivity",
    "prefrontal_BA9_cortex_sensitivity",
    "cerebellum_only_sensitivity"
  ),
  expression_model = c(
    "brain_expression_grouped_primary_public_R",
    "brain_expression_cortex_only_sensitivity_R",
    "brain_expression_prefrontal_cortex_sensitivity_R",
    "brain_expression_cerebellum_only_sensitivity_R",
    "brain_expression_cortex_plus_GSE36315_prefrontal_custom_annotation_sensitivity_R",
    "brain_expression_prefrontal_plus_gse36315_custom_annotated_sensitivity_R",
    "brain_expression_cerebellum_plus_gse36315_custom_annotated_sensitivity_R"
  ),
  interpretation_caveat = c(
    "Umbrella post-mortem brain synthesis; not a claim of anatomical homogeneity.",
    "Cortical sensitivity model; still combines cortical regions and platforms.",
    "Prefrontal sensitivity model; methylation is BA9-focused while expression is broader prefrontal cortex.",
    "Cerebellum sensitivity model; treated as region-specific sensitivity, not whole-brain evidence.",
    "Adds custom-annotated GSE36315 prefrontal expression to cortex sensitivity; custom probe-remapping layer.",
    "Adds custom-annotated GSE36315 prefrontal expression; custom probe-remapping layer.",
    "Adds custom-annotated GSE36315 cerebellum expression; custom probe-remapping layer."
  )
)

fwrite(pair_defs, file.path(out_dir, "01_inputs/brain_subtissue_pair_definitions.csv"))

model_summary <- res[tissue == "brain", .(
  genes = .N,
  FDR_significant = sum(FDR < 0.05, na.rm = TRUE),
  mKH_supported = sum(mKH_supported == TRUE, na.rm = TRUE),
  DL_nonzero = sum(DL_nonzero == TRUE, na.rm = TRUE),
  median_I2 = as.numeric(median(as.numeric(I2), na.rm = TRUE)),
  median_k = as.numeric(median(as.numeric(k), na.rm = TRUE))
), by = .(omic, model, role)]
fwrite(model_summary, file.path(out_dir, "01_inputs/brain_model_summary_available_for_subtissue_check.csv"))

combine_pair <- function(row) {
  m <- res[tissue == "brain" & omic == "methylation" & model == row$methylation_model]
  e <- res[tissue == "brain" & omic == "expression" & model == row$expression_model]
  if (!nrow(m) || !nrow(e)) return(NULL)
  m <- m[!duplicated(gene)]
  e <- e[!duplicated(gene)]
  joined <- merge(
    m[, .(gene, methylation_model = model, methylation_effect = effect,
          methylation_p = p_value, methylation_FDR = FDR,
          methylation_DL = DL_nonzero, methylation_mKH = mKH_supported,
          methylation_k = k, methylation_I2 = I2)],
    e[, .(gene, expression_model = model, expression_effect = effect,
          expression_p = p_value, expression_FDR = FDR,
          expression_DL = DL_nonzero, expression_mKH = mKH_supported,
          expression_k = k, expression_I2 = I2)],
    by = "gene"
  )
  if (!nrow(joined)) return(NULL)
  joined[, pair_id := row$pair_id]
  joined[, evidence_layer := row$evidence_layer]
  joined[, interpretation_caveat := row$interpretation_caveat]
  joined[, methylation_z := signed_z(methylation_effect, methylation_p)]
  joined[, expression_z := signed_z(expression_effect, expression_p)]
  joined[, fisher_chisq := -2 * (log(pmax(methylation_p, 1e-300)) + log(pmax(expression_p, 1e-300)))]
  joined[, fisher_p := pchisq(fisher_chisq, df = 4, lower.tail = FALSE)]
  joined[, fisher_FDR := p.adjust(fisher_p, method = "BH")]
  joined[, combined_score := -log10(pmax(fisher_p, 1e-300))]
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

all_convergence <- rbindlist(lapply(seq_len(nrow(pair_defs)), function(i) combine_pair(pair_defs[i])), fill = TRUE)
fwrite(all_convergence, file.path(out_dir, "02_gene_level_convergence/brain_subtissue_convergence_gene_scores.csv"))

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
), by = .(pair_id, evidence_layer)]
fwrite(gene_summary, file.path(out_dir, "02_gene_level_convergence/brain_subtissue_convergence_summary.csv"))

top_genes <- all_convergence[, head(.SD, 100), by = pair_id]
fwrite(top_genes, file.path(out_dir, "02_gene_level_convergence/top_100_brain_subtissue_convergence_genes.csv"))

run_ora <- function(genes, background, label, direction) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  background <- unique(background[!is.na(background) & background != ""])
  t2g <- unique(term2gene[symbol %in% background])
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
      pair_id = label,
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

run_fgsea_one <- function(stats, label, direction) {
  stats <- stats[is.finite(stats)]
  stats <- stats[!duplicated(names(stats))]
  pathways <- split(term2gene$symbol, paste(term2gene$database, term2gene$term, sep = "::"))
  pathways <- lapply(pathways, unique)
  fg <- suppressWarnings(fgsea::fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500, eps = 0))
  fg <- as.data.table(fg)
  if (!nrow(fg)) return(fg)
  fg[, c("database", "term") := tstrsplit(pathway, "::", fixed = TRUE)]
  fg[, pair_id := label]
  fg[, direction_set := direction]
  fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  fg[, theme := theme_for_term(term)]
  setorder(fg, padj, pval)
  fg
}

ora_all <- list()
fgsea_all <- list()
for (pid in unique(all_convergence$pair_id)) {
  dt <- all_convergence[pair_id == pid]
  bg <- dt$gene
  n5 <- ceiling(nrow(dt) * 0.05)
  n1 <- ceiling(nrow(dt) * 0.01)
  dt_order <- dt[order(fisher_p, -abs_direction_product)]
  same_order <- dt[direction_class == "same_direction"][order(fisher_p)]
  opp_order <- dt[direction_class == "opposite_direction"][order(fisher_p)]
  ora_all[[paste0(pid, "_top5")]] <- run_ora(head(dt_order$gene, n5), bg, pid, "top_5pct_any_direction")
  ora_all[[paste0(pid, "_top1")]] <- run_ora(head(dt_order$gene, n1), bg, pid, "top_1pct_any_direction")
  ora_all[[paste0(pid, "_same5")]] <- run_ora(head(same_order$gene, min(n5, nrow(same_order))), bg, pid, "same_direction_top_5pct")
  ora_all[[paste0(pid, "_opp5")]] <- run_ora(head(opp_order$gene, min(n5, nrow(opp_order))), bg, pid, "opposite_direction_top_5pct")

  comb_stats <- dt$combined_score
  names(comb_stats) <- dt$gene
  same_stats <- dt$same_direction_score
  names(same_stats) <- dt$gene
  opp_stats <- dt$opposite_direction_score
  names(opp_stats) <- dt$gene
  fgsea_all[[paste0(pid, "_combined")]] <- run_fgsea_one(comb_stats, pid, "combined_evidence")
  fgsea_all[[paste0(pid, "_same")]] <- run_fgsea_one(same_stats, pid, "same_direction_signed")
  fgsea_all[[paste0(pid, "_opposite")]] <- run_fgsea_one(opp_stats, pid, "opposite_direction_signed")
}

ora <- rbindlist(ora_all, fill = TRUE)
fgsea_res <- rbindlist(fgsea_all, fill = TRUE)
fwrite(ora, file.path(out_dir, "03_pathway_enrichment/brain_subtissue_ORA_results.csv"))
fwrite(fgsea_res, file.path(out_dir, "03_pathway_enrichment/brain_subtissue_fgsea_results.csv"))

theme_summary_ora <- if (nrow(ora)) {
  ora[adjusted_p < 0.05, .(
    significant_terms = .N,
    best_adjusted_p = min(adjusted_p),
    representative_terms = paste(head(term[order(adjusted_p)], 8), collapse = " | ")
  ), by = .(pair_id, direction_set, theme)]
} else data.table()

theme_summary_fgsea <- if (nrow(fgsea_res)) {
  fgsea_res[padj < 0.05, .(
    significant_terms = .N,
    best_adjusted_p = min(padj),
    max_abs_NES = max(abs(NES), na.rm = TRUE),
    representative_terms = paste(head(term[order(padj)], 8), collapse = " | ")
  ), by = .(pair_id, direction_set, theme)]
} else data.table()

theme_summary_ora[, method := "top-percentile ORA"]
theme_summary_fgsea[, method := "preranked fgsea"]
theme_summary <- rbindlist(list(theme_summary_ora, theme_summary_fgsea), fill = TRUE)
theme_summary <- merge(theme_summary, pair_defs[, .(pair_id, evidence_layer, interpretation_caveat)], by = "pair_id", all.x = TRUE)
setorder(theme_summary, pair_id, best_adjusted_p)
fwrite(theme_summary, file.path(out_dir, "03_pathway_enrichment/brain_subtissue_theme_summary.csv"))

top_terms <- rbindlist(list(
  ora[adjusted_p < 0.05, .(method = "top-percentile ORA", pair_id, direction_set, database, term,
                           adjusted_p, theme, overlap, NES = NA_real_)],
  fgsea_res[padj < 0.05, .(method = "preranked fgsea", pair_id, direction_set, database, term,
                            adjusted_p = padj, theme, overlap = size, NES)]
), fill = TRUE)
top_terms <- merge(top_terms, pair_defs[, .(pair_id, evidence_layer, interpretation_caveat)], by = "pair_id", all.x = TRUE)
fwrite(top_terms, file.path(out_dir, "03_pathway_enrichment/brain_subtissue_significant_terms.csv"))

theme_presence <- theme_summary[theme != "other", .(
  significant_methods = paste(sort(unique(method)), collapse = "; "),
  significant_direction_sets = paste(sort(unique(direction_set)), collapse = "; "),
  term_count = sum(significant_terms, na.rm = TRUE),
  best_adjusted_p = min(best_adjusted_p, na.rm = TRUE),
  representative_terms = paste(unique(unlist(strsplit(paste(representative_terms, collapse = " | "), " \\| "))), collapse = " | ")
), by = .(pair_id, evidence_layer, theme)]
theme_presence[, score := pmin(-log10(pmax(best_adjusted_p, 1e-300)), 80)]
setorder(theme_presence, pair_id, best_adjusted_p)
fwrite(theme_presence, file.path(out_dir, "03_pathway_enrichment/brain_subtissue_theme_presence_matrix.csv"))

pair_levels <- c(
  "grouped_postmortem_brain_reference",
  "cortex_sensitivity",
  "prefrontal_BA9_sensitivity",
  "cerebellum_sensitivity",
  "cortex_plus_gse36315_custom_annotated",
  "prefrontal_plus_gse36315_custom_annotated",
  "cerebellum_plus_gse36315_custom_annotated"
)
plot_theme <- theme_presence[theme != "other"]
plot_theme[, pair_id := factor(pair_id, levels = pair_levels)]
plot_theme[, theme := fct_reorder(theme, score, .fun = max)]

p_heat <- ggplot(plot_theme, aes(x = pair_id, y = theme)) +
  geom_tile(aes(fill = score), colour = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(score >= 2, sprintf("%.1f", score), "")), size = 2.8, colour = "white") +
  scale_fill_viridis_c(option = "magma", name = expression(-log[10]("best adjusted p")), limits = c(0, max(plot_theme$score, na.rm = TRUE))) +
  scale_x_discrete(labels = c(
    "Grouped\nreference", "Cortex", "Prefrontal\nBA9", "Cerebellum",
    "Cortex +\nGSE36315", "Prefrontal +\nGSE36315", "Cerebellum +\nGSE36315"
  )) +
  labs(title = "Brain subtissue methylation-expression pathway convergence",
       subtitle = "Cells show significant pathway themes in paired methylation and expression sensitivity models.",
       x = NULL, y = NULL,
       caption = "GSE36315 models use custom probe remapping and should be treated as sensitivity only.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = "#4b5563"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "04_figures/Figure_brain_subtissue_convergence_theme_heatmap.png"), p_heat, width = 11.5, height = 6.5, dpi = 320)
ggsave(file.path(out_dir, "04_figures/Figure_brain_subtissue_convergence_theme_heatmap.pdf"), p_heat, width = 11.5, height = 6.5, device = cairo_pdf)

top_plot <- top_terms[theme != "other"][order(pair_id, adjusted_p), head(.SD, 10), by = pair_id]
if (nrow(top_plot)) {
  top_plot[, score := pmin(-log10(pmax(adjusted_p, 1e-300)), 80)]
  top_plot[, term_short := str_to_sentence(gsub("_", " ", gsub("REACTOME_|WP_|HALLMARK_", "", term)))]
  top_plot[, term_short := str_trunc(term_short, 58)]
  top_plot[, pair_id := factor(pair_id, levels = pair_levels)]
  top_plot[, term_short := fct_reorder(term_short, score)]

  p_terms <- ggplot(top_plot, aes(x = score, y = term_short)) +
    geom_segment(aes(x = 0, xend = score, yend = term_short, colour = theme), linewidth = 0.65, alpha = 0.65) +
    geom_point(aes(size = overlap, fill = theme), shape = 21, colour = "white", stroke = 0.3) +
    facet_wrap(~ pair_id, scales = "free_y", ncol = 2,
               labeller = labeller(pair_id = c(
                 grouped_postmortem_brain_reference = "Grouped reference",
                 cortex_sensitivity = "Cortex",
                 prefrontal_BA9_sensitivity = "Prefrontal BA9",
                 cerebellum_sensitivity = "Cerebellum",
                 cortex_plus_gse36315_custom_annotated = "Cortex + GSE36315",
                 prefrontal_plus_gse36315_custom_annotated = "Prefrontal + GSE36315",
                 cerebellum_plus_gse36315_custom_annotated = "Cerebellum + GSE36315"
               ))) +
    scale_fill_viridis_d(option = "turbo", end = 0.86, name = "Theme") +
    scale_colour_viridis_d(option = "turbo", end = 0.86, guide = "none") +
    scale_size_continuous(range = c(2, 6), name = "Overlap / size") +
    labs(title = "Top significant brain subtissue convergence terms",
         subtitle = "Terms are ranked within each paired model by adjusted p-value.",
         x = expression(-log[10]("adjusted p")), y = NULL) +
    theme_minimal(base_size = 9.5) +
    theme(panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(colour = "#4b5563"),
          legend.position = "bottom")
  ggsave(file.path(out_dir, "04_figures/Figure_brain_subtissue_top_terms.png"), p_terms, width = 13, height = 11, dpi = 320)
  ggsave(file.path(out_dir, "04_figures/Figure_brain_subtissue_top_terms.pdf"), p_terms, width = 13, height = 11, device = cairo_pdf)
}

openxlsx::write.xlsx(
  list(
    pair_definitions = as.data.frame(pair_defs),
    model_summary = as.data.frame(model_summary),
    convergence_summary = as.data.frame(gene_summary),
    top_100_genes = as.data.frame(top_genes),
    theme_presence = as.data.frame(theme_presence),
    significant_terms = as.data.frame(top_terms),
    ORA = as.data.frame(ora),
    fgsea = as.data.frame(fgsea_res)
  ),
  file.path(out_dir, "brain_subtissue_cross_omic_convergence_results.xlsx"),
  overwrite = TRUE
)

key_theme_compact <- theme_presence[theme %in% c(
  "mitochondrial respiration / oxidative phosphorylation",
  "immune / antigen processing / MHC",
  "synaptic / neuronal signalling",
  "ribosomal / translation",
  "proteostasis / protein processing",
  "cell cycle / proliferation"
), .(pair_id, theme, term_count, best_adjusted_p, significant_methods, representative_terms)]
setorder(key_theme_compact, pair_id, best_adjusted_p)

report_lines <- c(
  "# Brain Subtissue Cross-Omic Convergence Report",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "This sensitivity analysis asks whether brain methylation-expression pathway convergence is visible only in the grouped post-mortem brain model, or whether comparable themes appear within region/subtissue models. It uses final audited model-gene outputs and does not rerun methylation or expression extraction.",
  "",
  "## Paired Models",
  "",
  paste(capture.output(print(pair_defs)), collapse = "\n"),
  "",
  "## Gene-Level Convergence Summary",
  "",
  paste(capture.output(print(gene_summary)), collapse = "\n"),
  "",
  "## Key Significant Theme Summary",
  "",
  paste(capture.output(print(key_theme_compact)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "The grouped post-mortem brain model is defensible only as an umbrella synthesis designed to maximise power across sparse post-mortem brain datasets. It should not be described as evidence that all brain regions share the same ASD-associated molecular profile.",
  "",
  "Region/subtissue sensitivity analyses are therefore important. If the same pathway themes recur in cortex, prefrontal, and/or cerebellum sensitivity analyses, the grouped-brain pathway findings are less likely to depend on a single anatomical grouping. If themes differ by subtissue, that should be reported as biological and technical heterogeneity rather than treated as a problem to hide.",
  "",
  "The custom-annotated GSE36315 models are useful sensitivity checks but should remain supplementary because they rely on custom probe remapping of an incompletely annotated platform.",
  "",
  "## Suggested Reporting",
  "",
  "Use 'grouped post-mortem brain' rather than simply 'brain'. Report it as a primary umbrella model with region/subtissue sensitivity analyses. Avoid wording such as 'brain-wide' or 'region-general' unless directly supported by the sensitivity results.",
  "",
  "Suggested sentence: 'Because available post-mortem brain datasets varied by region, the grouped post-mortem brain model was interpreted as an umbrella synthesis rather than as evidence of anatomically homogeneous effects; cortex, prefrontal and cerebellum sensitivity analyses were used to assess whether the main pathway themes were region-dependent.'"
)
writeLines(report_lines, file.path(out_dir, "05_reports/brain_subtissue_cross_omic_convergence_report.md"))

recommendation_lines <- c(
  "# Grouped-Brain Reporting Recommendation",
  "",
  "## Bottom Line",
  "",
  "The grouped post-mortem brain models can be reported, but only with careful wording. They are defensible as power-preserving umbrella syntheses across sparse post-mortem brain datasets, not as claims of homogeneous ASD effects across all brain regions.",
  "",
  "## What To Report In Main Text",
  "",
  "- Use the grouped post-mortem brain model for the main broad brain estimate.",
  "- State that no grouped-brain primary methylation or expression model produced FDR-significant individual genes, while interval-supported and pathway-level signals were present.",
  "- Include a short sentence that region/subtissue sensitivity analyses were run and that detailed results are provided in supplementary tables.",
  "- Mention recurring pathway themes only where they persist beyond the grouped model.",
  "",
  "## What To Keep In Supplement",
  "",
  "- Full cortex, prefrontal, cerebellum, BA19, platform, and custom-annotated GSE36315 sensitivity tables.",
  "- Full subtissue pathway enrichment outputs.",
  "- Custom-annotated GSE36315 remapping details.",
  "",
  "## Wording To Avoid",
  "",
  "- Avoid 'brain-wide effects'.",
  "- Avoid 'region-general methylation/expression changes'.",
  "- Avoid implying participant-level methylation-expression coupling.",
  "- Avoid presenting GSE36315 custom-annotated results as primary evidence."
)
writeLines(recommendation_lines, file.path(out_dir, "05_reports/grouped_brain_reporting_recommendation.md"))

qc_lines <- c(
  "# Brain Subtissue Convergence QC",
  "",
  "- Used final standardised model-gene outputs from results/pathway_enrichment.",
  "- Used final MSigDB term-to-gene table from the same pathway package.",
  "- Official paired subtissue checks included cortex, prefrontal/BA9, and cerebellum.",
  "- Custom-annotated GSE36315 expression sensitivities were included only as labelled sensitivity layers.",
  "- Shared gene backgrounds were computed separately for each methylation-expression pair.",
  "- FDR, modified Knapp-Hartung, and DerSimonian-Laird support categories were not conflated.",
  "- Outputs were written only to the configured brain-subtissue convergence results directory."
)
writeLines(qc_lines, file.path(out_dir, "06_quality_control/brain_subtissue_convergence_QC_report.md"))

pkg_versions <- data.table(
  package = c("data.table", "fgsea", "ggplot2", "forcats", "stringr", "scales", "viridis", "openxlsx"),
  version = vapply(c("data.table", "fgsea", "ggplot2", "forcats", "stringr", "scales", "viridis", "openxlsx"),
                   function(p) as.character(packageVersion(p)), character(1))
)
fwrite(pkg_versions, file.path(out_dir, "06_quality_control/R_package_versions.csv"))

file_index <- data.table(file_path = normalizePath(list.files(out_dir, recursive = TRUE, full.names = TRUE),
                                                   winslash = "/", mustWork = FALSE))
file_index[, bytes := file.info(file_path)$size]
fwrite(file_index, file.path(out_dir, "06_quality_control/file_index.csv"))

readme <- c(
  "# Brain Subtissue Cross-Omic Convergence",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript scripts/05_convergence_pathway/run_brain_subtissue_cross_omic_convergence.R",
  "```",
  "",
  "This package tests whether grouped post-mortem brain convergence themes are also visible in region/subtissue sensitivity models."
)
writeLines(readme, file.path(out_dir, "README_reproducibility.md"))

message("Finished brain subtissue convergence package: ", out_dir)

