#!/usr/bin/env Rscript

# Enhanced pathway figures for ASD methylation/expression/convergence synthesis.
#
# Figures are generated from final pathway-enrichment output tables. The goal is
# to create publication-facing visuals that emphasise the tissue-by-omic biology
# without changing the underlying enrichment results.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(ggalluvial)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(patchwork)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = FALSE)
env_root <- Sys.getenv("ASD_REPO_ROOT", unset = "")
package_root <- if (nzchar(env_root)) normalizePath(env_root, winslash = "/", mustWork = TRUE) else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source_pkg <- file.path(package_root, "results", "pathway_enrichment")
source_report_dir <- file.path(source_pkg, "08_reports")
out_dir <- file.path(package_root, "results", "figures", "pathway")

for (d in c("01_source_tables", "02_figure_source_data", "03_figures", "04_reports", "05_quality_control")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

copy_if <- function(src, dest) if (file.exists(src) && !file.exists(dest)) file.copy(src, dest, overwrite = FALSE)
copy_if(file.path(source_report_dir, "all_enrichment_results_normalised.csv"),
        file.path(out_dir, "01_source_tables/all_enrichment_results_normalised.csv"))
copy_if(file.path(source_report_dir, "pathway_theme_summary.csv"),
        file.path(out_dir, "01_source_tables/pathway_theme_summary.csv"))
copy_if(file.path(source_report_dir, "top_pathway_themes_overall.csv"),
        file.path(out_dir, "01_source_tables/top_pathway_themes_overall.csv"))

all <- fread(file.path(out_dir, "01_source_tables/all_enrichment_results_normalised.csv"))
theme <- fread(file.path(out_dir, "01_source_tables/pathway_theme_summary.csv"))

theme_order <- c(
  "cell cycle / proliferation",
  "immune / antigen processing / MHC",
  "mitochondrial respiration / oxidative phosphorylation",
  "synaptic / neuronal signalling",
  "ribosomal / translation",
  "proteostasis / protein processing",
  "WNT / adhesion / extracellular matrix",
  "chromatin / epigenetic regulation",
  "metabolism"
)

theme_labels <- c(
  "cell cycle / proliferation" = "Cell cycle /\nproliferation",
  "immune / antigen processing / MHC" = "Immune /\nMHC",
  "mitochondrial respiration / oxidative phosphorylation" = "Mitochondrial\nrespiration",
  "synaptic / neuronal signalling" = "Synaptic /\nneuronal",
  "ribosomal / translation" = "Ribosomal /\ntranslation",
  "proteostasis / protein processing" = "Proteostasis",
  "WNT / adhesion / extracellular matrix" = "WNT / adhesion /\nECM",
  "chromatin / epigenetic regulation" = "Chromatin /\nepigenetic",
  "metabolism" = "Metabolism"
)

omic_labels <- c(
  methylation = "Methylation",
  expression = "Expression",
  methylation_expression_convergence = "Cross-omic\nconvergence"
)

clean_model_label <- function(x) {
  x <- gsub("_R$", "", x)
  x <- gsub("brain_expression_grouped_primary_public", "Brain expression primary", x)
  x <- gsub("brain_expression_cortex_only_sensitivity", "Brain cortex sensitivity", x)
  x <- gsub("brain_expression_cerebellum_only_sensitivity", "Brain cerebellum sensitivity", x)
  x <- gsub("brain_expression_microarray_only_sensitivity", "Brain microarray sensitivity", x)
  x <- gsub("blood_expression_primary_public", "Blood expression primary", x)
  x <- gsub("blood_methylation_primary", "Blood methylation primary", x)
  x <- gsub("_", " ", x)
  str_squish(x)
}

cap_score <- function(p) pmin(-log10(pmax(p, 1e-300)), 60)

# -------------------------------------------------------------------------
# Figure 1: tissue-by-omic pathway evidence map
# -------------------------------------------------------------------------

focus <- theme[
  theme %in% theme_order &
    evidence_tier == "primary_or_threshold_robust" &
    tissue %in% c("blood", "brain", "LCL", "placenta")
]

evidence_map <- focus[, .(
  significant_terms = sum(significant_terms, na.rm = TRUE),
  best_adjusted_p = min(best_adjusted_p, na.rm = TRUE),
  methods = paste(sort(unique(unlist(strsplit(paste(methods, collapse = ";"), ";")))), collapse = ";"),
  representative_terms = paste(head(unique(unlist(strsplit(paste(representative_terms, collapse = " | "), "\\s+\\|\\s+"))), 6), collapse = " | ")
), by = .(tissue, omic, theme)]

evidence_map[, evidence_score := cap_score(best_adjusted_p)]
evidence_map[, omic_label := factor(omic_labels[omic], levels = omic_labels[c("methylation", "expression", "methylation_expression_convergence")])]
evidence_map[, tissue := factor(tissue, levels = c("blood", "brain", "placenta", "LCL"), labels = c("Blood", "Post-mortem brain", "Placenta", "LCL"))]
evidence_map[, theme := factor(theme, levels = rev(theme_order), labels = theme_labels[rev(theme_order)])]

all_combos <- CJ(
  tissue = factor(c("Blood", "Post-mortem brain", "Placenta", "LCL"), levels = c("Blood", "Post-mortem brain", "Placenta", "LCL")),
  omic_label = factor(omic_labels[c("methylation", "expression", "methylation_expression_convergence")],
                      levels = omic_labels[c("methylation", "expression", "methylation_expression_convergence")]),
  theme = factor(theme_labels[rev(theme_order)], levels = theme_labels[rev(theme_order)])
)
evidence_map_full <- merge(all_combos, evidence_map, by = c("tissue", "omic_label", "theme"), all.x = TRUE)
evidence_map_full[is.na(significant_terms), `:=`(significant_terms = 0, evidence_score = 0, best_adjusted_p = NA_real_)]
evidence_map_full[, evidence_band := cut(evidence_score, breaks = c(-Inf, 0, 2, 5, 10, Inf),
                                         labels = c("No FDR terms", "Nominal FDR", "Moderate", "Strong", "Very strong"))]
fwrite(evidence_map_full, file.path(out_dir, "02_figure_source_data/Figure1_pathway_evidence_map_source.csv"))

p1 <- ggplot(evidence_map_full, aes(x = omic_label, y = theme)) +
  geom_tile(aes(fill = evidence_score), colour = "white", linewidth = 0.45, width = 0.92, height = 0.92) +
  geom_point(aes(size = pmax(significant_terms, 0)), shape = 21, fill = "white", colour = "#1f2933", alpha = 0.82, stroke = 0.25) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, max(evidence_map_full$evidence_score, na.rm = TRUE)),
                       name = expression(-log[10]("adjusted p")), breaks = c(0, 5, 15, 30, 60)) +
  scale_size_area(max_size = 8, name = "FDR terms", breaks = c(1, 10, 50, 150, 300)) +
  labs(title = "Pathway-level evidence separates peripheral and post-mortem brain signals",
       subtitle = "Tiles show strongest adjusted enrichment per tissue, omic layer and biological theme; point size shows the number of significant enriched terms.",
       x = NULL, y = NULL,
       caption = "Only primary/threshold-supported enrichment inputs are shown; exploratory DL-screening results are excluded.") +
  theme_minimal(base_family = "sans", base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "#4b5563"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "03_figures/Figure1_pathway_evidence_map.png"), p1, width = 13.5, height = 6.8, dpi = 320)
ggsave(file.path(out_dir, "03_figures/Figure1_pathway_evidence_map.pdf"), p1, width = 13.5, height = 6.8, device = cairo_pdf)

# -------------------------------------------------------------------------
# Figure 2: focused brain biology dot plot
# -------------------------------------------------------------------------

term_keep <- all[
  tissue == "brain" &
    evidence_tier == "primary_or_threshold_robust" &
    significant_05 == TRUE &
    theme %in% c("mitochondrial respiration / oxidative phosphorylation",
                 "immune / antigen processing / MHC",
                 "synaptic / neuronal signalling",
                 "ribosomal / translation") &
    omic %in% c("expression", "methylation_expression_convergence")
]

term_keep <- term_keep[
  !grepl("molecular_function|cellular_component|biological_process|binding$|metabolic process$|cellular anatomical structure",
         term_name, ignore.case = TRUE)
]
term_keep[, term_short := term_name]
term_keep[, term_short := gsub("REACTOME_|WP_|HALLMARK_", "", term_short)]
term_keep[, term_short := str_replace_all(term_short, "_", " ")]
term_keep[, term_short := str_to_sentence(term_short)]
term_keep[, term_short := str_trunc(term_short, 62)]
term_keep[, score := cap_score(adjusted_p)]
term_keep[, layer := fifelse(omic == "expression", "Brain expression", "Brain methylation-expression convergence")]
term_keep[, evidence_source := paste(layer, gene_set_type, sep = " | ")]

brain_dot <- term_keep[
  order(theme, adjusted_p),
  head(.SD, 8),
  by = .(theme, layer)
]
brain_dot[, theme := factor(theme, levels = rev(c("mitochondrial respiration / oxidative phosphorylation",
                                                  "immune / antigen processing / MHC",
                                                  "synaptic / neuronal signalling",
                                                  "ribosomal / translation")),
                            labels = theme_labels[rev(c("mitochondrial respiration / oxidative phosphorylation",
                                                        "immune / antigen processing / MHC",
                                                        "synaptic / neuronal signalling",
                                                        "ribosomal / translation"))])]
brain_dot[, layer := factor(layer, levels = c("Brain expression", "Brain methylation-expression convergence"))]
brain_dot[, term_short := fct_reorder(term_short, score)]
fwrite(brain_dot, file.path(out_dir, "02_figure_source_data/Figure2_brain_focused_dotplot_source.csv"))

p2 <- ggplot(brain_dot, aes(x = score, y = term_short)) +
  geom_segment(aes(x = 0, xend = score, yend = term_short, colour = theme), linewidth = 0.55, alpha = 0.55) +
  geom_point(aes(size = input_gene_count, fill = theme), shape = 21, colour = "white", stroke = 0.3, alpha = 0.95) +
  facet_grid(theme ~ layer, scales = "free_y", space = "free_y") +
  scale_fill_viridis_d(option = "turbo", end = 0.86, name = "Theme") +
  scale_colour_viridis_d(option = "turbo", end = 0.86, guide = "none") +
  scale_size_continuous(range = c(2, 7), name = "Input genes") +
  labs(title = "Post-mortem brain enrichment centres on mitochondrial, immune and synaptic systems",
       subtitle = "Focused top terms from robust brain expression and methylation-expression convergence inputs.",
       x = expression(-log[10]("adjusted p")), y = NULL,
       caption = "Disease-labelled pathway terms were retained only where they represent gene-set containers; interpretation should focus on the underlying biological processes.") +
  theme_minimal(base_family = "sans", base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
    strip.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "#4b5563"),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "03_figures/Figure2_brain_focused_pathway_dotplot.png"), p2, width = 12.8, height = 11, dpi = 320)
ggsave(file.path(out_dir, "03_figures/Figure2_brain_focused_pathway_dotplot.pdf"), p2, width = 12.8, height = 11, device = cairo_pdf)

# -------------------------------------------------------------------------
# Figure 3: convergence/alluvial-style theme flow
# -------------------------------------------------------------------------

flow_dat <- all[
  evidence_tier == "primary_or_threshold_robust" &
    significant_05 == TRUE &
    tissue %in% c("blood", "brain") &
    theme %in% theme_order &
    omic %in% c("methylation", "expression", "methylation_expression_convergence")
]
flow <- flow_dat[, .(
  terms = .N,
  best_adjusted_p = min(adjusted_p, na.rm = TRUE)
), by = .(tissue, omic, theme)]
flow[, weight := pmin(cap_score(best_adjusted_p), 40) * pmax(1, log10(terms + 1))]
flow[, tissue := str_to_title(tissue)]
flow[, omic := factor(omic_labels[omic], levels = omic_labels[c("methylation", "expression", "methylation_expression_convergence")])]
flow[, theme := factor(theme, levels = theme_order, labels = theme_labels[theme_order])]
fwrite(flow, file.path(out_dir, "02_figure_source_data/Figure3_cross_omic_theme_flow_source.csv"))

p3 <- ggplot(flow, aes(axis1 = tissue, axis2 = omic, axis3 = theme, y = weight)) +
  geom_alluvium(aes(fill = theme), alpha = 0.78, width = 1/12, knot.pos = 0.45) +
  geom_stratum(width = 1/8, fill = "#f8fafc", colour = "#334155", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.1, colour = "#111827") +
  scale_x_discrete(limits = c("Tissue", "Omic layer", "Pathway theme"), expand = c(.08, .08)) +
  scale_fill_viridis_d(option = "turbo", end = 0.86, name = "Theme") +
  labs(title = "Pathway evidence flows from tissue-specific omic signals into shared biological themes",
       subtitle = "Flow width combines enrichment strength and number of significant terms for blood and brain robust inputs.",
       x = NULL, y = "Evidence-weighted pathway support",
       caption = "This is a visual summary of enriched pathway themes, not a causal network or participant-level regulatory model.") +
  theme_minimal(base_family = "sans", base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "#4b5563"),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "03_figures/Figure3_cross_omic_theme_flow.png"), p3, width = 13, height = 7.3, dpi = 320)
ggsave(file.path(out_dir, "03_figures/Figure3_cross_omic_theme_flow.pdf"), p3, width = 13, height = 7.3, device = cairo_pdf)

# -------------------------------------------------------------------------
# Figure 4: theme network
# -------------------------------------------------------------------------

node_themes <- theme_order
node_tissue <- c("blood", "brain")
node_omics <- c("methylation", "expression", "methylation_expression_convergence")

theme_nodes <- data.table(name = paste0("theme:", node_themes), label = theme_labels[node_themes], type = "theme")
tissue_nodes <- data.table(name = paste0("tissue:", node_tissue), label = str_to_title(node_tissue), type = "tissue")
omic_nodes <- data.table(name = paste0("omic:", node_omics), label = omic_labels[node_omics], type = "omic")
nodes <- rbindlist(list(tissue_nodes, omic_nodes, theme_nodes), fill = TRUE)

edge_base <- flow_dat[, .(
  terms = .N,
  best_adjusted_p = min(adjusted_p, na.rm = TRUE)
), by = .(tissue, omic, theme)]
edge_base[, weight := pmin(cap_score(best_adjusted_p), 50) * log10(terms + 1)]
tissue_edges <- edge_base[, .(from = paste0("tissue:", tissue), to = paste0("omic:", omic), weight = sum(weight), edge_type = "tissue_to_omic"), by = .(tissue, omic)]
omic_theme_edges <- edge_base[, .(from = paste0("omic:", omic), to = paste0("theme:", theme), weight = sum(weight), edge_type = "omic_to_theme"), by = .(omic, theme)]
edges <- rbindlist(list(tissue_edges[, .(from, to, weight, edge_type)], omic_theme_edges[, .(from, to, weight, edge_type)]), fill = TRUE)
fwrite(nodes, file.path(out_dir, "02_figure_source_data/Figure4_theme_network_nodes.csv"))
fwrite(edges, file.path(out_dir, "02_figure_source_data/Figure4_theme_network_edges.csv"))

graph <- tbl_graph(nodes = as.data.frame(nodes), edges = as.data.frame(edges), directed = FALSE)

p4 <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), colour = "#64748b", show.legend = FALSE) +
  geom_node_point(aes(size = type, fill = type), shape = 21, colour = "white", stroke = 0.5) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3.4, lineheight = 0.9, colour = "#111827", max.overlaps = Inf) +
  scale_edge_width(range = c(0.3, 4.2)) +
  scale_size_manual(values = c(theme = 9, tissue = 11, omic = 10), guide = "none") +
  scale_fill_manual(values = c(theme = "#22c55e", tissue = "#0ea5e9", omic = "#f97316"), guide = "none") +
  labs(title = "Compact network view of tissue, omic layer and pathway-theme support",
       subtitle = "Edges are weighted by enrichment strength and term counts across robust blood and brain inputs.",
       caption = "Visual grouping is descriptive; it should not be interpreted as a mechanistic network model.") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "#4b5563"),
    plot.caption = element_text(colour = "#4b5563")
  )

ggsave(file.path(out_dir, "03_figures/Figure4_tissue_omic_pathway_theme_network.png"), p4, width = 11.5, height = 8.5, dpi = 320)
ggsave(file.path(out_dir, "03_figures/Figure4_tissue_omic_pathway_theme_network.pdf"), p4, width = 11.5, height = 8.5, device = cairo_pdf)

# -------------------------------------------------------------------------
# Figure 5: compact top themes bar/ridge-style panel
# -------------------------------------------------------------------------

bar_dat <- theme[
  tissue %in% c("blood", "brain") &
    evidence_tier == "primary_or_threshold_robust" &
    theme %in% theme_order
]
bar_dat[, score := cap_score(best_adjusted_p)]
bar_summ <- bar_dat[, .(
  score = max(score, na.rm = TRUE),
  significant_terms = sum(significant_terms, na.rm = TRUE),
  best_adjusted_p = min(best_adjusted_p, na.rm = TRUE)
), by = .(tissue, omic, theme)]
bar_summ[, tissue := factor(str_to_title(tissue), levels = c("Blood", "Brain"))]
bar_summ[, omic := factor(omic_labels[omic], levels = omic_labels[c("methylation", "expression", "methylation_expression_convergence")])]
bar_summ[, theme := factor(theme, levels = rev(theme_order), labels = theme_labels[rev(theme_order)])]
fwrite(bar_summ, file.path(out_dir, "02_figure_source_data/Figure5_top_theme_strength_panel_source.csv"))

p5 <- ggplot(bar_summ, aes(x = score, y = theme, fill = omic)) +
  geom_col(position = position_dodge2(width = 0.78, preserve = "single"), width = 0.68, colour = "white", linewidth = 0.2) +
  facet_wrap(~ tissue, ncol = 1) +
  scale_fill_manual(values = c("Methylation" = "#2563eb", "Expression" = "#dc2626", "Cross-omic\nconvergence" = "#7c3aed"),
                    name = "Omic layer") +
  labs(title = "Blood and brain show different pathway priorities",
       subtitle = "Blood is dominated by cell-state/proliferation and immune terms; brain is dominated by mitochondrial, immune and synaptic systems.",
       x = expression("Best theme enrichment strength ("*-log[10]("adjusted p")*")"), y = NULL) +
  theme_minimal(base_family = "sans", base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "#4b5563"),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "03_figures/Figure5_blood_brain_theme_strength_panel.png"), p5, width = 10.5, height = 8.3, dpi = 320)
ggsave(file.path(out_dir, "03_figures/Figure5_blood_brain_theme_strength_panel.pdf"), p5, width = 10.5, height = 8.3, device = cairo_pdf)

source_files <- list.files(file.path(out_dir, "02_figure_source_data"), full.names = TRUE)
fig_files <- list.files(file.path(out_dir, "03_figures"), full.names = TRUE)

openxlsx::write.xlsx(
  list(
    figure1_evidence_map = as.data.frame(evidence_map_full),
    figure2_brain_dotplot = as.data.frame(brain_dot),
    figure3_theme_flow = as.data.frame(flow),
    figure4_nodes = as.data.frame(nodes),
    figure4_edges = as.data.frame(edges),
    figure5_theme_strength = as.data.frame(bar_summ)
  ),
  file.path(out_dir, "02_figure_source_data/enhanced_pathway_figure_source_data.xlsx"),
  overwrite = TRUE
)

script_index <- data.table(
  script = "build_enhanced_pathway_figures.R",
  role = "production figure-generation script",
  packages = paste(c("data.table", "ggplot2", "dplyr", "tidyr", "forcats", "stringr", "scales",
                     "viridis", "ggrepel", "ggalluvial", "igraph", "ggraph", "tidygraph",
                     "patchwork", "openxlsx"), collapse = "; ")
)
fwrite(script_index, file.path(out_dir, "05_quality_control/script_role_manifest.csv"))

pkg_versions <- data.table(
  package = c("data.table", "ggplot2", "dplyr", "tidyr", "forcats", "stringr", "scales",
              "viridis", "ggrepel", "ggalluvial", "igraph", "ggraph", "tidygraph",
              "patchwork", "openxlsx"),
  installed = vapply(c("data.table", "ggplot2", "dplyr", "tidyr", "forcats", "stringr", "scales",
                       "viridis", "ggrepel", "ggalluvial", "igraph", "ggraph", "tidygraph",
                       "patchwork", "openxlsx"), requireNamespace, logical(1), quietly = TRUE),
  version = vapply(c("data.table", "ggplot2", "dplyr", "tidyr", "forcats", "stringr", "scales",
                     "viridis", "ggrepel", "ggalluvial", "igraph", "ggraph", "tidygraph",
                     "patchwork", "openxlsx"),
                   function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_,
                   character(1))
)
fwrite(pkg_versions, file.path(out_dir, "05_quality_control/R_package_versions.csv"))

file_index <- data.table(
  file_path = normalizePath(list.files(out_dir, recursive = TRUE, full.names = TRUE), winslash = "/", mustWork = FALSE)
)
file_index[, bytes := file.info(file_path)$size]
fwrite(file_index, file.path(out_dir, "05_quality_control/enhanced_pathway_figure_file_index.csv"))

recommendations <- c(
  "# Enhanced Pathway Figure Recommendations",
  "",
  "## Best Main-Text Figure",
  "",
  "`Figure1_pathway_evidence_map` is the strongest main-text figure. It shows the whole story compactly: blood is more cell-state/proliferation/immune-weighted, while post-mortem brain has stronger mitochondrial, immune/MHC and synaptic pathway support.",
  "",
  "## Best Biological Detail Figure",
  "",
  "`Figure2_brain_focused_pathway_dotplot` is the best figure for showing the actual brain pathway terms. It should be used either as a main-text figure panel or as a high-value supplement.",
  "",
  "## Best Conceptual Figure",
  "",
  "`Figure3_cross_omic_theme_flow` communicates the convergence idea most clearly, but it should be described as an evidence-flow visual summary rather than a formal causal model.",
  "",
  "## Optional Supplementary Figure",
  "",
  "`Figure4_tissue_omic_pathway_theme_network` is visually interesting but more abstract. It is useful for talks or graphical supplement material, not essential for the main paper.",
  "",
  "## Compact Alternative",
  "",
  "`Figure5_blood_brain_theme_strength_panel` is the simplest figure if journal space is tight. It gives the blood/brain contrast cleanly without requiring readers to parse many terms.",
  "",
  "## Recommended Figure Set",
  "",
  "- Main text: Figure 1 as the primary systems-level summary.",
  "- Main text or supplement: Figure 2 for detailed brain pathway biology.",
  "- Supplement: Figure 3 or Figure 5 depending on journal space.",
  "",
  "## Interpretation Caveats",
  "",
  "- These figures summarise enrichment evidence, not causal pathway activation.",
  "- DL-screening-only results are excluded from the main evidence-map figures.",
  "- Disease-labelled gene sets should be interpreted by their component biology, not as evidence for those diseases.",
  "- Cross-omic convergence is summary-level and tissue-matched, not participant-level regulatory coupling."
)
writeLines(recommendations, file.path(out_dir, "04_reports/enhanced_pathway_figure_recommendations.md"))

qc <- c(
  "# Enhanced Pathway Figure QC",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Source Tables",
  "",
  paste0("- all_enrichment_results_normalised rows: ", nrow(all)),
  paste0("- pathway_theme_summary rows: ", nrow(theme)),
  "",
  "## Figures Created",
  "",
  paste0("- ", basename(fig_files), collapse = "\n"),
  "",
  "## QC Checks",
  "",
  "- Figures use final pathway enrichment outputs only.",
  "- Figure source-data CSV files were written for every figure.",
  "- Main evidence-map figures restrict headline interpretation to primary/threshold-robust enrichment inputs.",
  "- Exploratory DL-screening results are not used for headline visual evidence.",
  "- Standard R packages were used for plotting: ggplot2, ggalluvial, igraph/ggraph/tidygraph, patchwork and supporting tidyverse packages.",
  "- Outputs were written only to the configured figure output directory."
)
writeLines(qc, file.path(out_dir, "05_quality_control/enhanced_pathway_figure_QC.md"))

readme <- c(
  "# Enhanced Pathway Figures R Package",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript scripts/06_figures/build_enhanced_pathway_figures.R",
  "```",
  "",
  "The script reads final pathway enrichment tables from `results/pathway_enrichment/08_reports/` and writes figure source data, PNG/PDF figures, recommendations and QC outputs."
)
writeLines(readme, file.path(out_dir, "README_reproducibility.md"))

message("Enhanced pathway figures written to: ", out_dir)

