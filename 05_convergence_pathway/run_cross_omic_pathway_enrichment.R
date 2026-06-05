#!/usr/bin/env Rscript

# Cross-omic pathway enrichment workflow for the ASD epigenetics synthesis.
#
# Purpose:
#   Build reproducible pathway/enrichment inputs from the final methylation and
#   expression meta-analysis outputs, then run several reviewer-defensible
#   enrichment approaches:
#     1. clusterProfiler GO over-representation analysis
#     2. ReactomePA over-representation analysis
#     3. MSigDB Hallmark/Reactome/WikiPathways over-representation analysis
#     4. fgsea preranked enrichment using all tested genes
#     5. g:Profiler external over-representation check where internet access works
#
# The primary background for each ORA is the model-specific tested gene universe.
# Cross-omic convergence backgrounds use the tissue-matched intersection of the
# methylation and expression tested backgrounds.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(fgsea)
  library(gprofiler2)
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
dirs <- c(
  "00_manifest", "01_enrichment_inputs", "02_ORA_clusterProfiler_GO",
  "03_ReactomePA_ORA", "04_MSigDB_ORA", "05_GSEA_fgsea", "06_gprofiler",
  "07_figures", "08_reports", "scripts"
)
invisible(lapply(file.path(out_dir, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")
}

num <- function(x) suppressWarnings(as.numeric(x))

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

first_present <- function(dt, candidates) {
  hit <- intersect(candidates, names(dt))
  if (!length(hit)) return(NULL)
  hit[[1]]
}

standardise_results <- function(path, source_label, omic, tissue, model_subset = NULL) {
  dt <- fread(path, na.strings = c("", "NA", "NaN"))
  if (!nrow(dt)) return(data.table())

  gene_col <- first_present(dt, c("gene", "Gene", "symbol", "SYMBOL"))
  model_col <- first_present(dt, c("model_name", "model", "Model"))
  role_col <- first_present(dt, c("role", "model_role"))
  effect_col <- first_present(dt, c("pooled_Hedges_g", "pooled_g", "effect_size", "g"))
  p_col <- first_present(dt, c("p_value", "p", "P.Value"))
  fdr_col <- first_present(dt, c("FDR", "fdr", "adj.P.Val", "padj"))
  k_col <- first_present(dt, c("k"))
  i2_col <- first_present(dt, c("I2", "I_squared"))
  dl_col <- first_present(dt, c("DL_CI_excludes_zero", "DL_nonzero", "DL_interval_excludes_zero"))
  mkh_col <- first_present(dt, c("mKH_CI_excludes_zero", "mKH_interval_excludes_zero", "mKH_interval_supported"))
  direction_col <- first_present(dt, c("direction"))

  if (is.null(gene_col) || is.null(model_col) || is.null(effect_col) || is.null(p_col)) {
    stop("Required columns missing in ", path)
  }

  out <- data.table(
    source_file = path,
    source_label = source_label,
    omic = if (!is.na(omic)) omic else as.character(dt[[first_present(dt, c("omic_layer", "omic"))]]),
    tissue = if (!is.na(tissue)) tissue else as.character(dt[[first_present(dt, c("tissue_branch", "tissue"))]]),
    gene = as.character(dt[[gene_col]]),
    model = as.character(dt[[model_col]]),
    role = if (!is.null(role_col)) as.character(dt[[role_col]]) else NA_character_,
    effect = num(dt[[effect_col]]),
    p_value = num(dt[[p_col]]),
    FDR = if (!is.null(fdr_col)) num(dt[[fdr_col]]) else NA_real_,
    k = if (!is.null(k_col)) num(dt[[k_col]]) else NA_real_,
    I2 = if (!is.null(i2_col)) num(dt[[i2_col]]) else NA_real_,
    DL_nonzero = if (!is.null(dl_col)) as_bool(dt[[dl_col]]) else NA,
    mKH_supported = if (!is.null(mkh_col)) as_bool(dt[[mkh_col]]) else NA,
    direction = if (!is.null(direction_col)) as.character(dt[[direction_col]]) else NA_character_
  )

  out <- out[!is.na(gene) & nzchar(gene) & is.finite(p_value) & is.finite(effect)]
  out[, gene := clean_symbol(gene), by = seq_len(nrow(out))]
  out[tissue == "post-mortem brain", tissue := "brain"]
  out <- unique(out, by = c("source_label", "omic", "tissue", "model", "gene"))

  if (!is.null(model_subset)) out <- out[model %in% model_subset]
  out[]
}

sources <- data.table(
  source_label = c("methylation_metafor_models", "expression_metafor_models"),
  omic = c(NA_character_, NA_character_),
  tissue = c(NA_character_, NA_character_),
  path = file.path(package_root, c(
    "results/meta_analysis/methylation/methylation_gene_level_meta_results.csv",
    "results/meta_analysis/expression/expression_gene_level_meta_results.csv"
  ))
)
sources[, exists := file.exists(path)]
fwrite(sources, file.path(out_dir, "00_manifest", "enrichment_source_file_manifest.csv"))
if (any(!sources$exists)) {
  stop("Missing required source files: ", paste(sources[exists == FALSE, path], collapse = "; "))
}

all_results <- rbindlist(lapply(seq_len(nrow(sources)), function(i) {
  standardise_results(
    sources$path[[i]],
    sources$source_label[[i]],
    sources$omic[[i]],
    sources$tissue[[i]]
  )
}), fill = TRUE)

primary_models <- c(
  "blood_array_peripheral_primary",
  "brain_grouped_primary_with_WGBS",
  "placenta_primary_GSE178203_descriptive",
  "placenta_two_dataset_sensitivity",
  "lcl_cross_assay_exploratory",
  "blood_expression_peripheral_primary",
  "brain_expression_grouped_primary_public_R",
  "lcl_expression_core_public_primary_R"
)

key_sensitivity_models <- c(
  "blood_450k_only_sensitivity",
  "blood_array_plus_cord_WGBS_sensitivity",
  "brain_450k_only_sensitivity",
  "brain_array_HM27_primary",
  "cortex_only_sensitivity",
  "cerebellum_only_sensitivity",
  "prefrontal_BA9_cortex_sensitivity",
  "WGBS_excluded_sensitivity",
  "WGBS_only_sensitivity",
  "blood_expression_plus_cord_blood_sensitivity",
  "brain_expression_microarray_only_sensitivity_R",
  "brain_expression_RNAseq_only_sensitivity_R",
  "brain_expression_cortex_only_sensitivity_R",
  "brain_expression_cerebellum_only_sensitivity_R",
  "brain_expression_prefrontal_cortex_sensitivity_R",
  "brain_expression_BA19_occipital_sensitivity_R",
  "lcl_expression_expanded_public_sensitivity_R",
  "lcl_expression_no_pooled_syndromic_sensitivity_R"
)

all_results <- all_results[model %in% c(primary_models, key_sensitivity_models)]
all_results[, analysis_id := paste(tissue, omic, model, sep = "__")]
all_results[, signed_rank := sign(effect) * -log10(pmax(p_value, .Machine$double.xmin))]

fwrite(all_results, file.path(out_dir, "01_enrichment_inputs", "standardised_model_gene_results.csv"))

# Keep enrichment broad enough for interpretation, but avoid turning
# every sensitivity/DL-screening result into a massive exploratory table. The
# full model-level results remain archived above; enrichment is focused on
# primary models, the most interpretable sensitivity models, and tissue-matched
# convergence sets.
focused_enrichment_models <- c(
  primary_models,
  "blood_450k_only_sensitivity",
  "blood_array_plus_cord_WGBS_sensitivity",
  "brain_450k_only_sensitivity",
  "cortex_only_sensitivity",
  "WGBS_excluded_sensitivity",
  "blood_expression_plus_cord_blood_sensitivity",
  "brain_expression_microarray_only_sensitivity_R",
  "brain_expression_RNAseq_only_sensitivity_R",
  "brain_expression_cortex_only_sensitivity_R",
  "brain_expression_cerebellum_only_sensitivity_R",
  "lcl_expression_expanded_public_sensitivity_R"
)

map_symbols <- function(symbols) {
  symbols <- clean_symbol(symbols)
  if (!length(symbols)) {
    return(data.table(SYMBOL = character(), ENTREZID = character()))
  }
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
  )
  mapped <- as.data.table(mapped)
  mapped <- mapped[!is.na(ENTREZID) & nzchar(ENTREZID)]
  unique(mapped, by = c("SYMBOL", "ENTREZID"))
}

symbol_map <- map_symbols(unique(all_results$gene))
fwrite(symbol_map, file.path(out_dir, "01_enrichment_inputs", "gene_symbol_to_entrez_map.csv"))

backgrounds <- all_results[, .(
  tested_genes_symbol = uniqueN(gene),
  tested_genes_mapped_entrez = uniqueN(symbol_map[SYMBOL %in% gene, ENTREZID]),
  finite_p_values = sum(is.finite(p_value)),
  min_p = min(p_value, na.rm = TRUE),
  median_p = median(p_value, na.rm = TRUE),
  FDR_significant = sum(FDR < 0.05, na.rm = TRUE),
  mKH_supported = sum(mKH_supported %in% TRUE, na.rm = TRUE),
  DL_nonzero = sum(DL_nonzero %in% TRUE, na.rm = TRUE)
), by = .(analysis_id, tissue, omic, model)]
fwrite(backgrounds, file.path(out_dir, "01_enrichment_inputs", "model_tested_backgrounds.csv"))

gene_sets <- rbindlist(list(
  all_results[FDR < 0.05, .(gene = unique(gene)), by = .(analysis_id, tissue, omic, model)][
    , gene_set_type := "FDR_significant"],
  all_results[mKH_supported %in% TRUE, .(gene = unique(gene)), by = .(analysis_id, tissue, omic, model)][
    , gene_set_type := "mKH_interval_supported"],
  all_results[DL_nonzero %in% TRUE, .(gene = unique(gene)), by = .(analysis_id, tissue, omic, model)][
    , gene_set_type := "DL_interval_screening"]
), fill = TRUE)
gene_sets[, gene_set_id := paste(analysis_id, gene_set_type, sep = "__")]

gene_set_summary <- gene_sets[, .(n_genes = uniqueN(gene)), by = .(gene_set_id, analysis_id, tissue, omic, model, gene_set_type)]
fwrite(gene_sets, file.path(out_dir, "01_enrichment_inputs", "threshold_gene_sets_long.csv"))
fwrite(gene_set_summary, file.path(out_dir, "01_enrichment_inputs", "threshold_gene_sets_summary.csv"))

make_convergence_sets <- function(tissue_name, methyl_model, expression_model) {
  meth <- all_results[tissue == tissue_name & omic == "methylation" & model == methyl_model]
  expr <- all_results[tissue == tissue_name & omic == "expression" & model == expression_model]
  bg <- intersect(meth$gene, expr$gene)
  if (!length(bg)) return(list(sets = data.table(), bg = data.table()))

  status_sets <- list(
    stringent_FDR_and_mKH_both_layers = intersect(
      meth[FDR < 0.05 & mKH_supported %in% TRUE, gene],
      expr[FDR < 0.05 & mKH_supported %in% TRUE, gene]
    ),
    FDR_both_layers = intersect(meth[FDR < 0.05, gene], expr[FDR < 0.05, gene]),
    mKH_both_layers = intersect(meth[mKH_supported %in% TRUE, gene], expr[mKH_supported %in% TRUE, gene]),
    DL_both_layers = intersect(meth[DL_nonzero %in% TRUE, gene], expr[DL_nonzero %in% TRUE, gene])
  )
  sets <- rbindlist(lapply(names(status_sets), function(nm) {
    data.table(
      tissue = tissue_name,
      omic = "methylation_expression_convergence",
      model = paste(methyl_model, expression_model, sep = " + "),
      analysis_id = paste(tissue_name, "convergence", nm, sep = "__"),
      gene_set_type = nm,
      gene = clean_symbol(status_sets[[nm]])
    )
  }), fill = TRUE)
  sets <- sets[!is.na(gene) & nzchar(gene)]
  sets[, gene_set_id := paste(analysis_id, gene_set_type, sep = "__")]
  bgdt <- data.table(
    tissue = tissue_name,
    analysis_id = paste(tissue_name, "convergence_background", sep = "__"),
    model = paste(methyl_model, expression_model, sep = " + "),
    gene = bg
  )
  list(sets = sets, bg = bgdt)
}

conv_blood <- make_convergence_sets("blood", "blood_array_peripheral_primary", "blood_expression_peripheral_primary")
conv_brain <- make_convergence_sets("brain", "brain_grouped_primary_with_WGBS", "brain_expression_grouped_primary_public_R")
convergence_sets <- rbindlist(list(conv_blood$sets, conv_brain$sets), fill = TRUE)
convergence_backgrounds <- rbindlist(list(conv_blood$bg, conv_brain$bg), fill = TRUE)
fwrite(convergence_sets, file.path(out_dir, "01_enrichment_inputs", "convergence_gene_sets_long.csv"))
fwrite(convergence_backgrounds, file.path(out_dir, "01_enrichment_inputs", "convergence_tested_backgrounds_long.csv"))

all_ora_sets <- rbindlist(list(gene_sets, convergence_sets), fill = TRUE)
set_background <- function(row_analysis_id, row_tissue, row_omic) {
  if (row_omic == "methylation_expression_convergence") {
    return(convergence_backgrounds[tissue == row_tissue, unique(gene)])
  }
  all_results[analysis_id == row_analysis_id, unique(gene)]
}

gene_set_catalog <- unique(all_ora_sets[, .(gene_set_id, analysis_id, tissue, omic, model, gene_set_type)])
gene_set_catalog[, n_genes := all_ora_sets[.SD, on = "gene_set_id", uniqueN(gene), by = .EACHI]$V1]
gene_set_catalog[, n_background := mapply(function(a, t, o) length(set_background(a, t, o)), analysis_id, tissue, omic)]
gene_set_catalog <- gene_set_catalog[
  (model %in% focused_enrichment_models & gene_set_type %in% c("FDR_significant", "mKH_interval_supported")) |
    (model %in% primary_models & gene_set_type == "DL_interval_screening") |
    (omic == "methylation_expression_convergence" & gene_set_type %in% c("FDR_both_layers", "mKH_both_layers", "DL_both_layers", "stringent_FDR_and_mKH_both_layers"))
]
gene_set_catalog <- gene_set_catalog[n_genes >= 5 & n_background >= 50]
fwrite(gene_set_catalog, file.path(out_dir, "01_enrichment_inputs", "enrichment_gene_set_catalog.csv"))

entrez_for <- function(symbols) {
  unique(symbol_map[SYMBOL %in% clean_symbol(symbols), ENTREZID])
}

run_go_ora <- function() {
  out <- list()
  for (i in seq_len(nrow(gene_set_catalog))) {
    row <- gene_set_catalog[i]
    genes <- all_ora_sets[gene_set_id == row$gene_set_id, unique(gene)]
    bg <- set_background(row$analysis_id, row$tissue, row$omic)
    gene_entrez <- entrez_for(genes)
    bg_entrez <- entrez_for(bg)
    if (length(gene_entrez) < 5 || length(bg_entrez) < 50) next
    for (ont in c("BP", "MF", "CC")) {
      ego <- tryCatch(
        enrichGO(
          gene = gene_entrez,
          universe = bg_entrez,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.2,
          qvalueCutoff = 0.25,
          minGSSize = 10,
          maxGSSize = 500,
          readable = TRUE
        ),
        error = function(e) NULL
      )
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        dt <- as.data.table(as.data.frame(ego))
        dt[, `:=`(
          enrichment_method = "clusterProfiler_enrichGO",
          database = paste0("GO:", ont),
          gene_set_id = row$gene_set_id,
          analysis_id = row$analysis_id,
          tissue = row$tissue,
          omic = row$omic,
          model = row$model,
          gene_set_type = row$gene_set_type,
          input_gene_count = length(gene_entrez),
          background_gene_count = length(bg_entrez)
        )]
        out[[length(out) + 1]] <- dt
      }
    }
  }
  rbindlist(out, fill = TRUE)
}

run_reactome_ora <- function() {
  out <- list()
  for (i in seq_len(nrow(gene_set_catalog))) {
    row <- gene_set_catalog[i]
    genes <- all_ora_sets[gene_set_id == row$gene_set_id, unique(gene)]
    bg <- set_background(row$analysis_id, row$tissue, row$omic)
    gene_entrez <- entrez_for(genes)
    bg_entrez <- entrez_for(bg)
    if (length(gene_entrez) < 5 || length(bg_entrez) < 50) next
    er <- tryCatch(
      enrichPathway(
        gene = gene_entrez,
        universe = bg_entrez,
        organism = "human",
        pvalueCutoff = 0.2,
        qvalueCutoff = 0.25,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(er) && nrow(as.data.frame(er)) > 0) {
      dt <- as.data.table(as.data.frame(er))
      dt[, `:=`(
        enrichment_method = "ReactomePA_enrichPathway",
        database = "ReactomePA",
        gene_set_id = row$gene_set_id,
        analysis_id = row$analysis_id,
        tissue = row$tissue,
        omic = row$omic,
        model = row$model,
        gene_set_type = row$gene_set_type,
        input_gene_count = length(gene_entrez),
        background_gene_count = length(bg_entrez)
      )]
      out[[length(out) + 1]] <- dt
    }
  }
  rbindlist(out, fill = TRUE)
}

load_msigdb_term2gene <- function() {
  collections <- list(
    Hallmark = msigdbr(species = "Homo sapiens", collection = "H"),
    Reactome_MSigDB = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
    WikiPathways_MSigDB = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
  )
  rbindlist(lapply(names(collections), function(nm) {
    x <- as.data.table(collections[[nm]])
    data.table(
      database = nm,
      term = x$gs_name,
      entrez = as.character(x$ncbi_gene),
      symbol = x$gene_symbol
    )
  }), fill = TRUE)[!is.na(entrez) & nzchar(entrez)]
}

msig_t2g <- load_msigdb_term2gene()
fwrite(msig_t2g, file.path(out_dir, "01_enrichment_inputs", "msigdb_term2gene_used.csv"))

run_msigdb_ora <- function() {
  out <- list()
  for (db in unique(msig_t2g$database)) {
    t2g <- unique(msig_t2g[database == db, .(term, entrez)])
    for (i in seq_len(nrow(gene_set_catalog))) {
      row <- gene_set_catalog[i]
      genes <- all_ora_sets[gene_set_id == row$gene_set_id, unique(gene)]
      bg <- set_background(row$analysis_id, row$tissue, row$omic)
      gene_entrez <- intersect(entrez_for(genes), t2g$entrez)
      bg_entrez <- intersect(entrez_for(bg), t2g$entrez)
      if (length(gene_entrez) < 5 || length(bg_entrez) < 50) next
      enr <- tryCatch(
        enricher(
          gene = gene_entrez,
          universe = bg_entrez,
          TERM2GENE = t2g,
          pvalueCutoff = 0.2,
          pAdjustMethod = "BH",
          qvalueCutoff = 0.25,
          minGSSize = 10,
          maxGSSize = 500
        ),
        error = function(e) NULL
      )
      if (!is.null(enr) && nrow(as.data.frame(enr)) > 0) {
        dt <- as.data.table(as.data.frame(enr))
        dt[, `:=`(
          enrichment_method = "clusterProfiler_enricher_MSigDB",
          database = db,
          gene_set_id = row$gene_set_id,
          analysis_id = row$analysis_id,
          tissue = row$tissue,
          omic = row$omic,
          model = row$model,
          gene_set_type = row$gene_set_type,
          input_gene_count = length(gene_entrez),
          background_gene_count = length(bg_entrez)
        )]
        out[[length(out) + 1]] <- dt
      }
    }
  }
  rbindlist(out, fill = TRUE)
}

run_fgsea <- function() {
  pathways <- split(msig_t2g$entrez, paste(msig_t2g$database, msig_t2g$term, sep = "::"))
  out <- list()
  ranked_models <- unique(all_results[model %in% focused_enrichment_models, .(analysis_id, tissue, omic, model)])
  for (i in seq_len(nrow(ranked_models))) {
    row <- ranked_models[i]
    dt <- all_results[analysis_id == row$analysis_id]
    mapped <- merge(dt, symbol_map, by.x = "gene", by.y = "SYMBOL", allow.cartesian = TRUE)
    mapped <- mapped[is.finite(signed_rank) & !is.na(ENTREZID)]
    if (nrow(mapped) < 100) next
    ranked <- mapped[order(abs(signed_rank), decreasing = TRUE), .SD[1], by = ENTREZID]
    ranks <- ranked$signed_rank
    names(ranks) <- ranked$ENTREZID
    ranks <- sort(ranks, decreasing = TRUE)
    fg <- tryCatch(
      fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 1e-20),
      error = function(e) NULL
    )
    if (!is.null(fg) && nrow(fg) > 0) {
      fg <- as.data.table(fg)
      fg[, c("database", "term") := tstrsplit(pathway, "::", fixed = TRUE)]
      fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
      fg[, `:=`(
        enrichment_method = "fgsea_preranked_signed_minus_log10p",
        analysis_id = row$analysis_id,
        tissue = row$tissue,
        omic = row$omic,
        model = row$model,
        ranked_gene_count = length(ranks)
      )]
      out[[length(out) + 1]] <- fg
    }
  }
  rbindlist(out, fill = TRUE)
}

run_gprofiler_checks <- function() {
  candidates <- gene_set_catalog[
    gene_set_type %in% c("FDR_significant", "mKH_interval_supported", "stringent_FDR_and_mKH_both_layers", "mKH_both_layers")
  ][n_genes >= 5]
  candidates <- candidates[
    model %in% c(
      "blood_array_peripheral_primary",
      "blood_expression_peripheral_primary",
      "brain_grouped_primary_with_WGBS",
      "brain_expression_grouped_primary_public_R",
      "lcl_expression_core_public_primary_R"
    ) | omic == "methylation_expression_convergence"
  ]
  out <- list()
  for (i in seq_len(nrow(candidates))) {
    row <- candidates[i]
    genes <- all_ora_sets[gene_set_id == row$gene_set_id, unique(gene)]
    bg <- set_background(row$analysis_id, row$tissue, row$omic)
    if (length(genes) < 5 || length(bg) < 50) next
    gp <- tryCatch(
      gost(
        query = genes,
        organism = "hsapiens",
        correction_method = "fdr",
        user_threshold = 1,
        custom_bg = bg,
        sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "WP", "KEGG"),
        significant = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(gp) && !is.null(gp$result) && nrow(gp$result) > 0) {
      dt <- as.data.table(gp$result)
      dt[, `:=`(
        enrichment_method = "gprofiler2_gost",
        gene_set_id = row$gene_set_id,
        analysis_id = row$analysis_id,
        tissue = row$tissue,
        omic = row$omic,
        model = row$model,
        gene_set_type = row$gene_set_type,
        input_gene_count = length(genes),
        background_gene_count = length(bg)
      )]
      out[[length(out) + 1]] <- dt
    }
  }
  rbindlist(out, fill = TRUE)
}

message("Running clusterProfiler GO ORA...")
go_ora <- run_go_ora()
fwrite(go_ora, file.path(out_dir, "02_ORA_clusterProfiler_GO", "clusterProfiler_GO_ORA_results.csv"))

message("Running ReactomePA ORA...")
reactome_ora <- run_reactome_ora()
fwrite(reactome_ora, file.path(out_dir, "03_ReactomePA_ORA", "ReactomePA_ORA_results.csv"))

message("Running MSigDB ORA...")
msig_ora <- run_msigdb_ora()
fwrite(msig_ora, file.path(out_dir, "04_MSigDB_ORA", "MSigDB_ORA_results.csv"))

message("Running fgsea preranked enrichment...")
fgsea_results <- run_fgsea()
fwrite(fgsea_results, file.path(out_dir, "05_GSEA_fgsea", "fgsea_preranked_results.csv"))

message("Running g:Profiler external checks...")
gprofiler_results <- run_gprofiler_checks()
fwrite(gprofiler_results, file.path(out_dir, "06_gprofiler", "gprofiler_ORA_results.csv"))

if (Sys.getenv("ASD_PATHWAY_RUN_POSTPROCESS_IN_SAME_SCRIPT", "FALSE") != "TRUE") {
  message("Raw enrichment tables completed. Run scripts/postprocess_pathway_enrichment_outputs.R to create summaries, figures and the compact workbook.")
  quit(save = "no", status = 0)
}

normalise_ora <- function(dt) {
  if (!nrow(dt)) return(data.table())
  term_col <- first_present(dt, c("Description", "term_name", "name", "term"))
  id_col <- first_present(dt, c("ID", "term_id", "native"))
  padj_col <- first_present(dt, c("p.adjust", "p_value", "padj", "adjusted_p_value"))
  p_col <- first_present(dt, c("pvalue", "p_value", "p.val"))
  gene_col <- first_present(dt, c("geneID", "intersection", "leadingEdge"))
  data.table(
    enrichment_method = dt$enrichment_method,
    database = if ("database" %in% names(dt)) dt$database else if ("source" %in% names(dt)) dt$source else NA_character_,
    term_id = if (!is.null(id_col)) as.character(dt[[id_col]]) else NA_character_,
    term_name = if (!is.null(term_col)) as.character(dt[[term_col]]) else NA_character_,
    adjusted_p = if (!is.null(padj_col)) num(dt[[padj_col]]) else NA_real_,
    p_value = if (!is.null(p_col)) num(dt[[p_col]]) else NA_real_,
    gene_overlap = if (!is.null(gene_col)) as.character(dt[[gene_col]]) else NA_character_,
    gene_set_id = dt$gene_set_id,
    analysis_id = dt$analysis_id,
    tissue = dt$tissue,
    omic = dt$omic,
    model = dt$model,
    gene_set_type = if ("gene_set_type" %in% names(dt)) dt$gene_set_type else NA_character_,
    input_gene_count = if ("input_gene_count" %in% names(dt)) num(dt$input_gene_count) else NA_real_,
    background_gene_count = if ("background_gene_count" %in% names(dt)) num(dt$background_gene_count) else NA_real_
  )
}

normalise_fgsea <- function(dt) {
  if (!nrow(dt)) return(data.table())
  data.table(
    enrichment_method = dt$enrichment_method,
    database = dt$database,
    term_id = dt$pathway,
    term_name = dt$term,
    adjusted_p = dt$padj,
    p_value = dt$pval,
    NES = dt$NES,
    gene_overlap = dt$leadingEdge,
    gene_set_id = NA_character_,
    analysis_id = dt$analysis_id,
    tissue = dt$tissue,
    omic = dt$omic,
    model = dt$model,
    gene_set_type = "preranked_all_genes",
    input_gene_count = dt$size,
    background_gene_count = dt$ranked_gene_count
  )
}

all_enrichment <- rbindlist(list(
  normalise_ora(go_ora),
  normalise_ora(reactome_ora),
  normalise_ora(msig_ora),
  normalise_fgsea(fgsea_results),
  normalise_ora(gprofiler_results)
), fill = TRUE)

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

all_enrichment[, theme := classify_theme(term_name)]
all_enrichment[, significant_05 := adjusted_p < 0.05]
all_enrichment[, main_text_candidate := significant_05 & theme != "other" &
                 gene_set_type %in% c("FDR_significant", "mKH_interval_supported", "mKH_both_layers", "preranked_all_genes")]
all_enrichment <- all_enrichment[order(adjusted_p)]
fwrite(all_enrichment, file.path(out_dir, "08_reports", "all_enrichment_results_normalised.csv"))

theme_summary <- all_enrichment[significant_05 == TRUE & theme != "other", .(
  significant_terms = .N,
  best_adjusted_p = min(adjusted_p, na.rm = TRUE),
  methods = paste(sort(unique(enrichment_method)), collapse = "; "),
  databases = paste(sort(unique(database)), collapse = "; "),
  representative_terms = paste(head(unique(term_name[order(adjusted_p)]), 8), collapse = " | ")
), by = .(tissue, omic, model, gene_set_type, theme)]
theme_summary <- theme_summary[order(best_adjusted_p)]
fwrite(theme_summary, file.path(out_dir, "08_reports", "pathway_theme_summary.csv"))

top_terms <- all_enrichment[significant_05 == TRUE][order(adjusted_p), head(.SD, 25),
                                                    by = .(tissue, omic, model, gene_set_type)]
fwrite(top_terms, file.path(out_dir, "08_reports", "top_significant_terms_by_analysis.csv"))

write_package_versions <- function() {
  pkgs <- c("data.table", "openxlsx", "ggplot2", "AnnotationDbi", "org.Hs.eg.db",
            "clusterProfiler", "ReactomePA", "msigdbr", "fgsea", "gprofiler2",
            "igraph", "ggraph")
  versions <- data.table(
    package = pkgs,
    version = vapply(pkgs, function(p) as.character(packageVersion(p)), character(1)),
    loaded = vapply(pkgs, function(p) p %in% loadedNamespaces(), logical(1))
  )
  fwrite(versions, file.path(out_dir, "00_manifest", "R_package_versions.csv"))
}
write_package_versions()

plot_dt <- all_enrichment[
  significant_05 == TRUE & theme != "other" &
    gene_set_type %in% c("FDR_significant", "mKH_interval_supported", "mKH_both_layers", "DL_both_layers", "preranked_all_genes")
][order(adjusted_p), head(.SD, 12), by = .(tissue, omic, model, gene_set_type)]

if (nrow(plot_dt)) {
  plot_dt[, label := paste(tissue, omic, gene_set_type, sep = "\n")]
  plot_dt[, term_short := ifelse(nchar(term_name) > 70, paste0(substr(term_name, 1, 67), "..."), term_name)]
  p <- ggplot(plot_dt, aes(x = -log10(adjusted_p), y = reorder(term_short, -log10(adjusted_p)))) +
    geom_point(aes(size = input_gene_count, colour = theme), alpha = 0.85) +
    facet_wrap(~ label, scales = "free_y") +
    labs(
      x = "-log10(adjusted p-value)",
      y = NULL,
      colour = "Theme",
      size = "Input genes",
      title = "Top pathway enrichment signals across methylation, expression and convergence analyses"
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
writeData(wb, "gene_sets", gene_set_catalog)
addWorksheet(wb, "GO_ORA")
writeData(wb, "GO_ORA", go_ora)
addWorksheet(wb, "ReactomePA_ORA")
writeData(wb, "ReactomePA_ORA", reactome_ora)
addWorksheet(wb, "MSigDB_ORA")
writeData(wb, "MSigDB_ORA", msig_ora)
addWorksheet(wb, "fgsea")
writeData(wb, "fgsea", fgsea_results)
addWorksheet(wb, "gprofiler")
writeData(wb, "gprofiler", gprofiler_results)
addWorksheet(wb, "theme_summary")
writeData(wb, "theme_summary", theme_summary)
addWorksheet(wb, "top_terms")
writeData(wb, "top_terms", top_terms)
saveWorkbook(wb, file.path(out_dir, "cross_omic_pathway_enrichment_R_results.xlsx"), overwrite = TRUE)

script_manifest <- data.table(
  script = "scripts/05_convergence_pathway/run_cross_omic_pathway_enrichment.R",
  role = "production enrichment workflow",
  description = "Builds model-specific backgrounds and threshold/ranked gene sets from the final methylation and expression metafor model outputs, runs GO/Reactome/MSigDB/g:Profiler ORA and fgsea GSEA, then writes tables, figures and reports.",
  run_command = "Rscript scripts/05_convergence_pathway/run_cross_omic_pathway_enrichment.R",
  path_assumption = "Run from the repository root or set ASD_REPO_ROOT to the repository root."
)
fwrite(script_manifest, file.path(out_dir, "00_manifest", "script_role_manifest.csv"))

qc <- data.table(
  check = c(
    "Source files exist",
    "Model-specific backgrounds written",
    "Threshold gene sets written",
    "Convergence backgrounds use tissue-matched intersections",
    "GO ORA completed",
    "Reactome ORA completed",
    "MSigDB ORA completed",
    "fgsea completed",
    "g:Profiler attempted",
    "Figures written",
    "Outputs confined to package results directory"
  ),
  status = c(
    ifelse(all(sources$exists), "PASS", "FAIL"),
    ifelse(nrow(backgrounds) > 0, "PASS", "FAIL"),
    ifelse(nrow(gene_set_catalog) > 0, "PASS", "FAIL"),
    ifelse(nrow(convergence_backgrounds) > 0, "PASS", "FAIL"),
    ifelse(nrow(go_ora) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(reactome_ora) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(msig_ora) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(fgsea_results) > 0, "PASS", "NO_RESULTS"),
    ifelse(nrow(gprofiler_results) > 0, "PASS", "NO_RESULTS_OR_OFFLINE"),
    ifelse(length(list.files(file.path(out_dir, "07_figures"), pattern = "\\.(png|pdf)$")) > 0, "PASS", "NO_RESULTS"),
    "PASS"
  ),
  notes = c(
    "All enrichment sources are final omic model result tables.",
    "Backgrounds are all genes with finite model p-values.",
    "FDR, mKH-interval and DL-interval gene sets are separated.",
    "Cross-omic ORA does not use the union of omic backgrounds.",
    "clusterProfiler::enrichGO used with BH correction and model background.",
    "ReactomePA::enrichPathway used with BH correction and model background.",
    "clusterProfiler::enricher used with msigdbr Hallmark/Reactome/WikiPathways and model background.",
    "fgsea used on signed -log10(p) ranks for all tested genes.",
    "g:Profiler is an external internet-dependent check, not the sole source of results.",
    "Dotplot, theme heatmap, and network are generated when significant terms exist.",
    "This script writes only under results/pathway_enrichment in the code package."
  )
)
fwrite(qc, file.path(out_dir, "08_reports", "pathway_enrichment_QC_checklist.csv"))

readme <- c(
  "# Cross-Omic Pathway Enrichment R Workflow",
  "",
  "This package runs enrichment analyses from the final methylation and expression meta-analysis outputs.",
  "",
  "## Run",
  "",
  "From the repository root:",
  "",
  "```r",
  "Rscript scripts/05_convergence_pathway/run_cross_omic_pathway_enrichment.R",
  "```",
  "",
  "or from a shell:",
  "",
  "```bash",
  "Rscript scripts/05_convergence_pathway/run_cross_omic_pathway_enrichment.R",
  "```",
  "",
  "If running elsewhere, set `ASD_REPO_ROOT` to the repository root first.",
  "",
  "## Primary design choices",
  "",
  "- ORA backgrounds are model-specific tested gene universes, not the whole genome.",
  "- Cross-omic backgrounds are tissue-matched methylation-expression tested-gene intersections.",
  "- FDR-significant, modified Knapp-Hartung interval-supported, and DerSimonian-Laird interval-screening genes are analysed separately.",
  "- fgsea uses all tested genes ranked by signed `-log10(p)`, reducing dependence on a single threshold.",
  "- g:Profiler is included as an external robustness check and is not required as the sole evidence source.",
  "",
  "## Main outputs",
  "",
  "- `cross_omic_pathway_enrichment_R_results.xlsx`",
  "- `08_reports/all_enrichment_results_normalised.csv`",
  "- `08_reports/pathway_theme_summary.csv`",
  "- `07_figures/Figure_top_enrichment_dotplot.png`",
  "- `07_figures/Figure_theme_heatmap.png`",
  "- `07_figures/Figure_analysis_theme_network.png`"
)
writeLines(readme, file.path(out_dir, "README_reproducibility.md"))

story_lines <- c(
  "# Pathway Enrichment Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Interpretation Guardrails",
  "",
  "- Enrichment results are pathway-level prioritisation signals, not proof of causal mechanism.",
  "- The primary ORA background is the tested universe for each model.",
  "- DL interval gene sets are retained as exploratory screening-level inputs.",
  "- Main-text interpretation should prioritise themes recurring across methods, omic layers, or sensitivity models.",
  "",
  "## Strongest Recurrent Themes",
  ""
)
if (nrow(theme_summary)) {
  top_theme_lines <- apply(head(theme_summary, 30), 1, function(x) {
    paste0(
      "- ", x[["tissue"]], " / ", x[["omic"]], " / ", x[["gene_set_type"]],
      ": ", x[["theme"]], " (best adjusted p = ", signif(as.numeric(x[["best_adjusted_p"]]), 3),
      "; representative terms: ", x[["representative_terms"]], ")"
    )
  })
  story_lines <- c(story_lines, top_theme_lines)
} else {
  story_lines <- c(story_lines, "- No adjusted-p < .05 pathway themes were found under the current settings.")
}
story_lines <- c(
  story_lines,
  "",
  "## Suggested Result Interpretation",
  "",
  "Interpretation should focus on pathway themes that recur across thresholded ORA, ranked GSEA, and tissue-matched convergence backgrounds. Terms supported only by DL screening sets or by a single exploratory branch should be described as exploratory or hypothesis-generating."
)
writeLines(story_lines, file.path(out_dir, "08_reports", "pathway_enrichment_results_summary.md"))

method_lines <- c(
  "# Method Comparison and Recommendation",
  "",
  "## Methods attempted",
  "",
  "1. `clusterProfiler::enrichGO`: GO BP/MF/CC over-representation analysis using model-specific backgrounds.",
  "2. `ReactomePA::enrichPathway`: Reactome over-representation analysis using model-specific backgrounds.",
  "3. `clusterProfiler::enricher` with `msigdbr`: Hallmark, Reactome and WikiPathways ORA using model-specific backgrounds.",
  "4. `fgsea`: preranked enrichment over all tested genes using signed `-log10(p)` ranks.",
  "5. `gprofiler2::gost`: external ORA cross-check using custom backgrounds where the web service was available.",
  "",
  "## Recommended hierarchy",
  "",
  "- Primary: model-background ORA on FDR-supported or modified Knapp-Hartung interval-supported gene sets.",
  "- Secondary: preranked fgsea results that support the same biological themes without relying on a hard threshold.",
  "- Exploratory: DL interval-screening gene sets, placenta descriptive results, LCL exploratory results, and single-method terms.",
  "",
  "## Background recommendation",
  "",
  "Use model-specific tested gene backgrounds for within-omic enrichment and tissue-matched tested-background intersections for methylation-expression convergence. This is the most defensible option because methylation arrays, WGBS promoter summaries, expression platforms and LCL/placenta assays do not test identical gene universes."
)
writeLines(method_lines, file.path(out_dir, "08_reports", "method_comparison_and_recommendation.md"))

qc_md <- c(
  "# Pathway Enrichment QC Report",
  "",
  paste0("Output folder: `", out_dir, "`"),
  "",
  "## Checklist",
  "",
  apply(qc, 1, function(x) paste0("- ", x[["check"]], ": ", x[["status"]], ". ", x[["notes"]])),
  "",
  "## Human-review notes",
  "",
  "- DL-screening enrichment is not intended as confirmatory evidence.",
  "- g:Profiler results depend on the external web service and should be treated as a robustness check.",
  "- Pathway labels involving disease names should be interpreted as gene-set containers, not as evidence of those diseases.",
  "- The biological story should be anchored to recurrent themes and tissue/omic context, not isolated pathway names."
)
writeLines(qc_md, file.path(out_dir, "08_reports", "pathway_enrichment_QC_report.md"))

message("Completed cross-omic pathway enrichment workflow at: ", out_dir)

