#!/usr/bin/env Rscript

# PRISMA-style dataset search and validation flow diagram.
#
# This script redraws the final dataset search counts in a matched two-lane
# PRISMA-style format for paper/supplement use. It does not rerun searches.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
env_root <- Sys.getenv("ASD_REPO_ROOT", unset = "")
package_root <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

out_root <- file.path(package_root, "results", "figures", "dataset_search_flow")
source_dir <- file.path(out_root, "figure_source_data")
figure_dir <- file.path(out_root, "figures")
report_dir <- file.path(out_root, "reports")
qc_dir <- file.path(out_root, "qc")
for (d in c(source_dir, figure_dir, report_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

wrap_text <- function(x, width = 32) {
  vapply(x, function(s) {
    parts <- unlist(strsplit(s, "\n", fixed = TRUE))
    paste(vapply(parts, function(part) paste(strwrap(part, width = width), collapse = "\n"), character(1)),
          collapse = "\n")
  }, character(1))
}

nodes <- data.table(
  id = c(
    "m_id", "m_val", "m_eligible", "m_included", "m_screen_excl", "m_elig_excl", "m_control",
    "e_id", "e_val", "e_eligible", "e_included", "e_screen_excl", "e_elig_excl", "e_control"
  ),
  lane = c(rep("DNA methylation", 7), rep("Gene expression", 7)),
  stage = c(
    "Identification", "Screening", "Eligibility", "Included", "Screening", "Eligibility", "Included",
    "Identification", "Screening", "Eligibility", "Included", "Screening", "Eligibility", "Included"
  ),
  box_type = c(
    "main", "main", "main", "included", "excluded", "excluded", "controlled",
    "main", "main", "main", "included", "excluded", "excluded", "controlled"
  ),
  x = c(2.45, 2.45, 2.45, 2.45, 5.10, 5.10, 5.10, 8.25, 8.25, 8.25, 8.25, 10.90, 10.90, 10.90),
  y = c(7.55, 6.05, 4.55, 2.90, 6.05, 4.55, 2.90, 7.55, 6.05, 4.55, 2.90, 6.05, 4.55, 2.90),
  w = c(2.45, 2.45, 2.45, 2.45, 2.40, 2.40, 2.40, 2.45, 2.45, 2.45, 2.45, 2.40, 2.40, 2.40),
  h = c(0.80, 0.80, 0.80, 0.98, 0.88, 1.08, 0.88, 0.80, 0.80, 0.80, 0.98, 0.88, 1.08, 0.88),
  label = c(
    "Raw/accession records screened\nn = 97\nrepository/accession molecular records",
    "Validated methylation-related routes\nn = 53",
    "Eligible broad ASD-control methylation routes\nn = 22",
    "Open datasets included in completed methylation analyses\nn = 19\nBlood 6 | Brain 9 | Placenta 2 | LCL 2",
    "Records not retained after screening\nn = 44\nnon-methylation or non-relevant molecular records",
    "Routes excluded or held\nn = 31\nnot eligible/insufficient broad data, duplicate/superseries, context-only or non-comparable assay",
    "Controlled-access methylation routes retained for application/check\nn = 3",
    "Raw expression hits screened\nn = 408\n173 initial GEO hits + 235 supplemental broad GEO hits",
    "Validated expression routes\nn = 59",
    "Eligible broad expression routes\nn = 29",
    "Public expression routes included in completed analyses\nn = 26\nBlood-family 8 | Brain 12 | LCL 6 | Placenta 0",
    "Hits not retained after screening\nn = 349\nduplicates, irrelevant records or non-expression/non-ASD hits",
    "Routes excluded or held\nn = 30\nphenotype/matrix review, context/specialised, not usable or duplicate/overlap",
    "Controlled-access expression routes retained for application/check\nn = 3"
  )
)

nodes[, label_wrapped := wrap_text(label, width = 31)]
nodes[, label_wrapped := wrap_text(label, width = 28)]
nodes[box_type %in% c("excluded", "controlled"), label_wrapped := wrap_text(label, width = 31)]

edges <- data.table(
  from = c(
    "m_id", "m_val", "m_val", "m_eligible", "m_eligible", "m_included",
    "e_id", "e_val", "e_val", "e_eligible", "e_eligible", "e_included"
  ),
  to = c(
    "m_val", "m_screen_excl", "m_eligible", "m_elig_excl", "m_included", "m_control",
    "e_val", "e_screen_excl", "e_eligible", "e_elig_excl", "e_included", "e_control"
  )
)

edges <- merge(edges, nodes[, .(from = id, x_from = x, y_from = y, w_from = w, h_from = h)], by = "from")
edges <- merge(edges, nodes[, .(to = id, x_to = x, y_to = y, w_to = w, h_to = h)], by = "to")
edges[, `:=`(
  x_start = fifelse(abs(x_to - x_from) < 0.01, x_from, x_from + sign(x_to - x_from) * w_from / 2),
  y_start = fifelse(abs(x_to - x_from) < 0.01, y_from - h_from / 2, y_from),
  x_end = fifelse(abs(x_to - x_from) < 0.01, x_to, x_to - sign(x_to - x_from) * w_to / 2),
  y_end = fifelse(abs(x_to - x_from) < 0.01, y_to + h_to / 2, y_to)
)]

stage_labels <- data.table(
  stage = c("Identification", "Screening", "Eligibility", "Included"),
  x = 0.72,
  y = c(7.55, 6.05, 4.55, 2.90)
)

lane_labels <- data.table(
  x = c(2.45, 8.25),
  y = c(8.50, 8.50),
  label = c("DNA methylation", "Gene expression")
)

fills <- c(
  main = "#DCEEFF",
  included = "#BFE3FF",
  excluded = "#EEF6FD",
  controlled = "#E6F2FF"
)
strokes <- c(
  main = "#2F6F9F",
  included = "#1E6E9E",
  excluded = "#6B8AA3",
  controlled = "#2F6F9F"
)
linetypes <- c(main = "solid", included = "solid", excluded = "solid", controlled = "dashed")

fwrite(nodes[, .(id, lane, stage, box_type, label, x, y, w, h)], file.path(source_dir, "dataset_search_prisma_nodes.csv"))
fwrite(edges[, .(from, to)], file.path(source_dir, "dataset_search_prisma_edges.csv"))

p <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
    linewidth = 0.55,
    colour = "#60758A",
    lineend = "round"
  ) +
  geom_rect(
    data = nodes,
    aes(
      xmin = x - w / 2,
      xmax = x + w / 2,
      ymin = y - h / 2,
      ymax = y + h / 2,
      fill = box_type,
      colour = box_type,
      linetype = box_type
    ),
    linewidth = 0.8
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label_wrapped),
    size = 3.0,
    lineheight = 0.88,
    colour = "#102033",
    family = "sans"
  ) +
  geom_text(
    data = lane_labels,
    aes(x = x, y = y, label = label),
    size = 4.6,
    fontface = "bold",
    colour = "#102033"
  ) +
  geom_text(
    data = stage_labels,
    aes(x = x, y = y, label = stage),
    angle = 90,
    size = 3.2,
    fontface = "bold",
    colour = "#3C5D78"
  ) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_colour_manual(values = strokes, guide = "none") +
  scale_linetype_manual(values = linetypes, guide = "none") +
  coord_cartesian(xlim = c(0.2, 12.35), ylim = c(2.05, 8.45), expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 10, 8, 10)
  )

ggsave(file.path(figure_dir, "Figure_dataset_search_flow_PRISMA_style.png"), p, width = 12.4, height = 6.6, dpi = 360)
ggsave(file.path(figure_dir, "Figure_dataset_search_flow_PRISMA_style.pdf"), p, width = 12.4, height = 6.6)
ggsave(file.path(figure_dir, "Figure_dataset_search_flow_PRISMA_style.tiff"), p, width = 12.4, height = 6.6, dpi = 360, compression = "lzw")

qc <- data.table(
  item = c(
    "methylation_records_screened",
    "methylation_validated_routes",
    "methylation_not_retained_after_screening",
    "methylation_eligible_routes",
    "methylation_routes_excluded_or_held",
    "methylation_open_analysed",
    "methylation_controlled_access",
    "expression_raw_hits_screened",
    "expression_validated_routes",
    "expression_not_retained_after_screening",
    "expression_routes_retained_completed_or_future",
    "expression_routes_excluded_or_held",
    "expression_public_completed_analysis_routes",
    "expression_controlled_access"
  ),
  count = c(97, 53, 44, 22, 31, 19, 3, 408, 59, 349, 29, 30, 26, 3),
  note = c(
    "Final repository/accession molecular inventory; not a simple methylation-only GEO raw-hit count.",
    "Final methylation-related route records retained.",
    "Calculated as 97 - 53; records not retained in the methylation route flow.",
    "Includes 19 open analysed routes and 3 controlled-access application/check routes.",
    "Calculated as 53 - 22; reason categories summarised in figure.",
    "Completed broad methylation analyses.",
    "Controlled-access methylation application/check routes.",
    "173 initial GEO expression hits plus 235 supplemental broad GEO hits.",
    "Final validated expression route inventory.",
    "Calculated as 408 - 59; screening-level removal/deduplication/ineligible raw hits.",
    "Calculated as 26 completed public analysis routes plus 3 controlled-access future routes.",
    "Calculated as 59 - 29; reason categories summarised in figure.",
  "Completed public expression analysis routes used in the current study.",
    "Controlled-access expression application/check routes."
  )
)
fwrite(qc, file.path(qc_dir, "dataset_search_flow_PRISMA_style_QC_counts.csv"))

message("PRISMA-style dataset flow diagram written to: ", normalizePath(figure_dir, winslash = "/", mustWork = TRUE))
