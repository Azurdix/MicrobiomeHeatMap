# =============================================================================================
# Author: Mateusz Glenszczyk
# Email: mateusz.glenszczyk@gmail.com
# Date: 2026-04-01
# Description: Microbiome Heatmap - Top 30 Families Relative Abundance Visualization
# =============================================================================================
#
# **Purpose:**
# Generates a log-transformed heatmap of the top 30 bacterial families' relative abundance across samples.
# Part of a PhD project comparing microbiome composition across multiple samples.
# It visualizes patterns and clusters in microbiome composition for comparative analysis.
# My sample included two species (Parasteatoda vs Pardosa) + three spider egg sac conditions (Environment + Silk + Eggs).
#
# **Inputs:**
# - `top30_family_relative_short_clean.tsv`: Pre-filtered relative abundance table (top 30 families).
# - `metadata_plik.tsv`: Sample metadata with group assignments (used for column annotation).
#
# **Key Steps:**
# - Log10-transforms relative abundance data for better contrast.
# - Annotates columns by experimental group using metadata.
# - Generates interactive and high-resolution heatmaps (PNG/PDF).
#
# **Outputs:**
# - `top30_family_heatmap-color.png`: High-resolution heatmap image.
# - `top30_family_heatmap-color.pdf`: Vector-format heatmap for publications.
#
# **Dependencies:**
# - R (>= 4.0.0)
# - Packages: `pheatmap`
#
# **Notes:**
# - NA values are replaced with 0 before log-transformation.
# - Symmetric color breaks ensure balanced visualization of up-/down-regulated taxa.
# - For troubleshooting: PDF output is guaranteed; PNG may fail on some systems.
# =============================================================================================



#This to clean environment
rm(list = ls())

#This, so you can even do a heatmap
library(pheatmap)

# Read Family Table
tab <- read.table(
  "06_Exports/family_per_sample_export/top30_family_relative_short_clean.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# Read Metadata
meta <- read.table(
  "01_Metadata/metadata_plik.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Columns to Numeric
tab <- as.data.frame(
  lapply(tab, as.numeric),
  row.names = rownames(tab),
  check.names = FALSE
)

# Sample Order
meta <- meta[order(meta$class, meta$label), ]
sample_order <- meta$`sample-id`
sample_order <- sample_order[sample_order %in% colnames(tab)]

# Column Order
tab <- tab[, sample_order, drop = FALSE]

# NA to 0
tab[is.na(tab)] <- 0

# Deletion of 0
tab <- tab[rowSums(tab) > 0, , drop = FALSE]

# Column Adnotation
annotation_col <- meta[
  match(colnames(tab), meta$`sample-id`),
  "class",
  drop = FALSE
]
rownames(annotation_col) <- colnames(tab)

# Log-transform, so the heatmap has better contrast.
tab_log <- log10(as.matrix(tab) + 0.01)

# SLAYING with colors (Serve only with the contrasting ones!!! Otherwise prepare to flop)
hm_cols <- colorRampPalette(c("#31079C", "white", "#C80946"))(100)

# Symetric breaks
max_val <- max(abs(tab_log))
hm_breaks <- seq(-max_val, max_val, length.out = 101)

# Heatmap appears in R window. Might be useless if on Windows - I had some issues on Linux with it.
pheatmap(
  tab_log,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (log-transformed relative abundance)"
)

# PNG
png(
  "06_Exports/family_per_sample_export/top30_family_heatmap-color.png",
  width = 1800,
  height = 1400,
  res = 200
)

pheatmap(
  tab_log,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (log-transformed relative abundance)"
)

dev.off()

# PDF
pdf(
  "06_Exports/family_per_sample_export/top30_family_heatmap-color.pdf",
  width = 12,
  height = 9
)
# If PNG did not save properly, this thingie down here will definitley save it.
png("06_Exports/family_per_sample_export/top30_family_heatmap-color.png",
    width = 1800, height = 1400, res = 200)

pheatmap(
  tab_log,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (log-transformed relative abundance)"
)

dev.off()

