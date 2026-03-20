# =============================================================================================
# Title: Microbiome Heatmap - Top 30 Families Relative Abundance (Per-Sample Visualization)
# Author: Mateusz Glenszczyk
# Email: mateusz.glenszczyk@gmail.com
# Date: 2026-04-01
# Project: PhD project - spider-associated microbiome comparison
# =============================================================================================
#
# Description:
# This script generates a log-transformed heatmap of the top 30 bacterial families based on
# relative abundance data for individual samples. Unlike the aggregated version, each column
# in the heatmap represents a single sample rather than a grouped mean profile.
#
# Final figure settings:
# - bacterial families (rows) are sorted alphabetically
# - samples (columns) are displayed in fixed biological group order:
#     PTSD-EGGS, PTSD-SILK, PTSD-ENV, PRD-EGGS, PRD-SILK, PRD-ENV
# - samples within each group are sorted alphabetically by sample ID
#
# Group labels:
# - PTSD = Parasteatoda
# - PRD  = Pardosa
#
# Metadata requirements:
# - sample-id
# - class
# - and either:
#     * species
#   or
#     * label
#
# Expected condition labels:
# - eggs / silk / env
#   (the script also accepts variants such as "environment")
#
# Input files:
# - 06_Exports/family_per_sample_export/top30_family_relative_short_clean.tsv
# - 01_Metadata/metadata_plik.tsv
#
# Output files:
# - 06_Exports/family_per_sample_export/top30_family_heatmap-color.png
# - 06_Exports/family_per_sample_export/top30_family_heatmap-color.pdf
#
# Dependencies:
# - R (>= 4.0.0)
# - pheatmap
#
# Notes:
# - NA values are replaced with 0 before log-transformation.
# - Symmetric color breaks are used for balanced visualization.
# - Row clustering is disabled to preserve alphabetical family order.
# - Column clustering is disabled to preserve predefined sample order.
# =============================================================================================

# Clean environment
rm(list = ls())

# Load package
library(pheatmap)

# =========================
# 1. Read input files
# =========================

tab <- read.table(
  "06_Exports/family_per_sample_export/top30_family_relative_short_clean.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- read.table(
  "01_Metadata/metadata_plik.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# =========================
# 2. Prepare abundance table
# =========================

tab <- as.data.frame(
  lapply(tab, as.numeric),
  row.names = rownames(tab),
  check.names = FALSE
)

# Keep only shared samples
meta <- meta[meta$`sample-id` %in% colnames(tab), , drop = FALSE]
tab  <- tab[, meta$`sample-id`, drop = FALSE]

# Replace NA with 0
tab[is.na(tab)] <- 0

# Remove taxa absent across all samples
tab <- tab[rowSums(tab) > 0, , drop = FALSE]

# =========================
# 3. Detect metadata fields
# =========================

if ("species" %in% colnames(meta)) {
  species_raw <- meta$species
} else if ("label" %in% colnames(meta)) {
  species_raw <- meta$label
} else {
  stop("Metadata must contain either a 'species' column or a 'label' column.")
}

if (!"class" %in% colnames(meta)) {
  stop("Metadata must contain a 'class' column.")
}

class_raw <- meta$class

species_clean <- trimws(tolower(species_raw))
class_clean   <- trimws(tolower(class_raw))

# =========================
# 4. Build sample group labels
# =========================

species_code <- ifelse(
  grepl("parasteatoda", species_clean), "PTSD",
  ifelse(grepl("pardosa", species_clean), "PRD", NA)
)

class_code <- ifelse(
  grepl("egg", class_clean), "EGGS",
  ifelse(grepl("silk", class_clean), "SILK",
         ifelse(grepl("env|environment", class_clean), "ENV", NA))
)

meta$group <- paste(species_code, class_code, sep = "-")

# Remove unmapped rows
meta <- meta[!is.na(species_code) & !is.na(class_code), , drop = FALSE]
tab  <- tab[, meta$`sample-id`, drop = FALSE]

# =========================
# 5. Define fixed group order and sample order
# =========================

desired_order <- c(
  "PTSD-EGGS",
  "PTSD-SILK",
  "PTSD-ENV",
  "PRD-EGGS",
  "PRD-SILK",
  "PRD-ENV"
)

meta$group <- factor(meta$group, levels = desired_order)

# Sort by biological group first, then alphabetically by sample ID
meta <- meta[order(meta$group, meta$`sample-id`), , drop = FALSE]

sample_order <- meta$`sample-id`
sample_order <- sample_order[sample_order %in% colnames(tab)]

# Reorder abundance table
tab <- tab[, sample_order, drop = FALSE]

# =========================
# 6. Sort families alphabetically
# =========================

tab <- tab[sort(rownames(tab)), , drop = FALSE]

# =========================
# 7. Column annotation and colors
# =========================

annotation_col <- data.frame(
  Group = meta$group[match(colnames(tab), meta$`sample-id`)],
  row.names = colnames(tab),
  check.names = FALSE
)

annotation_colors <- list(
  Group = c(
    "PTSD-EGGS" = "#FFAB5C",
    "PTSD-SILK" = "#8DA3FC",
    "PTSD-ENV"  = "#5FFC6F",
    "PRD-EGGS"  = "#E66F00",
    "PRD-SILK"  = "#345DFA",
    "PRD-ENV"   = "#026B0C"
  )
)

# =========================
# 8. Log transformation
# =========================

tab_log <- log10(as.matrix(tab) + 0.01)

# =========================
# 9. Heatmap settings
# =========================

hm_cols <- colorRampPalette(c("#31079C", "white", "#C80946"))(100)

max_val <- max(abs(tab_log))
hm_breaks <- seq(-max_val, max_val, length.out = 101)

# =========================
# 10. Preview heatmap
# =========================

pheatmap(
  tab_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  cellwidth = 18,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (log-transformed relative abundance)"
)

# =========================
# 11. Save PNG
# =========================

png(
  "06_Exports/family_per_sample_export/top30_family_heatmap-color.png",
  width = 1800,
  height = 1800,
  res = 200
)

pheatmap(
  tab_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  cellwidth = 18,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (log-transformed relative abundance)"
)

dev.off()

# =========================
# 12. Save PDF
# =========================

pdf(
  "06_Exports/family_per_sample_export/top30_family_heatmap-color.pdf",
  width = 9,
  height = 9
)

pheatmap(
  tab_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  cellwidth = 18,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (log-transformed relative abundance)"
)

dev.off()

