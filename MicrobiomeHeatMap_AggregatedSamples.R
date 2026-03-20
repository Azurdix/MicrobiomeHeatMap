# =============================================================================================
# Title: Aggregated Microbiome Heatmap - Top 30 Families Relative Abundance
# Author: Mateusz Glenszczyk
# Email: mateusz.glenszczyk@gmail.com
# Date: 2026-04-01
# Project: PhD project - spider-associated microbiome comparison
# =============================================================================================
#
# Description:
# This script generates a log-transformed heatmap of the top 30 bacterial families based on
# relative abundance data. Individual samples are aggregated into predefined biological groups
# before visualization.
#
# Final figure settings:
# - bacterial families (rows) can be displayed either alphabetically or with row clustering
# - aggregated sample groups (columns) are displayed in a fixed order:
#     PTSD-EGGS, PTSD-SILK, PTSD-ENV, PRD-EGGS, PRD-SILK, PRD-ENV
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
# Output files:
# - 06_Exports/family_per_sample_export/top30_family_heatmap_aggregated.png
# - 06_Exports/family_per_sample_export/top30_family_heatmap_aggregated.pdf
#
# Notes:
# - If row dendrogram is turned off, families remain in alphabetical order.
# - If row dendrogram is turned on, families are reordered by clustering.
# - Column clustering is optional, but turning it on will override the fixed group order.
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

# Keep only samples present in both metadata and abundance table
meta <- meta[meta$`sample-id` %in% colnames(tab), , drop = FALSE]
tab  <- tab[, meta$`sample-id`, drop = FALSE]

# Replace NA with 0
tab[is.na(tab)] <- 0

# Remove taxa absent across all samples
tab <- tab[rowSums(tab) > 0, , drop = FALSE]

# =========================
# 3. Detect metadata columns
# =========================

# Species can be stored either in "species" or in "label"
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

# Clean text
species_clean <- trimws(tolower(species_raw))
class_clean   <- trimws(tolower(class_raw))

# =========================
# 4. Build final group labels
# =========================

# Species mapping
species_code <- ifelse(
  grepl("parasteatoda", species_clean), "PTSD",
  ifelse(grepl("pardosa", species_clean), "PRD", NA)
)

# Condition mapping
class_code <- ifelse(
  grepl("egg", class_clean), "EGGS",
  ifelse(grepl("silk", class_clean), "SILK",
         ifelse(grepl("env|environment", class_clean), "ENV", NA))
)

meta$group <- paste(species_code, class_code, sep = "-")

# Remove rows that could not be mapped
meta <- meta[!is.na(species_code) & !is.na(class_code), , drop = FALSE]
tab  <- tab[, meta$`sample-id`, drop = FALSE]

# Fixed order required in final heatmap
desired_order <- c(
  "PTSD-EGGS",
  "PTSD-SILK",
  "PTSD-ENV",
  "PRD-EGGS",
  "PRD-SILK",
  "PRD-ENV"
)

group_order <- desired_order[desired_order %in% meta$group]

if (length(group_order) == 0) {
  stop("No matching aggregated groups were found. Check metadata values in species/label and class.")
}

# =========================
# 5. Aggregate samples by group
# =========================

tab_matrix <- as.matrix(tab)

agg_list <- lapply(group_order, function(g) {
  cols <- which(meta$group == g)
  
  if (length(cols) == 1) {
    as.numeric(tab_matrix[, cols])
  } else {
    rowMeans(tab_matrix[, cols, drop = FALSE], na.rm = TRUE)
  }
})

tab_agg <- do.call(cbind, agg_list)
tab_agg <- as.matrix(tab_agg)

rownames(tab_agg) <- rownames(tab_matrix)
colnames(tab_agg) <- group_order

# Remove taxa absent across all aggregated groups
tab_agg <- tab_agg[rowSums(tab_agg) > 0, , drop = FALSE]

# Sort microbial families alphabetically
# This order will only remain visible if row clustering is turned off.
tab_agg <- tab_agg[sort(rownames(tab_agg)), , drop = FALSE]

# =========================
# 6. Column annotation and colors
# =========================

annotation_col <- data.frame(
  Group = colnames(tab_agg),
  row.names = colnames(tab_agg),
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
# 7. Log-transform
# =========================

tab_log <- log10(tab_agg + 0.01)

# =========================
# 8. Heatmap settings
# =========================

hm_cols <- colorRampPalette(c("#31079C", "white", "#C80946"))(100)

max_val <- max(abs(tab_log))
hm_breaks <- seq(-max_val, max_val, length.out = 101)

# =========================
# 8A. Dendrogram settings
# =========================

# TRUE  = show dendrogram
# FALSE = hide dendrogram

show_row_dendrogram <- TRUE
show_col_dendrogram <- FALSE

# Recommended setup for your project:
# show_row_dendrogram <- TRUE
# show_col_dendrogram <- FALSE

# =========================
# 9. Plot in R session
# =========================

pheatmap(
  tab_log,
  cluster_rows = show_row_dendrogram,
  cluster_cols = show_col_dendrogram,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 9,
  cellwidth = 18,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (aggregated relative abundance)"
)

# =========================
# 10. Save PNG
# =========================

png(
  "06_Exports/family_per_sample_export/top30_family_heatmap_aggregated.png",
  width = 1400,
  height = 1400,
  res = 200
)

pheatmap(
  tab_log,
  cluster_rows = show_row_dendrogram,
  cluster_cols = show_col_dendrogram,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 9,
  cellwidth = 18,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (aggregated relative abundance)"
)

dev.off()

# =========================
# 11. Save PDF
# =========================

pdf(
  "06_Exports/family_per_sample_export/top30_family_heatmap_aggregated.pdf",
  width = 9,
  height = 9
)

pheatmap(
  tab_log,
  cluster_rows = show_row_dendrogram,
  cluster_cols = show_col_dendrogram,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 9,
  cellwidth = 18,
  border_color = NA,
  color = hm_cols,
  breaks = hm_breaks,
  main = "Top 30 families (aggregated relative abundance)"
)

dev.off()
