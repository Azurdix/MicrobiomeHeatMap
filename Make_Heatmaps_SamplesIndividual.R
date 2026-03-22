# =============================================================================================
# Title: Microbiome Heatmap - Top 30 Taxa (Individual Samples)
# Author: Mateusz Glenszczyk
# Project: PhD project - spider-associated microbiome comparison
# =============================================================================================
#
# Description:
# This script generates heatmaps of the top n (30) bacterial taxa (family and genus level)
# based on QIIME2-derived feature tables at the level of individual samples.
#
# Workflow summary:
# - Taxonomic strings are simplified to the final rank (family/genus) and cleaned
# - Ambiguous taxa (e.g., "uncultured", "unknown", "subgroup", "soil/terrestrial group",
#   "Incertae Sedis", "chloroplast", "mitochondria") are removed
# - Raw counts are converted to relative abundance (%)
# - Top 30 taxa are selected based on mean abundance across samples
# - Data are log10-transformed (with pseudocount) for visualization
#
# Heatmap:
# - Rows: taxa
# - Columns: individual samples
# - Samples are ordered by biological groups:
#     PRD ENV → PRD SILK → PRD EGGS → PTSD ENV → PTSD SILK → PTSD EGGS
# - Row clustering enabled; column clustering disabled
# - Color scale represents % read abundance (log10 scale)
#
# Output:
# - PDF and high-resolution PNG heatmaps
# - Processed relative abundance tables
#
# Note:
# This script visualizes sample-level variation; for statistical analysis
# (e.g., DESeq2), raw count data should be used.
#
# =============================================================================================





rm(list = ls())

library(pheatmap)
library(grid)

# ============================================================================================
# SETTINGS
# ============================================================================================

workdir <- "/home/Azu/Pulpit/data_NGS/HeatMapy/qiime_heatmap_work"
setwd(workdir)

family_file <- "family_table.filtered.tsv"
genus_file  <- "genus_table.filtered.tsv"

top_n <- 30
pseudocount_percent <- 0.01

legend_percent_values <- c(0.01, 0.10, 1.00, 10.00)
legend_breaks <- log10(legend_percent_values)

scale_min <- log10(0.01)
scale_max <- log10(10.00)
scale_breaks <- seq(scale_min, scale_max, length.out = 101)

sample_type_colors <- c(
  "PRD ENV"   = "#A1D99B",
  "PRD SILK"  = "#9ECAE1",
  "PRD EGGS"  = "#FEE391",
  "PTSD ENV"  = "#31A354",
  "PTSD SILK" = "#3182BD",
  "PTSD EGGS" = "#E6AB02"
)

heat_colors <- colorRampPalette(
  c("#253494", "#2C7FB8", "#41B6C4", "#A1DAB4", "#FFFFBF", "#FDAE61", "#D7191C")
)(100)

# ============================================================================================
# FUNCTIONS
# ============================================================================================

read_qiime_tsv <- function(file) {
  lines <- readLines(file, warn = FALSE)
  
  if (length(lines) == 0) {
    stop(paste("Empty file:", file))
  }
  
  if (startsWith(lines[1], "#")) {
    df <- read.delim(file, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    df <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  }
  
  colnames(df)[1] <- "Taxon"
  df
}

clean_taxon_name <- function(x) {
  x <- trimws(x)
  parts <- strsplit(x, ";")[[1]]
  x <- tail(parts, 1)
  x <- sub("^[a-z]__", "", x, ignore.case = TRUE)
  x <- trimws(x)
  x <- gsub("^_+$", "", x)
  x <- trimws(x)
  
  if (is.na(x) || x == "") {
    return(NA_character_)
  }
  
  x
}

is_unwanted_taxon <- function(x) {
  if (is.na(x) || trimws(x) == "") {
    return(TRUE)
  }
  
  x2 <- tolower(trimws(x))
  
  patterns <- c(
    "^_+$",
    "uncultured",
    "unclassified",
    "unknown",
    "incertae[ _-]?sedis",
    "subgroup",
    "soil[ _-]?group",
    "terrestrial[ _-]?group",
    "chloroplast",
    "mitochondria",
    "norank",
    "metagenome",
    "environmental",
    "candidate",
    "ambiguous",
    "bacteriap[0-9]*",
    "unknown_family",
    "unknown_genus",
    "^[0-9]+-[0-9]+$",
    "^[a-z0-9]+-[a-z0-9-]+$"
  )
  
  any(sapply(patterns, function(p) grepl(p, x2, perl = TRUE)))
}

get_sample_group <- function(sample_name) {
  s <- toupper(trimws(sample_name))
  
  species <- if (grepl("^PRD", s)) {
    "PRD"
  } else if (grepl("^PTSD", s)) {
    "PTSD"
  } else {
    NA_character_
  }
  
  compartment <- if (grepl("EGGS", s)) {
    "EGGS"
  } else if (grepl("SILK", s)) {
    "SILK"
  } else if (grepl("(^|[-_])E[0-9]+$", s)) {
    "ENV"
  } else {
    NA_character_
  }
  
  if (is.na(species) || is.na(compartment)) {
    return(NA_character_)
  }
  
  paste(species, compartment)
}

get_group_rank <- function(group_name) {
  order_levels <- c(
    "PRD ENV",
    "PRD SILK",
    "PRD EGGS",
    "PTSD ENV",
    "PTSD SILK",
    "PTSD EGGS"
  )
  match(group_name, order_levels)
}

prepare_heatmap_object <- function(file, top_n = 30, pseudocount_percent = 0.01) {
  df <- read_qiime_tsv(file)
  df <- df[!is.na(df$Taxon), , drop = FALSE]
  
  sample_cols <- setdiff(colnames(df), "Taxon")
  
  for (col in sample_cols) {
    df[[col]] <- as.numeric(df[[col]])
    df[[col]][is.na(df[[col]])] <- 0
  }
  
  df$Taxon <- vapply(df$Taxon, clean_taxon_name, character(1))
  df <- df[!is.na(df$Taxon), , drop = FALSE]
  df <- df[!vapply(df$Taxon, is_unwanted_taxon, logical(1)), , drop = FALSE]
  
  if (nrow(df) == 0) {
    stop(paste("No taxa left after filtering in:", file))
  }
  
  agg <- aggregate(df[, sample_cols, drop = FALSE], by = list(Taxon = df$Taxon), FUN = sum)
  rownames(agg) <- agg$Taxon
  agg$Taxon <- NULL
  
  col_sums <- colSums(agg)
  rel <- sweep(agg, 2, col_sums, "/")
  rel[is.na(rel)] <- 0
  
  mean_abund <- rowMeans(rel)
  rel <- rel[order(mean_abund, decreasing = TRUE), , drop = FALSE]
  rel <- rel[seq_len(min(top_n, nrow(rel))), , drop = FALSE]
  
  rel_percent <- rel * 100
  rel_percent_out <- rel_percent
  
  rel_percent_plot <- rel_percent
  rel_percent_plot[rel_percent_plot < pseudocount_percent] <- pseudocount_percent
  rel_plot_log <- log10(rel_percent_plot)
  
  sample_info <- data.frame(
    Sample = colnames(rel_plot_log),
    SampleType = vapply(colnames(rel_plot_log), get_sample_group, character(1)),
    stringsAsFactors = FALSE
  )
  
  if (any(is.na(sample_info$SampleType))) {
    print(sample_info)
    stop("Some samples could not be assigned to Sample Type.")
  }
  
  sample_info$Rank <- vapply(sample_info$SampleType, get_group_rank, integer(1))
  sample_info <- sample_info[order(sample_info$Rank, sample_info$Sample), , drop = FALSE]
  
  rel_plot_log <- rel_plot_log[, sample_info$Sample, drop = FALSE]
  rel_percent_out <- rel_percent_out[, sample_info$Sample, drop = FALSE]
  
  annotation_col <- data.frame(
    "Sample Type" = factor(
      sample_info$SampleType,
      levels = c("PRD ENV", "PRD SILK", "PRD EGGS", "PTSD ENV", "PTSD SILK", "PTSD EGGS")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(annotation_col) <- sample_info$Sample
  
  group_sizes <- as.numeric(table(annotation_col[["Sample Type"]]))
  gaps_col <- cumsum(group_sizes)
  gaps_col <- gaps_col[-length(gaps_col)]
  
  list(
    mat = rel_plot_log,
    rel_percent = rel_percent_out,
    annotation_col = annotation_col,
    gaps_col = gaps_col
  )
}

draw_heatmap <- function(mat, annotation_col, gaps_col, main_title) {
  pheatmap(
    mat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_colors = list("Sample Type" = sample_type_colors),
    gaps_col = gaps_col,
    
    show_colnames = TRUE,
    show_rownames = TRUE,
    
    fontsize_row = 10,
    fontsize_col = 10,
    
    cellwidth = 22,
    cellheight = 16,
    
    border_color = "grey85",
    
    treeheight_row = 90,
    treeheight_col = 0,
    
    color = heat_colors,
    breaks = scale_breaks,
    legend_breaks = legend_breaks,
    legend_labels = sprintf("%.2f", legend_percent_values),
    
    main = main_title,
    angle_col = 45,
    silent = TRUE
  )
}


save_heatmap_files <- function(obj, title_text, pdf_name, png_name) {
  hm <- draw_heatmap(
    mat = obj$mat,
    annotation_col = obj$annotation_col,
    gaps_col = obj$gaps_col,
    main_title = title_text
  )
  
  # poszerzenie legendy
  legend_id <- which(hm$gtable$layout$name == "legend")
  if (length(legend_id) > 0) {
    legend_col <- hm$gtable$layout$l[legend_id]
    hm$gtable$widths[legend_col] <- unit(2.4, "cm")
  }
  
  # usuń obramowanie tylko z górnego paska adnotacji
  ann_id <- which(hm$gtable$layout$name %in% c("col_annotation", "annotation_col"))
  if (length(ann_id) > 0) {
    for (i in ann_id) {
      grob <- hm$gtable$grobs[[i]]
      
      if (!is.null(grob$children)) {
        for (j in seq_along(grob$children)) {
          child <- grob$children[[j]]
          if (!is.null(child$gp)) {
            child$gp$col <- NA
            grob$children[[j]] <- child
          }
        }
      }
      
      if (!is.null(grob$gp)) {
        grob$gp$col <- NA
      }
      
      hm$gtable$grobs[[i]] <- grob
    }
  }
  
  pdf(pdf_name, width = 16, height = 11)
  grid.newpage()
  grid.draw(hm$gtable)
  dev.off()
  
  png(png_name, width = 4800, height = 3300, res = 300)
  grid.newpage()
  grid.draw(hm$gtable)
  dev.off()
}

# ============================================================================================
# MAIN
# ============================================================================================

cat("Working directory:\n")
cat(getwd(), "\n\n")

family_obj <- prepare_heatmap_object(
  file = family_file,
  top_n = top_n,
  pseudocount_percent = pseudocount_percent
)

genus_obj <- prepare_heatmap_object(
  file = genus_file,
  top_n = top_n,
  pseudocount_percent = pseudocount_percent
)

cat("=== FAMILY sample grouping ===\n")
print(data.frame(
  Sample = colnames(family_obj$mat),
  SampleType = family_obj$annotation_col[["Sample Type"]]
))

cat("\n=== GENUS sample grouping ===\n")
print(data.frame(
  Sample = colnames(genus_obj$mat),
  SampleType = genus_obj$annotation_col[["Sample Type"]]
))

write.table(
  family_obj$rel_percent,
  file = "family_heatmap_input_percent.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  genus_obj$rel_percent,
  file = "genus_heatmap_input_percent.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

save_heatmap_files(
  obj = family_obj,
  title_text = "Top 30 bacterial families across samples",
  pdf_name = "FINAL_family_heatmap.pdf",
  png_name = "FINAL_family_heatmap.png"
)

save_heatmap_files(
  obj = genus_obj,
  title_text = "Top 30 bacterial genera across samples",
  pdf_name = "FINAL_genus_heatmap.pdf",
  png_name = "FINAL_genus_heatmap.png"
)

cat("\nSaved files:\n")
print(list.files(pattern = "^FINAL_.*\\.(pdf|png)$"))

