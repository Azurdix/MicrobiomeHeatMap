# Microbiome Heatmaps - PRD / PTSD / Compartments + Core Microbiome

Welcome!

This directory contains five scripts for generating microbiome heatmaps.

## Overview

The scripts in this folder reflect three stages of development.

### Early versions

The first two scripts (`MicrobiomeHeatMap_*`) are earlier versions developed during the initial stage of building and learning the metagenomic analysis workflow.  
They generate valid heatmaps and were useful for preliminary data exploration, visualization, and testing of the pipeline.

### Updated versions

The newer scripts (`Make_Heatmaps_*`) were developed after refining the analytical workflow and improving the biological interpretability of the results.  
These scripts are currently the **recommended versions** for generating publication-ready heatmaps.

### Core Microbiome Version

The very last script, **`CoreMicrobiome`**, introduces additional filtering levels beyond those used in the standard heatmap workflow.

## Why were the scripts updated?

During the analysis, it became clear that QIIME2-derived outputs may still contain a substantial number of ambiguous or low-informative taxa, such as:

- `uncultured`
- `unclassified`
- `incertae sedis`
- `subgroup`, `environmental group`, `terrestrial group`, etc.
- technical placeholder labels such as `d__Bacteria_X` or `d__Subgroup_Y`

Such entries introduce noise and may obscure biological interpretation.  
To address this, the updated scripts include an additional post-processing filtering step that removes ambiguous taxa still present after the standard QIIME2 workflow.

---

## Key Improvements in the Updated Scripts

- **Additional taxonomic filtering**  
  Removes ambiguous and low-informative taxa to improve clarity and biological interpretability.

- **Relative abundance scaling (% read abundance)**  
  Converts data to relative abundance values expressed as percentages and visualizes them on a log10 scale.

- **Improved legend design**  
  Uses more interpretable abundance levels:
  - `0.01`
  - `0.10`
  - `1.00`
  - `10.00` (%)

- **Non-symmetric color scale**  
  Better reflects the strongly skewed structure of microbiome data.

- **Flexible top-N selection**  
  Allows the user to define how many of the most abundant taxa (families or genera) are included in the final heatmap.

- **Improved visualization aesthetics**  
  Provides a more consistent layout, clearer annotation, and improved readability.

---

## Current Status

The updated scripts represent the currently recommended workflow for generating publication-ready microbiome heatmaps.

Further refinements may be introduced in the future, but the current versions are considered stable, complete, and suitable for final analyses.
