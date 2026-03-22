Hello and welcome!
This directory contains four scripts related to microbiome heatmap generation.

Early versions:
The first two scripts (MicrobiomeHeatMap_*) represent early versions developed during my learning phase of the metagenomic analysis pipeline.
They produce valid heatmaps and were useful for initial exploration and visualization of the data.

Updated scripts:
The newer scripts (Make_Heatmaps_*) were developed after refining the analytical approach.
Why the update?

During the analysis, it became clear that outputs from the QIIME2 workflow still contain a substantial amount of ambiguous or uninformative taxa, such as:

uncultured
unclassified
incertae sedis
subgroup / environmental group / terrestrial group etc.
technical labels (e.g. d__Bacteria_X, d_Subgroup_Y)

These entries introduce noise and can obscure biological interpretation.
The updated scripts address this by applying an additional post-processing filtering step, removing ambiguous taxa that remain after standard QIIME2 processing.

Key improvements
- Additional taxonomic filtering
   Removes ambiguous and low-informative taxa to improve biological clarity.
- Relative abundance scaling (% read abundance)
   Data are converted to percentages and visualized on a log10 scale.
- Improved legend design
   Uses interpretable abundance levels:
    0.01, 0.10, 1.00, 10.00 (%)
- Non-symmetric color scale
   Better reflects the structure of metagenomic data (which are inherently skewed).
- Flexible top-N selection
   Users can define how many top taxa (families or genera) are included.
- Improved visualization aesthetics
   More consistent layout, annotation, and readability.


Current status:
The updated scripts represent the recommended workflow for generating publication-ready heatmaps.
Further refinements may be introduced in the future, but the current version is considered stable and complete.
