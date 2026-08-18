# 12 — Figures and supplementary tables

Scripts to generate all manuscript figures and supplementary tables.

| Script | Output |
|--------|--------|
| `figure1_module_architecture.R` | Figure 1: seh-IS655-like module in rlmH-cstA interval |
| `figure2_CC1_distribution.R` | Figure 2: CC1 ST distribution of the module |
| `figure3_ST2187_interval.R` | Figure 3: ST2187 rlmH-cstA locus synteny |
| `figure4_recombination.R` | Figure 4: Mosaic ancestry sites in cepa_18 |
| `figure5_CC1_phylogeny.R` | Figure 5: CC1 phylogeny with module presence/absence |
| `figureS2_IS655_coding_capacity.py` | Figure S2: IS655-like ORF structure vs IS655 reference |
| `figureS3_mash_heatmap.R` | Figure S3: Mash distance heatmap (53 genomes) |
| `generate_supplementary_tables.py` | Tables S1-S16 → supplementary_tables_FEMSLE.xlsx |

Notes:
- Figures use R / ggtree (ggplot2 ≥ 4.0.0); gheatmap incompatible with ggplot2 4.0.0
  — replaced with scale_color_manual + hollow shapes.
- Figure S1 (ST1/ST2187 phylogram) was produced from IQ-TREE output (Panaroo core tree)
  using 79_plot_rooted_core_tree_ggtree.R (upstream of this section); see workflow.md.
