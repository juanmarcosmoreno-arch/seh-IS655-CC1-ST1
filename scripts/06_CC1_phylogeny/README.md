# 06 — CC1 core-genome phylogeny and ancestral state reconstruction

Pipeline: SKA2 → Gubbins → IQ-TREE → module-state mapping → Fitch parsimony

| Script | Software | Purpose |
|--------|----------|---------|
| `01_ska_panaroo_iqtree.sh` | SKA2 0.5.1, Gubbins 3.4.3, IQ-TREE 3.1.2 | Core-genome alignment, recombination masking, phylogeny |
| `02_sanitize_alignment.py` | | Fix SKA2 / dendropy header format issues |
| `03_resolve_CC1_phylogeny.py` | dendropy | Midpoint rooting and tree processing |
| `04_map_module_states.R` | R / ggtree | Map seh+IS655-like states onto phylogeny |
| `05_fitch_parsimony.py` | dendropy | Fitch parsimony reconstruction |
| `06_fitch_sensitivity.py` | dendropy | Sensitivity of the Fitch reconstruction to partial-tip coding and rooting (Table S17) |

IQ-TREE parameters: GTR+F+ASC, 1000 ultrafast bootstrap, 1000 SH-aLRT

Fitch parsimony: 5 partial assemblies coded as unknown (?).
Result: 1 ancestral gain, 0 losses.

Sensitivity analysis (`06_fitch_sensitivity.py`, Table S17): re-runs the reconstruction
under all 6 combinations of partial-tip coding (unknown/absent/present) × rooting
(midpoint/ST1660 outgroup). The single ancestral gain is identical in every case; the
loss count changes only under "absent" coding (5 spurious tip-level losses), which is
why partial tips are coded as unknown in the main analysis (`05_fitch_parsimony.py`).
Output: `results/summary_tables/fitch_sensitivity_summary.tsv`.
