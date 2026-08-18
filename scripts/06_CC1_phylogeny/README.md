# 06 — CC1 core-genome phylogeny and ancestral state reconstruction

Pipeline: SKA2 → Gubbins → IQ-TREE → module-state mapping → Fitch parsimony

| Script | Software | Purpose |
|--------|----------|---------|
| `01_ska_panaroo_iqtree.sh` | SKA2 0.5.1, Gubbins 3.4.3, IQ-TREE 3.1.2 | Core-genome alignment, recombination masking, phylogeny |
| `02_sanitize_alignment.py` | | Fix SKA2 / dendropy header format issues |
| `03_resolve_CC1_phylogeny.py` | dendropy | Midpoint rooting and tree processing |
| `04_map_module_states.R` | R / ggtree | Map seh+IS655-like states onto phylogeny |
| `05_fitch_parsimony.py` | dendropy | Fitch parsimony reconstruction |

IQ-TREE parameters: GTR+F+ASC, 1000 ultrafast bootstrap, 1000 SH-aLRT

Fitch parsimony: 5 partial assemblies coded as unknown (?).
Result: 1 ancestral gain, 0 losses.
