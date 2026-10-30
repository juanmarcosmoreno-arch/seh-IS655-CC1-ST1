# 09 — Recombination analysis

Formal tests for mosaic ancestry in cepa_18 (ST2187 genome with ST1-like rlmH-cstA locus).

Three-way alignment: cepa_18 (recipient) × 23EV612 (ST1 donor) × GCA_031889545.1 (ST2187 ref)
Alignment: ~47,937 bp (flanking_alignment.fasta in results/st2187_recombination_test/)

| Script | Software | Purpose |
|--------|----------|---------|
| `01_ancestry_sites_fisher_test.py` | scipy.stats | Ancestry-informative site classification + Fisher exact test |
| `02_3SEQ_run.py` | 3SEQ v1.9.1 beta | Mosaic test; all 6 oriented triplets; Dunn-Sidak correction |
| `03_GENECONV_parse.py` | GENECONV v1.81a | Parse GENECONV output; Bonferroni correction; 10,000 permutations |
| `04_recombination_robustness.py` | MAFFT, 3SEQ, GENECONV | Robustness of the recombination call to the choice of ST2187-negative reference (Table S18) |

Key results:
- Fisher: OR = 1778.7; P < 10^-300
- 3SEQ: Dunn-Sidak P = 8.58 × 10^-265
- GENECONV: 3 fragments; KA P = 9.93 × 10^-95
- Concordant breakpoints: left ~15,460 nt, right ~20,891 nt (cepa_18 coords)
- Recombinant tract ~5,432 bp, fully enclosing seh+IS655-like (17,424-19,432 nt)

RDP5 was NOT used (Windows-only; intentionally absent).

Robustness check (`04_recombination_robustness.py`, Table S18): repeats the same
3-way-alignment → 3SEQ → GENECONV pipeline with 3 alternative ST2187-negative
references (from Table S11, locus_state=ancestral_defense_block) in place of
GCA_031889545.1. P-values stay in the same order of magnitude and breakpoints/
recombinant length are unchanged (3SEQ) or differ by <1 bp (GENECONV) across all
4 references. GCA_031889225.1 was excluded: rlmH and cstA fall on two different
contigs in that assembly, so the flanking-window method cannot anchor it.
Output: `results/summary_tables/recombination_robustness_summary.tsv`.
