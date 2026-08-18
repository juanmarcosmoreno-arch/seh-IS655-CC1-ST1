# 02 — Typing and identity

Species confirmation, MLST sequence typing, and pairwise genomic distances.

| Script | Software | Purpose |
|--------|----------|---------|
| `01_mlst.R` | mlst 2.33.1, BLAST+ 2.16.0 | 7-locus S. aureus MLST (pubMLST scheme) |
| `03_mash_pairwise.py` | Mash 2.3 | Pairwise Mash distances for 53 ST1/ST2187 genomes |
| `04_screen_public_st1.py` | BLAST+ | Screen public FASTA files for ST1 candidates |

FastANI 1.34 (species confirmation ≥95% ANI) was part of the overnight pipeline
and is not separately wrapped here. FastANI output tables are in the local results.
