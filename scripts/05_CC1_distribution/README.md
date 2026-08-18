# 05 — CC1 distribution and ST1 near-fixation search

Module prevalence across CC1 panel and public ST1 genomes.

| Script | Purpose |
|--------|---------|
| `01_PRJNA1263585_read_screen.py` | k-mer screen of 75 Jahnen et al. 2025 CC1 Illumina read sets |
| `02_module_screen_CC1.py` | Survey seh+IS655-like module across CC1/ST multi-ST panel |
| `03_ST1_NCBI_search.py` | NCBI search + MLST verification of ST1 near-fixation |

Near-fixation search (NCBI Datasets CLI v18.29.1):
- 500 assemblies retrieved
- 94 ST1 candidates downloaded
- 91/94 confirmed ST1 by MLST
- 91/91 positive for complete seh+IS655-like module (identity ≥99.9%)
