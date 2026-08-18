# 03 — IS655-like element characterization

Define boundaries, coding capacity, TIR/TSD, and ORF integrity of the IS655-like
element adjacent to seh in all 46 complete loci.

| Script | Purpose |
|--------|---------|
| `01_scan_TIR_TSD.py` | Extract 100-bp flanks; check inverted repeats and 2–14 bp direct repeats |
| `02_characterize_IS3_architecture.py` | tBLASTn vs ISfinder IS3-family proteins; coding span projections |
| `03_delimit_IS655_boundaries.py` | Two-anchor boundary delimitation using ORF1/ORF2 protein projections |
| `04_assess_ORF_integrity.py` | ORF1/ORF2 stops, 15-aa deletion, DDE→DDN, A7G motif |
| `05_independent_IS_verification.py` | Cross-check using IS-finding tools (verification) |

Key results:
- Element length: 1,211 bp (identical across 46 loci)
- TIR: 9 nt exact match; TSD: not detected
- ORF1: 2 internal stops (TAA¹⁶¹, TGA²⁷²)
- ORF2: 1 internal stop (TAA¹⁰⁶⁹), 15-aa deletion, E¹⁹²→N (DDE→DDN)
- No A7G frameshift motif; autonomous transposition not supported
