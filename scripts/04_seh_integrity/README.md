# 04 — seh coding and regulatory integrity

Audit of seh CDS and immediately upstream regulatory region across all 46 loci.

| Script | Purpose |
|--------|---------|
| `01_seh_integrity.py` | CDS length (726 bp), internal stop codons, nucleotide/AA alleles |
| `02_extract_upstream_region.py` | Extract 653-nt upstream region; SD and sigma-A motif search |

Important: Promoter motif predictions (SD, -10/-35 sigma-A) are computational and
exploratory. They do not constitute experimental evidence of seh expression.

Results: 46/46 intact ORFs; 0 internal stops; 2 nucleotide alleles, 2 protein
alleles; 45/46 identical to 23EV612.
