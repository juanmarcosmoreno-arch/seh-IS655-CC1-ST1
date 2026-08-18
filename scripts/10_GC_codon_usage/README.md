# 10 — GC content and codon usage

GC analysis of the seh+IS655-like module relative to the genome background.

| Script | Purpose |
|--------|---------|
| `01_GC_codon_analysis.py` | Per-CDS GC, GC1, GC2, GC3; seh and IS655-like transposase GC3 |
| `02_GC3_window_test.py` | Random 8,571-bp window test; module GC vs background |

Reference: 23EV612 (Bakta contig_1:1880276-1888846 = rlmH-cstA module, 8,571 bp)

Window test parameters:
- 10,000 windows; window size = 8,571 bp; SEED = 2187
- Sampling: random.sample (unique start positions; overlapping windows allowed)
  NOTE: manuscript Methods says "non-overlapping" — this is incorrect;
  the actual implementation uses overlapping windows.

Results:
- seh GC3: 15.70% → 3.4th percentile (< 96.6% of 2,752 core CDS)
- Module GC: 29.18% vs background 32.90% ± 0.80%
- Module more AT-rich than 97.0% of windows (P = 0.030)
