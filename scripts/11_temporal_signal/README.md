# 11 — Temporal signal audit

Root-to-tip temporal signal evaluation of the ST1 CC1 phylogeny.

| Script | Software | Purpose |
|--------|----------|---------|
| `01_temporal_signal.py` | dendropy, scipy, numpy | Root-to-tip regression; 10,000 date randomizations |

This analysis evaluates whether the ST1 phylogeny shows sufficient temporal
structure for molecular clock analyses. It does NOT estimate substitution rates
or dates of divergence. The conclusion of the manuscript is that the signal
is present but exploratory root-date estimates are not interpreted as precise
dates of module acquisition.
