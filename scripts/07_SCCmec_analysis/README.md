# 07 — SCCmec analysis

SCCmec typing of ST1 MRSA isolates using a static BLAST database.

| Script | Purpose |
|--------|---------|
| `01_sccmec_subtyping.py` | BLAST-based SCCmec typing; ccr + mec class |
| `02_seh_sccmec_coexistence.py` | Test for seh module coexistence with SCCmec in ST1 MRSA |

SCCmecFinder static databases (2017-11-17):
- single_genes_database_20171117.fasta (thresholds: ≥90% id, ≥60% query coverage)
- whole_cassette_SCCmec_database_REFERENCE_20171117.fasta (supporting evidence only)

Note: SCCmecFinder was not used as a versioned web service; the static database
files were used directly with BLAST+ 2.16.0.
