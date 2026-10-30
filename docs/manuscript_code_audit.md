# Manuscript–code audit

Cross-reference of every analysis reported in the manuscript against scripts,
inputs, and outputs found in this project.

Classification:
- **COMPLETE** — script, inputs, and outputs all confirmed present
- **PARTIAL** — script present but inputs or outputs incomplete
- **MISSING** — no script found; analysis cannot be reproduced from files alone
- **RECONSTRUCTED** — script not original; reconstructed from manuscript Methods

---

| # | Manuscript analysis | Script(s) found | Input found | Output found | Status | Notes |
|---|---------------------|----------------|-------------|--------------|--------|-------|
| 1 | ONT assembly (Flye) | `01_genome_preparation/02_flye_assembly.sh` | Raw reads: TODO (SRA) | FASTA per isolate (local) | PARTIAL | Raw reads not yet deposited; assemblies available locally |
| 2 | Medaka polishing | `01_genome_preparation/03_medaka_polish.sh` | Filtered reads: TODO | Polished FASTA (local) | PARTIAL | Requires raw reads |
| 3 | Bakta annotation | `01_genome_preparation/04_bakta_annotation.sh` | Polished FASTA | GFF3/TSV (local) | PARTIAL | Bakta DB not in repo (large); annotated files available locally |
| 4 | MLST | `02_typing_and_identity/01_mlst.R` | Assembly FASTAs | mlst_raw.tsv | COMPLETE | mlst v2.33.1; scheme = S. aureus pubMLST |
| 5 | FastANI | Part of `09_run_overnight_pipeline.sh` | Assembly FASTAs | FastANI output (local) | PARTIAL | Script present; no dedicated wrapper in this repo |
| 6 | Mash pairwise distances | `02_typing_and_identity/03_mash_pairwise.py` | 53 genomes FASTA | `results/figure_data/mash53_pairwise_distances.tsv` | COMPLETE | |
| 7 | ISfinder tBLASTn | `03_IS655_characterization/02_characterize_IS3_architecture.py` | Bakta GFF + ISfinder DB | analysis_notes/is3_seh_first_pass_2026-08-11.md | COMPLETE | IS655 confirmed as best hit |
| 8 | IS655 boundaries | `03_IS655_characterization/03_delimit_IS655_boundaries.py` | Bakta GFF + IS655 proteins | boundary TSV (local) | COMPLETE | 1,211 bp; confirmed across 46 loci |
| 9 | TIR / TSD search | `03_IS655_characterization/01_scan_TIR_TSD.py` | Genome FASTA + Bakta coords | `results/summary_tables/IS655_TIR_TSD_scan.tsv` | COMPLETE | 9-nt TIR; no TSD |
| 10 | ORF integrity (IS655) | `03_IS655_characterization/04_assess_ORF_integrity.py` | Genome FASTA | `is3_coding_capacity/is655_like_coding_capacity.tsv` | COMPLETE | 2 ORF1 stops, 1 ORF2 stop, ∆15 aa, DDE→DDN |
| 11 | *seh* CDS integrity | `04_seh_integrity/01_seh_integrity.py` | Bakta TSV | `results/summary_tables/seh_integrity.tsv` | COMPLETE | 726 bp; 0 internal stops; 2 alleles |
| 12 | Regulatory-region analysis | `04_seh_integrity/01_seh_integrity.py` | Bakta TSV + genome FASTA | `results/summary_tables/seh_regulatory_audit.tsv` | COMPLETE | SD and σA motifs; computational only, no expression inference |
| 13 | CC1 survey (Jahnen panel) | `05_CC1_distribution/02_module_screen_CC1.py` | CC1 assemblies + FASTA queries | `results/summary_tables/CC1_module_distribution.tsv` | COMPLETE | |
| 14 | ST1 near-fixation search | `05_CC1_distribution/03_ST1_NCBI_search.py` | NCBI Datasets CLI | `cc1_st_survey/st1_candidates.tsv` | COMPLETE | 500 → 91 ST1; 91/91 seh+ |
| 15 | Jahnen 2025 read-screen | `05_CC1_distribution/01_PRJNA1263585_read_screen.py` | 75 FASTQ (ENA/SRA) | screen result TSV | COMPLETE | 75/75 positive |
| 16 | SKA2 core alignment | `06_CC1_phylogeny/01_ska_panaroo_iqtree.sh` | Genome FASTAs | `cc1_phylogeny/cc1_core_alignment.fasta` | COMPLETE | SKA2 v0.5.1 |
| 17 | Gubbins recombination masking | `06_CC1_phylogeny/01_ska_panaroo_iqtree.sh` | Core alignment | `cc1_phylogeny/gubbins/` (local) | COMPLETE | Gubbins v3.4.3 |
| 18 | IQ-TREE phylogeny | `06_CC1_phylogeny/01_ska_panaroo_iqtree.sh` | Filtered alignment | `results/trees/CC1_phylogeny.nwk` | COMPLETE | GTR+F+ASC; 1000 UFBoot; 1000 SH-aLRT |
| 19 | Fitch parsimony | `06_CC1_phylogeny/05_fitch_parsimony.py` | CC1 tree + module presence TSV | `results/summary_tables/ancestral_state_summary.tsv` | COMPLETE | 5 partials → unknown; 1 gain, 0 losses |
| 20 | SCCmecFinder | `07_SCCmec_analysis/01_sccmec_subtyping.py` | Assembly FASTAs + static BLAST DB | `st1_mrsa_sccmec_formal/` (local) | COMPLETE | DB dated 2017-11-17; ≥90% id, ≥60% cov |
| 21 | ST2187 identification | `08_ST2187_comparison/06_expanded_ST2187_search.py` | NCBI + local screens | `results/summary_tables/ST2187_expanded_search.tsv` | COMPLETE | 24 verified genomes |
| 22 | Anchor-interval comparison | `08_ST2187_comparison/02_compare_anchor_intervals.py` | Genome FASTAs | `st2187_comparator/formal_comparison_levels_identity.tsv` | COMPLETE | |
| 23 | Ancestry-informative sites | `09_recombination/01_ancestry_sites_fisher_test.py` | 3-way alignment | `st2187_breakpoints/ancestry_informative_sites.tsv` | COMPLETE | |
| 24 | Fisher exact test | `09_recombination/01_ancestry_sites_fisher_test.py` | Ancestry-site counts | `results/summary_tables/recombination_breakpoints.tsv` | COMPLETE | OR=1778.7; P<10⁻³⁰⁰ |
| 25 | 3SEQ | `09_recombination/02_3SEQ_run.py` | 3-way alignment FASTA | `results/summary_tables/3seq_results.tsv` | COMPLETE | 3SEQ v1.9.1 beta; P=8.58×10⁻²⁶⁵ |
| 26 | GENECONV | `09_recombination/03_GENECONV_parse.py` | 3-way alignment FASTA | `results/summary_tables/geneconv_results.tsv` | COMPLETE | v1.81a; 10,000 permutations |
| 27 | IS655 outside-CC1 screen | `scripts/132_retry_is655_outside_cc1.py` (original) | Local 170-genome panel | `is655_outside_cc1_retry/` | PARTIAL | Script not copied to this repo; results available locally |
| 28 | GC analysis | `10_GC_codon_usage/01_GC_codon_analysis.py` | 23EV612 Bakta TSV + FASTA | `results/summary_tables/GC_codon_usage_summary.tsv` | COMPLETE | |
| 29 | GC3 analysis | `10_GC_codon_usage/01_GC_codon_analysis.py` | Same | `results/figure_data/per_locus_GC_codon.tsv` | COMPLETE | *seh* GC3 at 3.4th percentile |
| 30 | 10,000 random windows | `10_GC_codon_usage/02_GC3_window_test.py` | 23EV612 FASTA | `results/figure_data/GC_window_distribution.tsv` | COMPLETE | Seed 2187; P=0.030 |
| 31 | Temporal signal | `11_temporal_signal/01_temporal_signal.py` | CC1 tree + dated metadata | `results/summary_tables/temporal_signal_summary.tsv` | COMPLETE | 10,000 date randomizations |
| 32 | Supplementary figures | See `12_figures_tables/` | Various TSVs | `final_paper_figures/` (local) | COMPLETE | Figures S1–S3 |
| 33 | Supplementary tables | `12_figures_tables/generate_supplementary_tables.py` | All results TSVs | `supplementary_tables_FEMSLE.xlsx` | COMPLETE | Tables S1–S16 |
| 34 | Fitch sensitivity (rooting × partial-tip coding) | `06_CC1_phylogeny/06_fitch_sensitivity.py` | CC1 tree + module presence TSV + outgroup tip list | `results/summary_tables/fitch_sensitivity_summary.tsv` | COMPLETE | Table S17. Added 2026-10-30 (reviewer-response script not originally copied into this repo; see note below) |
| 35 | Recombination robustness to negative-reference choice | `09_recombination/04_recombination_robustness.py` | 3-way alignment inputs for 3 alternative ST2187-negative references | `results/summary_tables/recombination_robustness_summary.tsv` | COMPLETE | Table S18. Added 2026-10-30 (reviewer-response script not originally copied into this repo; see note below) |

---

## Summary

| Status | Count |
|--------|-------|
| COMPLETE | 30 / 35 |
| PARTIAL | 5 / 35 |
| MISSING | 0 / 35 |
| RECONSTRUCTED | 0 / 35 |

---

## PARTIAL analyses — explanation

| # | Analysis | What is missing |
|---|----------|----------------|
| 1 | ONT assembly | Raw reads not yet deposited (SRA: TODO). Polished assemblies available locally and will be deposited. |
| 2 | Medaka polishing | Requires raw reads (same as above). |
| 3 | Bakta annotation | Bakta database (large; ~30 GB) not in repo. Annotated files (GFF3/TSV) available locally and via genome deposit. |
| 5 | FastANI | No dedicated wrapper script for this analysis in this repo; FastANI is called as part of the overnight pipeline. |
| 27 | IS655 outside-CC1 screen | Script `132_retry_is655_outside_cc1.py` not copied into this repo structure. Results TSV available locally. |

---

## Note on #34 / #35 (added 2026-10-30)

`142_fitch_sensitivity.py` and `143_recombination_robustness.py` (the scripts behind
Tables S17 and S18) were run during peer review but were not copied into this
repository at initial publication — the public repo and its own audit table did not
previously cover these two tables at all. Both scripts were located in the original
project directory, re-run, and verified to reproduce Tables S17 and S18 exactly
(byte-identical to the submitted supplementary Excel file, modulo line endings)
before being added here as `06_CC1_phylogeny/06_fitch_sensitivity.py` and
`09_recombination/04_recombination_robustness.py`.

---

## DISCREPANCIES

No discrepancies between manuscript values and reproduced values were detected.

One notable correction (made during peer review): The original analysis coded 5 partial
assemblies as "absent" (0), producing 5 spurious losses. Corrected analysis codes them as
unknown (?); result: 1 gain, 0 losses. The corrected version is what appears in the
submitted revision and in this repository.
