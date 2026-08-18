# Analysis workflow

Step-by-step description of all major analyses. Each block includes purpose, inputs,
software, key parameters, outputs, and the manuscript result it supports.

---

## Step 01 — Genome preparation (local isolates)

**Purpose:** Process raw Oxford Nanopore reads into polished, annotated assemblies.

**Input:** Raw POD5/FAST5 reads per barcode from GridION run.

### 01a. Quality filtering

| Item | Detail |
|------|--------|
| Software | Filtlong 0.3.1 |
| Script | `scripts/01_genome_preparation/01_filtlong.sh` |
| Parameters | `--min_length 1000 --keep_percent 90 --target_bases 500000000` |
| Input | Raw FASTQ per barcode |
| Output | Filtered FASTQ per barcode |

### 01b. Assembly

| Item | Detail |
|------|--------|
| Software | Flye 2.9.6 |
| Script | `scripts/01_genome_preparation/02_flye_assembly.sh` |
| Parameters | `--nano-hq` (R10.4.1 SUP reads) |
| Output | `assembly.fasta` per isolate |

### 01c. Polishing

| Item | Detail |
|------|--------|
| Software | Medaka 2.1.1 |
| Script | `scripts/01_genome_preparation/03_medaka_polish.sh` |
| Model | `r1041_e82_400bps_sup_v5.2.0` |
| Output | Polished consensus FASTA per isolate |

### 01d. Annotation

| Item | Detail |
|------|--------|
| Software | Bakta 1.11.4 |
| Script | `scripts/01_genome_preparation/04_bakta_annotation.sh` |
| Database | Bakta full database (`databases/bakta_db`) |
| Output | GFF3, TSV, GBK per isolate |
| Result used | All downstream analyses |

---

## Step 02 — Typing and identity

**Purpose:** Confirm species (FastANI), assign MLST sequence type, compute pairwise genomic distances.

### 02a. MLST

| Item | Detail |
|------|--------|
| Software | mlst 2.33.1, BLAST+ 2.16.0 |
| Script | `scripts/02_typing_and_identity/01_mlst.R` |
| Scheme | MLST 7-locus *S. aureus* (pubMLST) |
| Output | `results/mlst/mlst_raw.tsv` |
| Result used | Genome inclusion criteria (ST1 or ST2187) |

### 02b. FastANI

| Item | Detail |
|------|--------|
| Software | FastANI 1.34 |
| Script | Part of `09_run_overnight_pipeline.sh` (upstream step) |
| Purpose | Species confirmation (≥95% ANI to *S. aureus* reference) |

### 02c. Mash pairwise distances

| Item | Detail |
|------|--------|
| Software | Mash 2.3 |
| Script | `scripts/02_typing_and_identity/03_mash_pairwise.py` |
| Input | 53 ST1 + ST2187 genomes (phylogenomics set) |
| Output | `results/figure_data/mash53_pairwise_distances.tsv` |
| Figure | Figure S3 (Mash heatmap) |

---

## Step 03 — IS655-like element characterization

**Purpose:** Define the boundaries, TIR, TSD, and coding capacity of the IS655-like element
adjacent to *seh* in each genome.

### 03a. TIR and TSD scan

| Item | Detail |
|------|--------|
| Script | `scripts/03_IS655_characterization/01_scan_TIR_TSD.py` |
| Method | Extract 100 bp flanks from IS3-family CDS coordinates; check inverted repeats; search 2–14 bp direct repeats flanking the element |
| Result | 9-nt exact TIR; no detectable TSD across 46 loci |
| Output | `results/summary_tables/IS655_TIR_TSD_scan.tsv` |

### 03b. IS3-family architecture

| Item | Detail |
|------|--------|
| Script | `scripts/03_IS655_characterization/02_characterize_IS3_architecture.py` |
| Method | tBLASTn of IS3-family transposase from ISfinder; reports coding span and element boundary projections |
| Result | IS655 is best hit for all 46 loci; Bakta CDS underpredicts element length |

### 03c. Element boundary delimitation

| Item | Detail |
|------|--------|
| Script | `scripts/03_IS655_characterization/03_delimit_IS655_boundaries.py` |
| Method | Two-anchor approach: ORF1 (94 aa) and ORF2 (288 aa) protein projections; ranks boundaries by TIR identity |
| Result | Element = 1,211 bp across all 46 complete loci |

### 03d. Coding capacity assessment

| Item | Detail |
|------|--------|
| Script | `scripts/03_IS655_characterization/04_assess_ORF_integrity.py` |
| Method | Align IS655 ORF1/ORF2 to each candidate; count in-frame stops and subject gaps; check DDE positions (aa 139, 157, 192 in ORF2); check A(n)G frameshift motif |
| Result | 2 internal stops in ORF1 (TAA¹⁶¹, TGA²⁷²); 1 stop in ORF2 (TAA¹⁰⁶⁹); 15-aa deletion in ORF2; DDE→DDN (E¹⁹²→N); no A₇G motif; autonomous transposition not supported |
| Output | `results/seh_st1_bovine/is3_coding_capacity/is655_like_coding_capacity.tsv` |
| Figure | Figure S2 |

---

## Step 04 — *seh* coding and regulatory integrity

**Purpose:** Verify that *seh* is intact and has recognizable regulatory elements in all genomes.

| Script | Purpose |
|--------|---------|
| `scripts/04_seh_integrity/01_seh_integrity.py` | Audit *seh* ORF (726 bp, 0 internal stops, nucleotide/AA alleles) |
| `scripts/04_seh_integrity/02_extract_upstream_region.py` | Extract 653-nt upstream region; find Shine-Dalgarno and σA-like −10/−35 motifs (computational, not experimental) |

**Results:** 46/46 intact ORFs; 2 nucleotide alleles, 2 protein alleles; 45/46 identical
to 23EV612 (GCA_001680935.1 differs at one position, one AA substitution).

**Output:** `results/summary_tables/seh_integrity.tsv`, `results/summary_tables/seh_regulatory_audit.tsv`

---

## Step 05 — CC1 distribution and near-fixation search

**Purpose:** Determine the prevalence of the *seh*–IS655-like module across CC1 and ST1.

### 05a. PRJNA1263585 CC1 panel screen

| Item | Detail |
|------|--------|
| Script | `scripts/05_CC1_distribution/01_PRJNA1263585_read_screen.py` |
| Method | k-mer screen of 75 Jahnen et al. 2025 Illumina read sets for diagnostic *seh* k-mers; streaming without storing reads |
| Result | 75/75 CC1/CC1660 read sets positive for *seh* k-mers |

### 05b. CC1 multi-ST survey

| Item | Detail |
|------|--------|
| Script | `scripts/05_CC1_distribution/02_module_screen_CC1.py` |
| Input | CC1 assemblies (Jahnen 2025 + historical ST772 + references) |
| Output | `results/summary_tables/CC1_module_distribution.tsv` |
| Figure | Figure 5 (CC1 phylogeny) |

### 05c. ST1 near-fixation search (NCBI)

| Item | Detail |
|------|--------|
| Script | `scripts/05_CC1_distribution/03_ST1_NCBI_search.py` |
| Tool | NCBI Datasets CLI v18.29.1 |
| Query | "Staphylococcus aureus ST1" |
| Retrieved | 500 assemblies |
| ST1 confirmed by MLST | 91/94 candidates |
| Result | 91/91 confirmed ST1 carried complete *seh*+IS655-like module (identity ≥99.9%) |

---

## Step 06 — CC1 core-genome phylogeny and ancestral reconstruction

**Purpose:** Build a midpoint-rooted CC1 phylogeny and infer the evolutionary history
of the *seh*–IS655-like module.

### 06a. Core-genome phylogeny

| Step | Software | Script | Parameters |
|------|----------|--------|------------|
| Alignment | SKA2 0.5.1 | `scripts/06_CC1_phylogeny/01_ska_panaroo_iqtree.sh` | |
| Recombination masking | Gubbins 3.4.3 | same | GTRGAMMA model |
| Phylogeny | IQ-TREE 3.1.2 | same | GTR+F+ASC, 1000 ultrafast bootstrap, 1000 SH-aLRT |
| Output tree | | | `results/trees/CC1_phylogeny.nwk` |

### 06b. Module-state mapping and ancestral reconstruction

| Item | Detail |
|------|--------|
| Script | `scripts/06_CC1_phylogeny/05_fitch_parsimony.py` |
| Method | Fitch parsimony; 5 partial loci coded as unknown (?); 1 ancestral gain, 0 losses |
| 5 partial loci coded as unknown | GCF_051474285.1, GCF_051474845.1, GCF_051475345.1, GCF_051475525.1, GCF_051475795.1 |
| Output | `results/summary_tables/ancestral_state_summary.tsv` |
| Figure | Figure 5 |
| Manuscript result | Single ancestral gain in ST1/ST81 ancestor; no parsimonious losses |

---

## Step 07 — SCCmec analysis

**Purpose:** Determine SCCmec type in ST1 MRSA isolates and verify coexistence with the module.

| Item | Detail |
|------|--------|
| Script | `scripts/07_SCCmec_analysis/01_sccmec_subtyping.py` |
| Database | Static BLAST DB, dated 2017-11-17 |
| Thresholds | ≥90% identity, ≥60% query coverage per HSP |
| Method | CCr + mec class combination; subtype only with subtype-specific markers |
| Output | `results/seh_st1_bovine/st1_mrsa_sccmec_formal/` |
| Manuscript | SCCmec IV predominant in ST1 MRSA; module coexists |

---

## Step 08 — ST2187 comparative genomics

**Purpose:** Identify all available ST2187 genomes, verify by MLST, and characterize the
*rlmH*–*cstA* interval in each.

| Script | Purpose |
|--------|---------|
| `01_build_locus_state_matrix.py` | Build presence/absence matrix for ST1 and ST2187 |
| `02_compare_anchor_intervals.py` | Compare *rlmH*–*cstA* intervals between ST1 and ST2187 |
| `03_formalize_comparison.py` | Formalize nucleotide identity levels at multiple regions |
| `04_map_breakpoints.py` | Map recombination breakpoints in cepa_18 |
| `05_characterize_ancestral_block.py` | Characterize the ancestral defense-island block in ST2187 |
| `06_expanded_ST2187_search.py` | Expanded search: 16 Souza 2023 + 2 Naushad 2020 + 5 Maine 2022 |

**24 verified ST2187 genomes:** see `results/summary_tables/ST2187_expanded_search.tsv`.

---

## Step 09 — Recombination analysis

**Purpose:** Formally test mosaic ancestry in cepa_18 (ST2187 genome with ST1-like locus).

Three-way alignment: cepa_18 (recipient), 23EV612 (ST1 donor), GCA_031889545.1 (ST2187 reference). Alignment ~47,937 bp.

### 09a. Ancestry-informative site analysis and Fisher test

| Item | Detail |
|------|--------|
| Script | `scripts/09_recombination/01_ancestry_sites_fisher_test.py` |
| Method | Classify sites inside/outside *rlmH*–*cstA* as ST1-like or ST2187-like; Fisher exact test |
| Result | OR = 1778.7; P < 10⁻³⁰⁰ |
| Output | `results/summary_tables/recombination_breakpoints.tsv` |

### 09b. 3SEQ

| Item | Detail |
|------|--------|
| Software | 3SEQ v1.9.1 beta build 260813 |
| Script | `scripts/09_recombination/02_3SEQ_run.py` |
| Method | All 6 oriented triplet comparisons; Dunn-Sidak multiple-testing correction |
| Result | cepa_18 as recombinant child: Dunn-Sidak P = 8.58×10⁻²⁶⁵ |
| Output | `results/summary_tables/3seq_results.tsv` |

### 09c. GENECONV

| Item | Detail |
|------|--------|
| Software | GENECONV v1.81a |
| Script | `scripts/09_recombination/03_GENECONV_parse.py` |
| Parameters | 10,000 permutations, Bonferroni correction |
| Result | 3 inner fragments; KA P = 9.93×10⁻⁹⁵; concordant breakpoints with 3SEQ |
| Output | `results/summary_tables/geneconv_results.tsv` |

**Note:** RDP5 was not used (Windows-only software; intentionally absent from this pipeline).

**Concordant breakpoints:** left ~15,460 nt, right ~20,891 nt (cepa_18 coords).
Recombinant tract ~5,432 bp; completely encloses *seh*+IS655-like (cepa_18: 17,424–19,432 nt).

---

## Step 10 — GC content and codon usage

**Purpose:** Quantify AT enrichment of the *seh*–IS655-like module relative to genome background.

| Item | Detail |
|------|--------|
| Scripts | `scripts/10_GC_codon_usage/01_GC_codon_analysis.py`, `02_GC3_window_test.py` |
| Reference | 23EV612 (complete ST1 reference) |
| Module interval | *rlmH*–*cstA* (Bakta contig_1:1880276-1888846, 8,571 bp) |
| *seh* GC3 | 15.70% → 3.4th percentile among 2,752 core CDS |
| Module GC | 29.18% vs. background mean 32.90% ± 0.80% |
| Window test | 10,000 random 8,571-bp windows (unique start positions; overlapping windows allowed); module more AT-rich than 97.0% (P = 0.030) |
| Outputs | `results/summary_tables/GC_codon_usage_summary.tsv`, `results/figure_data/GC_window_distribution.tsv` |

---

## Step 11 — Temporal signal audit

**Purpose:** Evaluate root-to-tip temporal signal before any molecular-clock analysis.

| Item | Detail |
|------|--------|
| Script | `scripts/11_temporal_signal/01_temporal_signal.py` |
| Tools | dendropy, scipy.stats.linregress, numpy |
| Input | CC1 phylogeny (Step 06) + collection dates from metadata |
| Method | Root-to-tip distance vs. collection year; 10,000 date randomizations |
| Output | `results/summary_tables/temporal_signal_summary.tsv`, `results/figure_data/root_to_tip_points.tsv` |
| Note | Analysis evaluates temporal signal only; no molecular clock rate estimation was performed |

---

## Step 12 — Figures and supplementary tables

| Figure | Script | Key data |
|--------|--------|---------|
| Figure 1 (module architecture) | `scripts/12_figures_tables/figure1_module_architecture.R` | Bakta annotation + rlmH–cstA coordinates |
| Figure 2 (CC1 distribution) | `scripts/12_figures_tables/figure2_CC1_distribution.R` | CC1 module presence/absence |
| Figure 3 (ST2187 interval) | `scripts/12_figures_tables/figure3_ST2187_interval.R` | Anchor interval comparison |
| Figure 4 (recombination) | `scripts/12_figures_tables/figure4_recombination.R` | Ancestry-informative sites |
| Figure 5 (CC1 phylogeny) | `scripts/12_figures_tables/figure5_CC1_phylogeny.R` | CC1 phylogeny + module states |
| Figure S1 (ST1/ST2187 phylogram) | original file; script not separately stored | Panaroo+IQ-TREE output |
| Figure S2 (IS655 coding capacity) | `scripts/12_figures_tables/figureS2_IS655_coding_capacity.py` | IS655-like TSV data |
| Figure S3 (Mash heatmap) | `scripts/12_figures_tables/figureS3_mash_heatmap.R` | mash53_pairwise_distances.tsv |
| Tables S1–S16 | `scripts/12_figures_tables/generate_supplementary_tables.py` | All results TSVs; output: `supplementary_tables_FEMSLE.xlsx` |
