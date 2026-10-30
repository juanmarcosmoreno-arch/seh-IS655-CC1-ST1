# Genomic analysis of the seh–IS655-like module in *Staphylococcus aureus* CC1/ST1

This repository contains scripts, analysis commands, metadata and processed data supporting the manuscript:

> **"The seh enterotoxin gene resides in a conserved IS655-like genomic module in *Staphylococcus aureus* CC1/ST1"**
>
> submitted to *FEMS Microbiology Letters*.

---

## Study overview

This study characterizes a conserved genomic module in *S. aureus* CC1/ST1 bovine mastitis
isolates in which the *seh* enterotoxin gene is flanked by an IS655-like element
(IS3 family) in the *rlmH*–*cstA* chromosomal interval. The analyses address:

1. **Architecture and coding capacity of the IS655-like element** adjacent to *seh*:
   element boundaries, TIR/TSD, ORF1/ORF2 structure, internal stop codons, 15-aa
   deletion in ORF2, DDE→DDN catalytic substitution, absence of A₇G frameshift motif.
2. **Conservation of the *seh*–IS655-like module** across ST1 and related CC1 lineages:
   near-fixation in ST1 (91/91 public genomes), presence in CC1 panel (Jahnen et al. 2025).
3. **Coexistence of the module with SCCmec** in ST1 MRSA isolates.
4. **Comparison of the ST1-like locus with the alternative ST2187 defense-island locus**:
   24 verified ST2187 genomes screened; mosaic genome cepa_18 carries both locus types.
5. **Evidence for homologous recombination** in ST2187 mosaic genome:
   ancestry-informative sites (Fisher OR=1778.7), 3SEQ (P=8.58×10⁻²⁶⁵),
   GENECONV (KA P=9.93×10⁻⁹⁵).
6. **GC content and codon-usage characteristics**: *seh* GC3 at 3.4th percentile
   among 2,752 core CDS; module GC (29.18%) lower than 97% of 10,000 random windows
   (P=0.030).
7. **Temporal structure of the ST1 population**: root-to-tip regression, date-randomization
   test (10,000 permutations). No molecular clock dating performed.

---

## Repository organization

```
seh-IS655-CC1-ST1/
├── metadata/                  Genome metadata tables
│   ├── genome_metadata.tsv    All 48 genomes (47 ST1 + 1 ST2187) with assembly info
│   ├── st1_genomes.tsv        ST1 genomes only
│   ├── st2187_genomes.tsv     24 verified ST2187 genomes with locus screening results
│   └── public_accessions.tsv  Public GenBank/SRA accessions
│
├── scripts/
│   ├── 01_genome_preparation/ Filtlong, Flye, Medaka, Bakta
│   ├── 02_typing_and_identity/ MLST, FastANI, Mash
│   ├── 03_IS655_characterization/ TIR/TSD scan, ORF integrity, element boundaries
│   ├── 04_seh_integrity/      seh CDS audit, upstream regulatory region
│   ├── 05_CC1_distribution/   CC1 panel screening, ST1 near-fixation search
│   ├── 06_CC1_phylogeny/      SKA2, Gubbins, IQ-TREE, Fitch parsimony
│   ├── 07_SCCmec_analysis/    SCCmecFinder typing
│   ├── 08_ST2187_comparison/  ST2187 identification, anchor-interval comparison
│   ├── 09_recombination/      Ancestry-informative sites, Fisher test, 3SEQ, GENECONV
│   ├── 10_GC_codon_usage/     GC/GC3 analysis, sliding-window test
│   ├── 11_temporal_signal/    Root-to-tip regression, date-randomization audit
│   └── 12_figures_tables/     Figure generation (R/Python) and supplementary tables
│
├── config/
│   ├── software_versions.txt  Verified tool versions from conda-meta files
│   ├── parameters.tsv         Key analysis parameters with sources
│   └── environment_notes.md   Computational environment details
│
├── results/
│   ├── summary_tables/        Key TSV outputs supporting figures and tables
│   ├── trees/                 Phylogenetic tree files (Newick)
│   ├── alignments/            Core-genome alignment files
│   └── figure_data/           Data files used to generate figures
│
└── docs/
    ├── workflow.md             Step-by-step analysis pipeline
    ├── reproducibility_notes.md  Fully reproducible vs. reconstructed vs. missing
    ├── manuscript_code_audit.md  33-analysis reproducibility audit
    ├── release_checklist.md    Pre-publication checklist
    └── security_audit.md       Sensitive-information audit
```

---

## Workflow

```
Sequencing / public genomes
        ↓
Assembly and annotation
(Filtlong → Flye → Medaka → Bakta)
        ↓
MLST / FastANI / Mash
        ↓
IS655-like characterization
(boundary delimitation, TIR/TSD, ORF integrity, coding capacity)
        ↓
seh coding and upstream-region analysis
        ↓
CC1 distribution (CC1 panel screen + ST1 near-fixation search)
        ↓
CC1 phylogeny and ancestral reconstruction
(SKA2 → Gubbins → IQ-TREE → Fitch parsimony)
        ↓
SCCmec analysis
        ↓
ST2187 comparative genomics
(identification, verification, anchor-interval comparison)
        ↓
Recombination analyses
(ancestry-informative sites → Fisher → 3SEQ → GENECONV)
        ↓
GC / codon usage
        ↓
Temporal signal audit
        ↓
Figures and supplementary tables
```

---

## Data availability

### Local Uruguay isolates (8 genomes)

Raw Nanopore reads and assembled/annotated genomes for 7 local ST1 isolates
(cepa_17, cepa_21, cepa_30, cepa_52, cepa_98, cepa_107, cepa_202) and 1 ST2187
isolate (cepa_18) will be deposited at:

```
BioProject:           PRJNA1525248 (registered; embargoed until publication — not yet public at NCBI)
SRA accessions:       TODO
Genome assemblies:    TODO
Zenodo archive:       TODO (pending — see note below)
```

> **Note on the Zenodo archive (2026-10-30):** a DOI (`10.5281/zenodo.21066330`) was
> previously added here and in `CITATION.cff`, but it in fact archives a different,
> earlier GitHub repository (`seh-st1-mobile-element-reproducibility`), not this one.
> It has been removed pending a correct Zenodo archive of this repository.

### Public genomes

Public genomes are identified by accession numbers listed in `metadata/public_accessions.tsv`.
Key public datasets used:

- **23EV612** (ST1 bovine reference, New Zealand): GenBank CP160024, SAMN42208606, PRJNA1130542
- **Jahnen et al. 2025 CC1/CC1660 panel**: BioProject PRJNA1263585 (75 Illumina read sets)
- **Souza et al. 2023 ST2187**: BioProject PRJNA741582 (16 genomes)
- **Naushad et al. 2020 ST2187**: BioProject PRJNA599195 (2 verified ST2187 genomes)
- **Maine 2022 ST2187**: BioProject PRJNA974108 (5 genomes)
- **Australia ST1 bovine genomes**: see `metadata/public_accessions.tsv`

---

## Citation

If you use this repository or the data within it, please cite the associated manuscript
(see `CITATION.cff` for machine-readable citation metadata).

A Zenodo DOI badge will be added here once this repository has its own correct
Zenodo archive (see note in Data availability above).

---

## License

Scripts are released under the MIT License (see `LICENSE`).
