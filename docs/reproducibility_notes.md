# Reproducibility notes

---

## Fully reproducible

The following analyses can be fully reproduced from publicly available data (deposited
assemblies, public accessions) using the scripts in this repository:

- **MLST typing** (`02_typing_and_identity/01_mlst.R`): mlst v2.33.1; scheme files from pubMLST.
- **Mash pairwise distances** (`02_typing_and_identity/03_mash_pairwise.py`): from deposited assemblies.
- **IS655-like element characterization** (`03_IS655_characterization/`): all scripts + local FASTA inputs.
- **seh coding and regulatory integrity** (`04_seh_integrity/`): from deposited assemblies + Bakta annotation.
- **CC1 survey module screen** (`05_CC1_distribution/02_module_screen_CC1.py`): from CC1 assemblies + query sequences in `results/seh_st1_bovine/cc1_st_survey/`.
- **Jahnen 2025 read-screen** (`05_CC1_distribution/01_PRJNA1263585_read_screen.py`): public FASTQ from PRJNA1263585.
- **CC1 phylogeny** (`06_CC1_phylogeny/`): from deposited assemblies; genome FASTA → SKA2 → Gubbins → IQ-TREE.
- **Fitch parsimony** (`06_CC1_phylogeny/05_fitch_parsimony.py`): from CC1 tree + `CC1_module_distribution.tsv`.
- **SCCmec subtyping** (`07_SCCmec_analysis/01_sccmec_subtyping.py`): from deposited assemblies + static BLAST DB.
- **ST2187 expanded search** (`08_ST2187_comparison/06_expanded_ST2187_search.py`): from deposited assemblies + public NCBI downloads.
- **Recombination analyses** (`09_recombination/`): from three-way alignment FASTA (`flanking_alignment.fasta`).
- **GC and codon usage** (`10_GC_codon_usage/`): from 23EV612 Bakta TSV + genome FASTA.
- **Temporal signal** (`11_temporal_signal/`): from CC1 tree + dated metadata.
- **All supplementary tables** (`12_figures_tables/generate_supplementary_tables.py`): from results TSVs.
- **Figures S2 and S3** (`12_figures_tables/figureS2_IS655_coding_capacity.py`, `figureS3_mash_heatmap.R`): from results TSVs.

---

## Reconstructed from manuscript

No scripts were reconstructed from scratch. All scripts in this repository are original
analysis scripts from the project.

---

## Missing information

The following items are not yet available in this repository and require manual action
before public release:

### Critical (blocking GitHub publication)

```
TODO: BioProject accession for local Uruguay isolates
TODO: SRA accessions for 7 Uruguay ST1 isolates (cepa_17/21/30/52/98/107/202)
TODO: SRA accession for Uruguay ST2187 isolate (cepa_18)
TODO: GenBank/ENA genome assembly accessions for all 8 local isolates
TODO: Final author list for LICENSE, CITATION.cff, and README
TODO: GitHub repository URL
TODO: Zenodo DOI (after GitHub release)
```

### Non-blocking

```
TODO: Exact R package versions for ggtree, ape, ggplot2 (≥4.0.0 confirmed), dplyr
TODO: Exact versions of numpy, scipy, dendropy (Python; gridion_mlst env)
TODO: Verify Bakta DB version used (database files at databases/bakta_db/)
TODO: BioSample accessions for Naushad 2020 ST2187 isolates
      (GCA_013404735.1 and GCA_013406015.1; metadata from publication filled in)
TODO: Script for IS655-outside-CC1 screen (132_retry_is655_outside_cc1.py)
      — add to 05_CC1_distribution/ if included in final repo
TODO: Figure S1 generation script — the phylogram exists as a PDF but the
      dedicated R script (79_plot_rooted_core_tree_ggtree.R) is the upstream
      candidate; verify and copy if needed
```

### Notes on local genome paths

All local file paths (e.g., `/Volumes/JM/BioInf/GridION/...`) appear in the original
scripts and have been preserved. Before GitHub publication, these should be parameterized
as relative paths or environment variables. The scripts in this repo still contain
local absolute paths — this is intentional for audit purposes; users running the scripts
must set paths to their own data locations.

### RDP5

RDP5 was not used in this analysis (Windows-only software). Files named `rdp5_results.tsv`
in `st2187_recombination_test/` are empty or contain placeholder entries. This is
documented intentionally.
