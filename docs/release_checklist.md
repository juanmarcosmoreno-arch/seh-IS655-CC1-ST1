# Release checklist

Steps required before creating the public GitHub release.

## Accessions and identifiers

- [ ] Deposit raw reads for 8 local Uruguay isolates (NCBI SRA / ENA)
- [ ] Obtain SRA accessions; update `metadata/genome_metadata.tsv`
- [ ] Deposit assembled/annotated genomes for 8 local isolates
- [ ] Obtain GenBank/ENA assembly accessions; update all relevant metadata files
- [ ] Create BioProject; add accession to README.md Data Availability section
- [ ] Verify all public genome accessions in `metadata/public_accessions.tsv` are correct
- [ ] Add GitHub repository URL to README.md, CITATION.cff
- [ ] Create Zenodo archive after GitHub release v1.0.0
- [ ] Add Zenodo DOI to README.md and CITATION.cff
- [ ] Add data availability statement with BioProject, SRA, and Zenodo DOI to manuscript

## Table S11

- [ ] Verify Table S11 contains all 24 ST2187 genomes (current: 24 rows confirmed)
- [ ] Add BioSample IDs for Naushad 2020 genomes (GCA_013404735.1, GCA_013406015.1)

## Authorship and citation

- [ ] Add final complete author list to LICENSE (replace "TODO")
- [ ] Add final author list, ORCIDs, and date-released to CITATION.cff
- [ ] Add ORCID URLs and affiliation notes to README.md if desired

## Code quality

- [ ] Remove all local absolute paths from scripts (parameterize inputs)
      Key pattern: `/Volumes/JM/BioInf/GridION/` → use `$PROJECT_ROOT` or `../..`
- [ ] Remove any sensitive information identified in security audit
- [ ] Verify all scripts run from a clean environment (test at minimum Steps 09–12)
- [ ] Add brief usage note at top of each README in script subdirectories
- [ ] Verify software_versions.txt R package versions (run `sessionInfo()`)
- [ ] Verify software_versions.txt Python package versions (run `pip freeze`)

## Reproducibility

- [ ] Run `scripts/12_figures_tables/generate_supplementary_tables.py` from clean env
- [ ] Verify supplementary Excel Tables S1–S16 match manuscript
- [ ] Run `scripts/12_figures_tables/figureS2_IS655_coding_capacity.py`
- [ ] Run `scripts/12_figures_tables/figureS3_mash_heatmap.R`
- [ ] Confirm `results/trees/CC1_phylogeny.nwk` is the final midpoint-rooted tree

## Security

- [ ] Run security audit (see `docs/security_audit.md`) on all files to be committed
- [ ] Confirm no local usernames, farm/producer names, or internal passwords are present
- [ ] Confirm no NCBI credentials or SSH keys are present

## Git and GitHub

- [ ] `git init` in `seh-IS655-CC1-ST1/`
- [ ] `git add` only tracked files (verify with `git status`; nothing in `.gitignore`)
- [ ] Create initial commit; verify commit message
- [ ] Create private GitHub repository `seh-IS655-CC1-ST1`
- [ ] Push; make public after manuscript acceptance
- [ ] Create GitHub release v1.0.0 tagged at final commit
- [ ] Archive release in Zenodo (auto-sync or manual upload)

## Manuscript

- [ ] Add Data Availability section citing BioProject, Zenodo DOI, and GitHub URL
- [ ] Confirm Methods references to accession numbers match deposited ones
- [ ] Correct "non-overlapping" → "randomly sampled" in GC window Methods text
      (actual implementation: unique start positions, overlapping windows allowed)
