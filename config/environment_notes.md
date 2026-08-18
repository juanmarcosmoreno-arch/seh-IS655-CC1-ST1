# Computational environment notes

## Platform

- Hardware: Oxford Nanopore GridION (sequencing); Apple Silicon Mac (analysis)
- OS: macOS (darwin, arm64)
- Shell: zsh
- Conda-compatible environment manager: micromamba (prefix `.micromamba/`)

## Conda environments

Each analysis stage uses a dedicated micromamba environment to avoid dependency conflicts.

| Environment | Primary tools |
|-------------|--------------|
| `gridion_nanopore` | Filtlong 0.3.1, Flye 2.9.6, FastANI 1.34, Mash 2.3 |
| `gridion_medaka` | Medaka 2.1.1, Python 3.9 |
| `gridion_bakta` | Bakta 1.11.4, BLAST+ 2.16.0 |
| `gridion_mlst` | mlst 2.33.1, BLAST+ 2.16.0, Python 3.10, dendropy, scipy, numpy, openpyxl |
| `gridion_recombination` | SKA2 0.5.1, Gubbins 3.4.3, IQ-TREE 2.4.0, R 4.5.3 |
| `gridion_pangenome_pip` | Panaroo 1.6.1, IQ-TREE 3.1.2, Python 3.10 |

## External tools (not conda-managed)

- 3SEQ v1.9.1 beta build 260813: binary at `tools/3seq/3seq`
- GENECONV v1.81a: binary at `tools/geneconv/geneconv`

## R packages

System R 4.5.1 was used for ggtree-based figure generation. Key R packages:
- `ggtree`, `ape`, `ggplot2` 4.0.0, `dplyr`, `patchwork`

Note: ggplot2 4.0.0 introduced breaking changes affecting `gheatmap`; figure scripts
were updated to use `scale_color_manual` / `scale_fill_manual` with hollow shapes instead.

## Reproducibility notes for local genomes

Raw Nanopore reads for local Uruguay isolates are deposited as:
- BioProject: TODO
- SRA accessions: TODO

The assembled and annotated genomes (Flye + Medaka + Bakta) are deposited as:
- Genome assembly accessions: TODO

Without raw reads, the bioinformatic analysis starting from Step 01 cannot be fully
reproduced for local genomes. Steps 03 onward (IS655 characterization and beyond)
can be reproduced from the deposited assemblies.
