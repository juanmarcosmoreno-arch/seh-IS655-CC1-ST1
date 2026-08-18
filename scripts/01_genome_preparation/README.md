# 01 — Genome preparation

Read quality filtering, de novo assembly, consensus polishing, and annotation
of local Oxford Nanopore sequencing data.

| Script | Software | Purpose |
|--------|----------|---------|
| `01_filtlong.sh` | Filtlong 0.3.1 | Filter reads: min_length=1000, keep_percent=90, target_bases=500M |
| `02_flye_assembly.sh` | Flye 2.9.6 | De novo assembly from R10.4.1 SUP reads (--nano-hq) |
| `03_medaka_polish.sh` | Medaka 2.1.1 | Consensus polishing; model r1041_e82_400bps_sup_v5.2.0 |
| `04_bakta_annotation.sh` | Bakta 1.11.4 | Genome annotation; Bakta full database |

Scripts read barcode-to-sample mapping from `metadata/barcode_muestras.tsv`.
Public genomes (NCBI) were downloaded with NCBI Datasets CLI v18.29.1 and
annotated with the same Bakta version.

**Note:** This step applies only to local Uruguay isolates. Public genome
assemblies were used directly after download.
