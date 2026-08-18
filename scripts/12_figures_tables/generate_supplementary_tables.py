#!/usr/bin/env python3
"""Compile all 16 supplementary tables into a single Excel workbook.

One sheet per table (Table S1 … Table S16).
Each sheet: title row, description row, blank row, then the data with bold headers.
"""
from __future__ import annotations
import csv
from pathlib import Path

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

ROOT = Path("/Volumes/JM/BioInf/GridION")
RES  = ROOT / "results/seh_st1_bovine"
OUT  = RES  / "supplementary_tables_FEMSLE.xlsx"

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------
BLUE_FILL   = PatternFill("solid", fgColor="1565C0")
GRAY_FILL   = PatternFill("solid", fgColor="ECEFF1")
TITLE_FONT  = Font(name="Arial", bold=True, size=10, color="FFFFFF")
HEADER_FONT = Font(name="Arial", bold=True, size=9)
DATA_FONT   = Font(name="Arial", size=9)
DESC_FONT   = Font(name="Arial", italic=True, size=9, color="546E7A")
WRAP        = Alignment(wrap_text=True, vertical="top")
TOP         = Alignment(vertical="top")

def set_col_width(ws) -> None:
    """Auto-fit column widths (capped at 50)."""
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            try:
                v = str(cell.value or "")
                max_len = max(max_len, min(len(v), 60))
            except Exception:
                pass
        ws.column_dimensions[col_letter].width = min(max_len + 2, 52)

def write_table(ws, title: str, description: str, rows: list[list]) -> None:
    """Write title, description, blank row, header, data rows."""
    # Row 1: title
    ws.append([title])
    title_cell = ws.cell(1, 1)
    title_cell.font  = Font(name="Arial", bold=True, size=11)

    # Row 2: description
    ws.append([description])
    ws.cell(2, 1).font      = DESC_FONT
    ws.cell(2, 1).alignment = Alignment(wrap_text=True)
    ws.row_dimensions[2].height = 40

    # Row 3: blank
    ws.append([])

    if not rows:
        ws.append(["[No data]"])
        return

    headers = rows[0]
    n_cols  = len(headers)

    # Row 4: header
    ws.append(headers)
    hrow = ws.max_row
    for c in range(1, n_cols + 1):
        cell = ws.cell(hrow, c)
        cell.fill      = BLUE_FILL
        cell.font      = Font(name="Arial", bold=True, size=9, color="FFFFFF")
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[hrow].height = 28

    # Data rows
    for row in rows[1:]:
        ws.append(row)
        r = ws.max_row
        fill = GRAY_FILL if (r - hrow) % 2 == 0 else None
        for c in range(1, n_cols + 1):
            cell = ws.cell(r, c)
            cell.font      = DATA_FONT
            cell.alignment = TOP
            if fill:
                cell.fill = fill

    ws.freeze_panes = ws.cell(hrow + 1, 1)
    set_col_width(ws)


def read_tsv(path: Path) -> list[list]:
    if not path.exists():
        return [["FILE NOT FOUND", str(path)]]
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        return [row for row in reader]


# ---------------------------------------------------------------------------
# Table definitions  (sheet_name, title, description, tsv_path or list[tsv_path])
# ---------------------------------------------------------------------------
TABLES = [
    (
        "Table S1",
        "Table S1. Genome metadata and assembly characteristics",
        "Genome metadata, assembly technology and sequencing coverage for all genomes analyzed. "
        "Includes complete accession list, assembly source, ST, host, country, year, and "
        "near-fixation search reproducibility details (query tool: NCBI Datasets CLI v18.29.1; "
        "500 broad records retrieved; 94 ST1 candidates downloaded; 91/94 ST1 confirmed by MLST; "
        "91/91 seh-positive).",
        [RES / "final_tables/supplementary_table_genomes_loci_tree.tsv",
         RES / "final_tables/TableS1_assembly_technology_evidence.tsv"],
    ),
    (
        "Table S2",
        "Table S2. IS655-like element nucleotide identity vs 23EV612",
        "IS655-like element nucleotide identity vs 23EV612 per genome (pident, alignment length, "
        "query coverage).",
        [RES / "final_tables/final_spacer_seh_identity_vs_23EV612.tsv"],
    ),
    (
        "Table S3",
        "Table S3. seh coding integrity per genome",
        "seh coding integrity per genome (726 bp CDS, internal stop codons, allele classification, "
        "amino-acid variants relative to 23EV612).",
        [RES / "seh_integrity/seh_coding_integrity.tsv"],
    ),
    (
        "Table S4",
        "Table S4. Upstream regulatory audit",
        "Upstream regulatory audit (653-nt interval): Shine-Dalgarno motif, sigma-A-like −10/−35 "
        "promoter pair, per genome.",
        [RES / "seh_integrity/seh_upstream_regulatory_audit.tsv"],
    ),
    (
        "Table S5",
        "Table S5. ISfinder tBLASTn best hits for IS655-like element",
        "ISfinder tBLASTn best hits for the IS655-like element (IS3 family): subject IS name, "
        "pident, alignment length, E-value, IS family assignment.",
        [RES / "is3_coding_capacity/is655_orf_tblastn_alignments.tsv"],
    ),
    (
        "Table S6",
        "Table S6. CC1 ST survey — module presence/absence",
        "CC1 ST survey: seh+IS655-like module presence/absence per genome (Jahnen 2025 "
        "PRJNA1263585 panel + ST772 historical + references; 99 genomes total). "
        "Includes module_status coding: Present / Partial (unknown) / Absent.",
        [RES / "cc1_st_survey/cc1_st_survey_seh_is655_presence.tsv"],
    ),
    (
        "Table S7",
        "Table S7. ST1 MRSA SCCmec subtyping",
        "ST1 MRSA SCCmec subtyping details per genome: SCCmecFinder best hit, type, identity, "
        "coverage, coexistence with seh+IS655-like module.",
        [RES / "st1_mrsa_sccmec_formal/st1_mrsa_sccmec_subtyping.tsv"],
    ),
    (
        "Table S8",
        "Table S8. ST2187 anchor-interval comparison and recombination breakpoints",
        "ST2187 anchor-interval comparison: identity levels at multiple genomic regions (core, "
        "module, flanking) vs ST1 reference 23EV612 and ST2187-negative background. "
        "Includes recombination breakpoint mapping.",
        [RES / "st2187_comparator/formal_comparison_levels_identity.tsv"],
    ),
    (
        "Table S9",
        "Table S9. Fisher mosaic test — ancestry-informative sites",
        "Fisher mosaic test data: counts of ancestry-informative sites inside and outside the "
        "rlmH–cstA interval, odds ratio and P-value for enrichment of ST1-like sites within "
        "the recombinant tract.",
        [RES / "st2187_breakpoints/ancestry_informative_sites.tsv"],
    ),
    (
        "Table S10",
        "Table S10. IS655-like distribution — 170-genome local panel",
        "IS655-like element distribution in the 170-genome local panel: assembly ID, ST, "
        "panel source, IS655-like screen result.",
        [RES / "is655_outside_cc1_retry/local_genome_panel_inventory.tsv"],
    ),
    (
        "Table S11",
        "Table S11. Expanded ST2187 search results",
        "Expanded ST2187 search results: genome metadata, defense-block coverage and identity, "
        "and seh/IS655-like screening for all 24 verified ST2187 genomes. "
        "Cohorts: Souza_2023_multihost (16 genomes, PRJNA741582, northern New England USA), "
        "Naushad_2020_Canada (2 genomes, PRJNA599195, Canada), "
        "Maine_2022 (5 genomes, PRJNA974108, USA: Maine), "
        "Uruguay_local_positive (cepa_18, local sequencing). "
        "Columns: seh_present/identity/coverage — BLAST vs 23EV612 seh (726 bp); "
        "is655_complete/identity/coverage — BLAST vs IS655-like element (1211 bp); "
        "defense_block_complete/identity/coverage — BLAST vs ancestral rlmH–cstA interval.",
        [RES / "st2187_expanded_panel/st2187_full_table_s11.tsv"],
    ),
    (
        "Table S12",
        "Table S12. 3SEQ recombination test results",
        "3SEQ recombination test results: triplet comparisons, corrected Dunn-Sidak probabilities, "
        "recombinant segment length, and breakpoint intervals for the cepa18 ST2187 mosaic genome.",
        [RES / "st2187_recombination_test/3seq_results.tsv"],
    ),
    (
        "Table S13",
        "Table S13. Temporal-signal audit",
        "Temporal-signal audit: root-to-tip regression statistics (slope, R², permutation P), "
        "date-randomization probabilities (10,000 permutations), and exploratory root-date "
        "bootstrap 95% CI (10,000 replicates). Tool: custom Python (dendropy + scipy.stats + numpy; "
        "scripts/135_audit_temporal_signal_st1.py). Tree: Gubbins recombination-filtered final tree.",
        [RES / "temporal_signal_audit/temporal_signal_summary.tsv",
         RES / "temporal_signal_audit/root_to_tip_points.tsv"],
    ),
    (
        "Table S14",
        "Table S14. GENECONV recombination test results",
        "GENECONV recombination test results: inner fragments for the cepa18;23EV612 pair, "
        "fragment coordinates, lengths, polymorphic site counts, and corrected KA P-values. "
        "Tool: GENECONV v1.81a.",
        [RES / "st2187_recombination_test/geneconv_results.tsv"],
    ),
    (
        "Table S15",
        "Table S15. CC1 ancestral state reconstruction",
        "CC1 ancestral state reconstruction: inferred node states, gain/loss events, and "
        "bootstrap support values. Fitch parsimony with five partial loci coded as unknown (?) "
        "— reconstructed one ancestral gain and zero losses. Script: 141_fitch_partials_unknown.py.",
        [RES / "cc1_phylogeny/cc1_ancestral_state_reconstruction_v2.tsv"],
    ),
    (
        "Table S16",
        "Table S16. GC content and codon usage analysis",
        "GC content and codon usage: per-locus GC, GC1, GC2, GC3 for all 23EV612 CDS. "
        "Includes GC3 percentile ranks for seh (3.4th percentile; lower than 96.6% of core CDS) "
        "and IS655-like transposase (27.7th percentile), and the 10,000 random 8,571-bp window "
        "comparison (module more AT-rich than 97.0% of windows; empirical P = 0.030). "
        "Script: 139_gc_codon_analysis.py + 142_gc3_percentile_window.py.",
        [RES / "gc_codon_analysis/per_locus_gc_codon.tsv",
         RES / "gc_codon_analysis/gc3_percentile_summary.tsv",
         RES / "gc_codon_analysis/module_window_gc.tsv"],
    ),
]


# ---------------------------------------------------------------------------
# Build workbook
# ---------------------------------------------------------------------------
wb = openpyxl.Workbook()
wb.remove(wb.active)  # remove default sheet

for sheet_name, title, description, paths in TABLES:
    ws = wb.create_sheet(title=sheet_name)

    # For tables with multiple TSV sources, concatenate with a separator row
    if isinstance(paths, list) and len(paths) > 1:
        combined: list[list] = []
        for i, p in enumerate(paths):
            rows = read_tsv(p)
            if not rows:
                continue
            if combined:
                # separator
                combined.append([f"— {p.name} —"])
                combined.extend(rows)  # include header of each sub-table
            else:
                combined.extend(rows)
        write_table(ws, title, description, combined)
    else:
        p = paths[0] if isinstance(paths, list) else paths
        write_table(ws, title, description, read_tsv(p))

    # Merge title cell across all columns
    n_cols = ws.max_column or 1
    if n_cols > 1:
        ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=n_cols)
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=n_cols)

wb.save(OUT)
print(f"Saved: {OUT}")
print(f"Sheets: {[s.title for s in wb.worksheets]}")
