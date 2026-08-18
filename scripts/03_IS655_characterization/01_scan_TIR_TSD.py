#!/usr/bin/env python3
"""Scan seh-locus IS3-family transposase boundaries for TIR/TSD signals.

This script uses Bakta CDS boundaries for the nearest IS3-family transposase in
the final locus table. These coordinates mark the transposase CDS, not
necessarily the complete insertion-sequence element, so negative TIR/TSD calls
should be interpreted as "not detected at the annotated CDS boundaries".

For each representative genome, the script:

* extracts 100 bp flanks around the annotated IS3 CDS;
* orients the sequence relative to the IS3 CDS strand;
* compares terminal IS3 CDS sequences as inverted repeats (candidate TIRs);
* searches for exact 2-14 bp direct repeats immediately flanking the CDS
  (candidate TSDs);
* writes a TSV table and an annotated FASTA of left flank | IS3 CDS | right
  flank in oriented coordinates.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FINAL_TABLE = ROOT / "results/seh_st1_bovine/final_tables/final_genome_locus_tree_table.tsv"
DEFAULT_OUT_TSV = ROOT / "results/seh_st1_bovine/final_tables/is3_tir_tsd_boundary_scan.tsv"
DEFAULT_OUT_FASTA = ROOT / "results/seh_st1_bovine/final_tables/is3_tir_tsd_boundary_sequences.fasta"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

FIELDS = [
    "genome_id",
    "group",
    "tree_phylo_group",
    "IS_locus_tag",
    "contig",
    "is3_start",
    "is3_end",
    "is3_strand",
    "oriented_left_flank_start",
    "oriented_left_flank_end",
    "oriented_right_flank_start",
    "oriented_right_flank_end",
    "flank_length_bp",
    "best_terminal_TIR_length_bp",
    "best_terminal_TIR_identity",
    "best_terminal_TIR_mismatches",
    "left_terminal_sequence",
    "right_terminal_sequence",
    "right_terminal_reverse_complement",
    "TIR_call",
    "TSD_length_bp",
    "TSD_sequence",
    "left_flanking_TSD_sequence",
    "right_flanking_TSD_sequence",
    "TSD_call",
    "left_flank_100bp_oriented",
    "right_flank_100bp_oriented",
    "interpretation",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def select_panel(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    sample_ids = [
        "23EV612",
        "cepa_17_barcode01",
        "cepa_21_barcode03",
        "cepa_30_barcode06",
        "B1822_S7_L001",
        "B916_S394",
        "B1180_S15",
        "GCF_051990345.1",
        "GCF_051990445.1",
        "GCF_051990565.1",
        "GCA_001680935.1",
        "GCA_001681155.1",
        "GCA_024360105.1",
        "GCA_030734435.1",
        "GCA_029415515.1",
        "cepa_18_barcode02",
    ]
    by_id = {row["sample_id"]: row for row in rows if row.get("nearest_is3_locus_tag")}
    return [by_id[sample_id] for sample_id in sample_ids if sample_id in by_id]


def best_terminal_tir(element: str, min_len: int, max_len: int) -> dict[str, str]:
    best: tuple[float, int, int, str, str, str] | None = None
    max_len = min(max_len, len(element))
    for length in range(min_len, max_len + 1):
        left = element[:length]
        right = element[-length:]
        right_rc = rc(right)
        mismatches = sum(a != b for a, b in zip(left, right_rc))
        identity = (length - mismatches) / length
        score = (identity, length, -mismatches, left, right, right_rc)
        if best is None or score > best:
            best = score
    if best is None:
        return {
            "length": "NA",
            "identity": "NA",
            "mismatches": "NA",
            "left": "",
            "right": "",
            "right_rc": "",
            "call": "not_tested",
        }
    identity, length, neg_mismatches, left, right, right_rc = best
    mismatches = -neg_mismatches
    call = "candidate_TIR" if length >= 10 and identity >= 0.80 else "not_detected"
    return {
        "length": str(length),
        "identity": f"{identity:.3f}",
        "mismatches": str(mismatches),
        "left": left,
        "right": right,
        "right_rc": right_rc,
        "call": call,
    }


def find_tsd(left_flank: str, right_flank: str, min_len: int, max_len: int) -> dict[str, str]:
    for length in range(max_len, min_len - 1, -1):
        left = left_flank[-length:]
        right = right_flank[:length]
        if len(left) == length and left == right:
            return {
                "length": str(length),
                "sequence": left,
                "left": left,
                "right": right,
                "call": "candidate_TSD",
            }
    return {
        "length": "0",
        "sequence": "",
        "left": left_flank[-max_len:] if left_flank else "",
        "right": right_flank[:max_len] if right_flank else "",
        "call": "not_detected",
    }


def fasta_path(row: dict[str, str], root: Path) -> Path:
    return (root / row["bakta_tsv"]).with_suffix(".fna")


def oriented_context(row: dict[str, str], root: Path, flank: int) -> dict[str, str]:
    seqs = read_fasta(fasta_path(row, root))
    contig = row["seh_contig"]
    seq = seqs[contig]
    start = int(row["nearest_is3_start"])
    end = int(row["nearest_is3_end"])
    strand = row["nearest_is3_strand"]
    lo, hi = min(start, end), max(start, end)

    left_start = max(1, lo - flank)
    left_end = lo - 1
    right_start = hi + 1
    right_end = min(len(seq), hi + flank)

    left = seq[left_start - 1 : left_end].upper()
    element = seq[lo - 1 : hi].upper()
    right = seq[right_start - 1 : right_end].upper()

    if strand == "-":
        return {
            "left": rc(right),
            "element": rc(element),
            "right": rc(left),
            "left_start": str(right_end),
            "left_end": str(right_start),
            "right_start": str(left_end),
            "right_end": str(left_start),
        }
    return {
        "left": left,
        "element": element,
        "right": right,
        "left_start": str(left_start),
        "left_end": str(left_end),
        "right_start": str(right_start),
        "right_end": str(right_end),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final-table", type=Path, default=DEFAULT_FINAL_TABLE)
    parser.add_argument("--out-tsv", type=Path, default=DEFAULT_OUT_TSV)
    parser.add_argument("--out-fasta", type=Path, default=DEFAULT_OUT_FASTA)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--tir-min", type=int, default=10)
    parser.add_argument("--tir-max", type=int, default=40)
    parser.add_argument("--tsd-min", type=int, default=2)
    parser.add_argument("--tsd-max", type=int, default=14)
    args = parser.parse_args()

    rows = select_panel(read_tsv(args.final_table))
    out_rows: list[dict[str, str]] = []
    fasta_records: list[tuple[str, str]] = []

    for row in rows:
        ctx = oriented_context(row, args.root, args.flank)
        tir = best_terminal_tir(ctx["element"], args.tir_min, args.tir_max)
        tsd = find_tsd(ctx["left"], ctx["right"], args.tsd_min, args.tsd_max)
        interpretation = (
            "No clear TIR/TSD at annotated IS3 CDS boundaries"
            if tir["call"] == "not_detected" and tsd["call"] == "not_detected"
            else "Candidate boundary signal detected; inspect sequence context"
        )
        out_rows.append(
            {
                "genome_id": row["sample_id"],
                "group": row["analysis_group"],
                "tree_phylo_group": row["tree_phylo_group"],
                "IS_locus_tag": row["nearest_is3_locus_tag"],
                "contig": row["seh_contig"],
                "is3_start": row["nearest_is3_start"],
                "is3_end": row["nearest_is3_end"],
                "is3_strand": row["nearest_is3_strand"],
                "oriented_left_flank_start": ctx["left_start"],
                "oriented_left_flank_end": ctx["left_end"],
                "oriented_right_flank_start": ctx["right_start"],
                "oriented_right_flank_end": ctx["right_end"],
                "flank_length_bp": str(args.flank),
                "best_terminal_TIR_length_bp": tir["length"],
                "best_terminal_TIR_identity": tir["identity"],
                "best_terminal_TIR_mismatches": tir["mismatches"],
                "left_terminal_sequence": tir["left"],
                "right_terminal_sequence": tir["right"],
                "right_terminal_reverse_complement": tir["right_rc"],
                "TIR_call": tir["call"],
                "TSD_length_bp": tsd["length"],
                "TSD_sequence": tsd["sequence"],
                "left_flanking_TSD_sequence": tsd["left"],
                "right_flanking_TSD_sequence": tsd["right"],
                "TSD_call": tsd["call"],
                "left_flank_100bp_oriented": ctx["left"],
                "right_flank_100bp_oriented": ctx["right"],
                "interpretation": interpretation,
            }
        )
        header = (
            f"{row['sample_id']}|{row['nearest_is3_locus_tag']}|"
            f"{row['tree_phylo_group']}|left100_lower_is3_upper_right100_lower"
        )
        fasta_records.append((header, ctx["left"].lower() + ctx["element"].upper() + ctx["right"].lower()))

    write_tsv(args.out_tsv, out_rows)
    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.out_fasta.open("w", encoding="utf-8") as handle:
        for header, seq in fasta_records:
            handle.write(f">{header}\n{wrap(seq)}\n")

    print(f"Wrote {len(out_rows)} rows to {args.out_tsv}")
    print(f"Wrote annotated boundary FASTA to {args.out_fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
