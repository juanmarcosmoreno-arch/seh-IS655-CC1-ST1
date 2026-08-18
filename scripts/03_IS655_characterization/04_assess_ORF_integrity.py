#!/usr/bin/env python3
"""Assess autonomous coding capacity of the IS655-like candidate element.

The analysis aligns IS655 ORF1 and ORF2 proteins to each candidate nucleotide
sequence, counts in-frame stops/deletions, checks the reference DDE catalytic
positions (ORF2 residues 139, 157 and 192), searches the ORF overlap
for a canonical A(n)G frameshift motif, and reports the longest conventional
ORF. Calls are predictions about autonomous coding capacity, not experimental
transposition assays.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BOUNDARY_DIR = ROOT / "results/seh_st1_bovine/is3_boundaries"
DEFAULT_QUERY = ROOT / "results/seh_st1_bovine/is3_architecture/selected_is3_reference_proteins.fasta"
DEFAULT_OUTDIR = ROOT / "results/seh_st1_bovine/is3_coding_capacity"
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}
CATALYTIC_QUERY_POSITIONS = (139, 157, 192)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            current = raw[1:].split("|")[0]
            records[current] = []
        elif current is not None:
            records[current].append(raw.strip())
    return {name: "".join(parts).upper() for name, parts in records.items()}


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_tblastn(query: Path, subject: Path, output: Path) -> None:
    fields = "qseqid qlen sseqid sstart send pident length qstart qend evalue bitscore qseq sseq"
    subprocess.run(
        ["tblastn", "-query", str(query), "-subject", str(subject), "-evalue", "1e-5", "-seg", "no", "-outfmt", f"6 {fields}", "-out", str(output)],
        check=True,
    )


def parse_best_alignments(path: Path) -> dict[str, dict[int, dict[str, str]]]:
    fields = ["query_id", "query_length_aa", "subject_id", "subject_start", "subject_end", "pident", "alignment_length_aa", "query_start", "query_end", "evalue", "bitscore", "query_alignment", "subject_alignment"]
    best: dict[tuple[str, int], dict[str, str]] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        row = dict(zip(fields, raw.split("\t")))
        match = re.match(r"IS655_ORF([12])_", row["query_id"])
        if not match:
            continue
        sample = row["subject_id"].split("|")[0]
        orf = int(match.group(1))
        key = (sample, orf)
        if key not in best or float(row["bitscore"]) > float(best[key]["bitscore"]):
            best[key] = row
    output: dict[str, dict[int, dict[str, str]]] = {}
    for (sample, orf), row in best.items():
        output.setdefault(sample, {})[orf] = row
    return output


def aligned_query_map(row: dict[str, str]) -> dict[int, str]:
    query_pos = int(row["query_start"])
    mapping: dict[int, str] = {}
    for query_aa, subject_aa in zip(row["query_alignment"], row["subject_alignment"]):
        if query_aa != "-":
            mapping[query_pos] = subject_aa
            query_pos += 1
    return mapping


def aligned_stop_sites(row: dict[str, str], seq: str) -> tuple[str, str]:
    subject_pos = int(row["subject_start"])
    positions: list[str] = []
    codons: list[str] = []
    for subject_aa in row["subject_alignment"]:
        if subject_aa == "*":
            positions.append(str(subject_pos))
            codons.append(seq[subject_pos - 1 : subject_pos + 2])
        if subject_aa != "-":
            subject_pos += 3
    return ",".join(positions), ",".join(codons)


def longest_orf(seq: str) -> tuple[int, int, int, str, str]:
    best = (0, 0, 0, "", "")
    for frame in range(3):
        possible_starts: list[tuple[int, str]] = []
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i : i + 3]
            if codon in STARTS:
                possible_starts.append((i, codon))
            if codon in STOPS:
                for start, start_codon in possible_starts:
                    aa_length = (i + 3 - start) // 3 - 1
                    candidate = (aa_length, start + 1, i + 3, start_codon, codon)
                    if candidate > best:
                        best = candidate
                possible_starts = []
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundary-dir", type=Path, default=DEFAULT_BOUNDARY_DIR)
    parser.add_argument("--query", type=Path, default=DEFAULT_QUERY)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    subject = args.boundary_dir / "is3_candidate_elements.fasta"
    sequences = read_fasta(subject)
    raw = args.outdir / "is655_orf_tblastn_alignments.tsv"
    run_tblastn(args.query, subject, raw)
    alignments = parse_best_alignments(raw)

    rows: list[dict[str, object]] = []
    for sample, seq in sequences.items():
        orf1, orf2 = alignments[sample][1], alignments[sample][2]
        overlap_start = max(min(int(orf1["subject_start"]), int(orf1["subject_end"])), min(int(orf2["subject_start"]), int(orf2["subject_end"])))
        overlap_end = min(max(int(orf1["subject_start"]), int(orf1["subject_end"])), max(int(orf2["subject_start"]), int(orf2["subject_end"])))
        overlap_seq = seq[overlap_start - 1 : overlap_end]
        frameshift_matches = [match.group(0) for match in re.finditer(r"A{4,}G", overlap_seq)]
        catalytic_map = aligned_query_map(orf2)
        catalytic_state = "".join(catalytic_map.get(position, "?") for position in CATALYTIC_QUERY_POSITIONS)
        aa_len, start, end, start_codon, stop_codon = longest_orf(seq)
        orf1_stops = orf1["subject_alignment"].count("*")
        orf2_stops = orf2["subject_alignment"].count("*")
        orf1_stop_positions, orf1_stop_codons = aligned_stop_sites(orf1, seq)
        orf2_stop_positions, orf2_stop_codons = aligned_stop_sites(orf2, seq)
        intact_canonical_orfab = orf1_stops == 0 and orf2_stops == 0 and catalytic_state == "DDE" and bool(frameshift_matches)
        rows.append(
            {
                "sample_id": sample,
                "element_length_bp": len(seq),
                "orf1_reference_coverage": f"{(int(orf1['query_end'])-int(orf1['query_start'])+1)/int(orf1['query_length_aa']):.4f}",
                "orf1_identity": orf1["pident"],
                "orf1_internal_stops_in_alignment": orf1_stops,
                "orf1_stop_positions_in_candidate": orf1_stop_positions,
                "orf1_stop_codons": orf1_stop_codons,
                "orf1_subject_gaps": orf1["subject_alignment"].count("-"),
                "orf2_reference_coverage": f"{(int(orf2['query_end'])-int(orf2['query_start'])+1)/int(orf2['query_length_aa']):.4f}",
                "orf2_identity": orf2["pident"],
                "orf2_internal_stops_in_alignment": orf2_stops,
                "orf2_stop_positions_in_candidate": orf2_stop_positions,
                "orf2_stop_codons": orf2_stop_codons,
                "orf2_subject_gaps": orf2["subject_alignment"].count("-"),
                "orf_overlap_start": overlap_start,
                "orf_overlap_end": overlap_end,
                "orf_overlap_length_bp": max(0, overlap_end - overlap_start + 1),
                "canonical_AnG_frameshift_motif": ",".join(frameshift_matches),
                "canonical_AnG_frameshift_motif_present": bool(frameshift_matches),
                "catalytic_reference_positions": "D139,D157,E192",
                "catalytic_residues_in_candidate": catalytic_state,
                "catalytic_DDE_intact": catalytic_state == "DDE",
                "longest_conventional_orf_start": start,
                "longest_conventional_orf_end": end,
                "longest_conventional_orf_length_aa": aa_len,
                "longest_conventional_orf_start_codon": start_codon,
                "longest_conventional_orf_stop_codon": stop_codon,
                "intact_canonical_orfab_call": intact_canonical_orfab,
                "predicted_autonomous_coding_capacity": "not_supported_disrupted_orfA_orfB_and_DDE" if not intact_canonical_orfab else "canonical_orfAB_architecture_retained",
                "interpretation": "Autonomous transposition is unlikely; mobilization in trans is not tested" if not intact_canonical_orfab else "Coding architecture compatible with autonomous transposition; activity not demonstrated",
            }
        )
    write_tsv(args.outdir / "is655_like_coding_capacity.tsv", rows)
    overview = [
        {"metric": "loci_analyzed", "value": len(rows)},
        {"metric": "orf1_with_internal_stops", "value": sum(int(row["orf1_internal_stops_in_alignment"]) > 0 for row in rows)},
        {"metric": "orf2_with_internal_stops", "value": sum(int(row["orf2_internal_stops_in_alignment"]) > 0 for row in rows)},
        {"metric": "canonical_frameshift_motif_present", "value": sum(row["canonical_AnG_frameshift_motif_present"] is True for row in rows)},
        {"metric": "catalytic_DDE_intact", "value": sum(row["catalytic_DDE_intact"] is True for row in rows)},
        {"metric": "intact_canonical_orfab", "value": sum(row["intact_canonical_orfab_call"] is True for row in rows)},
        {"metric": "shared_catalytic_state", "value": ",".join(sorted({str(row["catalytic_residues_in_candidate"]) for row in rows}))},
        {"metric": "shared_longest_conventional_orf_length_aa", "value": ",".join(sorted({str(row["longest_conventional_orf_length_aa"]) for row in rows}))},
    ]
    write_tsv(args.outdir / "is655_like_coding_capacity_overview.tsv", overview)
    print(f"Assessed coding capacity for {len(rows)} loci; wrote {args.outdir}")


if __name__ == "__main__":
    main()
