#!/usr/bin/env python3
"""Project and evaluate boundaries of the IS655-like element adjacent to seh.

Two independent IS655 protein anchors are used:

* ORF1: 94 aa, coordinates 80-364 in the 1,221 bp IS655 reference
* ORF2: 288 aa, coordinates 325-1191 in the reference

Candidate boundaries are searched near the projected ends, must contain the
entire translated-homology span, and are ranked using inverted-repeat identity,
repeat length, exact flanking target-site duplication, and distance from the
protein-anchor projection. Short AT-rich repeats are explicitly retained as
provisional rather than treated as definitive element boundaries.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCH = ROOT / "results/seh_st1_bovine/is3_architecture"
DEFAULT_OUTDIR = ROOT / "results/seh_st1_bovine/is3_boundaries"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
IS655_LENGTH = 1221
ORF_COORDS = {1: (80, 364), 2: (325, 1191)}


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line)
    return {key: "".join(value).upper() for key, value in records.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in rows:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def select_anchors(hits: list[dict[str, str]]) -> dict[str, dict[int, dict[str, str]]]:
    anchors: dict[str, dict[int, dict[str, str]]] = {}
    for hit in hits:
        if hit["reference_element"] != "IS655":
            continue
        orf = int(hit["reference_orf_index"])
        if orf not in ORF_COORDS or float(hit["query_coverage"]) < 0.80 or hit["hit_strand"] != "+":
            continue
        anchors.setdefault(hit["sample_id"], {})[orf] = hit
    return anchors


def project_bounds(sample_anchors: dict[int, dict[str, str]]) -> tuple[int, int, str]:
    if 1 not in sample_anchors or 2 not in sample_anchors:
        raise ValueError("Both IS655 ORF anchors are required")
    left_estimates: list[int] = []
    right_estimates: list[int] = []
    for orf, hit in sample_anchors.items():
        ref_start, ref_end = ORF_COORDS[orf]
        subject_start = min(int(hit["subject_start"]), int(hit["subject_end"]))
        subject_end = max(int(hit["subject_start"]), int(hit["subject_end"]))
        if int(hit["query_start"]) == 1:
            left_estimates.append(subject_start - (ref_start - 1))
        if int(hit["query_end"]) >= int(hit["query_length_aa"]) - 1:
            right_estimates.append(subject_end + (IS655_LENGTH - ref_end))
    # ORF1 is close to the left end and ORF2 to the right end. Use the closest
    # anchor at each end so an internal indel is not incorrectly propagated
    # across the entire reference element. Keep all estimates as diagnostics.
    left = min(left_estimates)
    right = min(right_estimates)
    evidence = f"left={','.join(map(str,left_estimates))};right={','.join(map(str,right_estimates))}"
    return left, right, evidence


def best_tsd(seq: str, left: int, right: int) -> tuple[int, str]:
    for length in range(14, 1, -1):
        upstream = seq[left - 1 - length : left - 1]
        downstream = seq[right : right + length]
        if len(upstream) == length and upstream == downstream:
            return length, upstream
    return 0, ""


def boundary_candidates(
    seq: str,
    projected_left: int,
    projected_right: int,
    homology_left: int,
    homology_right: int,
    radius: int,
) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    for left in range(max(1, projected_left - radius), projected_left + radius + 1):
        if left > homology_left:
            continue
        for right in range(projected_right - radius, min(len(seq), projected_right + radius) + 1):
            if right < homology_right or right <= left:
                continue
            for repeat_len in range(9, 41):
                left_repeat = seq[left - 1 : left - 1 + repeat_len]
                right_repeat = rc(seq[right - repeat_len : right])
                if len(left_repeat) != repeat_len or len(right_repeat) != repeat_len:
                    continue
                mismatches = sum(a != b for a, b in zip(left_repeat, right_repeat))
                identity = (repeat_len - mismatches) / repeat_len
                if identity < 0.75:
                    continue
                tsd_len, tsd = best_tsd(seq, left, right)
                deviation = abs(left - projected_left) + abs(right - projected_right)
                score = identity * 100 + repeat_len * 0.25 + tsd_len * 2 - deviation * 0.10
                candidates.append(
                    {
                        "candidate_left": left,
                        "candidate_right": right,
                        "candidate_length_bp": right - left + 1,
                        "tir_length_bp": repeat_len,
                        "tir_identity": f"{identity:.4f}",
                        "tir_mismatches": mismatches,
                        "left_terminal_sequence": left_repeat,
                        "right_terminal_reverse_complement": right_repeat,
                        "tsd_length_bp": tsd_len,
                        "tsd_sequence": tsd,
                        "total_boundary_deviation_bp": deviation,
                        "candidate_score": f"{score:.4f}",
                    }
                )
    return sorted(candidates, key=lambda row: float(row["candidate_score"]), reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--architecture-dir", type=Path, default=DEFAULT_ARCH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--radius", type=int, default=60)
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    windows = read_fasta(args.architecture_dir / "seh_adjacent_is3_oriented_windows.fasta")
    summary = {row["sample_id"]: row for row in read_tsv(args.architecture_dir / "is3_architecture_summary.tsv")}
    anchors = select_anchors(read_tsv(args.architecture_dir / "selected_is3_tblastn_best_hits.tsv"))
    ranked_rows: list[dict[str, object]] = []
    best_rows: list[dict[str, object]] = []
    element_records: list[tuple[str, str]] = []

    for sample_id, seq in windows.items():
        projected_left, projected_right, projection_evidence = project_bounds(anchors[sample_id])
        homology_left = int(summary[sample_id]["translated_homology_start"])
        homology_right = int(summary[sample_id]["translated_homology_end"])
        candidates = boundary_candidates(seq, projected_left, projected_right, homology_left, homology_right, args.radius)
        for rank, candidate in enumerate(candidates[: args.top], 1):
            ranked_rows.append(
                {
                    "sample_id": sample_id,
                    "rank": rank,
                    "projected_left": projected_left,
                    "projected_right": projected_right,
                    "projected_length_bp": projected_right - projected_left + 1,
                    "projection_evidence": projection_evidence,
                    "translated_homology_left": homology_left,
                    "translated_homology_right": homology_right,
                    **candidate,
                }
            )
        if not candidates:
            raise RuntimeError(f"No boundary candidate found for {sample_id}")
        best = candidates[0]
        quality = "provisional_short_TIR_no_TSD"
        if int(best["tir_length_bp"]) >= 15 and float(best["tir_identity"]) >= 0.80 and int(best["tsd_length_bp"]) >= 2:
            quality = "moderate_candidate_TIR_plus_TSD"
        elif int(best["tir_length_bp"]) >= 15 and float(best["tir_identity"]) >= 0.80:
            quality = "provisional_TIR_without_TSD"
        best_row = {
            "sample_id": sample_id,
            "projected_left": projected_left,
            "projected_right": projected_right,
            "projected_length_bp": projected_right - projected_left + 1,
            "projection_evidence": projection_evidence,
            "translated_homology_left": homology_left,
            "translated_homology_right": homology_right,
            **best,
            "boundary_quality": quality,
            "interpretation": "Protein-anchor projection is strong; exact nucleotide boundaries remain provisional",
        }
        best_rows.append(best_row)
        left, right = int(best["candidate_left"]), int(best["candidate_right"])
        element_records.append((f"{sample_id}|candidate_IS655_like|{left}-{right}|length={right-left+1}", seq[left - 1 : right]))

    write_tsv(args.outdir / "is3_boundary_candidates_ranked.tsv", ranked_rows)
    write_tsv(args.outdir / "is3_boundary_best_candidate.tsv", best_rows)
    write_fasta(args.outdir / "is3_candidate_elements.fasta", element_records)
    element_sequences = [seq for _, seq in element_records]
    sequence_counts = Counter(element_sequences)
    element_by_sample = {name.split("|")[0]: seq for name, seq in element_records}
    reference_sample = "23EV612"
    reference_seq = element_by_sample[reference_sample]
    comparison_rows: list[dict[str, object]] = []
    for sample_id, seq in element_by_sample.items():
        differences = [i + 1 for i, (a, b) in enumerate(zip(reference_seq, seq)) if a != b]
        comparison_rows.append(
            {
                "sample_id": sample_id,
                "candidate_length_bp": len(seq),
                "differences_vs_23EV612": len(differences),
                "variant_positions_vs_23EV612": ",".join(map(str, differences)),
                "identity_vs_23EV612": f"{100 * (len(seq) - len(differences)) / len(seq):.4f}",
                "identical_to_23EV612": seq == reference_seq,
            }
        )
    write_tsv(args.outdir / "is3_candidate_element_identity_vs_23EV612.tsv", comparison_rows)
    overview = [
        {"metric": "loci_analyzed", "value": len(best_rows)},
        {"metric": "unique_projected_boundaries", "value": len({(r["projected_left"], r["projected_right"]) for r in best_rows})},
        {"metric": "unique_best_candidate_boundaries", "value": len({(r["candidate_left"], r["candidate_right"]) for r in best_rows})},
        {"metric": "unique_candidate_element_sequences", "value": len(sequence_counts)},
        {"metric": "most_common_candidate_sequence_count", "value": sequence_counts.most_common(1)[0][1]},
        {"metric": "boundary_quality_counts", "value": ";".join(f"{k}:{v}" for k, v in Counter(r["boundary_quality"] for r in best_rows).items())},
    ]
    write_tsv(args.outdir / "is3_boundary_overview.tsv", overview)
    print(f"Evaluated projected boundaries for {len(best_rows)} loci; wrote {args.outdir}")


if __name__ == "__main__":
    main()
