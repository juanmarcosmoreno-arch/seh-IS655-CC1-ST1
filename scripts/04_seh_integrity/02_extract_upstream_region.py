#!/usr/bin/env python3
"""Extract seh flanking windows and the nearest IS3-seh spacer sequences."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
SUMMARY_TSV = PROJECT_DIR / "results" / "seh_st1_bovine" / "seh_locus_summary.tsv"
FEATURES_TSV = PROJECT_DIR / "results" / "seh_st1_bovine" / "seh_locus_features.tsv"
OUT_DIR = PROJECT_DIR / "results" / "seh_st1_bovine" / "flanks"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def bakta_fasta_from_tsv(raw_path: str) -> Path:
    tsv = PROJECT_DIR / raw_path
    return tsv.with_suffix(".fna")


def coordinate_gap(a: dict[str, str], b: dict[str, str]) -> int:
    a_start, a_end = int(a["feature_start"]), int(a["feature_end"])
    b_start, b_end = int(b["feature_start"]), int(b["feature_end"])
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def inter_feature_bounds(a: dict[str, str], b: dict[str, str]) -> tuple[int, int]:
    a_start, a_end = int(a["feature_start"]), int(a["feature_end"])
    b_start, b_end = int(b["feature_start"]), int(b["feature_end"])
    if a_end < b_start:
        return a_end + 1, b_start - 1
    if b_end < a_start:
        return b_end + 1, a_start - 1
    return 1, 0


def signed_gap(seh: dict[str, str], is3: dict[str, str]) -> int:
    seh_mid = (int(seh["feature_start"]) + int(seh["feature_end"])) / 2
    is3_mid = (int(is3["feature_start"]) + int(is3["feature_end"])) / 2
    sign = -1 if is3_mid < seh_mid else 1
    return sign * coordinate_gap(seh, is3)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = [r for r in read_tsv(SUMMARY_TSV) if r.get("status") == "ok" and r.get("seh_present") == "TRUE"]
    features = read_tsv(FEATURES_TSV)
    features_by_sample: dict[str, list[dict[str, str]]] = {}
    for row in features:
        features_by_sample.setdefault(row["sample_id"], []).append(row)

    fasta_cache: dict[Path, dict[str, str]] = {}
    flank_records: dict[int, list[tuple[str, str]]] = {10_000: [], 20_000: []}
    spacer_records: list[tuple[str, str]] = []
    spacer_rows: list[dict[str, str | int]] = []

    for row in summary:
        sample = row["sample_id"]
        bakta_tsv = row.get("bakta_tsv", "")
        if not bakta_tsv:
            continue
        fasta_path = bakta_fasta_from_tsv(bakta_tsv)
        if fasta_path not in fasta_cache:
            if not fasta_path.exists():
                spacer_rows.append({"sample_id": sample, "status": "missing_bakta_fna", "bakta_fna": str(fasta_path.relative_to(PROJECT_DIR))})
                continue
            fasta_cache[fasta_path] = read_fasta(fasta_path)
        contigs = fasta_cache[fasta_path]
        contig = row["seh_contig"]
        seq = contigs.get(contig)
        if seq is None:
            spacer_rows.append({"sample_id": sample, "status": "missing_contig_in_bakta_fna", "seh_contig": contig})
            continue

        seh_start = int(row["seh_start"])
        seh_end = int(row["seh_end"])
        for flank_bp in flank_records:
            start = max(1, seh_start - flank_bp)
            end = min(len(seq), seh_end + flank_bp)
            flank = seq[start - 1:end]
            header = (
                f"{sample}|{contig}:{start}-{end}|seh={seh_start}-{seh_end}"
                f"|flank={flank_bp}|truncated_left={start == 1}|truncated_right={end == len(seq)}"
            )
            flank_records[flank_bp].append((header, flank))

        sample_features = features_by_sample.get(sample, [])
        seh_features = [f for f in sample_features if f.get("category") == "seh" and f.get("contig") == contig]
        is3_features = [f for f in sample_features if f.get("category") == "is3_transposase" and f.get("contig") == contig]
        if not seh_features or not is3_features:
            spacer_rows.append({"sample_id": sample, "status": "missing_seh_or_is3_feature", "seh_contig": contig})
            continue
        seh_feature = seh_features[0]
        nearest_is3 = min(is3_features, key=lambda f: abs(signed_gap(seh_feature, f)))
        spacer_start, spacer_end = inter_feature_bounds(seh_feature, nearest_is3)
        if spacer_end < spacer_start:
            spacer = ""
        else:
            spacer = seq[spacer_start - 1:spacer_end]
        coord_gap = coordinate_gap(seh_feature, nearest_is3)
        signed_coord_gap = signed_gap(seh_feature, nearest_is3)
        header = (
            f"{sample}|{contig}:{spacer_start}-{spacer_end}|coordinate_gap={coord_gap}"
            f"|sequence_length={len(spacer)}|signed_coordinate_gap={signed_coord_gap}"
        )
        spacer_records.append((header, spacer))
        spacer_rows.append({
            "sample_id": sample,
            "dataset": row["dataset"],
            "country": row["country"],
            "seh_contig": contig,
            "seh_start": seh_feature["feature_start"],
            "seh_end": seh_feature["feature_end"],
            "seh_strand": seh_feature["strand"],
            "is3_start": nearest_is3["feature_start"],
            "is3_end": nearest_is3["feature_end"],
            "is3_strand": nearest_is3["strand"],
            "coordinate_gap_bp": coord_gap,
            "signed_coordinate_gap_bp": signed_coord_gap,
            "between_features_sequence_start": spacer_start,
            "between_features_sequence_end": spacer_end,
            "between_features_sequence_length_bp": len(spacer),
            "status": "ok",
        })

    for flank_bp, records in flank_records.items():
        out = OUT_DIR / f"seh_flanks_{flank_bp // 1000}kb.fasta"
        with out.open("w", encoding="utf-8") as handle:
            for header, seq in records:
                handle.write(f">{header}\n{wrap(seq)}\n")

    with (OUT_DIR / "is3_seh_between_features.fasta").open("w", encoding="utf-8") as handle:
        for header, seq in spacer_records:
            handle.write(f">{header}\n{wrap(seq)}\n")

    fields = [
        "sample_id", "dataset", "country", "seh_contig", "seh_start", "seh_end", "seh_strand",
        "is3_start", "is3_end", "is3_strand", "coordinate_gap_bp", "signed_coordinate_gap_bp",
        "between_features_sequence_start", "between_features_sequence_end",
        "between_features_sequence_length_bp", "status",
    ]
    with (OUT_DIR / "is3_seh_intergenic_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(spacer_rows)

    print(f"Wrote {OUT_DIR / 'seh_flanks_10kb.fasta'}")
    print(f"Wrote {OUT_DIR / 'seh_flanks_20kb.fasta'}")
    print(f"Wrote {OUT_DIR / 'is3_seh_between_features.fasta'}")
    print(f"Wrote {OUT_DIR / 'is3_seh_intergenic_summary.tsv'}")


if __name__ == "__main__":
    main()
