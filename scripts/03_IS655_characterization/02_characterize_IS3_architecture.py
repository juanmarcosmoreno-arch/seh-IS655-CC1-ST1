#!/usr/bin/env python3
"""Characterize the IS3-family region adjacent to seh.

The Bakta feature called "IS3 family transposase" is a CDS annotation, not an
insertion-sequence boundary.  This workflow therefore extracts an oriented
window around that CDS and searches it in translated space against selected
ISfinder IS3-family proteins.  It reports the coding span supported by the
reference ORFs and a conservative projection of the complete element bounds.

The projection is intended to guide boundary/TIR/TSD inspection; it is not, by
itself, a formal ISfinder element call.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TABLE = ROOT / "results/seh_st1_bovine/final_tables/final_genome_locus_tree_table.tsv"
DEFAULT_IS_FAA = ROOT / "resources/isfinder_thanhleviet/IS.faa"
DEFAULT_OUTDIR = ROOT / "results/seh_st1_bovine/is3_architecture"
TARGETS = ("IS655", "ISClsp1", "ISSau2", "ISSau4")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


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
                records[current].append(line)
    return {key: "".join(value).upper() for key, value in records.items()}


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def genome_fasta(row: dict[str, str]) -> Path:
    return (ROOT / row["bakta_tsv"]).with_suffix(".fna")


def extract_windows(rows: list[dict[str, str]], flank: int) -> tuple[list[dict[str, str]], list[tuple[str, str]]]:
    metadata: list[dict[str, str]] = []
    records: list[tuple[str, str]] = []
    for row in rows:
        if row.get("formal_locus_status") != "ok" or not row.get("nearest_is3_start"):
            continue
        seqs = read_fasta(genome_fasta(row))
        contig = row["seh_contig"]
        genome_seq = seqs[contig]
        start = int(row["nearest_is3_start"])
        end = int(row["nearest_is3_end"])
        lo, hi = sorted((start, end))
        window_start = max(1, lo - flank)
        window_end = min(len(genome_seq), hi + flank)
        window = genome_seq[window_start - 1 : window_end]
        cds_start = lo - window_start + 1
        cds_end = hi - window_start + 1
        if row["nearest_is3_strand"] == "-":
            window = rc(window)
            cds_start, cds_end = window_end - hi + 1, window_end - lo + 1
        sample = row["sample_id"]
        records.append((sample, window))
        metadata.append(
            {
                "sample_id": sample,
                "group": row["analysis_group"],
                "tree_phylo_group": row["tree_phylo_group"],
                "contig": contig,
                "genomic_window_start": str(window_start),
                "genomic_window_end": str(window_end),
                "original_is3_strand": row["nearest_is3_strand"],
                "oriented_bakta_cds_start": str(cds_start),
                "oriented_bakta_cds_end": str(cds_end),
                "bakta_cds_length_bp": str(hi - lo + 1),
                "window_length_bp": str(len(window)),
            }
        )
    return metadata, records


def select_reference_proteins(path: Path) -> tuple[list[dict[str, str]], list[tuple[str, str]]]:
    selected: list[tuple[str, str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    selected.append((header.split()[0], header, "".join(parts)))
                header, parts = line[1:], []
            elif header is not None:
                parts.append(line)
        if header is not None:
            selected.append((header.split()[0], header, "".join(parts)))

    counts: dict[str, int] = defaultdict(int)
    metadata: list[dict[str, str]] = []
    records: list[tuple[str, str]] = []
    for token, full_header, seq in selected:
        element = token.split("//")[-1]
        if element not in TARGETS:
            continue
        counts[element] += 1
        query_id = f"{element}_ORF{counts[element]}_aa{len(seq)}"
        records.append((query_id, seq))
        metadata.append(
            {
                "query_id": query_id,
                "reference_element": element,
                "reference_orf_index": str(counts[element]),
                "query_length_aa": str(len(seq)),
                "original_header": full_header,
            }
        )
    return metadata, records


def run_tblastn(query: Path, subject: Path, output: Path) -> None:
    fields = "qseqid qlen sseqid sstart send pident length qstart qend evalue bitscore"
    subprocess.run(
        [
            "tblastn",
            "-query",
            str(query),
            "-subject",
            str(subject),
            "-evalue",
            "1e-5",
            "-seg",
            "no",
            "-outfmt",
            f"6 {fields}",
            "-out",
            str(output),
        ],
        check=True,
    )


def parse_hits(path: Path, ref_meta: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    fields = ["query_id", "query_length_aa", "sample_id", "subject_start", "subject_end", "pident", "alignment_length_aa", "query_start", "query_end", "evalue", "bitscore"]
    hits: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            values = raw.rstrip("\n").split("\t")
            row = dict(zip(fields, values))
            qlen = int(row["query_length_aa"])
            qstart, qend = int(row["query_start"]), int(row["query_end"])
            row["query_coverage"] = f"{(qend - qstart + 1) / qlen:.4f}"
            row["hit_strand"] = "+" if int(row["subject_start"]) <= int(row["subject_end"]) else "-"
            row.update({k: v for k, v in ref_meta[row["query_id"]].items() if k != "query_id"})
            hits.append(row)
    return hits


def best_hits(hits: list[dict[str, str]]) -> list[dict[str, str]]:
    best: dict[tuple[str, str], dict[str, str]] = {}
    for hit in hits:
        key = (hit["sample_id"], hit["query_id"])
        if key not in best or float(hit["bitscore"]) > float(best[key]["bitscore"]):
            best[key] = hit
    return list(best.values())


def summarize(samples: list[dict[str, str]], hits: list[dict[str, str]]) -> list[dict[str, str]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in hits:
        if float(hit["query_coverage"]) >= 0.80 and hit["hit_strand"] == "+":
            by_sample[hit["sample_id"]].append(hit)

    output: list[dict[str, str]] = []
    for sample in samples:
        sample_id = sample["sample_id"]
        candidates = by_sample.get(sample_id, [])
        ranked: list[tuple[float, str, list[dict[str, str]]]] = []
        for element in TARGETS:
            element_hits = [hit for hit in candidates if hit["reference_element"] == element]
            if not element_hits:
                continue
            score = sum(float(hit["bitscore"]) for hit in element_hits)
            ranked.append((score, element, element_hits))
        ranked.sort(reverse=True)
        if not ranked:
            output.append({**sample, "best_reference_element": "none", "complete_reference_orfs_detected": "0", "translated_homology_start": "NA", "translated_homology_end": "NA", "translated_homology_span_bp": "NA", "span_vs_bakta_cds_bp": "NA", "interpretation": "No >=80% query-coverage selected-reference hit"})
            continue
        _, element, element_hits = ranked[0]
        starts = [min(int(hit["subject_start"]), int(hit["subject_end"])) for hit in element_hits]
        ends = [max(int(hit["subject_start"]), int(hit["subject_end"])) for hit in element_hits]
        span_start, span_end = min(starts), max(ends)
        span = span_end - span_start + 1
        bakta_len = int(sample["bakta_cds_length_bp"])
        output.append(
            {
                **sample,
                "best_reference_element": element,
                "complete_reference_orfs_detected": str(len(element_hits)),
                "translated_homology_start": str(span_start),
                "translated_homology_end": str(span_end),
                "translated_homology_span_bp": str(span),
                "span_vs_bakta_cds_bp": str(span - bakta_len),
                "interpretation": "Translated homology extends beyond Bakta CDS" if span > bakta_len + 30 else "Translated homology approximately matches Bakta CDS",
            }
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--is-faa", type=Path, default=DEFAULT_IS_FAA)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--flank", type=int, default=5000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    samples, windows = extract_windows(read_table(args.table), args.flank)
    ref_rows, ref_records = select_reference_proteins(args.is_faa)
    windows_fasta = args.outdir / "seh_adjacent_is3_oriented_windows.fasta"
    refs_fasta = args.outdir / "selected_is3_reference_proteins.fasta"
    raw_hits = args.outdir / "selected_is3_tblastn_raw.tsv"
    write_fasta(windows_fasta, windows)
    write_fasta(refs_fasta, ref_records)
    write_tsv(args.outdir / "selected_is3_reference_metadata.tsv", ref_rows)
    run_tblastn(refs_fasta, windows_fasta, raw_hits)

    ref_meta = {row["query_id"]: row for row in ref_rows}
    hits = best_hits(parse_hits(raw_hits, ref_meta))
    hits.sort(key=lambda row: (row["sample_id"], row["reference_element"], int(row["reference_orf_index"])))
    write_tsv(args.outdir / "selected_is3_tblastn_best_hits.tsv", hits)
    summary = summarize(samples, hits)
    write_tsv(args.outdir / "is3_architecture_summary.tsv", summary)
    print(f"Analyzed {len(samples)} complete loci; wrote {args.outdir}")


if __name__ == "__main__":
    main()
