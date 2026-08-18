#!/usr/bin/env python3
"""Compare the rlmH--cstA interval in ST1 and seh-positive/negative ST2187."""

from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/seh_st1_bovine/st1_st2187_phylogenomics/anchor_interval"

SAMPLES = {
    "23EV612": (
        "ST1_seh_positive",
        ROOT / "results/voss_2025_msphere_st1/bakta/23EV612/23EV612.fna",
        ROOT / "results/voss_2025_msphere_st1/bakta/23EV612/23EV612.gff3",
    ),
    "cepa_18_barcode02": (
        "ST2187_seh_positive",
        ROOT / "results/bakta/cepa_18_barcode02/cepa_18_barcode02.fna",
        ROOT / "results/bakta/cepa_18_barcode02/cepa_18_barcode02.gff3",
    ),
    **{
        accession: (
            "ST2187_seh_negative",
            ROOT / f"results/public_context_regional/bakta/{accession}/{accession}.fna",
            ROOT / f"results/public_context_regional/bakta/{accession}/{accession}.gff3",
        )
        for accession in (
            "GCA_031889225.1",
            "GCA_031889265.1",
            "GCA_031889285.1",
            "GCA_031889365.1",
            "GCA_031889545.1",
        )
    },
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name:
                records[name].append(line)
    return {name: "".join(parts).upper() for name, parts in records.items()}


def attr_value(attrs: str, key: str) -> str:
    match = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", attrs)
    return match.group(1) if match else ""


def read_gff(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) != 9:
                continue
            rows.append({
                "contig": parts[0], "type": parts[2], "start": int(parts[3]),
                "end": int(parts[4]), "strand": parts[6], "attrs": parts[8],
                "gene": attr_value(parts[8], "gene"),
                "name": attr_value(parts[8], "Name"),
                "product": attr_value(parts[8], "product"),
            })
    return rows


def label(row: dict[str, object]) -> str:
    return str(row["gene"] or row["name"] or row["product"] or row["type"])


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fasta_out = OUTDIR / "rlmH_cstA_intervals.fasta"
    summary_out = OUTDIR / "rlmH_cstA_architecture.tsv"
    summaries = []
    fasta_records = []

    for sample, (group, fasta, gff) in SAMPLES.items():
        features = read_gff(gff)
        anchors = {gene: [r for r in features if r["gene"] == gene] for gene in ("rlmH", "cstA")}
        if len(anchors["rlmH"]) != 1 or len(anchors["cstA"]) != 1:
            summaries.append({"sample": sample, "group": group, "status": "anchors_not_resolved"})
            continue
        rlmh, csta = anchors["rlmH"][0], anchors["cstA"][0]
        if rlmh["contig"] != csta["contig"]:
            summaries.append({"sample": sample, "group": group, "status": "anchors_on_different_contigs"})
            continue
        contig = str(rlmh["contig"])
        lo = min(int(rlmh["start"]), int(csta["start"]))
        hi = max(int(rlmh["end"]), int(csta["end"]))
        seq = read_fasta(fasta)[contig][lo - 1:hi]
        inner = [r for r in features if r["contig"] == contig and int(r["start"]) >= lo and int(r["end"]) <= hi and r["type"] == "CDS"]
        inner.sort(key=lambda r: int(r["start"]))
        architecture = " | ".join(f"{label(r)}({r['strand']})" for r in inner)
        fasta_records.append((f"{sample}|{group}|{contig}:{lo}-{hi}|len={len(seq)}", seq))
        summaries.append({
            "sample": sample, "group": group, "status": "resolved", "contig": contig,
            "start": lo, "end": hi, "interval_length_bp": len(seq),
            "cds_count": len(inner), "architecture": architecture,
        })

    with fasta_out.open("w") as handle:
        for header, seq in fasta_records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")
    fields = ["sample", "group", "status", "contig", "start", "end", "interval_length_bp", "cds_count", "architecture"]
    with summary_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)

    blast_out = OUTDIR / "all_vs_all_anchor_intervals.tsv"
    cmd = [
        "blastn", "-query", str(fasta_out), "-subject", str(fasta_out), "-task", "blastn",
        "-dust", "no", "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-out", str(blast_out),
    ]
    subprocess.run(cmd, check=True)
    print(summary_out)
    print(blast_out)


if __name__ == "__main__":
    main()
