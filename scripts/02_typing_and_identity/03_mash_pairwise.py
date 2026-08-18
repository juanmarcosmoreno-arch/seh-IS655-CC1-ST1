#!/usr/bin/env python3
"""Compute and summarize Mash distances for the 53-genome ST1/ST2187 panel."""

from __future__ import annotations

import csv
import statistics
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MASH = ROOT / ".micromamba/envs/gridion_nanopore/bin/mash"
SKETCH = ROOT / "results/seh_st1_bovine/st1_st2187_phylogenomics/mash53.msh"
OUT = ROOT / "results/seh_st1_bovine/st1_st2187_phylogenomics/mash53_pairwise_distances.tsv"
SUMMARY = ROOT / "results/seh_st1_bovine/st1_st2187_phylogenomics/mash53_st2187_summary.tsv"
POSITIVE = "cepa_18_barcode02"
NEGATIVES = {
    "GCA_031889225.1", "GCA_031889265.1", "GCA_031889285.1",
    "GCA_031889365.1", "GCA_031889545.1",
}


def sample_name(value: str) -> str:
    return Path(value.split()[0]).stem


def main() -> None:
    completed = subprocess.run(
        [str(MASH), "dist", str(SKETCH), str(SKETCH)],
        check=True, text=True, capture_output=True,
    )
    rows = []
    lookup = {}
    for line in completed.stdout.splitlines():
        ref, query, distance, pvalue, shared = line.split("\t")
        ref, query = sample_name(ref), sample_name(query)
        row = {
            "sample_a": ref, "sample_b": query, "mash_distance": distance,
            "p_value": pvalue, "shared_hashes": shared,
        }
        rows.append(row)
        lookup[(ref, query)] = float(distance)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    samples = sorted({row["sample_a"] for row in rows})
    st1 = [sample for sample in samples if sample != POSITIVE and sample not in NEGATIVES]
    neg_distances = [lookup[(POSITIVE, sample)] for sample in sorted(NEGATIVES)]
    st1_distances = [lookup[(POSITIVE, sample)] for sample in st1]
    summary_rows = [
        {
            "comparison": "positive_ST2187_to_negative_ST2187",
            "n": len(neg_distances), "minimum": min(neg_distances),
            "median": statistics.median(neg_distances), "maximum": max(neg_distances),
        },
        {
            "comparison": "positive_ST2187_to_ST1",
            "n": len(st1_distances), "minimum": min(st1_distances),
            "median": statistics.median(st1_distances), "maximum": max(st1_distances),
        },
    ]
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    for row in summary_rows:
        print(row)
    nearest = sorted((lookup[(POSITIVE, sample)], sample) for sample in samples if sample != POSITIVE)
    print("Nearest to positive ST2187:", nearest[:10])


if __name__ == "__main__":
    main()
