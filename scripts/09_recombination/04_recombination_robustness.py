#!/usr/bin/env python3
"""Robustness of the ST1->ST2187 recombination call to the choice of
ST2187-negative "parent" reference (Major comments 1, 6, 7).

The submitted manuscript used a single ST2187-negative comparator
(GCA_031889545.1) as the third sequence in the 3SEQ/GENECONV triplet test
with cepa18 (recipient, strain_18_barcode02) and 23EV612 (ST1-like donor).
This script repeats the exact same pipeline (MAFFT three-way window
alignment -> 3SEQ -full -> GENECONV -Gscale=0) with alternative
ST2187-negative references drawn from Table S11
(locus_state=ancestral_defense_block), and compiles a small comparison
table of P-value / minimum recombinant length / breakpoints across runs.

Does not touch any file used by the submitted manuscript; everything is
written to results/seh_st1_bovine/st2187_recombination_test/robustness/.
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
BAKTA_DIR = ROOT / "results/public_context_regional/bakta"
OUT = ROOT / "results/seh_st1_bovine/st2187_recombination_test/robustness"
MAFFT = ROOT / ".micromamba/envs/gridion_pangenome_pip/bin/mafft"
THREESEQ = ROOT / "tools/3seq/3seq"
PTABLE = ROOT / "tools/3seq/3seq_ptable_500"
GENECONV = ROOT / "tools/geneconv/geneconv"
FLANK = 15000

FIXED_SAMPLES = {
    "cepa18": (ROOT / "results/bakta/cepa_18_barcode02/cepa_18_barcode02.fna",
               ROOT / "results/bakta/cepa_18_barcode02/cepa_18_barcode02.gff3"),
    "23EV612": (ROOT / "results/voss_2025_msphere_st1/bakta/23EV612/23EV612.fna",
                ROOT / "results/voss_2025_msphere_st1/bakta/23EV612/23EV612.gff3"),
}

# Table S11 ancestral_defense_block genomes with Bakta annotation available
# locally AND rlmH/cstA on the same contig (required for the flanking-window
# method). GCA_031889225.1 was excluded: rlmH and cstA fall on two different
# contigs in that assembly (contig_402 vs contig_390), i.e. the assembly is
# itself fragmented at this locus. Souza_2023 cohort genomes (e.g.
# GCA_029101905.1, GCA_032282325.1) are not Bakta-annotated in this project
# and were not re-annotated for this sensitivity check.
ALT_NEGATIVES = [
    "GCA_031889265.1",
    "GCA_031889285.1",
    "GCA_031889365.1",
]
ORIGINAL_NEGATIVE = "GCA_031889545.1"


def anchors(gff: Path):
    d = {}
    for line in open(gff):
        if line.startswith("#"):
            continue
        x = line.rstrip().split("\t")
        if len(x) != 9:
            continue
        a = dict(p.split("=", 1) for p in x[8].split(";") if "=" in p)
        gene = a.get("gene", "")
        if x[2] == "CDS" and gene in {"rlmH", "cstA"}:
            d[gene] = {"contig": x[0], "start": int(x[3]), "end": int(x[4])}
    return d


def genome(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def bakta_paths(acc: str) -> tuple[Path, Path]:
    d = BAKTA_DIR / acc
    return d / f"{acc}.fna", d / f"{acc}.gff3"


def extract_window(name: str, fa: Path, gff: Path) -> str:
    a = anchors(gff)
    assert a["rlmH"]["contig"] == a["cstA"]["contig"], f"{name}: rlmH/cstA on different contigs"
    contig = a["rlmH"]["contig"]
    seq = genome(fa)[contig]
    lo = max(1, min(a["rlmH"]["start"], a["cstA"]["start"]) - FLANK)
    hi = min(len(seq), max(a["rlmH"]["end"], a["cstA"]["end"]) + FLANK)
    raw = seq[lo - 1:hi]
    forward = a["rlmH"]["start"] < a["cstA"]["start"]
    return raw if forward else str(Seq(raw).reverse_complement())


def align_to_recipient(seq: str, alignment_coordinate: int) -> int:
    return sum(base != "-" for base in seq[:alignment_coordinate])


def run_alignment(run_dir: Path, negative_acc: str) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    windows = {}
    windows["cepa18"] = extract_window("cepa18", *FIXED_SAMPLES["cepa18"])
    windows["23EV612"] = extract_window("23EV612", *FIXED_SAMPLES["23EV612"])
    windows[negative_acc] = extract_window(negative_acc, *bakta_paths(negative_acc))

    combined = run_dir / "three_wide_windows.fasta"
    with combined.open("w") as h:
        for n, s in windows.items():
            h.write(f">{n}\n{s}\n")

    env = os.environ.copy()
    env["MAFFT_BINARIES"] = str(MAFFT.parent.parent / "libexec/mafft")
    result = subprocess.run([str(MAFFT), "--auto", str(combined)],
                             text=True, capture_output=True, check=True, env=env)
    aligned = run_dir / "three_way_alignment.fasta"
    aligned.write_text(result.stdout)
    return aligned


def run_3seq(run_dir: Path, alignment: Path, negative_acc: str) -> dict:
    sub = run_dir / "3seq_run"
    sub.mkdir(exist_ok=True)
    flanking = sub.parent / "flanking_alignment.fasta"
    flanking.write_text(alignment.read_text())
    prefix = "robustness"
    cmd = [str(THREESEQ), "-full", "../flanking_alignment.fasta",
           "-ptable", str(PTABLE), "-#", "-id", prefix, "-noprompt"]
    completed = subprocess.run(cmd, cwd=sub, text=True, capture_output=True, check=True)
    (sub / f"{prefix}.stdout.txt").write_text(completed.stdout + completed.stderr)

    rec_csv = sub / f"{prefix}.3s.rec.csv"
    if not rec_csv.is_file() or rec_csv.stat().st_size == 0:
        return {"negative_reference": negative_acc, "method": "3SEQ",
                "result": "no_significant_triplet", "p_value": "NA",
                "min_recombinant_length_bp": "NA", "breakpoints_cepa18_window": "NA"}

    records = {}
    name = ""
    for raw in alignment.read_text().splitlines():
        if raw.startswith(">"):
            name = raw[1:].split()[0]
            records[name] = ""
        elif name:
            records[name] += raw.strip()
    recipient = records["cepa18"]

    rows = list(csv.DictReader(rec_csv.open()))
    if not rows:
        return {"negative_reference": negative_acc, "method": "3SEQ",
                "result": "no_significant_triplet", "p_value": "NA",
                "min_recombinant_length_bp": "NA", "breakpoints_cepa18_window": "NA"}
    row = rows[0]
    recipient_intervals = []
    for item in row["breakpoints"].split(" & "):
        left, right = (int(x) for x in item.split("-"))
        recipient_intervals.append(f"{align_to_recipient(recipient, left)}-{align_to_recipient(recipient, right)}")
    return {
        "negative_reference": negative_acc,
        "method": "3SEQ",
        "result": "significant",
        "p_value": row["DS(p)"],
        "min_recombinant_length_bp": row["min_rec_length"],
        "breakpoints_cepa18_window": ";".join(recipient_intervals),
    }


def run_geneconv(run_dir: Path, alignment: Path, negative_acc: str) -> dict:
    sub = run_dir / "geneconv_run"
    sub.mkdir(exist_ok=True)
    fasta = sub / "st2187_flanks.fasta"
    fasta.write_text(alignment.read_text())
    cmd = [str(GENECONV), "st2187_flanks.fasta", "geneconv_g0",
           "-Gscale=0", "-Numsim=10000", "-Startseed=2187",
           "-Dumptab", "-ExpFormat", "-Annotate", "-ListPairs"]
    subprocess.run(cmd, cwd=sub, text=True, capture_output=True, check=True)
    tab = sub / "geneconv_g0.tab"
    if not tab.is_file():
        return {"negative_reference": negative_acc, "method": "GENECONV",
                "result": "no_output", "p_value": "NA",
                "min_recombinant_length_bp": "NA", "breakpoints_cepa18_window": "NA"}

    header = None
    rows = []
    for line in tab.read_text().splitlines():
        if line.startswith("# Frag\tSequence"):
            header = ["fragment_type", "sequence_names", "sim_p", "bonferroni_KA_p",
                      "aligned_begin", "aligned_end", "aligned_length",
                      "seq1_begin", "seq1_end", "seq1_length",
                      "seq2_begin", "seq2_end", "seq2_length",
                      "num_poly", "num_dif", "total_difs", "mismatch_penalty"]
        elif line.startswith("GI\t") and header:
            x = line.rstrip().split("\t")
            rows.append(dict(zip(header, x)))
    focal = [r for r in rows if r["sequence_names"] == "cepa18;23EV612"]
    if not focal:
        return {"negative_reference": negative_acc, "method": "GENECONV",
                "result": "no_significant_fragment", "p_value": "NA",
                "min_recombinant_length_bp": "NA", "breakpoints_cepa18_window": "NA"}

    begins = [int(r["seq1_begin"]) for r in focal]
    ends = [int(r["seq1_end"]) for r in focal]
    best_p = min(float(r["bonferroni_KA_p"]) for r in focal)
    envelope_len = max(ends) - min(begins) + 1
    return {
        "negative_reference": negative_acc,
        "method": "GENECONV",
        "result": "significant",
        "p_value": f"{best_p:.3e}",
        "min_recombinant_length_bp": str(envelope_len),
        "breakpoints_cepa18_window": f"{min(begins)}-{max(ends)}",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_refs = [ORIGINAL_NEGATIVE] + ALT_NEGATIVES
    results = []
    for acc in all_refs:
        tag = "original" if acc == ORIGINAL_NEGATIVE else "alternative"
        print(f"\n=== {acc} ({tag}) ===")
        run_dir = OUT / acc
        aligned = run_alignment(run_dir, acc)
        n_bp = len(next(iter(
            {r.id: str(r.seq) for r in SeqIO.parse(aligned, "fasta")}.values()
        )))
        r3 = run_3seq(run_dir, aligned, acc)
        rg = run_geneconv(run_dir, aligned, acc)
        for r in (r3, rg):
            r["reference_role"] = tag
            r["alignment_length_bp"] = n_bp
            results.append(r)
            print(f"  {r['method']:9s} result={r['result']:22s} p={r['p_value']:>14s}  "
                  f"min_len={r['min_recombinant_length_bp']:>6s}  bp={r['breakpoints_cepa18_window']}")

    out_path = OUT / "recombination_robustness_summary.tsv"
    fieldnames = ["negative_reference", "reference_role", "method", "result",
                  "p_value", "min_recombinant_length_bp",
                  "breakpoints_cepa18_window", "alignment_length_bp"]
    with out_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print("\nWritten:", out_path)


if __name__ == "__main__":
    main()
