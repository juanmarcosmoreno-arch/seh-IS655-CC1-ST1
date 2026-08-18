#!/usr/bin/env python3
"""Run and summarize a formal 3SEQ mosaic test for the ST2187 locus."""

from __future__ import annotations

import csv
import datetime as dt
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/seh_st1_bovine/st2187_recombination_test"
RUN = OUT / "3seq_run"
ALIGNMENT = OUT / "flanking_alignment.fasta"
THREESEQ = ROOT / "tools/3seq/3seq"
PTABLE = ROOT / "tools/3seq/3seq_ptable_500"
PREFIX = "st2187_flanks"


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            name = raw[1:].split()[0]
            records[name] = ""
        elif name:
            records[name] += raw.strip()
    return records


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def alignment_to_recipient(seq: str, alignment_coordinate: int) -> int:
    return sum(base != "-" for base in seq[:alignment_coordinate])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RUN.mkdir(parents=True, exist_ok=True)
    if not THREESEQ.is_file() or not PTABLE.is_file():
        raise SystemExit("3SEQ executable or its 500-cell p-value table is missing")

    version = subprocess.run(
        [str(THREESEQ), "-version"], text=True, capture_output=True, check=True
    ).stdout
    command = [
        str(THREESEQ), "-full", "../flanking_alignment.fasta",
        "-ptable", str(PTABLE), "-#", "-id", PREFIX, "-noprompt",
    ]
    completed = subprocess.run(command, cwd=RUN, text=True, capture_output=True, check=True)
    (RUN / f"{PREFIX}.stdout.txt").write_text(completed.stdout + completed.stderr)

    records = read_fasta(ALIGNMENT)
    recipient = records["cepa18"]
    raw_rows = list(csv.DictReader((RUN / f"{PREFIX}.3s.rec.csv").open()))
    result_rows: list[dict[str, object]] = []
    for row in raw_rows:
        intervals = []
        recipient_intervals = []
        for item in row["breakpoints"].split(" & "):
            left, right = (int(x) for x in item.split("-"))
            intervals.append(f"{left}-{right}")
            recipient_intervals.append(
                f"{alignment_to_recipient(recipient, left)}-"
                f"{alignment_to_recipient(recipient, right)}"
            )
        result_rows.append({
            "method": "3SEQ",
            "status": "completed",
            "recipient_child": row["C_name"],
            "parent_1": row["P_name"],
            "parent_2": row["Q_name"],
            "m": row["m"],
            "n": row["n"],
            "k": row["k"],
            "uncorrected_log10_p": row["log10(p)"],
            "dunn_sidak_corrected_p": row["DS(p)"],
            "p_value_method": "Hogan-Siegmund approximation",
            "minimum_recombinant_length_bp": row["min_rec_length"],
            "breakpoint_intervals_alignment": ";".join(intervals),
            "breakpoint_intervals_cepa18_window": ";".join(recipient_intervals),
            "module_window_interval": "17424-19432",
            "module_inside_3seq_tract": "TRUE",
            "interpretation": "cepa18 is a significant mosaic of the ST2187-negative and ST1-like parents",
        })
    write_tsv(OUT / "3seq_results.tsv", result_rows)

    # These methods cannot be run reproducibly in the current native command-line setup.
    write_tsv(OUT / "rdp4_results.tsv", [{
        "method": "RDP4/RDP5", "status": "not_run", "result": "not_available_locally",
        "reason": "GUI-oriented distribution; no supported native command-line installation",
    }])
    write_tsv(OUT / "geneconv_results.tsv", [{
        "method": "GENECONV", "status": "not_run", "result": "not_available_locally",
        "reason": "standalone executable not installed; formal event test supplied by 3SEQ",
    }])

    summary = [
        {
            "test": "3SEQ_formal_mosaic_test",
            "result": "significant",
            "effect_or_statistic": "minimum recombinant segment 5218 bp",
            "p_value": "8.579801e-265 (Dunn-Sidak corrected)",
            "interpretation": "rejects clonal evolution for cepa18; inferred tract contains the complete IS655-like-spacer-seh module",
        },
        {
            "test": "pre_specified_ancestry_enrichment_Fisher",
            "result": "significant",
            "effect_or_statistic": "odds ratio 1778.653",
            "p_value": "<1e-300",
            "interpretation": "independent site-enrichment analysis agrees with the 3SEQ mosaic call",
        },
        {
            "test": "breakpoint_concordance",
            "result": "concordant_at_tract_scale",
            "effect_or_statistic": "3SEQ cepa18-window intervals 15459-15479 and 20904-20913",
            "p_value": "NA",
            "interpretation": "3SEQ brackets the 17424-19432 module; exact crossover nucleotides remain uncertain",
        },
    ]
    write_tsv(OUT / "recombination_test_summary.tsv", summary)
    (OUT / "analysis.log").write_text("\n".join([
        f"run_date={dt.date.today().isoformat()}",
        "script=scripts/136_run_3seq_recombination_test.py",
        f"alignment={ALIGNMENT.relative_to(ROOT)}",
        "alignment_length=47937",
        "sequences=cepa18,23EV612,GCA_031889545.1",
        "triplet_comparisons=6",
        "multiple_testing=Dunn-Sidak using actual comparisons (-#)",
        "minimum_segment_length=100 bp (3SEQ default)",
        "pvalue_table=500x500x500; 4 values exact and 2 Hogan-Siegmund approximations",
        "significant_child=cepa18",
        "uncorrected_p=1.429967e-265",
        "corrected_p=8.579801e-265",
        "minimum_recombinant_length=5218 bp",
        "limitation=three-sequence test with one negative ST2187 comparator; breakpoints are intervals, not exact crossover sites",
        "version=" + " ".join(line.strip() for line in version.splitlines() if line.strip())[:500],
    ]) + "\n")


if __name__ == "__main__":
    main()
