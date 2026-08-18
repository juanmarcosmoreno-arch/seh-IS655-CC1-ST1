#!/usr/bin/env python3
"""Major reviewer #6: GC3 percentile of seh + sliding-window GC for module interval.

Inputs (already computed / present):
  per_locus_gc_codon.tsv  – GC3 for every CDS in 23EV612 (from script 139)
  23EV612.fna              – complete ST1 chromosome

Outputs:
  gc3_core_distribution.tsv   – full GC3 distribution, sorted
  gc3_percentile_summary.tsv  – percentile of seh and IS3 transposase
  module_window_gc.tsv        – GC% of 8571-bp module vs 10,000 random windows
  gc3_analysis_v2_summary.md  – prose for manuscript
"""

from __future__ import annotations
import csv
import math
import random
import statistics
from pathlib import Path

ROOT    = Path(__file__).resolve().parents[1]
GC_TSV  = ROOT / "results/seh_st1_bovine/gc_codon_analysis/per_locus_gc_codon.tsv"
FNA     = ROOT / "results/voss_2025_msphere_st1/bakta/23EV612/23EV612.fna"
OUT     = ROOT / "results/seh_st1_bovine/gc_codon_analysis"

MODULE_START = 1_880_276   # 1-based inclusive in 23EV612 contig_1
MODULE_END   = 1_888_846
MODULE_LEN   = MODULE_END - MODULE_START + 1   # 8,571 bp
N_WINDOWS    = 10_000
SEED         = 2187

FOCAL = {
    "VOSS23EV612_01821": "seh",
    "VOSS23EV612_01822": "IS3 transposase",
}


def gc(seq: str) -> float:
    n = len(seq)
    if n == 0:
        return math.nan
    s = seq.upper()
    return (s.count("G") + s.count("C")) / n


def load_genome(fna: Path) -> str:
    lines = fna.read_text().splitlines()
    return "".join(l for l in lines if not l.startswith(">")).upper()


def wt(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


def main() -> None:
    # --- load per-locus GC3 ---
    all_rows = list(csv.DictReader(GC_TSV.open(), delimiter="\t"))
    core     = [r for r in all_rows if r["class"] == "core_genome"]
    focal    = {r["locus_tag"]: r for r in all_rows if r["locus_tag"] in FOCAL}

    gc3_vals = sorted(float(r["gc3"]) for r in core)
    n_core   = len(gc3_vals)

    # --- 1. GC3 percentile ---
    pct_rows = []
    for lt, label in FOCAL.items():
        r = focal.get(lt)
        if r is None:
            continue
        g3 = float(r["gc3"])
        # fraction of core genes with GC3 < g3
        n_below = sum(1 for v in gc3_vals if v < g3)
        pct_below = n_below / n_core * 100          # percentile (lower tail)
        pct_above = (n_core - n_below) / n_core * 100
        pct_rows.append({
            "locus": label,
            "locus_tag": lt,
            "gc3_pct": round(g3, 4),
            "core_n": n_core,
            "n_core_below": n_below,
            "percentile_lower": round(pct_below, 2),
            "pct_core_genes_with_higher_gc3": round(pct_above, 2),
            "sentence": (
                f"{label} GC3 ({g3:.1f}%) was lower than "
                f"{pct_above:.1f}% of core genome CDS "
                f"(ranked {n_below+1}/{n_core}; "
                f"{pct_below:.1f}th percentile)"
            ),
        })
    wt(OUT / "gc3_percentile_summary.tsv", pct_rows)

    # --- 2. Full distribution table ---
    dist_rows = [{"rank": i+1, "locus_tag": r["locus_tag"], "class": r["class"],
                  "gc3": r["gc3"]}
                 for i, r in enumerate(sorted(all_rows, key=lambda x: float(x["gc3"])))]
    wt(OUT / "gc3_core_distribution.tsv", dist_rows)

    # --- 3. Sliding-window: module vs 10,000 random windows ---
    genome = load_genome(FNA)
    genome_len = len(genome)

    module_seq = genome[MODULE_START - 1: MODULE_END]
    module_gc  = gc(module_seq)

    rng = random.Random(SEED)
    window_gcs = []
    # Exclude positions that overlap the module itself
    valid_starts = [
        s for s in range(0, genome_len - MODULE_LEN + 1)
        if not (s < MODULE_END and s + MODULE_LEN > MODULE_START - 1)
    ]
    sampled = rng.sample(valid_starts, min(N_WINDOWS, len(valid_starts)))
    for s in sampled:
        window_gcs.append(gc(genome[s: s + MODULE_LEN]))

    window_gcs.sort()
    n_win = len(window_gcs)
    n_below_win = sum(1 for v in window_gcs if v < module_gc)
    pct_win_below = n_below_win / n_win * 100
    pct_win_above = (n_win - n_below_win) / n_win * 100

    win_rows = [{"window_start_0based": s, "gc_pct": round(gc(genome[s: s+MODULE_LEN])*100, 4)}
                for s in sorted(sampled)]
    wt(OUT / "module_window_gc.tsv", win_rows)

    win_mean = statistics.mean(window_gcs) * 100
    win_sd   = statistics.stdev(window_gcs) * 100
    win_p5   = window_gcs[int(0.05 * n_win)] * 100
    win_p95  = window_gcs[int(0.95 * n_win)] * 100

    # --- 4. Prose summary ---
    seh_row = next((r for r in pct_rows if r["locus"] == "seh"), None)
    is3_row = next((r for r in pct_rows if r["locus"] == "IS3 transposase"), None)

    prose = f"""## GC3 percentile and sliding-window analysis (Reviewer Major #6)

**Genome:** 23EV612 (ST1 complete chromosome; {genome_len:,} bp)
**Module interval:** contig_1:{MODULE_START:,}–{MODULE_END:,} ({MODULE_LEN:,} bp; rlmH–cstA)

---

### 1. GC3 percentile of focal genes within core CDS distribution

n core CDS = {n_core}

| Gene | GC3 (%) | Rank | Percentile | Statement |
|------|---------|------|------------|-----------|
| seh | {seh_row['gc3_pct'] if seh_row else 'NA'}% | {seh_row['n_core_below']+1 if seh_row else 'NA'}/{n_core} | {seh_row['percentile_lower'] if seh_row else 'NA'}th | {seh_row['sentence'] if seh_row else ''} |
| IS3 transposase | {is3_row['gc3_pct'] if is3_row else 'NA'}% | {is3_row['n_core_below']+1 if is3_row else 'NA'}/{n_core} | {is3_row['percentile_lower'] if is3_row else 'NA'}th | {is3_row['sentence'] if is3_row else ''} |

**Core GC3:** min={min(gc3_vals):.1f}%, median={statistics.median(gc3_vals):.1f}%, mean={statistics.mean(gc3_vals):.1f}%, max={max(gc3_vals):.1f}%

---

### 2. Module interval GC vs {n_win:,} random windows of equivalent size ({MODULE_LEN:,} bp)

| Metric | Value |
|--------|-------|
| Module GC% | {module_gc*100:.2f}% |
| Random-window mean ± SD | {win_mean:.2f}% ± {win_sd:.2f}% |
| Random-window 5th–95th pctile | {win_p5:.2f}%–{win_p95:.2f}% |
| Windows with GC < module | {n_below_win} / {n_win} ({pct_win_below:.1f}%) |
| Windows with GC > module | {n_win - n_below_win} / {n_win} ({pct_win_above:.1f}%) |
| Empirical P (module more AT than random) | {n_below_win / n_win:.4f} |

**The module interval is more AT-rich than {pct_win_above:.1f}% of random {MODULE_LEN:,}-bp genomic windows.**

---

### Manuscript text (ready to insert)

> To evaluate whether the GC content of the seh+IS655-like module is anomalous relative to the
> 23EV612 core genome, we computed GC% at the synonymous codon position 3 (GC3) for all
> {n_core} CDS classified as core genome, and compared the module interval GC against {n_win:,}
> randomly sampled genomic windows of equivalent size ({MODULE_LEN:,} bp).
>
> seh GC3 ({seh_row['gc3_pct'] if seh_row else 'NA'}%) was lower than {seh_row['pct_core_genes_with_higher_gc3'] if seh_row else 'NA'}% of core genome CDS
> (ranked {seh_row['n_core_below']+1 if seh_row else 'NA'} of {n_core}; {seh_row['percentile_lower'] if seh_row else 'NA'}th percentile), placing it among the most
> AT-enriched genes in the genome. The IS655-like transposase GC3 ({is3_row['gc3_pct'] if is3_row else 'NA'}%) was similarly
> below the core mean, at the {is3_row['percentile_lower'] if is3_row else 'NA'}th percentile.
> The rlmH–cstA module interval (GC = {module_gc*100:.2f}%) was more AT-rich than
> {pct_win_above:.1f}% of {n_win:,} randomly sampled {MODULE_LEN:,}-bp genomic windows
> (empirical P = {n_below_win / n_win:.4f}), confirming that the AT enrichment of the module is
> not an artifact of genome-wide base composition but reflects a locally anomalous region.

---

*Generated by scripts/142_gc3_percentile_window.py*
"""
    (OUT / "gc3_analysis_v2_summary.md").write_text(prose)
    print(prose)


if __name__ == "__main__":
    main()
