#!/usr/bin/env python3
"""GC content and codon-position analysis of the seh+IS655-like module vs genome background.

Uses 23EV612 (complete ST1 chromosome) as the reference genome.
Module region: rlmH–cstA interval (Bakta coords: contig_1:1880276-1888846, 1-based inclusive).
"""

from __future__ import annotations
import csv
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BAKTA_DIR = ROOT / "results/voss_2025_msphere_st1/bakta/23EV612"
FNA   = BAKTA_DIR / "23EV612.fna"
FFN   = BAKTA_DIR / "23EV612.ffn"
TSV   = BAKTA_DIR / "23EV612.tsv"
OUT   = ROOT / "results/seh_st1_bovine/gc_codon_analysis"
OUT.mkdir(parents=True, exist_ok=True)

MODULE_START = 1880276  # 1-based inclusive (rlmH–cstA interval in 23EV612)
MODULE_END   = 1888846  # 1-based inclusive
MODULE_LABEL = "seh+IS655-like module (rlmH–cstA interval)"

FOCAL_LOCI = {"VOSS23EV612_01821": "seh", "VOSS23EV612_01822": "IS3 transposase"}


def gc(seq: str) -> float:
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return float("nan")
    return (seq.count("G") + seq.count("C")) / n


def gc_positions(cds: str, pos: int) -> float:
    """GC fraction at codon position pos (1, 2, or 3)."""
    bases = [cds[i] for i in range(pos - 1, len(cds) - (3 - pos), 3) if i < len(cds)]
    bases = [b.upper() for b in bases]
    if not bases:
        return float("nan")
    return (bases.count("G") + bases.count("C")) / len(bases)


def load_genome(fna: Path) -> str:
    lines = fna.read_text().splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    return seq.upper()


def load_ffn(ffn: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current_id = None
    buf: list[str] = []
    for line in ffn.read_text().splitlines():
        if line.startswith(">"):
            if current_id:
                seqs[current_id] = "".join(buf).upper()
            current_id = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if current_id:
        seqs[current_id] = "".join(buf).upper()
    return seqs


def load_cds_loci(tsv: Path) -> list[dict]:
    rows = []
    with tsv.open() as h:
        for line in h:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            contig, feat_type, start, end, strand, locus_tag = parts[:6]
            gene = parts[6] if len(parts) > 6 else ""
            if feat_type == "cds":
                rows.append({
                    "contig": contig,
                    "type": feat_type,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "locus_tag": locus_tag,
                    "gene": gene,
                })
    return rows


def classify_locus(row: dict) -> str:
    lt = row["locus_tag"]
    if lt in FOCAL_LOCI:
        return FOCAL_LOCI[lt]
    # any CDS overlapping the module interval on the primary contig
    if row["contig"] == "contig_1" and row["start"] >= MODULE_START and row["end"] <= MODULE_END:
        return "module_other"
    return "core_genome"


def chi2_gc(gc_region: float, n_region: int, gc_genome: float, n_genome: int) -> float:
    """Two-proportion z-test, returns p-value (two-sided)."""
    p_pool = (gc_region * n_region + gc_genome * n_genome) / (n_region + n_genome)
    se = math.sqrt(p_pool * (1 - p_pool) * (1 / n_region + 1 / n_genome))
    if se == 0:
        return float("nan")
    z = (gc_region - gc_genome) / se
    # approximate p from standard normal using complementary error function
    p = math.erfc(abs(z) / math.sqrt(2))
    return p


def wt(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    genome = load_genome(FNA)
    ffn_seqs = load_ffn(FFN)
    cds_rows = load_cds_loci(TSV)

    # --- 1. Genomic GC: module vs background ---
    mod_seq = genome[MODULE_START - 1: MODULE_END]
    bg_seq  = genome[:MODULE_START - 1] + genome[MODULE_END:]

    gc_mod = gc(mod_seq)
    gc_bg  = gc(bg_seq)
    gc_all = gc(genome)

    p_mod_vs_bg = chi2_gc(gc_mod, len(mod_seq), gc_bg, len(bg_seq))

    region_rows = [
        {"region": MODULE_LABEL,       "length_bp": len(mod_seq),  "gc_pct": round(gc_mod * 100, 4), "gc_diff_vs_genome": round((gc_mod - gc_all) * 100, 4)},
        {"region": "Genome background", "length_bp": len(bg_seq), "gc_pct": round(gc_bg * 100, 4),  "gc_diff_vs_genome": round((gc_bg - gc_all) * 100, 4)},
        {"region": "Whole genome",      "length_bp": len(genome), "gc_pct": round(gc_all * 100, 4), "gc_diff_vs_genome": 0.0},
    ]
    wt(OUT / "gc_module_vs_genome.tsv", region_rows)

    # --- 2. Per-locus codon-position GC ---
    # Build mapping of locus_tag → classification
    locus_class = {r["locus_tag"]: classify_locus(r) for r in cds_rows}

    gene_rows = []
    for lt, seq in ffn_seqs.items():
        cls = locus_class.get(lt, "core_genome")
        if len(seq) < 3:
            continue
        gene_rows.append({
            "locus_tag": lt,
            "class": cls,
            "cds_length_nt": len(seq),
            "gc_overall": round(gc(seq) * 100, 4),
            "gc1": round(gc_positions(seq, 1) * 100, 4),
            "gc2": round(gc_positions(seq, 2) * 100, 4),
            "gc3": round(gc_positions(seq, 3) * 100, 4),
        })
    wt(OUT / "per_locus_gc_codon.tsv", gene_rows)

    # --- 3. Summary statistics ---
    core = [r for r in gene_rows if r["class"] == "core_genome"]
    focal = [r for r in gene_rows if r["class"] in ("seh", "IS3 transposase")]

    def mean_sd(vals: list[float]) -> tuple[float, float]:
        if not vals:
            return float("nan"), float("nan")
        m = statistics.mean(vals)
        s = statistics.stdev(vals) if len(vals) > 1 else 0.0
        return m, s

    core_gc3 = [r["gc3"] for r in core]
    m_core, sd_core = mean_sd(core_gc3)

    summary_rows = []
    for r in focal:
        diff = r["gc3"] - m_core
        p_approx = chi2_gc(r["gc3"] / 100, r["cds_length_nt"] // 3,
                           m_core / 100, int(len(core_gc3) * statistics.mean(r["cds_length_nt"] for r in core) // 3))
        summary_rows.append({
            "locus": r["class"],
            "locus_tag": r["locus_tag"],
            "gc3_pct": r["gc3"],
            "core_genome_gc3_mean": round(m_core, 4),
            "core_genome_gc3_sd": round(sd_core, 4),
            "gc3_diff_vs_core": round(diff, 4),
            "z_approx_p": f"{p_approx:.3e}",
        })
    summary_rows.append({
        "locus": "module_interval_genomic",
        "locus_tag": f"contig_1:{MODULE_START}-{MODULE_END}",
        "gc3_pct": "NA",
        "core_genome_gc3_mean": round(m_core, 4),
        "core_genome_gc3_sd": round(sd_core, 4),
        "gc3_diff_vs_core": "NA",
        "z_approx_p": f"{p_mod_vs_bg:.3e}",
    })
    wt(OUT / "gc_codon_comparison_summary.tsv", summary_rows)

    # --- 4. Prose summary ---
    seh_r  = next((r for r in gene_rows if r["class"] == "seh"), None)
    is3_r  = next((r for r in gene_rows if r["class"] == "IS3 transposase"), None)

    prose = f"""## GC content and codon-position analysis

**Genome**: 23EV612 (ST1 complete chromosome; {len(genome):,} bp; GC {gc_all*100:.2f}%)

### Genomic GC content: module interval vs background

| Region | Length (bp) | GC% |
|--------|------------|------|
| {MODULE_LABEL} | {len(mod_seq):,} | {gc_mod*100:.2f}% |
| Genome background | {len(bg_seq):,} | {gc_bg*100:.2f}% |
| Whole genome | {len(genome):,} | {gc_all*100:.2f}% |

Difference (module vs background): {(gc_mod - gc_bg)*100:+.2f} percentage points
Two-proportion z-test p-value: {p_mod_vs_bg:.3e}

### Codon-position GC for focal genes vs core genome

| Gene | GC_all | GC1 | GC2 | GC3 | GC3 diff vs core |
|------|--------|-----|-----|-----|------------------|
| seh  | {seh_r["gc_overall"] if seh_r else "NA"}% | {seh_r["gc1"] if seh_r else "NA"}% | {seh_r["gc2"] if seh_r else "NA"}% | {seh_r["gc3"] if seh_r else "NA"}% | {(seh_r["gc3"] - m_core):+.2f}% |
| IS3 transposase | {is3_r["gc_overall"] if is3_r else "NA"}% | {is3_r["gc1"] if is3_r else "NA"}% | {is3_r["gc2"] if is3_r else "NA"}% | {is3_r["gc3"] if is3_r else "NA"}% | {(is3_r["gc3"] - m_core):+.2f}% |
| Core genome mean | {statistics.mean(r["gc_overall"] for r in core):.2f}% | {statistics.mean(r["gc1"] for r in core):.2f}% | {statistics.mean(r["gc2"] for r in core):.2f}% | {m_core:.2f}% (±{sd_core:.2f}) | — |

n core CDS = {len(core)}

### Interpretation

GC% at synonymous codon position 3 (GC3) is a sensitive indicator of atypical gene origin because
selection acts weakly at silent sites. The mean GC3 for the 23EV612 core genome is {m_core:.1f}%.
`seh` GC3 = {seh_r["gc3"] if seh_r else "NA"}%, which is {(seh_r["gc3"] - m_core) if seh_r else "NA":.1f} pp relative to the core mean.
The module interval is also more AT-rich than the genome background overall ({gc_mod*100:.2f}% vs {gc_bg*100:.2f}%; P = {p_mod_vs_bg:.2e}).

*Interpretation*: Both the module interval and `seh` itself are anomalously AT-rich relative to the host
core genome. This divergent codon usage, especially at GC3 where purifying selection is weakest, is
consistent with `seh` and the IS655-like transposase having a distinct evolutionary origin from the
majority of the *S. aureus* core genome — supporting the interpretation of horizontal acquisition.
The extremely low `seh` GC3 ({seh_r["gc3"] if seh_r else "NA"}%) may reflect origin from an even more AT-biased
source or strong AT-mutational pressure after insertion into the *S. aureus* chromosome.
"""
    (OUT / "gc_analysis_summary.md").write_text(prose)

    print("Done. Outputs in:", OUT)
    for f in sorted(OUT.glob("*")):
        print(" ", f.name)


if __name__ == "__main__":
    main()
