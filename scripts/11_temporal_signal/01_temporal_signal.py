#!/usr/bin/env python3
"""Audit root-to-tip temporal signal before attempting molecular dating."""

from __future__ import annotations

import csv
import datetime as dt
import re
from pathlib import Path

import dendropy
import numpy as np
from scipy.stats import linregress


ROOT = Path(__file__).resolve().parents[1]
TREE = ROOT / "results/seh_st1_bovine/phylogenomics_set/recombination_gubbins/run4_escalated_10iter/seh_st1_gubbins_10iter.final_tree.tre"
META = ROOT / "results/seh_st1_bovine/phylogenomics_set/seh_st1_phylogenomics_set.tsv"
OUT = ROOT / "results/seh_st1_bovine/temporal_signal_audit"
N_RANDOM = 10_000
N_BOOT = 10_000
SEED = 2187


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    years: dict[str, int] = {}
    missing: list[str] = []
    for row in csv.DictReader(META.open(), delimiter="\t"):
        match = re.search(r"(?:19|20)\d{2}", row["year"])
        if row["sample_id"] == "cepa_18_barcode02":
            continue  # ST2187 outgroup/recipient, not part of the ST1 clock regression.
        if match:
            years[row["sample_id"]] = int(match.group())
        else:
            missing.append(row["sample_id"])

    tree = dendropy.Tree.get(path=str(TREE), schema="newick", preserve_underscores=True)
    tree.is_rooted = True  # Gubbins topology is seeded on the ST2187 outgroup split.
    st1_root = tree.mrca(taxon_labels=list(years))
    root_distance = st1_root.distance_from_root()
    points = []
    for sample, year in years.items():
        node = tree.find_node_with_taxon_label(sample)
        if node is None:
            raise SystemExit(f"dated tip missing from tree: {sample}")
        points.append((sample, year, node.distance_from_root() - root_distance))

    x = np.asarray([p[1] for p in points], dtype=float)
    y = np.asarray([p[2] for p in points], dtype=float)
    fit = linregress(x, y)
    root_year = -fit.intercept / fit.slope
    rng = np.random.default_rng(SEED)

    random_slopes = np.empty(N_RANDOM)
    random_r2 = np.empty(N_RANDOM)
    for i in range(N_RANDOM):
        shuffled = rng.permutation(x)
        trial = linregress(shuffled, y)
        random_slopes[i] = trial.slope
        random_r2[i] = trial.rvalue**2
    permutation_p = (1 + np.count_nonzero(random_slopes >= fit.slope)) / (N_RANDOM + 1)
    permutation_r2_p = (1 + np.count_nonzero(random_r2 >= fit.rvalue**2)) / (N_RANDOM + 1)

    boot_roots = []
    boot_slopes = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(x[idx])) < 2:
            continue
        trial = linregress(x[idx], y[idx])
        if trial.slope > 0:
            boot_slopes.append(trial.slope)
            boot_roots.append(-trial.intercept / trial.slope)
    slope_ci = np.quantile(boot_slopes, [0.025, 0.975])
    root_ci = np.quantile(boot_roots, [0.025, 0.975])

    tip_rows = []
    for sample, year, distance in points:
        predicted = fit.intercept + fit.slope * year
        tip_rows.append({
            "sample_id": sample, "collection_year": year,
            "root_to_tip_branch_length": f"{distance:.6f}",
            "fitted_branch_length": f"{predicted:.6f}",
            "residual": f"{distance - predicted:.6f}",
        })
    write_tsv(OUT / "root_to_tip_points.tsv", tip_rows)
    write_tsv(OUT / "temporal_signal_summary.tsv", [{
        "dated_ST1_tips": len(points),
        "undated_ST1_tips": len(missing),
        "sampling_year_min": int(x.min()),
        "sampling_year_max": int(x.max()),
        "sampling_span_years": int(x.max() - x.min()),
        "slope_branch_units_per_year": f"{fit.slope:.6f}",
        "slope_bootstrap_95pct_CI": f"{slope_ci[0]:.6f}-{slope_ci[1]:.6f}",
        "pearson_r": f"{fit.rvalue:.6f}",
        "r_squared": f"{fit.rvalue**2:.6f}",
        "linear_regression_p": f"{fit.pvalue:.6g}",
        "date_randomization_slope_p": f"{permutation_p:.6g}",
        "date_randomization_r2_p": f"{permutation_r2_p:.6g}",
        "exploratory_ST1_root_year": f"{root_year:.1f}",
        "tip_bootstrap_root_year_95pct_CI": f"{root_ci[0]:.1f}-{root_ci[1]:.1f}",
        "dating_decision": "temporal signal detectable but insertion age not identifiable",
        "reason": "all sampled ST1 are module-positive; clock can at most date the sampled ST1 ancestor, not distinguish acquisition before versus within its unsampled stem",
    }])
    (OUT / "analysis.log").write_text("\n".join([
        f"run_date={dt.date.today().isoformat()}",
        "script=scripts/135_audit_temporal_signal_st1.py",
        f"tree={TREE.relative_to(ROOT)}",
        "tree_type=Gubbins recombination-filtered final tree; branch lengths in inferred substitutions/SNP units",
        "rooting=ST2187 outgroup split; regression root is MRCA of dated ST1 tips",
        f"dated_tips={len(points)}",
        f"undated_tips={','.join(missing)}",
        f"date_randomizations={N_RANDOM}",
        f"tip_bootstrap_replicates={N_BOOT}",
        "caveat=tip bootstrap does not model topology, clock, recombination-filter, or sampling-structure uncertainty",
        "interpretation=exploratory root date is not an insertion date and must not be reported as one",
    ]) + "\n")


if __name__ == "__main__":
    main()
