#!/usr/bin/env python3
"""Build the locus-state/haplotype matrix for the joint ST1-ST2187 tree."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "results/seh_st1_bovine/phylogenomics_set/assemblies"
OUT = ROOT / "results/seh_st1_bovine/st1_st2187_phylogenomics/locus_state_matrix.tsv"


def rows_by_id(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        return {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def main() -> None:
    coding = rows_by_id(ROOT / "results/seh_st1_bovine/seh_integrity/seh_coding_integrity.tsv")
    upstream = rows_by_id(ROOT / "results/seh_st1_bovine/seh_integrity/seh_upstream_regulatory_audit.tsv")
    element = rows_by_id(ROOT / "results/seh_st1_bovine/is3_boundaries/is3_candidate_element_identity_vs_23EV612.tsv")
    locus = rows_by_id(ROOT / "results/seh_st1_bovine/seh_locus_summary.tsv")
    public_negatives = [
        "GCA_031889225.1", "GCA_031889265.1", "GCA_031889285.1",
        "GCA_031889365.1", "GCA_031889545.1",
    ]
    samples = sorted({path.stem for path in PANEL.glob("*.fna")} | set(public_negatives))

    raw_haplotypes = {}
    for sample in samples:
        if sample in element:
            raw_haplotypes[sample] = (
                element[sample]["variant_positions_vs_23EV612"] or "ref",
                coding[sample]["nt_variant_positions"] or "ref",
                upstream[sample]["upstream_variant_positions"] or "ref",
            )
    haplotype_counts = Counter(raw_haplotypes.values())
    ordered = sorted(haplotype_counts, key=lambda h: (-haplotype_counts[h], h))
    haplotype_names = {haplotype: f"H{index}" for index, haplotype in enumerate(ordered, start=1)}

    rows = []
    for sample in samples:
        st = "ST2187" if sample == "cepa_18_barcode02" or sample in public_negatives else "ST1"
        if sample in public_negatives:
            state = "module_absent"
            seh_state = "absent"
            element_state = "no_full_element"
            completeness = "resolved_negative"
        elif sample in element:
            state = "complete_module"
            seh_state = "intact"
            element_state = "complete_degenerated_remnant"
            completeness = "complete"
        else:
            state = "seh_present_locus_at_contig_edge"
            seh_state = "intact" if locus.get(sample, {}).get("seh_intact_orf") == "TRUE" else "present"
            element_state = "unresolved"
            completeness = "unresolved_edge"
        haplotype = raw_haplotypes.get(sample)
        rows.append({
            "sample_id": sample,
            "st": st,
            "seh_state": seh_state,
            "is655_like_state": element_state,
            "locus_state": state,
            "locus_completeness": completeness,
            "module_haplotype": haplotype_names.get(haplotype, "NA"),
            "element_variants_vs_23EV612": haplotype[0] if haplotype else "NA",
            "seh_variants_vs_23EV612": haplotype[1] if haplotype else "NA",
            "upstream_variants_vs_23EV612": haplotype[2] if haplotype else "NA",
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT} ({len(rows)} samples; {len(haplotype_names)} complete-module haplotypes)")
    for haplotype, name in haplotype_names.items():
        print(name, haplotype_counts[haplotype], haplotype)


if __name__ == "__main__":
    main()
