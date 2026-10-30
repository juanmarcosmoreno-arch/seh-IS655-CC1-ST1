#!/usr/bin/env python3
"""Sensitivity analysis for the Fitch ancestral-state reconstruction (Major comment 3).

Re-runs the CC1 module presence/absence Fitch reconstruction under every
combination of:
  - partial-tip coding: unknown (?), absent (0), present (1)
  - rooting: midpoint (as in the submitted manuscript), or outgroup
    (ST1660 clade, results/seh_st1_bovine/cc1_phylogeny/outgroup_tips.txt)

The 5 "partial" tips are complete_module=False AND is655_like_present=True
(IS655-like detected but seh/spacer unrecoverable at a contig edge).

Does not overwrite any existing manuscript output; everything is written to
results/seh_st1_bovine/cc1_phylogeny/sensitivity/.
"""

from __future__ import annotations
import csv
from pathlib import Path
import dendropy

ROOT = Path(__file__).resolve().parents[1]
PHYLO = ROOT / "results/seh_st1_bovine/cc1_phylogeny"
OUT = PHYLO / "sensitivity"
OUT.mkdir(exist_ok=True)

CODINGS = ["unknown", "absent", "present"]
ROOTINGS = ["midpoint", "outgroup"]


def load_states(coding: str) -> tuple[dict[str, bool | None], list[str]]:
    pres_rows = list(csv.DictReader((PHYLO / "cc1_module_presence.tsv").open(), delimiter="\t"))
    states: dict[str, bool | None] = {}
    partial_tips: list[str] = []
    for r in pres_rows:
        if r["complete_module"] == "True":
            states[r["sample_id"]] = True
        elif r["is655_like_present"] == "True":
            partial_tips.append(r["sample_id"])
            if coding == "unknown":
                states[r["sample_id"]] = None
            elif coding == "absent":
                states[r["sample_id"]] = False
            elif coding == "present":
                states[r["sample_id"]] = True
        else:
            states[r["sample_id"]] = False
    return states, partial_tips


def load_tree(rooting: str) -> dendropy.Tree:
    tree = dendropy.Tree.get(path=str(PHYLO / "cc1_phylogeny.nwk"),
                              schema="newick", preserve_underscores=True)
    if rooting == "midpoint":
        tree.reroot_at_midpoint(update_bipartitions=False)
    elif rooting == "outgroup":
        outgroup_labels = (PHYLO / "outgroup_tips.txt").read_text().strip().split(",")
        outgroup_taxa = [t for t in tree.taxon_namespace if t.label in outgroup_labels]
        missing = set(outgroup_labels) - {t.label for t in outgroup_taxa}
        if missing:
            raise SystemExit(f"Outgroup tips not found in tree: {missing}")
        # The raw IQ-TREE newick has a multifurcating (unrooted) seed node, which
        # breaks dendropy's mrca/to_outgroup_position when the outgroup spans more
        # than one of its immediate children. Root arbitrarily first (midpoint) so
        # the tree is a proper bifurcating rooted structure, then re-root on the
        # outgroup MRCA within that structure.
        tree.reroot_at_midpoint(update_bipartitions=False)
        mrca = tree.mrca(taxa=outgroup_taxa)
        if mrca is tree.seed_node:
            raise SystemExit("Outgroup (ST1660) is not monophyletic under this topology")
        tree.to_outgroup_position(mrca, update_bipartitions=False)
    return tree


def run_fitch(tree: dendropy.Tree, states: dict[str, bool | None]) -> tuple[int, int]:
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            s = states.get(node.taxon.label)
            node._fitch = {True, False} if s is None else {s}
        else:
            kids = [c._fitch for c in node.child_node_iter()]
            inter = set.intersection(*kids)
            node._fitch = inter if inter else set.union(*kids)

    root = tree.seed_node
    root._state = False if False in root._fitch else next(iter(root._fitch))
    gains = losses = 0
    for node in tree.preorder_node_iter():
        if node is root:
            continue
        parent = node.parent_node._state
        node._state = parent if parent in node._fitch else next(iter(node._fitch))
        gains += int(not parent and node._state)
        losses += int(parent and not node._state)
    return gains, losses


def main() -> None:
    results = []
    for rooting in ROOTINGS:
        for coding in CODINGS:
            states, partial_tips = load_states(coding)
            tree = load_tree(rooting)
            gains, losses = run_fitch(tree, states)
            results.append({
                "rooting": rooting,
                "partial_coding": coding,
                "n_partial_tips": len(partial_tips),
                "parsimonious_gains": gains,
                "parsimonious_losses": losses,
            })
            print(f"rooting={rooting:9s} coding={coding:8s} -> gains={gains} losses={losses}")

    out_path = OUT / "fitch_sensitivity_summary.tsv"
    with out_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(results[0]), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print("\nWritten:", out_path)


if __name__ == "__main__":
    main()
