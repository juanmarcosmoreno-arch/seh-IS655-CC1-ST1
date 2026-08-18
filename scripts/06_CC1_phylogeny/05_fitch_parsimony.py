#!/usr/bin/env python3
"""Re-run Fitch parsimony coding the 5 partial loci as unknown (?).

Partials: complete_module=False AND is655_like_present=True.
These tips have IS655-like detectable but seh or spacer unrecoverable at a
contig edge; their state is ambiguous (assembly fragmentation, not confirmed
biological absence).  Coding them as {0,1} prevents the Fitch algorithm from
inferring spurious tip-level loss events.

True absences (complete_module=False AND is655_like_present=False) remain
coded as 0 (absent), as expected for ST1660/ST772 genomes.
"""

from __future__ import annotations
import csv
from pathlib import Path
import dendropy

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "results/seh_st1_bovine/cc1_phylogeny"


def wt(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


def main() -> None:
    pres_rows = list(csv.DictReader((OUT / "cc1_module_presence.tsv").open(), delimiter="\t"))

    # Build character state per tip
    # True → 1 (present)
    # False + is655_like=False → 0 (absent, real)
    # False + is655_like=True  → None (unknown, partial — assembly fragmentation)
    states: dict[str, bool | None] = {}
    partial_tips: list[str] = []
    for r in pres_rows:
        if r["complete_module"] == "True":
            states[r["sample_id"]] = True
        elif r["is655_like_present"] == "True":
            states[r["sample_id"]] = None   # unknown
            partial_tips.append(r["sample_id"])
        else:
            states[r["sample_id"]] = False  # true absence

    print(f"Partial tips coded as unknown (?): {len(partial_tips)}")
    for t in partial_tips:
        print(f"  {t}")

    tree = dendropy.Tree.get(path=str(OUT / "cc1_phylogeny.nwk"),
                             schema="newick", preserve_underscores=True)
    tree.reroot_at_midpoint(update_bipartitions=False)

    # Fitch pass 1: bottom-up Fitch sets
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            s = states.get(node.taxon.label)
            if s is None:
                node._fitch = {True, False}  # unknown: both states allowed
            else:
                node._fitch = {s}
        else:
            kids = [c._fitch for c in node.child_node_iter()]
            inter = set.intersection(*kids)
            node._fitch = inter if inter else set.union(*kids)

    # Fitch pass 2: top-down state assignment + event counting
    root = tree.seed_node
    root._state = False if False in root._fitch else next(iter(root._fitch))
    rows = []; gains = losses = 0
    for i, node in enumerate(tree.preorder_node_iter()):
        if node is not root:
            parent = node.parent_node._state
            node._state = parent if parent in node._fitch else next(iter(node._fitch))
            gains  += int(not parent and node._state)
            losses += int(parent and not node._state)
        label   = node.taxon.label if node.is_leaf() else f"node_{i}"
        support = str(node.label or "")
        partial = label in partial_tips
        event   = ("root" if node is root
                   else "gain" if not node.parent_node._state and node._state
                   else "loss" if node.parent_node._state and not node._state
                   else "none")
        rows.append({
            "node":                   label,
            "is_tip":                 node.is_leaf(),
            "coded_as_unknown":       partial,
            "inferred_module_state":  int(node._state),
            "fitch_set":              "|".join(str(int(x)) for x in sorted(node._fitch)),
            "node_support":           support,
            "event_from_parent":      event,
        })

    wt(OUT / "cc1_ancestral_state_reconstruction_v2.tsv", rows)
    tree.write(path=str(OUT / "cc1_phylogeny_rooted_v2.nwk"),
               schema="newick", suppress_rooting=False)
    with (OUT / "ancestral_state_summary_v2.tsv").open("w") as h:
        h.write(
            "metric\tvalue\n"
            "rooting\tmidpoint (CC1660 non-monophyletic)\n"
            f"partial_tips_coded_unknown\t{len(partial_tips)}\n"
            f"parsimonious_gains\t{gains}\n"
            f"parsimonious_losses\t{losses}\n"
            "note\tPartials (is655_like present but module incomplete) coded as unknown; "
            "losses from v1 reconstruction were entirely tip-level events in partial genomes\n"
        )

    print(f"\nResult: {gains} gain(s), {losses} loss(es)  [partials coded as unknown]")
    print("Outputs written to", OUT)


if __name__ == "__main__":
    main()
