#!/usr/bin/env python3
"""Prepare/finalize a resolved CC1 phylogeny panel and map the seh module."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/seh_st1_bovine/cc1_phylogeny"
JAHNEN = ROOT / "resources/jahnen_2025_cc1_cc1660/extracted/ncbi_dataset/data"
ST772 = ROOT / "resources/st772_cc1_comparators/ncbi_dataset/ncbi_dataset/data"
SURVEY = ROOT / "results/seh_st1_bovine/cc1_st_survey/cc1_st_survey_seh_is655_presence.tsv"
MLST = ROOT / ".micromamba/envs/gridion_mlst/bin/mlst"
BLASTN = ROOT / ".micromamba/envs/gridion_mlst/bin/blastn"
ANY2FASTA_BIN = Path("/opt/homebrew/Cellar/micromamba/2.6.2/pkgs/https/conda.anaconda.org/bioconda/noarch/any2fasta-0.8.1-hdfd78af_0/bin")
REFERENCE = ROOT / "results/voss_2025_msphere_st1/reference/CP160024_23EV612.fna"
MW2 = ROOT / "results/st1_negative_search/complete_chromosome_candidate_screen/extracted/ncbi_dataset/data/GCF_050908105.1/GCF_050908105.1_MW2_P0_assembly_genomic.fna"


def wt(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


def accession(path: Path) -> str:
    m = re.search(r"GC[AF]_\d+\.\d+", str(path))
    return m.group() if m else path.stem


def best_hit(query: Path, subject: Path) -> tuple[float, float]:
    run = subprocess.run([str(BLASTN), "-query", str(query), "-subject", str(subject),
        "-outfmt", "6 pident length qlen bitscore", "-max_hsps", "20"],
        text=True, capture_output=True, check=True)
    best = (0.0, 0.0, 0.0)
    for line in run.stdout.splitlines():
        pid, length, qlen, bits = map(float, line.split())
        candidate = (bits, pid, 100 * length / qlen)
        if candidate > best: best = candidate
    return best[1], best[2]


def prepare() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prior = {r["assembly"]: r for r in csv.DictReader(SURVEY.open(), delimiter="\t")}
    paths = list(JAHNEN.glob("GC*/*genomic.fna")) + list(ST772.glob("GC*/*genomic.fna"))
    paths += [REFERENCE, MW2]
    seen: set[str] = set(); rows = []
    qseh = ROOT / "results/seh_st1_bovine/cc1_st_survey/query_seh.fasta"
    qis = ROOT / "results/seh_st1_bovine/cc1_st_survey/query_is655_like.fasta"
    qsp = ROOT / "results/seh_st1_bovine/cc1_st_survey/query_spacer.fasta"
    for path in sorted(paths, key=lambda x: str(x)):
        acc = "23EV612" if path == REFERENCE else accession(path)
        if acc in seen: continue
        seen.add(acc)
        if acc in prior:
            st = prior[acc]["st"]
            seh = prior[acc]["seh_present"] == "True"
            is655 = prior[acc]["is655_like_present"] == "True"
            spacer = prior[acc]["spacer_present"] == "True"
            module = prior[acc]["complete_adjacent_module"] == "True"
            cohort = "Jahnen_2025"
        else:
            env = os.environ.copy()
            env["PATH"] = f"{ANY2FASTA_BIN}:{MLST.parent}:{env.get('PATH','')}"
            call = subprocess.run([str(MLST), str(path)], text=True, capture_output=True, check=True, env=env).stdout.rstrip().split("\t")
            st = call[2] if len(call) > 2 else "-"
            sh = best_hit(qseh, path); ih = best_hit(qis, path); ph = best_hit(qsp, path)
            seh = sh[0] >= 95 and sh[1] >= 90
            is655 = ih[0] >= 95 and ih[1] >= 95
            spacer = ph[0] >= 95 and ph[1] >= 90
            module = seh and is655 and spacer
            cohort = "ST772_historical" if "st772_cc1" in str(path) else "reference"
        contigs = sum(1 for line in path.open(errors="ignore") if line.startswith(">"))
        rows.append({"sample_id": acc, "cohort": cohort, "st": st,
            "assembly_level_proxy": "complete_chromosome" if contigs == 1 else "draft",
            "contigs": contigs, "seh_present": seh, "is655_like_present": is655,
            "spacer_present": spacer, "complete_module": module, "fasta": str(path)})
    wt(OUT / "cc1_genome_panel.tsv", rows)
    wt(OUT / "cc1_module_presence.tsv", [{k:r[k] for k in ("sample_id","cohort","st","seh_present","is655_like_present","spacer_present","complete_module")} for r in rows])
    with (OUT / "ska_input.tsv").open("w") as h:
        for r in rows: h.write(f"{r['sample_id']}\t{r['fasta']}\n")
    (OUT / "outgroup_tips.txt").write_text(",".join(r["sample_id"] for r in rows if r["st"] == "1660") + "\n")
    counts = Counter(r["st"] for r in rows)
    (OUT / "panel_summary.tsv").write_text("metric\tvalue\n" +
        f"genomes_total\t{len(rows)}\ncomplete_chromosomes\t{sum(r['contigs']==1 for r in rows)}\n" +
        f"drafts\t{sum(r['contigs']!=1 for r in rows)}\nST_distribution\t" + ";".join(f"{k}:{v}" for k,v in sorted(counts.items())) + "\n")


def finalize() -> None:
    import dendropy
    tree_path = OUT / "cc1_phylogeny.nwk"
    tree = dendropy.Tree.get(path=str(tree_path), schema="newick", preserve_underscores=True)
    states = {r["sample_id"]: r["complete_module"] == "True" for r in csv.DictReader((OUT/"cc1_module_presence.tsv").open(), delimiter="\t")}
    # CC1660 tips are not monophyletic in the ML topology, so collective
    # outgroup rooting would force a nonexistent split. Use midpoint rooting
    # for the explicit parsimony audit and report that limitation.
    tree.reroot_at_midpoint(update_bipartitions=False)
    # Fitch sets, then deterministic minimum-change assignment preferring parent state.
    for node in tree.postorder_node_iter():
        if node.is_leaf(): node._fitch = {states.get(node.taxon.label, False)}
        else:
            kids = [c._fitch for c in node.child_node_iter()]
            inter = set.intersection(*kids); node._fitch = inter if inter else set.union(*kids)
    root = tree.seed_node; root._state = False if False in root._fitch else next(iter(root._fitch))
    rows=[]; gains=losses=0
    for i,node in enumerate(tree.preorder_node_iter()):
        if node is not root:
            parent=node.parent_node._state; node._state=parent if parent in node._fitch else next(iter(node._fitch))
            gains += (not parent and node._state); losses += (parent and not node._state)
        label=node.taxon.label if node.is_leaf() else f"node_{i}"
        support=str(node.label or "")
        rows.append({"node":label,"is_tip":node.is_leaf(),"inferred_module_state":int(node._state),"fitch_set":"|".join(str(int(x)) for x in sorted(node._fitch)),"node_support":support,"event_from_parent":"root" if node is root else ("gain" if not node.parent_node._state and node._state else "loss" if node.parent_node._state and not node._state else "none")})
    wt(OUT/"cc1_ancestral_state_reconstruction.tsv",rows)
    tree.write(path=str(OUT/"cc1_phylogeny_rooted.nwk"),schema="newick",suppress_rooting=False)
    with (OUT/"ancestral_state_summary.tsv").open("w") as h:
        h.write(f"metric\tvalue\nrooting\tmidpoint (CC1660 non-monophyletic)\nparsimonious_gains\t{gains}\nparsimonious_losses\t{losses}\n")


if __name__ == "__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("mode",choices=["prepare","finalize"]); args=ap.parse_args()
    prepare() if args.mode=="prepare" else finalize()
