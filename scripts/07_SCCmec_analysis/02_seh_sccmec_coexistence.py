#!/usr/bin/env python3
"""Screen the local ST1 panel for coexistence of SCCmec markers and seh--IS655."""

from __future__ import annotations

import csv
import datetime as dt
import os
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/seh_st1_bovine/st1_mrsa_seh_sccmec"
ASSEMBLIES = ROOT / "results/seh_st1_bovine/phylogenomics_set/assemblies"
GFFS = ROOT / "results/seh_st1_bovine/phylogenomics_set/gff3"
MATRIX = ROOT / "results/seh_st1_bovine/st1_st2187_phylogenomics/locus_state_matrix.tsv"
MW2_GBFF = ROOT / "resources/mw2_ba000033/BA000033.2.gbff"
MW2_FASTA = ROOT / "resources/mw2_ba000033/BA000033.2.fna"
ELEMENTS = ROOT / "results/seh_st1_bovine/is3_boundaries/is3_candidate_elements.fasta"
MLST = ROOT / ".micromamba/envs/gridion_mlst/bin/mlst"


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def mw2_queries() -> dict[str, str]:
    record = SeqIO.read(MW2_GBFF, "genbank")
    queries = {}
    for f in record.features:
        if f.type != "CDS": continue
        q = f.qualifiers
        gene = (q.get("gene") or [""])[0]
        product = (q.get("product") or [""])[0].lower()
        key = gene if gene in {"mecA", "ccrA", "ccrB"} else ("seh" if "enterotoxin h" in product else "")
        if key:
            queries[key] = str(f.extract(record.seq)).upper()
    element = next(r for r in SeqIO.parse(ELEMENTS, "fasta") if r.id.startswith("23EV612|"))
    queries["is655_like"] = str(element.seq).upper()
    return queries


def blast(query: Path, subject: Path) -> dict[str, object] | None:
    fmt = "6 sseqid pident length qlen qstart qend sstart send evalue bitscore"
    run = subprocess.run(["blastn", "-query", str(query), "-subject", str(subject), "-dust", "no", "-outfmt", fmt], text=True, capture_output=True, check=True)
    hits = []
    for line in run.stdout.splitlines():
        x = line.split("\t")
        hits.append({"contig": x[0], "identity": float(x[1]), "length": int(x[2]), "qlen": int(x[3]), "qstart": int(x[4]), "qend": int(x[5]), "sstart": int(x[6]), "send": int(x[7]), "bitscore": float(x[9])})
    return max(hits, key=lambda h: (h["length"], h["identity"], h["bitscore"])) if hits else None


def present(hit: dict[str, object] | None, cov: float = 90, ident: float = 95) -> bool:
    return bool(hit and 100 * int(hit["length"]) / int(hit["qlen"]) >= cov and float(hit["identity"]) >= ident)


def gff_anchors(path: Path) -> tuple[str, int, str]:
    anchors = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"): continue
            x = line.rstrip().split("\t")
            if len(x) != 9 or x[2] != "CDS": continue
            attrs = dict(p.split("=", 1) for p in x[8].split(";") if "=" in p)
            gene = attrs.get("gene", "")
            if gene in {"cstA", "rlmH"}: anchors[gene] = (x[0], int(x[3]), int(x[4]))
    if set(anchors) != {"cstA", "rlmH"}: return "unresolved", 0, ""
    if anchors["cstA"][0] != anchors["rlmH"][0]: return "different_contigs", 0, ""
    lo = min(anchors["cstA"][1], anchors["rlmH"][1]); hi = max(anchors["cstA"][2], anchors["rlmH"][2])
    return "resolved", hi - lo + 1, anchors["cstA"][0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    query_seqs = mw2_queries()
    query_paths = {}
    for key, seq in query_seqs.items():
        path = OUT / f"query_{key}.fasta"; path.write_text(f">MW2_{key}\n{seq}\n"); query_paths[key] = path

    states = {r["sample_id"]: r for r in tsv(MATRIX) if r["st"] == "ST1"}
    genomes = [(p.stem, p, GFFS / f"{p.stem}.gff3", "local_ST1_panel") for p in sorted(ASSEMBLIES.glob("*.fna")) if p.stem in states]
    genomes.append(("MW2_BA000033.2", MW2_FASTA, Path(""), "reference_ST1_MRSA_IVa"))
    rows = []
    for sample, fasta, gff, source in genomes:
        hits = {key: blast(path, fasta) for key, path in query_paths.items()}
        calls = {key: present(hit, 95 if key == "is655_like" else 90, 95) for key, hit in hits.items()}
        seh_hit, is_hit = hits["seh"], hits["is655_like"]
        adjacent = bool(calls["seh"] and calls["is655_like"] and seh_hit["contig"] == is_hit["contig"] and max(min(int(seh_hit["sstart"]), int(seh_hit["send"])), min(int(is_hit["sstart"]), int(is_hit["send"]))) - min(max(int(seh_hit["sstart"]), int(seh_hit["send"])), max(int(is_hit["sstart"]), int(is_hit["send"]))) - 1 <= 1000)
        if sample.startswith("MW2_"):
            interval_status, interval_len, interval_contig = "not_measured_legacy_annotation", 0, "BA000033.2"
            mlst_call = "ST1 (published reference)"
        else:
            interval_status, interval_len, interval_contig = gff_anchors(gff)
            mlst_env = os.environ.copy()
            mlst_env["PATH"] = f"{MLST.parent}:{mlst_env.get('PATH', '')}"
            mlst_run = subprocess.run([str(MLST), str(fasta)], text=True, capture_output=True, check=True, env=mlst_env)
            mlst_call = mlst_run.stdout.rstrip().split("\t")[2] if mlst_run.stdout.strip() else "unresolved"
        mrsa = calls["mecA"]
        ccr = calls["ccrA"] and calls["ccrB"]
        scc = "MW2-like ccrAB+mec complex; subtype not formally assigned" if mrsa and ccr else ("mecA-positive; ccrAB not MW2-like" if mrsa else "no SCCmec inferred")
        locus = "intact_adjacent_seh_IS655" if calls["seh"] and calls["is655_like"] and adjacent else ("seh_only_or_disrupted" if calls["seh"] else "module_absent")
        row = {"sample": sample, "source": source, "mlst_call": mlst_call, "mecA_present": mrsa, "ccrA_present": calls["ccrA"], "ccrB_present": calls["ccrB"], "sccmec_marker_interpretation": scc, "seh_present": calls["seh"], "is655_like_present": calls["is655_like"], "seh_is655_adjacent": adjacent, "locus_state": locus, "anchor_status": interval_status, "cstA_rlmH_interval_bp": interval_len or "", "anchor_contig": interval_contig}
        for key, hit in hits.items():
            row[f"{key}_identity"] = f"{hit['identity']:.3f}" if hit else ""
            row[f"{key}_coverage_pct"] = f"{min(100.0, 100*int(hit['length'])/int(hit['qlen'])):.2f}" if hit else ""
        rows.append(row)

    inventory_fields = ["sample", "source", "mlst_call", "mecA_present", "ccrA_present", "ccrB_present", "sccmec_marker_interpretation"]
    presence_fields = ["sample", "mecA_present", "seh_present", "seh_identity", "seh_coverage_pct", "is655_like_present", "is655_like_identity", "is655_like_coverage_pct", "seh_is655_adjacent", "locus_state"]
    interval_fields = ["sample", "mecA_present", "anchor_status", "anchor_contig", "cstA_rlmH_interval_bp", "locus_state"]
    write_tsv(OUT / "st1_mrsa_genomes_inventory.tsv", rows, inventory_fields)
    write_tsv(OUT / "seh_is655_presence_in_st1_mrsa.tsv", rows, presence_fields)
    write_tsv(OUT / "cstA_rlmH_intervals_st1_mrsa.tsv", rows, interval_fields)
    write_tsv(OUT / "st1_mrsa_locus_summary.tsv", [r for r in rows if r["mecA_present"]], list(rows[0]))

    n_mrsa = sum(bool(r["mecA_present"]) for r in rows); n_coexist = sum(bool(r["mecA_present"] and r["locus_state"] == "intact_adjacent_seh_IS655") for r in rows)
    versions = subprocess.run(["blastn", "-version"], text=True, capture_output=True, check=True).stdout.splitlines()[0]
    (OUT / "analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}", "script=scripts/122_test_seh_sccmec_coexistence_st1_mrsa.py", f"blastn={versions}", "mlst=2.33.1", "presence_thresholds=seh/mecA/ccr >=90% coverage and >=95% identity; IS655-like >=95% coverage and >=95% identity", "SCCmec typing limitation=MW2 ccrA/ccrB and mecA marker matching is not formal SCCmec subtyping", f"genomes_screened={len(rows)}", f"mecA_positive={n_mrsa}", f"mecA_positive_with_intact_module={n_coexist}"]) + "\n")


if __name__ == "__main__": main()
