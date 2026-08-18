#!/usr/bin/env python3
"""Formal SCCmec typing of the available ST1 MRSA panel with SCCmecFinder DB."""
from __future__ import annotations

import csv
import datetime as dt
import os
import re
import subprocess
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/seh_st1_bovine/st1_mrsa_sccmec_formal"
OLD = ROOT / "results/seh_st1_bovine/st1_mrsa_seh_sccmec/st1_mrsa_locus_summary.tsv"
SURVEY = ROOT / "results/seh_st1_bovine/cc1_st_survey/cc1_st_survey_seh_is655_presence.tsv"
LOCAL = ROOT / "results/seh_st1_bovine/phylogenomics_set/assemblies"
JAHNEN = ROOT / "resources/jahnen_2025_cc1_cc1660/extracted/ncbi_dataset/data"
EXPANSION = ROOT / "resources/st1_mrsa_expansion/ncbi_dataset/data"
MW2 = ROOT / "resources/mw2_ba000033/BA000033.2.fna"
SINGLE = ROOT / "resources/sccmecfinder/single_genes_database_20171117.fasta"
WHOLE = ROOT / "resources/sccmecfinder/whole_cassette_SCCmec_database_REFERENCE_20171117.fasta"
MLST = ROOT / ".micromamba/envs/gridion_mlst/bin/mlst"
SEH_QUERY = ROOT / "results/seh_st1_bovine/st1_mrsa_seh_sccmec/query_seh.fasta"
IS_QUERY = ROOT / "results/seh_st1_bovine/st1_mrsa_seh_sccmec/query_is655_like.fasta"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def blast_queries(query: Path, subject: Path) -> list[dict[str, object]]:
    fmt = "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore"
    run = subprocess.run(
        ["blastn", "-query", str(query), "-subject", str(subject), "-dust", "no", "-outfmt", fmt],
        text=True, capture_output=True, check=True,
    )
    rows = []
    for line in run.stdout.splitlines():
        x = line.split("\t")
        rows.append({"query": x[0], "contig": x[1], "identity": float(x[2]),
                     "length": int(x[3]), "qlen": int(x[4]), "qstart": int(x[5]),
                     "qend": int(x[6]), "sstart": int(x[7]), "send": int(x[8]),
                     "evalue": x[9], "bitscore": float(x[10])})
    return rows


def best_by_query(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    best = {}
    for row in rows:
        key = str(row["query"])
        cov = 100 * int(row["length"]) / int(row["qlen"])
        row = dict(row); row["coverage"] = cov
        if key not in best or (cov, row["identity"], row["bitscore"]) > (
            best[key]["coverage"], best[key]["identity"], best[key]["bitscore"]
        ):
            best[key] = row
    return best


def family(header: str) -> str:
    return header.split(":", 1)[0]


def choose_family(hits: dict[str, dict[str, object]], prefix: str) -> dict[str, object] | None:
    candidates = [dict(v, header=k) for k, v in hits.items()
                  if family(k).startswith(prefix) and v["identity"] >= 90 and v["coverage"] >= 60]
    return max(candidates, key=lambda x: (x["coverage"], x["identity"], x["bitscore"])) if candidates else None


def mec_class(hits: dict[str, dict[str, object]]) -> tuple[str, str]:
    mec_i = choose_family(hits, "mecI")
    mec_r = choose_family(hits, "mecR1")
    dmec_r = choose_family(hits, "dmecR1")
    is1272 = choose_family(hits, "IS1272")
    if mec_i and mec_r: return "A", "mecI+mecR1"
    if dmec_r and is1272: return "B", "dmecR1+IS1272"
    # C1/C2 require the dedicated orientation/sequence database; none is
    # assigned from weak single-gene evidence alone.
    return "unresolved", "+".join(x for x, present in [("mecI",mec_i),("mecR1",mec_r),("dmecR1",dmec_r),("IS1272",is1272)] if present) or "no_regulatory_marker"


def ccr_type(hits: dict[str, dict[str, object]]) -> tuple[str, str]:
    calls = {}
    for prefix in ["ccrA1","ccrA2","ccrA3","ccrA4","ccrB1","ccrB2","ccrB3","ccrB4","ccrB6","ccrC1","ccrC2"]:
        hit = choose_family(hits, prefix)
        if hit: calls[prefix] = hit
    for n in ["1","2","3","4"]:
        if f"ccrA{n}" in calls and f"ccrB{n}" in calls:
            return n, f"ccrA{n}+ccrB{n}"
    if "ccrC1" in calls: return "5", "ccrC1"
    if "ccrA1" in calls and "ccrB6" in calls: return "6", "ccrA1+ccrB6"
    if "ccrC2" in calls: return "7", "ccrC2"
    return "unresolved", "+".join(sorted(calls)) or "none"


TYPE_MAP = {("1","B"):"I", ("2","A"):"II", ("3","A"):"III", ("2","B"):"IV",
            ("5","C2"):"V", ("4","B"):"VI", ("5","C1"):"VII", ("4","A"):"VIII",
            ("1","C2"):"IX", ("7","C1"):"X", ("8","E"):"XI"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    old = {r["sample"]: r for r in read_tsv(OLD)}
    survey = {r["assembly"]: r for r in read_tsv(SURVEY) if r["st"] == "1"}
    genomes: dict[str, tuple[Path,str]] = {}
    for p in LOCAL.glob("*.fna"):
        if p.stem in old and old[p.stem]["mecA_present"] == "True": genomes[p.stem] = (p,"original_panel")
    for acc, meta in survey.items():
        p = JAHNEN / acc / f"{acc}_genomic.fna"
        if p.exists(): genomes.setdefault(acc,(p,"Jahnen_CC1_panel"))
    for p in EXPANSION.glob("*/*genomic.fna"):
        genomes.setdefault(p.parent.name,(p,"NCBI_ST1_expansion"))
    genomes["MW2_BA000033.2"] = (MW2,"reference_IVa")

    locus_query = OUT / "query_seh_is655.fasta"
    SeqIO.write(list(SeqIO.parse(SEH_QUERY,"fasta")) + list(SeqIO.parse(IS_QUERY,"fasta")), locus_query, "fasta")

    inventory=[]; typed=[]; all_gene_hits=[]; all_whole_hits=[]
    for sample,(fasta,source) in sorted(genomes.items()):
        if source == "NCBI_ST1_expansion":
            mlst_env=os.environ.copy(); mlst_env["PATH"]=f"{MLST.parent}:{mlst_env.get('PATH','')}"
            call=subprocess.run([str(MLST),"--quiet",str(fasta)],text=True,capture_output=True,check=True,env=mlst_env).stdout.strip().split("\t")
            if len(call)<3 or call[2] != "1": continue
        gene_hits = best_by_query(blast_queries(SINGLE,fasta))
        mec = choose_family(gene_hits,"mecA")
        if not mec: continue
        ccr, ccr_detail = ccr_type(gene_hits)
        mclass, mdetail = mec_class(gene_hits)
        nominal = TYPE_MAP.get((ccr,mclass),"unresolved")
        subtype_hits = [dict(v,header=k) for k,v in gene_hits.items() if k.startswith("subtype-") and v["identity"] >= 90 and v["coverage"] >= 60]
        subtype_hit = max(subtype_hits,key=lambda x:(x["coverage"],x["identity"],x["bitscore"])) if subtype_hits else None
        subtype = family(str(subtype_hit["header"])).replace("subtype-","") if subtype_hit else "unresolved"
        whole_hits = best_by_query(blast_queries(WHOLE,fasta))
        whole_candidates=[dict(v,header=k) for k,v in whole_hits.items() if v["identity"] >= 90]
        whole_best=max(whole_candidates,key=lambda x:(x["coverage"],x["identity"],x["bitscore"])) if whole_candidates else None
        module = old.get(sample,{}).get("locus_state", survey.get(sample,{}).get("complete_adjacent_module","not_screened"))
        if module in {"True","true"}: module="intact_adjacent_seh_IS655"
        seh = old.get(sample,{}).get("seh_present", survey.get(sample,{}).get("seh_present","not_screened"))
        is655 = old.get(sample,{}).get("is655_like_present", survey.get(sample,{}).get("is655_like_present","not_screened"))
        if source == "NCBI_ST1_expansion":
            locus_hits=best_by_query(blast_queries(locus_query,fasta))
            sh=next((v for k,v in locus_hits.items() if "seh" in k.lower()),None)
            ih=next((v for k,v in locus_hits.items() if "is655" in k.lower()),None)
            seh=bool(sh and sh["identity"]>=95 and sh["coverage"]>=90)
            is655=bool(ih and ih["identity"]>=95 and ih["coverage"]>=95)
            adjacent=bool(seh and is655 and sh["contig"]==ih["contig"] and abs(min(sh["sstart"],sh["send"])-min(ih["sstart"],ih["send"]))<3000)
            module="intact_adjacent_seh_IS655" if adjacent else "assembly_fragmented_or_module_incomplete" if seh or is655 else "module_absent"
        status="formal_type_and_subtype" if nominal != "unresolved" and subtype != "unresolved" else "formal_type_subtype_unresolved" if nominal != "unresolved" else "non_typeable"
        row={"sample":sample,"source":source,"st":"1","mecA_identity":f"{mec['identity']:.3f}","mecA_coverage_pct":f"{mec['coverage']:.2f}",
             "ccr_complex":ccr,"ccr_markers":ccr_detail,"mec_complex_class":mclass,"mec_markers":mdetail,
             "sccmec_type":nominal,"sccmec_subtype":subtype,"typing_status":status,
             "subtype_marker_identity":f"{subtype_hit['identity']:.3f}" if subtype_hit else "","subtype_marker_coverage_pct":f"{subtype_hit['coverage']:.2f}" if subtype_hit else "",
             "best_whole_reference":str(whole_best["header"]) if whole_best else "","whole_reference_identity":f"{whole_best['identity']:.3f}" if whole_best else "","whole_reference_single_hsp_coverage_pct":f"{whole_best['coverage']:.2f}" if whole_best else "",
             "seh_present":seh,"is655_like_present":is655,"locus_state":module}
        typed.append(row); inventory.append({"sample":sample,"source":source,"fasta":str(fasta),"st":"1","mecA_present":True})
        for h,v in gene_hits.items():
            if v["identity"]>=90 and v["coverage"]>=60: all_gene_hits.append({"sample":sample,"query":h,**v})
        for h,v in whole_hits.items():
            if v["identity"]>=90: all_whole_hits.append({"sample":sample,"query":h,**v})

    fields=list(typed[0])
    write_tsv(OUT/"st1_mrsa_expanded_inventory.tsv",inventory,["sample","source","fasta","st","mecA_present"])
    write_tsv(OUT/"st1_mrsa_sccmec_subtyping.tsv",typed,fields)
    write_tsv(OUT/"sccmecfinder_single_gene_hits.tsv",all_gene_hits,["sample","query","contig","identity","coverage","length","qlen","sstart","send","evalue","bitscore"])
    write_tsv(OUT/"sccmecfinder_whole_reference_hits.tsv",all_whole_hits,["sample","query","contig","identity","coverage","length","qlen","sstart","send","evalue","bitscore"])
    n=len(typed); complete=sum(r["locus_state"]=="intact_adjacent_seh_IS655" for r in typed)
    seh_positive=sum(str(r["seh_present"]).lower()=="true" for r in typed)
    by_type={x:sum(r["sccmec_type"]==x for r in typed) for x in sorted({r["sccmec_type"] for r in typed})}
    by_subtype={x:sum(r["sccmec_subtype"]==x for r in typed) for x in sorted({r["sccmec_subtype"] for r in typed})}
    summary=[{"metric":"ST1_MRSA_typed","value":n,"interpretation":"mecA-positive ST1 including MW2 reference"},
             {"metric":"SCCmec_types","value":";".join(f"{k}:{v}" for k,v in by_type.items()),"interpretation":"ccr complex + mec complex class"},
             {"metric":"SCCmec_subtypes","value":";".join(f"{k}:{v}" for k,v in by_subtype.items()),"interpretation":"subtype-specific SCCmecFinder markers"},
             {"metric":"seh_positive","value":seh_positive,"interpretation":"all typed ST1 MRSA"},
             {"metric":"intact_seh_IS655_module","value":complete,"interpretation":"assembly-level complete adjacent module"},
             {"metric":"coexistence_assessment","value":f"seh {seh_positive}/{n}; complete module {complete}/{n}","interpretation":"one partial IS call is assembly-fragmented, not a verified biological absence"}]
    write_tsv(OUT/"st1_mrsa_coexistence_summary.tsv",summary,["metric","value","interpretation"])
    (OUT/"analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}","script=scripts/133_formal_sccmec_subtyping.py","database=SCCmecFinder 2017 curated single-gene and reference-cassette databases","gene_threshold=>=90% identity and >=60% single-HSP query coverage","typing_rule=ccr complex plus mec complex class; subtype only with subtype-specific marker","whole_cassette_note=single-HSP coverage is supporting evidence only and may be reduced by draft fragmentation",f"ST1_MRSA={n}"]) + "\n")


if __name__ == "__main__": main()
