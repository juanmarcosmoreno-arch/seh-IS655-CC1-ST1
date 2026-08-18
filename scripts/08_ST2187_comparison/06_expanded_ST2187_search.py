#!/usr/bin/env python3
"""Expanded ST2187 survey using allelic typing, literature cohorts and assemblies."""
from __future__ import annotations

import csv
import datetime as dt
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/seh_st1_bovine/st2187_expanded_search"
MLST_DB = ROOT / ".micromamba/envs/gridion_mlst/db/pubmlst/saureus/saureus.txt"
MLST_INFO = ROOT / ".micromamba/envs/gridion_mlst/db/pubmlst/saureus/saureus_info.json"
MLST = ROOT / ".micromamba/envs/gridion_mlst/bin/mlst"
SEH = ROOT / "results/seh_st1_bovine/st1_mrsa_seh_sccmec/query_seh.fasta"
IS655 = ROOT / "results/seh_st1_bovine/st1_mrsa_seh_sccmec/query_is655_like.fasta"
DEFENSE = ROOT / "results/seh_st1_bovine/st2187_ancestral_block/st2187_ancestral_block_sequence.fasta"
MULTI_META = ROOT / "resources/st2187_expanded_search/multihost_st2187_samples.tsv"
MULTI = ROOT / "resources/st2187_expanded_search/multihost_2023/ncbi_dataset/data"
CANADA = ROOT / "resources/st2187_expanded_search/prjna599195/ncbi_dataset/data"
MAINE = ROOT / "results/public_context_regional/extracted/ncbi_dataset/data"
LOCAL_POS = ROOT / "results/seh_st1_bovine/phylogenomics_set/assemblies/cepa_18_barcode02.fna"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle: return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def blast(query: Path, subject: Path) -> dict[str, object] | None:
    fmt="6 sseqid pident length qlen qcovhsp evalue bitscore sstart send"
    p=subprocess.run(["blastn","-query",str(query),"-subject",str(subject),"-dust","no","-outfmt",fmt],text=True,capture_output=True,check=True)
    hits=[]
    for line in p.stdout.splitlines():
        x=line.split("\t"); hits.append({"contig":x[0],"identity":float(x[1]),"length":int(x[2]),"qlen":int(x[3]),"qcov":float(x[4]),"evalue":x[5],"bitscore":float(x[6]),"start":int(x[7]),"end":int(x[8])})
    return max(hits,key=lambda h:(h["qcov"],h["identity"],h["bitscore"])) if hits else None


def mlst(fasta: Path) -> tuple[str,str]:
    env=os.environ.copy(); env["PATH"]=f"{MLST.parent}:{env.get('PATH','')}"
    p=subprocess.run([str(MLST),"--quiet",str(fasta)],text=True,capture_output=True,env=env,check=True)
    x=p.stdout.strip().split("\t")
    return (x[2] if len(x)>2 else "unresolved", ";".join(x[3:]) if len(x)>3 else "")


def report_for(path: Path, assembly: str, cohort: str, meta: dict[str,str]) -> dict[str,object]:
    st,profile=mlst(path); sh=blast(SEH,path); ih=blast(IS655,path); dh=blast(DEFENSE,path)
    sp=bool(sh and sh["identity"]>=95 and sh["qcov"]>=90)
    ip=bool(ih and ih["identity"]>=95 and ih["qcov"]>=95)
    adjacent=bool(sp and ip and sh["contig"]==ih["contig"] and abs(min(sh["start"],sh["end"])-min(ih["start"],ih["end"]))<3000)
    defense_complete=bool(dh and dh["identity"]>=95 and dh["qcov"]>=90)
    state="complete_seh_IS655_module" if adjacent else "ancestral_defense_block" if defense_complete else "seh_or_IS_partial" if sp or ip else "unresolved_or_other"
    return {"assembly":assembly,"cohort":cohort,"isolate":meta.get("isolate",""),"biosample":meta.get("biosample",""),"collection_date":meta.get("collection_date",""),"geography":meta.get("geography",""),"host":meta.get("host",""),"isolation_source":meta.get("isolation_source",""),"mlst_call":st,"allelic_profile":profile,
            "seh_present":sp,"seh_identity":f"{sh['identity']:.3f}" if sh else "","seh_coverage_pct":f"{sh['qcov']:.2f}" if sh else "",
            "is655_complete":ip,"is655_identity":f"{ih['identity']:.3f}" if ih else "","is655_coverage_pct":f"{ih['qcov']:.2f}" if ih else "",
            "seh_is655_adjacent":adjacent,"defense_block_complete":defense_complete,"defense_block_identity":f"{dh['identity']:.3f}" if dh else "","defense_block_coverage_pct":f"{dh['qcov']:.2f}" if dh else "","locus_state":state,"fasta":str(path)}


def assembly_report(root: Path) -> dict[str,dict[str,str]]:
    path=root/"assembly_data_report.jsonl"; out={}
    if not path.exists(): return out
    for line in path.open():
        d=json.loads(line); acc=d.get("accession",""); ai=d.get("assemblyInfo",d.get("assembly_info",{})) or {}; b=ai.get("biosample",{}); attrs={a.get("name",""):a.get("value","") for a in b.get("attributes",[])}
        out[acc]={"biosample":b.get("accession",""),"isolate":b.get("isolate",b.get("strain","")),"collection_date":b.get("collectionDate",b.get("collection_date",attrs.get("collection_date",""))),"geography":b.get("geoLocName",b.get("geo_loc_name",attrs.get("geo_loc_name",""))),"host":b.get("host",attrs.get("host","")),"isolation_source":b.get("isolationSource",b.get("isolation_source",attrs.get("isolation_source","")))}
    return out


def main() -> None:
    OUT.mkdir(parents=True,exist_ok=True)
    with MLST_DB.open() as h:
        reader=csv.DictReader(h,delimiter="\t"); profile=next(r for r in reader if r["ST"]=="2187")
    info=json.loads(MLST_INFO.read_text())
    pubrows=[{"record_type":"profile_definition","st":"2187",**profile,"database":"PubMLST S. aureus sequence definitions","database_download_date":info.get("download_date",""),"source_url":info.get("API","")},
             {"record_type":"published_cohort","st":"2187","arcC":"3","aroE":"317","glpF":"1","gmk":"1","pta":"1","tpi":"5","yqiL":"3","database":"Souza et al. 2023 supplementary Table S1","database_download_date":"","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC10463932/"}]
    write_tsv(OUT/"pubmlst_st2187_isolates.tsv",pubrows)

    rows=[]; multi_meta={r["biosample"]:{"isolate":r["isolate"],"biosample":r["biosample"],"collection_date":r["collection_date"],"geography":"USA: northern New England","host":r["host"],"isolation_source":r["isolation_site"]} for r in read_tsv(MULTI_META)}
    for root,cohort in [(MULTI,"Souza_2023_multihost"),(CANADA,"Naushad_2020_Canada")]:
        reports=assembly_report(root)
        for p in sorted(root.glob("*/*genomic.fna")):
            acc=p.parent.name; st,_=mlst(p)
            if st != "2187": continue
            meta=reports.get(acc,{})
            if cohort.startswith("Souza"): meta=multi_meta.get(meta.get("biosample",""),meta)
            rows.append(report_for(p,acc,cohort,meta))
    oldmeta={r["assembly"]:r for r in read_tsv(ROOT/"results/seh_st1_bovine/st2187_expanded_panel/st2187_expanded_metadata.tsv")}
    for acc in ["GCA_031889225.1","GCA_031889265.1","GCA_031889285.1","GCA_031889365.1","GCA_031889545.1"]:
        fs=[p for p in (MAINE/acc).glob("*genomic.fna") if not p.name.startswith("cds_")]; m=oldmeta.get(acc,{})
        if fs: rows.append(report_for(fs[0],acc,"Maine_2022",{"isolate":m.get("isolate",""),"biosample":m.get("biosample",""),"collection_date":m.get("collection_date",""),"geography":m.get("geography",""),"host":m.get("host",""),"isolation_source":m.get("isolation_source","")}))
    rows.append(report_for(LOCAL_POS,"cepa_18_barcode02","Uruguay_local_positive",{"isolate":"cepa_18_barcode02","geography":"Uruguay"}))
    write_tsv(OUT/"st2187_expanded_inventory.tsv",rows)

    ena=[{"query":"ST2187","search_scope":"ENA textual metadata","result":"no direct metadata hits recovered","interpretation":"ST2187 accessions were instead recovered from published supplementary tables and verified by in-silico MLST"}]
    for r in rows:
        if r["cohort"]=="Souza_2023_multihost": ena.append({"query":r["biosample"],"search_scope":"published SRA BioProject PRJNA741582","result":r["assembly"],"interpretation":"assembly linked to published BioSample; ST2187 confirmed in silico"})
    write_tsv(OUT/"ena_st2187_search_results.tsv",ena)
    lit=[{"study":"Naushad et al. 2020","year":2020,"reported_ST2187":3,"genomic_access":"PRJNA599195","geography":"Canada","source":"bovine milk","verified_assemblies_in_current_analysis":sum(r["cohort"]=="Naushad_2020_Canada" for r in rows),"url":"https://doi.org/10.1128/msystems.00063-20"},
         {"study":"Souza et al. 2023","year":2023,"reported_ST2187":16,"genomic_access":"PRJNA741582","geography":"northern New England, USA","source":"cattle disease isolates; mostly milk, one lung","verified_assemblies_in_current_analysis":sum(r["cohort"]=="Souza_2023_multihost" for r in rows),"url":"https://doi.org/10.1186/s12866-023-02985-9"},
         {"study":"Pennsylvania milk survey","year":2021,"reported_ST2187":28,"genomic_access":"not resolved from article","geography":"Pennsylvania, USA","source":"mastitic and bulk tank milk","verified_assemblies_in_current_analysis":0,"url":"https://doi.org/10.1371/journal.pone.0248528"},
         {"study":"US bulk tank milk survey","year":2021,"reported_ST2187":26,"genomic_access":"not resolved from article","geography":"USA","source":"bulk tank milk","verified_assemblies_in_current_analysis":0,"url":"https://doi.org/10.1186/s12917-021-02940-0"}]
    write_tsv(OUT/"st2187_literature_search.tsv",lit)
    genuine=[r for r in rows if r["mlst_call"]=="2187"]
    summary=[{"metric":"ST2187_profile","value":"3-317-1-1-1-5-3","interpretation":"arcC-aroE-glpF-gmk-pta-tpi-yqiL"},{"metric":"verified_genomes_total","value":len(genuine),"interpretation":"in-silico MLST ST2187"},{"metric":"verified_additional_beyond_Maine_and_local","value":sum(r["cohort"] not in {"Maine_2022","Uruguay_local_positive"} for r in genuine),"interpretation":"newly recovered in script 134"},{"metric":"seh_positive","value":sum(bool(r["seh_present"]) for r in genuine),"interpretation":"assembly thresholds"},{"metric":"complete_defense_block","value":sum(bool(r["defense_block_complete"]) for r in genuine),"interpretation":">=95% identity and >=90% query coverage"},{"metric":"locus_states","value":";".join(f"{s}:{sum(r['locus_state']==s for r in genuine)}" for s in sorted({r['locus_state'] for r in genuine})),"interpretation":"expanded panel"}]
    write_tsv(OUT/"st2187_expanded_search_summary.tsv",summary)
    versions=subprocess.run(["blastn","-version"],text=True,capture_output=True,check=True).stdout.splitlines()[0]
    (OUT/"analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}","script=scripts/134_expanded_st2187_search.py",f"blastn={versions}","mlst=2.33.1; PubMLST database download 2026-03-11","seh_threshold=>=95% identity, >=90% coverage","IS655_threshold=>=95% identity, >=95% coverage","defense_block_threshold=>=95% identity, >=90% coverage","literature_sources=primary articles and supplementary Table S1","limitation=PubMLST isolate database count not directly accessible; published cohorts reported separately"]) + "\n")


if __name__ == "__main__": main()
