#!/usr/bin/env python3
"""Survey seh--IS655 and the cstA--rlmH interval in the Jahnen 2025 assemblies."""

from __future__ import annotations

import csv, datetime as dt, json, os, subprocess, zipfile
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "resources/jahnen_2025_cc1_cc1660/extracted/ncbi_dataset/data"
REPORT = DATA / "assembly_data_report.jsonl"
PACKAGE = ROOT / "resources/jahnen_2025_cc1_cc1660/ncbi_dataset_refseq.zip"
OUT = ROOT / "results/seh_st1_bovine/cc1_st_survey"
REF = ROOT / "results/voss_2025_msphere_st1/bakta/23EV612/23EV612.fna"
ELEMENTS = ROOT / "results/seh_st1_bovine/is3_boundaries/is3_candidate_elements.fasta"
MLST = ROOT / ".micromamba/envs/gridion_mlst/bin/mlst"

def write_tsv(path, rows, fields):
    with path.open("w", newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def refs():
    rec=next(SeqIO.parse(REF,"fasta")); s=str(rec.seq).upper()
    element=str(next(r for r in SeqIO.parse(ELEMENTS,"fasta") if r.id.startswith("23EV612|")).seq).upper()
    return {"seh":s[1884414:1885140],"spacer":s[1885140:1885352],"is655_like":element,"cstA":str(Seq(s[1880275:1881343]).reverse_complement()),"rlmH":str(Seq(s[1888366:1888846]).reverse_complement())}

def blast(query, subject):
    fmt="6 sseqid pident length qlen qstart qend sstart send bitscore"
    p=subprocess.run(["blastn","-query",str(query),"-subject",str(subject),"-dust","no","-outfmt",fmt],text=True,capture_output=True,check=True)
    hits=[]
    for line in p.stdout.splitlines():
        x=line.split("\t"); hits.append({"contig":x[0],"identity":float(x[1]),"length":int(x[2]),"qlen":int(x[3]),"sstart":int(x[6]),"send":int(x[7]),"bitscore":float(x[8])})
    return max(hits,key=lambda h:(h["length"],h["identity"],h["bitscore"])) if hits else None

def call(hit,cov=90,ident=95): return bool(hit and 100*hit["length"]/hit["qlen"]>=cov and hit["identity"]>=ident)

def read_genome(path): return {r.id:str(r.seq).upper() for r in SeqIO.parse(path,"fasta")}

def main():
    OUT.mkdir(parents=True,exist_ok=True); qpaths={}
    for k,s in refs().items():
        p=OUT/f"query_{k}.fasta"; p.write_text(f">23EV612_{k}\n{s}\n"); qpaths[k]=p
    metadata={}
    # The four-accession repair package overwrote the extracted data report;
    # recover the complete report directly from the original bulk archive.
    if PACKAGE.exists():
        with zipfile.ZipFile(PACKAGE) as z, z.open("ncbi_dataset/data/assembly_data_report.jsonl") as h:
            for raw in h:
                x=json.loads(raw); metadata[x.get("accession","")] = x
    elif REPORT.exists():
        for line in REPORT.open():
            x=json.loads(line); metadata[x.get("accession","")] = x
    env=os.environ.copy(); env["PATH"]=f"{MLST.parent}:{env.get('PATH','')}"
    rows=[]; pos_records=[]; neg_records=[]
    fastas=sorted(DATA.glob("GCF_*/*.fna"))
    for i,fasta in enumerate(fastas,1):
        acc=fasta.parent.name
        ml=subprocess.run([str(MLST),str(fasta)],text=True,capture_output=True,check=True,env=env).stdout.strip().split("\t")
        scheme=ml[1] if len(ml)>1 else ""; st=ml[2] if len(ml)>2 else "-"
        hits={k:blast(p,fasta) for k,p in qpaths.items()}
        present={k:call(h,95 if k=="is655_like" else 90,95) for k,h in hits.items()}
        central=[hits[k] for k in ("seh","spacer","is655_like") if present[k]]
        same_contig=len(central)==3 and len({h["contig"] for h in central})==1
        if same_contig:
            lo=min(min(h["sstart"],h["send"]) for h in central); hi=max(max(h["sstart"],h["send"]) for h in central)
            architecture=hi-lo+1<=3000
        else: architecture=False
        anchor_ok=call(hits["cstA"],80,90) and call(hits["rlmH"],80,90) and hits["cstA"]["contig"]==hits["rlmH"]["contig"]
        interval_len=""; anchor_state="unresolved"
        genome=None
        if anchor_ok:
            anchor_state="resolved"; lo=min(hits["cstA"]["sstart"],hits["cstA"]["send"],hits["rlmH"]["sstart"],hits["rlmH"]["send"]); hi=max(hits["cstA"]["sstart"],hits["cstA"]["send"],hits["rlmH"]["sstart"],hits["rlmH"]["send"]); interval_len=hi-lo+1
            if not present["seh"]:
                genome=read_genome(fasta); seq=genome[hits["cstA"]["contig"]][lo-1:hi]; neg_records.append((f"{acc}|ST{st}|{hits['cstA']['contig']}:{lo}-{hi}",seq))
        if architecture and st!="1":
            if genome is None: genome=read_genome(fasta)
            lo=min(min(h["sstart"],h["send"]) for h in central); hi=max(max(h["sstart"],h["send"]) for h in central); seq=genome[central[0]["contig"]][max(0,lo-501):hi+500]; pos_records.append((f"{acc}|ST{st}|module_plus_500bp",seq))
        meta=metadata.get(acc,{}); bi=(meta.get("assemblyInfo") or meta.get("assembly_info") or {}).get("biosample",{})
        row={"assembly":acc,"biosample":bi.get("accession",""),"strain":bi.get("strain",bi.get("infraspecificNames",{}).get("strain","")),"mlst_scheme":scheme,"st":st,"seh_present":present["seh"],"seh_identity":f"{hits['seh']['identity']:.3f}" if hits["seh"] else "","seh_coverage_pct":f"{min(100,100*hits['seh']['length']/hits['seh']['qlen']):.2f}" if hits["seh"] else "","is655_like_present":present["is655_like"],"is655_identity":f"{hits['is655_like']['identity']:.3f}" if hits["is655_like"] else "","is655_coverage_pct":f"{min(100,100*hits['is655_like']['length']/hits['is655_like']['qlen']):.2f}" if hits["is655_like"] else "","spacer_present":present["spacer"],"complete_adjacent_module":architecture,"anchor_state":anchor_state,"cstA_rlmH_interval_bp":interval_len}
        rows.append(row)
    fields=list(rows[0]); write_tsv(OUT/"cc1_st_survey_seh_is655_presence.tsv",rows,fields)
    groups=defaultdict(list)
    for r in rows: groups[r["st"]].append(r)
    sums=[]
    for st,x in sorted(groups.items(),key=lambda kv:(kv[0]=="-",int(kv[0]) if kv[0].isdigit() else 999999)):
        sums.append({"st":st,"n":len(x),"seh_positive":sum(r["seh_present"] for r in x),"is655_positive":sum(r["is655_like_present"] for r in x),"complete_adjacent_module":sum(r["complete_adjacent_module"] for r in x),"resolved_anchor_interval":sum(r["anchor_state"]=="resolved" for r in x)})
    write_tsv(OUT/"cc1_st_survey_by_st_summary.tsv",sums,list(sums[0]))
    for path,recs in [(OUT/"non_st1_seh_positive_loci.fasta",pos_records),(OUT/"non_st1_seh_negative_intervals.fasta",neg_records)]:
        with path.open("w") as h:
            for name,seq in recs: h.write(f">{name}\n{seq}\n")
    counts=Counter(r["st"] for r in rows)
    (OUT/"analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}","script=scripts/123_survey_cstA_rlmH_across_cc1_sts.py",f"assemblies_screened={len(rows)}",f"mlst_ST_counts={dict(counts)}","input=RefSeq assemblies downloaded from PRJNA1263585 on 2026-08-11","methodological_deviation=assembly-level BLAST presence was used; read breadth/depth was not calculated","presence_thresholds=seh >=90% cov/95% id; IS655 >=95% cov/95% id; spacer >=90% cov/95% id"]) + "\n")

if __name__=="__main__": main()
