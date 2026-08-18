#!/usr/bin/env python3
"""Characterize and compare the ST2187-negative rlmH--cstA defense block."""

from __future__ import annotations
import csv, datetime as dt, subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/seh_st1_bovine/st2187_ancestral_block"
NEG_FA=ROOT/"results/public_context_regional/bakta/GCA_031889545.1/GCA_031889545.1.fna"
NEG_GFF=ROOT/"results/public_context_regional/bakta/GCA_031889545.1/GCA_031889545.1.gff3"
ST1_FA=ROOT/"results/voss_2025_msphere_st1/bakta/23EV612/23EV612.fna"
ST1_GFF=ROOT/"results/voss_2025_msphere_st1/bakta/23EV612/23EV612.gff3"
ISFAA=ROOT/"resources/isfinder_thanhleviet/IS.faa"

def parse_gff(path):
 rows=[]
 for line in open(path):
  if line.startswith("#"): continue
  x=line.rstrip().split("\t")
  if len(x)!=9: continue
  a=dict(p.split("=",1) for p in x[8].split(";") if "=" in p)
  rows.append({"contig":x[0],"type":x[2],"start":int(x[3]),"end":int(x[4]),"strand":x[6],"phase":x[7],"id":a.get("ID",""),"gene":a.get("gene",""),"name":a.get("Name",""),"product":a.get("product",a.get("Name","")),"dbxref":a.get("Dbxref","")})
 return rows

def fasta(path): return {r.id:str(r.seq).upper() for r in SeqIO.parse(path,"fasta")}
def anchors(rows): return {r["gene"]:r for r in rows if r["type"]=="CDS" and r["gene"] in {"rlmH","cstA"}}
def category(r):
 t=(r["gene"]+" "+r["product"]).lower()
 if r["gene"] in {"rlmH","cstA","cstR","tauE","dus"}: return "shared_core_or_flank"
 if r["gene"] in {"pycC","pycTM"} or "pycsar" in t or "cytidylate cyclase" in t: return "antiphage_defense_Pycsar"
 if "restriction" in t or "methyltransferase" in t: return "restriction_modification_defense"
 if "membrane" in t or "hypothetical" in t: return "membrane_or_unknown"
 if "seh" in t or "enterotoxin" in t: return "virulence"
 if "is3" in t or "transposase" in t: return "mobile_element_remnant"
 return "other"
def wt(path,rows,fields=None):
 with path.open("w",newline="") as h:
  fields=fields or list(rows[0]); w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(rows)

def main():
 OUT.mkdir(parents=True,exist_ok=True); nr=parse_gff(NEG_GFF); na=anchors(nr); contig=na["rlmH"]["contig"]; lo=min(na["rlmH"]["start"],na["cstA"]["start"]); hi=max(na["rlmH"]["end"],na["cstA"]["end"]); seq=fasta(NEG_FA)[contig][lo-1:hi]
 (OUT/"st2187_ancestral_block_sequence.fasta").write_text(f">GCA_031889545.1|{contig}:{lo}-{hi}|len={len(seq)}|rlmH_to_cstA\n{seq}\n")
 cds=[r for r in nr if r["type"]=="CDS" and r["contig"]==contig and r["start"]>=lo and r["end"]<=hi]
 proteins=[]; table=[]
 for i,r in enumerate(cds,1):
  dna=seq[r["start"]-lo:r["end"]-lo+1]; dna=str(Seq(dna).reverse_complement()) if r["strand"]=="-" else dna; aa=str(Seq(dna).translate(to_stop=True)); pid=r["id"] or f"block_cds_{i:02d}"; proteins.append((pid,aa))
  table.append({"order":i,"id":pid,"start_in_block":r["start"]-lo+1,"end_in_block":r["end"]-lo+1,"strand":r["strand"],"gene":r["gene"],"product":r["product"],"functional_category":category(r),"dbxref":r["dbxref"],"protein_length_aa":len(aa)})
 wt(OUT/"st2187_ancestral_block_cds_functional_table.tsv",table)
 with (OUT/"st2187_ancestral_block_proteins.faa").open("w") as h:
  for n,s in proteins:h.write(f">{n}\n{s}\n")
 with (OUT/"st2187_ancestral_block_annotation.gff").open("w") as h:
  h.write("##gff-version 3\n##source Bakta whole-genome annotation, coordinates projected to extracted block\n")
  for r in cds: h.write("\t".join(["ST2187_ancestral_block","Bakta","CDS",str(r["start"]-lo+1),str(r["end"]-lo+1),".",r["strand"],r["phase"],f"ID={r['id']};gene={r['gene']};product={r['product']}"])+"\n")
 # ISfinder protein homology.
 fmt="6 qseqid sseqid pident length qlen slen evalue bitscore"
 b=subprocess.run(["blastp","-query",str(OUT/"st2187_ancestral_block_proteins.faa"),"-subject",str(ISFAA),"-evalue","1e-5","-outfmt",fmt],text=True,capture_output=True,check=True)
 best={}
 for line in b.stdout.splitlines():
  x=line.split("\t"); rec={"query":x[0],"database":"ISfinder","subject":x[1],"identity":x[2],"alignment_length":x[3],"query_length":x[4],"subject_length":x[5],"evalue":x[6],"bitscore":float(x[7])}
  if x[0] not in best or rec["bitscore"]>best[x[0]]["bitscore"]:best[x[0]]=rec
 mobile=[]
 for r in table:
  t=(r["gene"]+" "+r["product"]).lower(); flags=[]
  for key in ("transposase","integrase","recombinase","phage","ccr","meca"): 
   if key in t: flags.append(key)
  hit=best.get(r["id"]); qcov=100*int(hit["alignment_length"])/int(hit["query_length"]) if hit else 0
  supported=bool(flags or (hit and float(hit["identity"])>=30 and qcov>=70 and int(hit["alignment_length"])>=100))
  mobile.append({"id":r["id"],"gene":r["gene"],"product":r["product"],"annotation_mobile_flags":"|".join(flags) or "none","ISfinder_hit":hit["subject"] if hit else "none","ISfinder_identity":hit["identity"] if hit else "","ISfinder_query_coverage_pct":f"{qcov:.2f}" if hit else "","ISfinder_evalue":hit["evalue"] if hit else "","interpretation":"supported_or_review_mobile_hit" if supported else ("low_similarity_nonspecific_hit" if hit else "no_supported_mobile_element_gene")})
 wt(OUT/"st2187_ancestral_block_mobile_element_search.tsv",mobile)
 # Compare gene inventories, excluding the two anchors.
 sr=parse_gff(ST1_GFF); sa=anchors(sr); slo=min(sa["rlmH"]["start"],sa["cstA"]["start"]); shi=max(sa["rlmH"]["end"],sa["cstA"]["end"]); scds=[r for r in sr if r["type"]=="CDS" and r["contig"]==sa["rlmH"]["contig"] and r["start"]>=slo and r["end"]<=shi]
 def key(r): return r["gene"] or r["product"]
 nk={key(r):r for r in cds}; sk={key(r):r for r in scds}; comp=[]
 for k in sorted(set(nk)|set(sk)):
  source_feature=nk[k] if k in nk else sk[k]
  comp.append({"feature":k,"in_ST2187_ancestral":k in nk,"in_ST1_like":k in sk,"ST2187_product":nk[k]["product"] if k in nk else "","ST1_product":sk[k]["product"] if k in sk else "","category":category(source_feature)})
 wt(OUT/"ancestral_vs_st1like_block_comparison.tsv",comp)
 summary=[{"metric":"ancestral_block_length_bp","value":len(seq),"interpretation":"rlmH-to-cstA inclusive"},{"metric":"ST1_like_block_length_bp","value":shi-slo+1,"interpretation":"rlmH-to-cstA inclusive"},{"metric":"length_difference_bp","value":len(seq)-(shi-slo+1),"interpretation":"ancestral minus ST1-like"},{"metric":"ancestral_CDS","value":len(cds),"interpretation":"including anchors"},{"metric":"Pycsar_genes","value":sum(r["functional_category"]=="antiphage_defense_Pycsar" for r in table),"interpretation":"pycC + pycTM antiphage system"},{"metric":"restriction_modification_genes","value":sum(r["functional_category"]=="restriction_modification_defense" for r in table),"interpretation":"Type I R-M components"},{"metric":"supported_mobile_element_genes","value":sum(x["interpretation"]=="supported_or_review_mobile_hit" for x in mobile),"interpretation":"annotation or ISfinder >=30% identity, >=70% query coverage and >=100 aa"},{"metric":"classification","value":"defense island-like locus","interpretation":"Pycsar plus Type I restriction-modification; no SCC, prophage, or IS machinery supported"}]
 wt(OUT/"block_comparison_summary.tsv",summary)
 # Local whole-block search across the available ST2187 panel; remote nt is a separate network step.
 local=[]
 for fa in sorted((ROOT/"results/public_context_regional/fasta").glob("GCA_031889*.fna")):
  q=OUT/"st2187_ancestral_block_sequence.fasta"; fmt2="6 sseqid pident length qlen qcovhsp evalue bitscore"; p=subprocess.run(["blastn","-query",str(q),"-subject",str(fa),"-dust","no","-outfmt",fmt2],text=True,capture_output=True,check=True)
  if p.stdout.strip():
   x=max((z.split("\t") for z in p.stdout.splitlines()),key=lambda z:float(z[6]));local.append({"search_scope":"local_ST2187_panel_not_NCBI_nt","assembly":fa.stem,"contig":x[0],"identity":x[1],"alignment_length":x[2],"query_length":x[3],"query_coverage_pct":x[4],"evalue":x[5],"bitscore":x[6]})
 wt(OUT/"st2187_ancestral_block_blastn_nt.tsv",local)
 (OUT/"analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}","script=scripts/125_characterize_st2187_ancestral_block.py","reference_negative=GCA_031889545.1","annotation=Bakta whole-genome annotation projected to extracted block","ISfinder_search=blastp evalue 1e-5 against resources/isfinder_thanhleviet/IS.faa","VFDB/Victors=not locally available; not claimed","SCC_screen=annotation keyword screen; no formal SCC classifier","blastn_nt_status=local ST2187 screen completed; remote NCBI nt search pending","Pycsar_reference=Tal et al. Cell 2021 PMID 34644530"]) + "\n")

if __name__=="__main__":main()
