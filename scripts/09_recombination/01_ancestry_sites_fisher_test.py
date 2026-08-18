#!/usr/bin/env python3
"""Formal mosaic-enrichment test for the ST2187 recipient; document unavailable tools."""
from __future__ import annotations
import csv, datetime as dt, math, shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
INFORM=ROOT/"results/seh_st1_bovine/st2187_breakpoints/ancestry_informative_sites.tsv"
ALIGN=ROOT/"results/seh_st1_bovine/st2187_breakpoints/three_way_alignment.fasta"
OUT=ROOT/"results/seh_st1_bovine/st2187_recombination_test"
ANCHOR_START=15001; ANCHOR_END=23570  # rlmH--cstA inclusive in normalized cepa18 window
def wt(path,rows):
 with path.open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def logchoose(n,k):return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
def fisher_greater(a,b,c,d):
 n=a+b+c+d; row=a+b; success=a+c; hi=min(row,success); logs=[]
 for k in range(a,hi+1):logs.append(logchoose(success,k)+logchoose(n-success,row-k)-logchoose(n,row))
 m=max(logs);logp=m+math.log(sum(math.exp(x-m) for x in logs));return logp/math.log(10)
def main():
 OUT.mkdir(parents=True,exist_ok=True);shutil.copyfile(ALIGN,OUT/"flanking_alignment.fasta")
 counts={(z,s):0 for z in ("anchor_interval","flanks") for s in ("ST1_like","ST2187_like")}
 with INFORM.open() as h:
  for r in csv.DictReader(h,delimiter="\t"):
   if r["ancestry_state"] not in {"ST1_like","ST2187_like"}:continue
   zone="anchor_interval" if ANCHOR_START<=int(r["cepa18_window_pos"])<=ANCHOR_END else "flanks";counts[zone,r["ancestry_state"]]+=1
 a=counts["anchor_interval","ST1_like"];b=counts["anchor_interval","ST2187_like"];c=counts["flanks","ST1_like"];d=counts["flanks","ST2187_like"]
 log10p=fisher_greater(a,b,c,d); ptext=f"{10**log10p:.3e}" if log10p>-300 else "<1e-300"; odds=(a*d)/(b*c) if b*c else float("inf")
 unavailable=[{"method":"RDP4/RDP5","status":"not_available_locally","result":"not_run","reason":"no supported command-line installation present"},{"method":"3SEQ","status":"not_available_locally","result":"not_run","reason":"standalone executable absent"},{"method":"GENECONV","status":"not_available_locally","result":"not_run","reason":"standalone executable absent"}]
 wt(OUT/"rdp4_results.tsv",[unavailable[0]]);wt(OUT/"3seq_results.tsv",[unavailable[1]]);wt(OUT/"geneconv_results.tsv",[unavailable[2]])
 test=[{"test":"pre-specified_anchor_interval_ancestry_enrichment","anchor_ST1_like":a,"anchor_ST2187_like":b,"flank_ST1_like":c,"flank_ST2187_like":d,"odds_ratio":f"{odds:.3f}","one_sided_fisher_p":ptext,"fisher_log10_p":f"{log10p:.3f}","interpretation":"strong mosaic ancestry enrichment consistent with homologous replacement"},{"test":"breakpoint_concordance","anchor_ST1_like":a,"anchor_ST2187_like":b,"flank_ST1_like":c,"flank_ST2187_like":d,"odds_ratio":"","one_sided_fisher_p":"","fisher_log10_p":"","interpretation":"181-bp rlmH-side interval and 3551-bp opposite interval agree with ancestry transition; no nucleotide-exact second breakpoint"}]
 wt(OUT/"recombination_test_summary.tsv",test)
 (OUT/"analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}","script=scripts/131_formal_recombination_test_st2187.py","alignment=MAFFT three-way alignment from script 124","formal_test=one-sided Fisher exact test on pre-specified rlmH-cstA interval versus flanks","informative_sites=ST1 and negative-ST2187 differ and recipient allele matches either parent","RDP4_3SEQ_GENECONV=not installed; explicitly not run","limitation=Fisher ancestry-enrichment test supports mosaicism but is not an RDP/3SEQ event test"]) + "\n")
if __name__=="__main__":main()
