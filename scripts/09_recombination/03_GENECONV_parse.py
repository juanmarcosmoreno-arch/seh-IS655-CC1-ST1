#!/usr/bin/env python3
"""Parse GENECONV and consolidate formal ST2187 recombination evidence."""

from __future__ import annotations
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/seh_st1_bovine/st2187_recombination_test"
TAB=OUT/"geneconv_run/geneconv_g0.tab"

def wt(path,rows):
 with path.open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)

def main():
 rows=[]; header=None
 for line in TAB.read_text().splitlines():
  if line.startswith("# Frag\tSequence"):
   header=["fragment_type","sequence_names","sim_p","bonferroni_KA_p","aligned_begin","aligned_end","aligned_length","seq1_begin","seq1_end","seq1_length","seq2_begin","seq2_end","seq2_length","num_poly","num_dif","total_difs","mismatch_penalty"]
  elif line.startswith("GI\t") and header:
   x=line.rstrip().split("\t"); rows.append(dict(zip(header,x)))
 focal=[r for r in rows if r["sequence_names"]=="cepa18;23EV612"]
 out=[]
 for r in focal:
  out.append({"method":"GENECONV 1.81a","status":"significant","sequence_pair":r["sequence_names"],"gscale":0,"permutations":10000,"global_permutation_p":r["sim_p"],"permutation_p_interpretation":"<1/10001" if r["sim_p"]=="0.0000" else r["sim_p"],"bonferroni_corrected_KA_p":r["bonferroni_KA_p"],"aligned_begin":r["aligned_begin"],"aligned_end":r["aligned_end"],"cepa18_begin":r["seq1_begin"],"cepa18_end":r["seq1_end"],"aligned_length":r["aligned_length"],"polymorphic_sites":r["num_poly"],"internal_mismatches":r["num_dif"]})
 wt(OUT/"geneconv_results.tsv",out)
 wt(OUT/"rdp5_results.tsv",[{"method":"RDP5 Beta 5.93/RDP5CL","status":"not_run","result":"not_available_in_native_environment","reason":"official distribution is Windows/VB; no Wine/RDP5CL installation was available; no result inferred from other methods"}])
 comp=[
 {"method":"structural ancestry markers","left_breakpoint":"14935-15115","right_breakpoint":"21065-24615","coordinate_system":"cepa18 normalized window","corrected_p":"NA","minimum_or_envelope_segment_bp":"5951-9681","concordance":"broad structural bounds"},
 {"method":"3SEQ","left_breakpoint":"15459-15479","right_breakpoint":"20904-20913","coordinate_system":"cepa18 normalized window","corrected_p":"8.579801e-265 Dunn-Sidak","minimum_or_envelope_segment_bp":"5218","concordance":"brackets complete module"},
 {"method":"GENECONV gscale=0","left_breakpoint":"15460","right_breakpoint":"20891","coordinate_system":"cepa18 normalized window; envelope of three adjacent zero-mismatch fragments","corrected_p":"best BC KA 9.93e-95; permutation <1/10001","minimum_or_envelope_segment_bp":"5432","concordance":"within structural bounds and nearly identical to 3SEQ limits"},
 {"method":"Fisher ancestry enrichment","left_breakpoint":"pre-specified anchor interval","right_breakpoint":"pre-specified anchor interval","coordinate_system":"cepa18 normalized window","corrected_p":"<1e-300","minimum_or_envelope_segment_bp":"NA","concordance":"independent site enrichment"},
 {"method":"RDP5","left_breakpoint":"NA","right_breakpoint":"NA","coordinate_system":"NA","corrected_p":"not run","minimum_or_envelope_segment_bp":"NA","concordance":"not assessed"}]
 wt(OUT/"breakpoint_comparison_table.tsv",comp)
 wt(OUT/"recombination_test_summary_v2.tsv",[
 {"method":"3SEQ","formal_result":"significant","key_statistic":"Dunn-Sidak P=8.579801e-265","interpretation":"cepa18 mosaic of negative ST2187 and ST1-like parent"},
 {"method":"GENECONV","formal_result":"significant","key_statistic":"global permutation P<1/10001; best BC KA P=9.93e-95","interpretation":"central cepa18-ST1 tract 15460-20891, concordant with 3SEQ"},
 {"method":"Fisher","formal_result":"significant","key_statistic":"OR=1778.653; P<1e-300","interpretation":"central ST1-like ancestry enrichment"},
 {"method":"RDP5","formal_result":"not tested","key_statistic":"NA","interpretation":"native executable unavailable"}])
 (OUT/"geneconv_run/analysis.log").write_text("GENECONV 1.81a; gscale=0; minPoly=2 default; 10000 permutations; seed=2187; Dumptab; ExpFormat; Annotate; ListPairs\nThe central tract is split into three adjacent zero-mismatch fragments by indel blocks; their cepa18-coordinate envelope is 15460-20891.\n")
if __name__=="__main__":main()
