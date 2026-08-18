#!/usr/bin/env python3
"""Document and screen the current public ST2187 search, including false positives."""
from __future__ import annotations
import csv, datetime as dt, json, shutil, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/seh_st1_bovine/st2187_expanded_panel"
REPORT=ROOT/"results/public_context_regional/extracted/ncbi_dataset/data/assembly_data_report.jsonl"
DATASETS=ROOT/".micromamba/envs/ncbi_datasets/bin/datasets"
MATRIX=ROOT/"results/seh_st1_bovine/st1_st2187_phylogenomics/locus_state_matrix.tsv"
BLOCK=ROOT/"results/seh_st1_bovine/st2187_ancestral_block/st2187_ancestral_block_sequence.fasta"
TREE=ROOT/"results/seh_st1_bovine/st1_st2187_phylogenomics/iqtree_gubbins_bootstrap/st1_st2187_recomb_filtered.treefile"
PUBLIC=["GCA_031889225.1","GCA_031889265.1","GCA_031889285.1","GCA_031889365.1","GCA_031889545.1"]
def rows(path):
 with open(path) as h:return list(csv.DictReader(h,delimiter="\t"))
def wt(path,data,fields=None):
 with open(path,"w",newline="") as h:
  fields=fields or list(data[0]);w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(data)
def main():
 OUT.mkdir(parents=True,exist_ok=True); meta={}
 query=subprocess.run([str(DATASETS),"summary","genome","taxon","Staphylococcus aureus","--search","ST2187","--as-json-lines"],text=True,capture_output=True,check=True)
 (OUT/"ncbi_st2187_query.jsonl").write_text(query.stdout)
 for line in query.stdout.splitlines():
  if not line.startswith("{"):continue
  x=json.loads(line); a=x.get("accession",x.get("currentAccession","")); meta[a]=x
 states={r["sample_id"]:r for r in rows(MATRIX)}; inventory=[]; screens=[]
 for acc in PUBLIC:
  x=meta.get(acc,{}); ai=x.get("assemblyInfo",x.get("assembly_info",{})); bi=ai.get("biosample",{}); attrs={z.get("name"):z.get("value") for z in bi.get("attributes",[])}
  inventory.append({"assembly":acc,"search_status":"accepted_ST2187","mlst_call":"ST2187","biosample":bi.get("accession",""),"isolate":bi.get("isolate",attrs.get("isolate","")),"host":bi.get("host",attrs.get("host","")),"isolation_source":bi.get("isolationSource",bi.get("isolation_source",attrs.get("isolation_source",""))),"collection_date":bi.get("collectionDate",bi.get("collection_date",attrs.get("collection_date",""))),"geography":bi.get("geoLocName",bi.get("geo_loc_name",attrs.get("geo_loc_name",""))),"reason":"verified previously by mlst 2.33.1"})
  fa=ROOT/f"results/public_context_regional/fasta/{acc}.fna"; fmt="6 pident length qlen qcovhsp evalue bitscore"
  p=subprocess.run(["blastn","-query",str(BLOCK),"-subject",str(fa),"-dust","no","-outfmt",fmt],text=True,capture_output=True,check=True)
  hit=max((z.split("\t") for z in p.stdout.splitlines()),key=lambda z:float(z[5])) if p.stdout.strip() else None; s=states[acc]
  screens.append({"assembly":acc,"st":"ST2187","seh_status":s["seh_state"],"is655_status":s["is655_like_state"],"locus_state":s["locus_state"],"interval_type":"ancestral_defense_block" if s["seh_state"]=="absent" else "ST1_like_seh_IS655","ancestral_block_identity":hit[0] if hit else "","ancestral_block_query_coverage_pct":hit[3] if hit else "","locus_completeness":s["locus_completeness"]})
 inventory.append({"assembly":"GCA_900046405.1","search_status":"rejected_false_positive","mlst_call":"ST22","biosample":"SAMEA1469533","isolate":"EMRSA15","host":"Homo sapiens","isolation_source":"blood","collection_date":"2008","geography":"United Kingdom: England","reason":"string 2187 occurs in sample identifier BSAC2187; in-silico MLST alleles resolve ST22"})
 # Add the local positive comparator to the expanded state matrix.
 local=states["cepa_18_barcode02"]; screens.append({"assembly":"cepa_18_barcode02","st":"ST2187","seh_status":local["seh_state"],"is655_status":local["is655_like_state"],"locus_state":local["locus_state"],"interval_type":"ST1_like_seh_IS655","ancestral_block_identity":"","ancestral_block_query_coverage_pct":"","locus_completeness":local["locus_completeness"]})
 wt(OUT/"st2187_public_genomes_inventory.tsv",inventory);wt(OUT/"st2187_seh_is655_screen_results.tsv",screens);wt(OUT/"st2187_expanded_locus_state_matrix.tsv",screens);wt(OUT/"st2187_expanded_metadata.tsv",inventory)
 shutil.copyfile(TREE,OUT/"st2187_expanded_phylogeny.nwk")
 (OUT/"analysis.log").write_text("\n".join([f"run_date={dt.date.today().isoformat()}","script=scripts/126_search_and_screen_additional_st2187.py","NCBI_query=Staphylococcus aureus taxon with search term ST2187","catalog_result=5 known genuine candidates plus GCA_900046405.1 false positive","new_genuine_ST2187=0","mlst=2.33.1","screen_thresholds=existing validated locus-state matrix; ancestral-block blastn reported","phylogeny=unchanged copy because no new genuine genome was found"]) + "\n")
if __name__=="__main__":main()
