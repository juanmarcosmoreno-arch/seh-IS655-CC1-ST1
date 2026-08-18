#!/usr/bin/env python3
"""Verify seh-locus IS3 calls with independent IS-focused signals.

This wrapper documents and runs the feasible local verification path on macOS
arm64:

* Full ISEScan was downloaded from xiezhq/ISEScan, but its bundled SSW C library
  uses x86 SSE intrinsics and does not compile natively on Apple Silicon.
* The ISEScan pHMM transposase models are still usable with HMMER. To reduce
  dependence on Bakta gene calls, ORFs are predicted in local locus windows with
  FragGeneScan, then searched with the ISEScan pHMMs.
* MobileElementFinder was not available as a standalone Bioconda/GitHub tool in
  this environment. As a documented database-based proxy, predicted ORFs are
  searched with BLASTP against the local ISfinder-sequences protein FASTA.

The output table reports Bakta coordinates, ISEScan-HMM/FragGeneScan hits,
MobileElementFinder availability, and ISfinder proxy hits.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "results/seh_st1_bovine/independent_is_detection"
DEFAULT_TABLE = ROOT / "results/seh_st1_bovine/final_tables/independent_is_tool_comparison.tsv"
FINAL_TABLE = ROOT / "results/seh_st1_bovine/final_tables/final_genome_locus_tree_table.tsv"
ISESCAN_HMM = ROOT / "tools/ISEScan/pHMMs/clusters.faa.hmm"
ISFINDER_FAA = ROOT / "deposit/seh_st1_mobile_element_reproducibility/resources/isfinder_thanhleviet/IS.faa"
ISFINDER_CSV = ROOT / "deposit/seh_st1_mobile_element_reproducibility/resources/isfinder_thanhleviet/IS.csv"

PANEL = [
    "23EV612",
    "cepa_17_barcode01",
    "cepa_21_barcode03",
    "cepa_30_barcode06",
    "cepa_52_barcode09",
    "cepa_98_barcode16",
    "cepa_107_barcode19",
    "cepa_202_barcode21",
    "cepa_18_barcode02",
    "B1822_S7_L001",
    "B916_S394",
    "B1180_S15",
    "GCF_051990345.1",
    "GCF_051990445.1",
    "GCF_051990565.1",
    "GCA_001680935.1",
    "GCA_001681155.1",
    "GCA_024360105.1",
    "GCA_030734435.1",
    "GCA_029415515.1",
]

FIELDS = [
    "genome_id",
    "group",
    "tree_phylo_group",
    "Bakta_call",
    "Bakta_coordinates",
    "ISEScan_call",
    "ISEScan_coordinates",
    "ISEScan_evalue",
    "ISEScan_bit_score",
    "MobileElementFinder_call",
    "MobileElementFinder_coordinates",
    "ISfinder_proxy_call",
    "ISfinder_proxy_coordinates",
    "ISfinder_proxy_identity",
    "ISfinder_proxy_coverage",
    "ISfinder_proxy_evalue",
    "ISfinder_proxy_bit_score",
    "position_overlap_Bakta_ISEScan",
    "position_overlap_Bakta_MobileElementFinder",
    "position_overlap_Bakta_ISfinder_proxy",
    "family_match_Bakta_ISEScan",
    "family_match_Bakta_MobileElementFinder",
    "family_match_Bakta_ISfinder_proxy",
    "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def write_fasta(path: Path, name: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            handle.write(seq[i : i + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(cmd: list[str], stdout: Path | None = None) -> None:
    stdout_handle = stdout.open("w", encoding="utf-8") if stdout else subprocess.DEVNULL
    try:
        subprocess.run(cmd, check=True, text=True, stdout=stdout_handle, stderr=subprocess.STDOUT)
    finally:
        if stdout:
            stdout_handle.close()


def parse_fgs_id(target: str) -> tuple[str, int, int, str] | None:
    match = re.match(r"^(.+)_(\d+)_(\d+)_([+-])$", target)
    if not match:
        return None
    seqid, start, end, strand = match.groups()
    return seqid, int(start), int(end), strand


def parse_hmm_family(model: str) -> str:
    return model.split("_")[0].replace(".faa", "")


def normalize_is_name(value: str) -> str:
    value = value.replace("html.2020//", "")
    return value.split()[0].split("~~~")[0].split("|")[0].split("_")[0]


def read_isfinder_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["Name"]: row for row in csv.DictReader(handle)}


def fasta_by_header(path: Path) -> dict[str, str]:
    return read_fasta(path)


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return max(a_start, b_start) <= min(a_end, b_end)


def coords(contig: str, start: int | str, end: int | str, strand: str) -> str:
    return f"{contig}:{start}-{end}({strand})"


def best_hmm_overlap(tblout: Path, offset: int, bakta: dict[str, str]) -> dict[str, str] | None:
    b_start = int(bakta["nearest_is3_start"])
    b_end = int(bakta["nearest_is3_end"])
    best: dict[str, str] | None = None
    with tblout.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            parsed = parse_fgs_id(parts[0])
            if parsed is None:
                continue
            seqid, local_start, local_end, strand = parsed
            start = offset + local_start - 1
            end = offset + local_end - 1
            if not overlap(start, end, b_start, b_end):
                continue
            hit = {
                "target": parts[0],
                "model": parts[2],
                "family": parse_hmm_family(parts[2]),
                "start": str(start),
                "end": str(end),
                "strand": strand,
                "evalue": parts[4],
                "bitscore": parts[5],
            }
            if best is None or float(hit["bitscore"]) > float(best["bitscore"]):
                best = hit
    return best


def best_blast_for_target(blastout: Path, target: str, is_meta: dict[str, dict[str, str]]) -> dict[str, str] | None:
    best: dict[str, str] | None = None
    with blastout.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            qseqid, sseqid, pident, qcovhsp, evalue, bitscore = line.rstrip("\n").split("\t")
            if qseqid != target:
                continue
            is_name = normalize_is_name(sseqid)
            meta = is_meta.get(is_name, {})
            hit = {
                "hit": is_name,
                "family": meta.get("Family", "NA"),
                "group": meta.get("Group", "unknown") or "unknown",
                "identity": pident,
                "coverage": qcovhsp,
                "evalue": evalue,
                "bitscore": bitscore,
            }
            if best is None or float(hit["bitscore"]) > float(best["bitscore"]):
                best = hit
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--out-table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--window", type=int, default=5000)
    parser.add_argument("--threads", type=str, default="2")
    parser.add_argument(
        "--hmmer-bin",
        type=Path,
        default=ROOT / ".micromamba_is_tools/envs/is_hmmer/bin/hmmsearch",
    )
    parser.add_argument(
        "--fgs-bin",
        type=Path,
        default=ROOT / ".micromamba_is_tools/envs/is_hmmer/bin/run_FragGeneScan.pl",
    )
    parser.add_argument("--blastp", default="blastp")
    args = parser.parse_args()

    rows_by_id = {row["sample_id"]: row for row in read_tsv(FINAL_TABLE)}
    is_meta = read_isfinder_metadata(ISFINDER_CSV)
    output_rows: list[dict[str, str]] = []

    for sample_id in PANEL:
        row = rows_by_id[sample_id]
        sample_dir = args.outdir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        assembly = (ROOT / row["bakta_tsv"]).with_suffix(".fna")
        seqs = read_fasta(assembly)
        contig = row["seh_contig"]
        seq = seqs[contig]
        b_start = int(row["nearest_is3_start"])
        b_end = int(row["nearest_is3_end"])
        win_start = max(1, min(b_start, b_end) - args.window)
        win_end = min(len(seq), max(b_start, b_end) + args.window)
        window_seq = seq[win_start - 1 : win_end]
        window_fna = sample_dir / f"{sample_id}.is3_window.fna"
        write_fasta(window_fna, "window", window_seq)

        fgs_prefix = sample_dir / sample_id
        faa = Path(str(fgs_prefix) + ".faa")
        if not faa.exists():
            run(
                [
                    str(args.fgs_bin),
                    f"-genome={window_fna}",
                    f"-out={fgs_prefix}",
                    "-complete=1",
                    "-train=complete",
                    f"-thread={args.threads}",
                ],
                stdout=sample_dir / "fraggenescan.log",
            )

        hmm_tbl = sample_dir / f"{sample_id}.isescan_hmm.tbl"
        if not hmm_tbl.exists():
            run(
                [
                    str(args.hmmer_bin),
                    "--tblout",
                    str(hmm_tbl),
                    "--noali",
                    "-E",
                    "1e-5",
                    str(ISESCAN_HMM),
                    str(faa),
                ],
                stdout=sample_dir / "isescan_hmmsearch.log",
            )
        hmm_hit = best_hmm_overlap(hmm_tbl, win_start, row)

        blast_out = sample_dir / f"{sample_id}.isfinder_blastp.tsv"
        if not blast_out.exists():
            proc = subprocess.run(
                [
                    args.blastp,
                    "-query",
                    str(faa),
                    "-subject",
                    str(ISFINDER_FAA),
                    "-max_target_seqs",
                    "5",
                    "-outfmt",
                    "6 qseqid sseqid pident qcovhsp evalue bitscore",
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            blast_out.write_text(proc.stdout, encoding="utf-8")

        isfinder_hit = best_blast_for_target(blast_out, hmm_hit["target"], is_meta) if hmm_hit else None

        bakta_coords = coords(contig, row["nearest_is3_start"], row["nearest_is3_end"], row["nearest_is3_strand"])
        isescan_coords = (
            coords(contig, hmm_hit["start"], hmm_hit["end"], hmm_hit["strand"]) if hmm_hit else "no_overlap_hit"
        )
        proxy_coords = isescan_coords if isfinder_hit else "no_overlap_hit"
        isescan_family_match = "TRUE" if hmm_hit and hmm_hit["family"] == "IS3" else "FALSE"
        proxy_family_match = "TRUE" if isfinder_hit and isfinder_hit["family"] == "IS3" else "FALSE"
        output_rows.append(
            {
                "genome_id": sample_id,
                "group": row["analysis_group"],
                "tree_phylo_group": row["tree_phylo_group"],
                "Bakta_call": row["nearest_is3_product"],
                "Bakta_coordinates": bakta_coords,
                "ISEScan_call": f"{hmm_hit['model']} ({hmm_hit['family']})" if hmm_hit else "no_overlap_hit",
                "ISEScan_coordinates": isescan_coords,
                "ISEScan_evalue": hmm_hit["evalue"] if hmm_hit else "NA",
                "ISEScan_bit_score": hmm_hit["bitscore"] if hmm_hit else "NA",
                "MobileElementFinder_call": "not_run_standalone_tool_not_available",
                "MobileElementFinder_coordinates": "NA",
                "ISfinder_proxy_call": (
                    f"{isfinder_hit['hit']} ({isfinder_hit['family']}/{isfinder_hit['group']})"
                    if isfinder_hit
                    else "no_overlap_hit"
                ),
                "ISfinder_proxy_coordinates": proxy_coords,
                "ISfinder_proxy_identity": isfinder_hit["identity"] if isfinder_hit else "NA",
                "ISfinder_proxy_coverage": isfinder_hit["coverage"] if isfinder_hit else "NA",
                "ISfinder_proxy_evalue": isfinder_hit["evalue"] if isfinder_hit else "NA",
                "ISfinder_proxy_bit_score": isfinder_hit["bitscore"] if isfinder_hit else "NA",
                "position_overlap_Bakta_ISEScan": "TRUE" if hmm_hit else "FALSE",
                "position_overlap_Bakta_MobileElementFinder": "NA",
                "position_overlap_Bakta_ISfinder_proxy": "TRUE" if isfinder_hit else "FALSE",
                "family_match_Bakta_ISEScan": isescan_family_match,
                "family_match_Bakta_MobileElementFinder": "NA",
                "family_match_Bakta_ISfinder_proxy": proxy_family_match,
                "notes": (
                    "Full ISEScan unavailable on macOS arm64 because bundled SSW uses x86 SSE; "
                    "ISEScan pHMMs were run on FragGeneScan ORFs. MobileElementFinder standalone not found; "
                    "ISfinder BLASTP is reported as database proxy."
                ),
            }
        )

    write_tsv(args.out_table, output_rows)
    print(f"Wrote {len(output_rows)} rows to {args.out_table}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
