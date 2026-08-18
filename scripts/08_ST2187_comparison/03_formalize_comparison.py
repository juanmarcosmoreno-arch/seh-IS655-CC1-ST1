#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
BASE = PROJECT / "results" / "seh_st1_bovine" / "st2187_comparator"
SYNTENY = BASE / "panel_synteny_context.tsv"
MODULES = BASE / "panel_IS3_spacer_seh_modules.fasta"
FLANKS = BASE / "panel_seh_flanks_20kb_normalized.fasta"
OUT_LEVELS = BASE / "formal_comparison_levels_identity.tsv"
OUT_BOUNDARIES = BASE / "st2187_conserved_block_boundaries.tsv"
OUT_REPEATS = BASE / "module_edge_repeat_scan.tsv"
OUT_PUBLIC = BASE / "public_st2187_seh_screen.tsv"

PUBLIC_ST2187 = [
    "GCA_031889225.1",
    "GCA_031889265.1",
    "GCA_031889285.1",
    "GCA_031889365.1",
    "GCA_031889545.1",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    out.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.upper())
    if header is not None:
        out.append((header, "".join(seq)))
    return out


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
    return seq.translate(table)[::-1]


def sample_from_header(header: str) -> str:
    return header.split("|", 1)[0]


def best_blast(query: Path, subject: Path, extra: list[str] | None = None) -> list[dict[str, str]]:
    fields = [
        "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qcovhsp", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    cmd = [
        "blastn", "-query", str(query), "-subject", str(subject),
        "-task", "blastn", "-dust", "no",
        "-outfmt", "6 " + " ".join(fields),
    ]
    if extra:
        cmd.extend(extra)
    proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return [dict(zip(fields, line.split("\t"))) for line in proc.stdout.splitlines()]


def top_hit(hits: list[dict[str, str]]) -> dict[str, str] | None:
    if not hits:
        return None
    return max(hits, key=lambda r: (float(r["bitscore"]), float(r["qcovhsp"]), float(r["pident"])))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def find_direct_repeat(left: str, right: str, max_k: int = 30) -> tuple[str, int]:
    best = ""
    for k in range(4, max_k + 1):
        if len(left) >= k and len(right) >= k and left[-k:] == right[:k]:
            best = left[-k:]
    return best, len(best)


def find_terminal_inverted_repeat(seq: str, max_k: int = 40) -> tuple[str, int]:
    best = ""
    for k in range(6, max_k + 1):
        if len(seq) >= 2 * k and seq[:k] == revcomp(seq[-k:]):
            best = seq[:k]
    return best, len(best)


def feature_search(path: Path, patterns: list[re.Pattern[str]]) -> bool:
    if not path.exists():
        return False
    text = path.read_text(errors="ignore")
    return any(p.search(text) for p in patterns)


def main() -> None:
    synteny = {r["sample_id"]: r for r in read_tsv(SYNTENY)}
    modules = {sample_from_header(h): s for h, s in read_fasta(MODULES)}
    flanks = {sample_from_header(h): s for h, s in read_fasta(FLANKS)}

    # The normalized module is IS3 + 653 nt spacer + 726 nt seh.
    ref_module = modules["23EV612"]
    is3_len = 630
    nucleus_start = is3_len
    nucleus_records = [(sample, seq[nucleus_start:]) for sample, seq in modules.items()]
    mobile_records = list(modules.items())
    extended_records = list(flanks.items())

    levels = [
        ("nucleo_minimo", "654 bp coordinate gap + seh", "conservacion_directa_locus_inmediato", nucleus_records),
        ("modulo_movil", "IS3 + 654 bp coordinate gap + seh", "arquitectura_movil_conservada", mobile_records),
        ("contexto_extendido", "+/-20 kb alrededor de seh", "conservacion_parcial_y_bordes", extended_records),
    ]

    identity_rows: list[dict[str, str]] = []
    tmp_dir = BASE / "_tmp_levels"
    tmp_dir.mkdir(exist_ok=True)
    try:
        for level_id, region, question, records in levels:
            ref_seq = dict(records)["23EV612"]
            ref = tmp_dir / f"ref_{level_id}.fasta"
            write_fasta(ref, [(f"23EV612_{level_id}", ref_seq)])
            for sample, seq in records:
                subject = tmp_dir / f"{sample}_{level_id}.fasta"
                write_fasta(subject, [(f"{sample}_{level_id}", seq)])
                hit = top_hit(best_blast(ref, subject))
                identity_rows.append({
                    "level": level_id,
                    "region": region,
                    "question": question,
                    "sample_id": sample,
                    "group": synteny.get(sample, {}).get("group", ""),
                    "st": synteny.get(sample, {}).get("st", ""),
                    "pident_vs_23EV612": "" if hit is None else hit["pident"],
                    "aligned_length_bp": "" if hit is None else hit["length"],
                    "query_length_bp": "" if hit is None else hit["qlen"],
                    "subject_length_bp": "" if hit is None else hit["slen"],
                    "query_coverage_pct": "" if hit is None else hit["qcovhsp"],
                    "bitscore": "" if hit is None else hit["bitscore"],
                })
    finally:
        for path in tmp_dir.glob("*"):
            path.unlink()
        tmp_dir.rmdir()

    with OUT_LEVELS.open("w", newline="") as handle:
        fields = [
            "level", "region", "question", "sample_id", "group", "st",
            "pident_vs_23EV612", "aligned_length_bp", "query_length_bp",
            "subject_length_bp", "query_coverage_pct", "bitscore",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(identity_rows)

    # Boundaries of high-identity blocks in ST2187 vs 23EV612 extended context.
    tmp_query = BASE / "_ref_flank23EV612.tmp.fasta"
    tmp_subject = BASE / "_cepa18_flank.tmp.fasta"
    write_fasta(tmp_query, [("23EV612_flank20kb", flanks["23EV612"])])
    write_fasta(tmp_subject, [("cepa18_flank20kb", flanks["cepa_18_barcode02"])])
    hits = best_blast(tmp_query, tmp_subject, ["-max_hsps", "10000"])
    tmp_query.unlink()
    tmp_subject.unlink()
    high = [
        h for h in hits
        if float(h["pident"]) >= 95 and int(h["length"]) >= 500
    ]
    q_intervals = merge_intervals([(int(h["qstart"]), int(h["qend"])) for h in high])
    s_intervals = merge_intervals([(int(h["sstart"]), int(h["send"])) for h in high])
    total_q = sum(end - start + 1 for start, end in q_intervals)
    total_s = sum(end - start + 1 for start, end in s_intervals)
    boundary_rows = [{
        "comparison": "23EV612_vs_cepa_18_ST2187",
        "threshold": "pident>=95,length>=500",
        "ref_23EV612_query_intervals": ";".join(f"{a}-{b}" for a, b in q_intervals),
        "cepa18_subject_intervals": ";".join(f"{a}-{b}" for a, b in s_intervals),
        "ref_conserved_bp": str(total_q),
        "cepa18_conserved_bp": str(total_s),
        "ref_flank_length_bp": str(len(flanks["23EV612"])),
        "cepa18_flank_length_bp": str(len(flanks["cepa_18_barcode02"])),
        "ref_conserved_fraction_pct": f"{total_q / len(flanks['23EV612']) * 100:.1f}",
        "cepa18_conserved_fraction_pct": f"{total_s / len(flanks['cepa_18_barcode02']) * 100:.1f}",
        "hsp_count": str(len(high)),
    }]
    with OUT_BOUNDARIES.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(boundary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(boundary_rows)

    repeat_rows: list[dict[str, str]] = []
    for sample, module in modules.items():
        flank = flanks[sample]
        idx = flank.find(module)
        if idx < 0:
            left = ""
            right = ""
            status = "module_not_found_in_normalized_flank"
        else:
            left = flank[max(0, idx - 60):idx]
            right = flank[idx + len(module):idx + len(module) + 60]
            status = "ok"
        dr, dr_len = find_direct_repeat(left, right)
        tir, tir_len = find_terminal_inverted_repeat(module)
        repeat_rows.append({
            "sample_id": sample,
            "group": synteny.get(sample, {}).get("group", ""),
            "st": synteny.get(sample, {}).get("st", ""),
            "status": status,
            "left_60bp_before_module": left,
            "right_60bp_after_module": right,
            "best_direct_repeat_len_bp": str(dr_len),
            "best_direct_repeat_seq": dr,
            "best_terminal_inverted_repeat_len_bp": str(tir_len),
            "best_terminal_inverted_repeat_left_seq": tir,
        })
    with OUT_REPEATS.open("w", newline="") as handle:
        fields = [
            "sample_id", "group", "st", "status", "left_60bp_before_module",
            "right_60bp_after_module", "best_direct_repeat_len_bp", "best_direct_repeat_seq",
            "best_terminal_inverted_repeat_len_bp", "best_terminal_inverted_repeat_left_seq",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(repeat_rows)

    seh_patterns = [re.compile(r"\bseh\b", re.I), re.compile(r"enterotoxin type H", re.I)]
    ref_module_public = BASE / "_public_ref_module.tmp.fasta"
    write_fasta(ref_module_public, [("23EV612_IS3_spacer_seh", ref_module)])
    public_rows: list[dict[str, str]] = []
    for acc in PUBLIC_ST2187:
        bakta = PROJECT / "results" / "public_context_regional" / "bakta" / acc / f"{acc}.tsv"
        fasta = PROJECT / "results" / "public_context_regional" / "bakta" / acc / f"{acc}.fna"
        amr = PROJECT / "results" / "public_context_regional" / "amrfinder" / f"{acc}.amrfinder_plus.protein.tsv"
        bakta_seh = feature_search(bakta, seh_patterns)
        amr_seh = feature_search(amr, seh_patterns)
        blast_hit = top_hit(best_blast(ref_module_public, fasta)) if fasta.exists() else None
        blast_like_module = bool(
            blast_hit
            and float(blast_hit["pident"]) >= 95
            and float(blast_hit["qcovhsp"]) >= 80
        )
        public_rows.append({
            "assembly_accession": acc,
            "st": "2187",
            "country": "USA:Maine",
            "host": "Bos_taurus",
            "source": "milk",
            "seh_by_bakta": "TRUE" if bakta_seh else "FALSE",
            "seh_by_amrfinder": "TRUE" if amr_seh else "FALSE",
            "is3_seh_module_like_23EV612_by_blast": "TRUE" if blast_like_module else "FALSE",
            "best_module_blast_pident": "" if blast_hit is None else blast_hit["pident"],
            "best_module_blast_aligned_length_bp": "" if blast_hit is None else blast_hit["length"],
            "best_module_blast_query_coverage_pct": "" if blast_hit is None else blast_hit["qcovhsp"],
            "seh_screen_interpretation": "seh_present" if bakta_seh or amr_seh or blast_like_module else "seh_absent",
            "bakta_tsv": str(bakta.relative_to(PROJECT)),
            "amrfinder_tsv": str(amr.relative_to(PROJECT)),
        })
    ref_module_public.unlink(missing_ok=True)
    with OUT_PUBLIC.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(public_rows[0].keys()))
        writer.writeheader()
        writer.writerows(public_rows)

    print(f"Wrote {OUT_LEVELS.relative_to(PROJECT)}")
    print(f"Wrote {OUT_BOUNDARIES.relative_to(PROJECT)}")
    print(f"Wrote {OUT_REPEATS.relative_to(PROJECT)}")
    print(f"Wrote {OUT_PUBLIC.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
