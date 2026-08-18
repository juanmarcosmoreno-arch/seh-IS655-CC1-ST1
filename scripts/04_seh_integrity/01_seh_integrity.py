#!/usr/bin/env python3
"""Audit seh coding integrity and its immediately upstream regulatory region.

This is a sequence-conservation analysis. Candidate Shine-Dalgarno and sigma-A
motifs are computational annotations and must not be interpreted as evidence of
transcription or translation.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = ROOT / "results/seh_st1_bovine/final_tables/final_spacer_seh_sequences.fasta"
DEFAULT_OUTDIR = ROOT / "results/seh_st1_bovine/seh_integrity"
COMP = str.maketrans("ACGTN", "TGCAN")

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    header: str | None = None
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines() + [">"]:
        if raw.startswith(">"):
            if header is not None:
                sample, feature, *_ = header.split("|")
                records.setdefault(sample, {})[feature] = "".join(parts).upper()
            header, parts = raw[1:], []
        else:
            parts.append(raw.strip())
    return records


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3))


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def valid_orf(seq: str) -> bool:
    protein = translate(seq)
    return len(seq) % 3 == 0 and seq[:3] == "ATG" and seq[-3:] in {"TAA", "TAG", "TGA"} and protein[:-1].count("*") == 0


def orient_to_seh(features: dict[str, str]) -> dict[str, str]:
    """Orient both extracted features in the seh coding direction."""
    seh, spacer = features["seh"], features["spacer"]
    if valid_orf(seh):
        return {"seh": seh, "spacer": spacer}
    if valid_orf(rc(seh)):
        return {"seh": rc(seh), "spacer": rc(spacer)}
    return {"seh": seh, "spacer": spacer}


def mismatches(a: str, b: str) -> list[int]:
    return [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]


def motif_mismatches(seq: str, motif: str) -> int:
    return sum(a != b for a, b in zip(seq, motif))


def find_best_motif(seq: str, motif: str, max_mismatches: int) -> tuple[int | None, str, int | None]:
    best: tuple[int, int, str] | None = None
    for i in range(len(seq) - len(motif) + 1):
        observed = seq[i : i + len(motif)]
        mm = motif_mismatches(observed, motif)
        candidate = (mm, -i, observed)
        if mm <= max_mismatches and (best is None or candidate < best):
            best = candidate
    if best is None:
        return None, "", None
    mm, neg_i, observed = best
    return -neg_i + 1, observed, mm


def find_sigma_a_pair(seq: str) -> tuple[str, str, str, str]:
    # Search only the proximal 250 nt. Coordinates are 1-based within spacer.
    proximal = seq[-250:]
    offset = len(seq) - len(proximal)
    candidates: list[tuple[int, int, int, str, str]] = []
    for i in range(len(proximal) - 6):
        minus35 = proximal[i : i + 6]
        mm35 = motif_mismatches(minus35, "TTGACA")
        if mm35 > 2:
            continue
        for spacer_len in range(15, 20):
            j = i + 6 + spacer_len
            if j + 6 > len(proximal):
                continue
            minus10 = proximal[j : j + 6]
            mm10 = motif_mismatches(minus10, "TATAAT")
            if mm10 <= 2:
                candidates.append((mm35 + mm10, i, j, minus35, minus10))
    if not candidates:
        return "NA", "", "NA", ""
    _, i, j, minus35, minus10 = min(candidates)
    return str(offset + i + 1), minus35, str(offset + j + 1), minus10


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--reference", default="23EV612")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    records = {sample: orient_to_seh(features) for sample, features in read_fasta(args.fasta).items()}
    if args.reference not in records:
        raise SystemExit(f"Reference {args.reference!r} not present")
    ref_seh = records[args.reference]["seh"]
    ref_spacer = records[args.reference]["spacer"]
    ref_protein = translate(ref_seh)

    coding_rows: list[dict[str, object]] = []
    regulatory_rows: list[dict[str, object]] = []
    for sample, features in records.items():
        seh, spacer = features["seh"], features["spacer"]
        protein = translate(seh)
        nt_changes = mismatches(ref_seh, seh)
        aa_changes = mismatches(ref_protein, protein)
        internal_stops = protein[:-1].count("*")
        rbs_pos, rbs_seq, rbs_mm = find_best_motif(spacer[-40:], "AGGAGG", 1)
        if rbs_pos is not None:
            rbs_pos = len(spacer) - 40 + rbs_pos
        m35_pos, m35_seq, m10_pos, m10_seq = find_sigma_a_pair(spacer)
        coding_rows.append(
            {
                "sample_id": sample,
                "seh_length_bp": len(seh),
                "length_multiple_of_3": len(seh) % 3 == 0,
                "start_codon": seh[:3],
                "terminal_codon": seh[-3:],
                "protein_length_aa_including_stop": len(protein),
                "internal_stop_codons": internal_stops,
                "nt_differences_vs_23EV612": len(nt_changes),
                "nt_variant_positions": ",".join(map(str, nt_changes)),
                "aa_differences_vs_23EV612": len(aa_changes),
                "aa_variant_positions": ",".join(map(str, aa_changes)),
                "intact_orf_call": len(seh) % 3 == 0 and seh[:3] == "ATG" and seh[-3:] in {"TAA", "TAG", "TGA"} and internal_stops == 0,
            }
        )
        regulatory_rows.append(
            {
                "sample_id": sample,
                "upstream_length_bp": len(spacer),
                "nt_differences_vs_23EV612": len(mismatches(ref_spacer, spacer)),
                "upstream_variant_positions": ",".join(map(str, mismatches(ref_spacer, spacer))),
                "candidate_rbs_start_in_upstream": rbs_pos if rbs_pos is not None else "NA",
                "candidate_rbs_sequence": rbs_seq,
                "candidate_rbs_mismatches_to_AGGAGG": rbs_mm if rbs_mm is not None else "NA",
                "rbs_to_start_distance_bp": len(spacer) - (int(rbs_pos) + len(rbs_seq) - 1) if rbs_pos is not None else "NA",
                "candidate_minus35_start": m35_pos,
                "candidate_minus35_sequence": m35_seq,
                "candidate_minus10_start": m10_pos,
                "candidate_minus10_sequence": m10_seq,
                "interpretation": "Candidate motifs only; no expression inference",
            }
        )

    write_tsv(args.outdir / "seh_coding_integrity.tsv", coding_rows)
    write_tsv(args.outdir / "seh_upstream_regulatory_audit.tsv", regulatory_rows)
    summary = [
        {"metric": "sequences_audited", "value": len(coding_rows)},
        {"metric": "intact_seh_orfs", "value": sum(row["intact_orf_call"] is True for row in coding_rows)},
        {"metric": "unique_seh_nucleotide_alleles", "value": len({features["seh"] for features in records.values()})},
        {"metric": "unique_SEH_protein_alleles", "value": len({translate(features["seh"]) for features in records.values()})},
        {"metric": "unique_upstream_653nt_alleles", "value": len({features["spacer"] for features in records.values()})},
        {"metric": "candidate_RBS_present", "value": sum(row["candidate_rbs_start_in_upstream"] != "NA" for row in regulatory_rows)},
        {"metric": "most_common_terminal_40nt", "value": Counter(features["spacer"][-40:] for features in records.values()).most_common(1)[0][0]},
    ]
    write_tsv(args.outdir / "seh_integrity_summary.tsv", summary)
    print(f"Audited {len(coding_rows)} seh loci; wrote {args.outdir}")


if __name__ == "__main__":
    main()
