#!/usr/bin/env python3
"""Run MLST on arbitrary public FASTA files and flag ST1 candidates."""

from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import subprocess
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = PROJECT_DIR / "results" / "seh_st1_bovine" / "public_mlst_screen"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta-glob",
        action="append",
        required=True,
        help="FASTA glob. Can be used multiple times, e.g. 'results/public*/fasta/*.fna'.",
    )
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Output directory.")
    parser.add_argument("--mlst-bin", default="", help="Path to mlst executable.")
    parser.add_argument("--threads", default="4", help="Threads passed to mlst when supported.")
    return parser.parse_args()


def find_mlst(explicit: str) -> str:
    candidates = [
        explicit,
        str(PROJECT_DIR / ".micromamba" / "envs" / "gridion_mlst" / "bin" / "mlst"),
        "mlst",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        found = shutil.which(candidate)
        if found:
            return found
        path = Path(candidate)
        if path.exists() and path.is_file():
            return str(path)
    raise SystemExit("No encuentro el ejecutable mlst. Pasa --mlst-bin=/ruta/a/mlst")


def sample_id_from_fasta(path: Path) -> str:
    name = path.name
    for suffix in (".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def quality(scheme: str, st: str) -> str:
    scheme_l = scheme.strip().lower()
    st_l = st.strip().lower()
    if scheme_l in {"", "-", "na", "nan", "unknown", "none"}:
        return "sin_esquema"
    if st_l in {"", "-", "0", "na", "nan", "unknown", "none", "?"}:
        return "perfil_sin_ST"
    return "ok"


def run_mlst(mlst_bin: str, fastas: list[Path], threads: str, raw_path: Path, stderr_path: Path) -> None:
    env = os.environ.copy()
    env["PATH"] = f"{Path(mlst_bin).parent}{os.pathsep}{env.get('PATH', '')}"
    args = [mlst_bin, "--threads", threads, *map(str, fastas)]
    with raw_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        proc = subprocess.run(args, stdout=stdout, stderr=stderr, text=True, env=env)
    if proc.returncode == 0:
        return

    stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    if "thread" in stderr.lower() or "unknown option" in stderr.lower() or "invalid option" in stderr.lower():
        args = [mlst_bin, *map(str, fastas)]
        with raw_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as err:
            proc = subprocess.run(args, stdout=stdout, stderr=err, text=True, env=env)
    if proc.returncode != 0:
        raise SystemExit(f"mlst termino con codigo {proc.returncode}. Ver {stderr_path}")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mlst_bin = find_mlst(args.mlst_bin)

    fastas: list[Path] = []
    for pattern in args.fasta_glob:
        fastas.extend(Path(p) for p in glob.glob(pattern))
    fastas = sorted({p.resolve() for p in fastas if p.exists() and p.is_file()})
    if not fastas:
        raise SystemExit("No encontre FASTA con los glob indicados.")

    raw_path = out_dir / "mlst_raw.tsv"
    stderr_path = out_dir / "mlst.stderr.log"
    summary_path = out_dir / "mlst_summary.tsv"
    st1_path = out_dir / "st1_candidates.tsv"

    run_mlst(mlst_bin, fastas, args.threads, raw_path, stderr_path)

    rows: list[dict[str, str]] = []
    with raw_path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            fields += [""] * max(0, 3 - len(fields))
            fasta, scheme, st = fields[:3]
            alleles = fields[3:]
            sample = sample_id_from_fasta(Path(fasta))
            q = quality(scheme, st)
            rows.append({
                "sample_id": sample,
                "fasta": fasta,
                "scheme": scheme,
                "st": st,
                "mlst_quality": q,
                "is_st1_candidate": "TRUE" if q == "ok" and st == "1" else "FALSE",
                "allele_profile": ";".join(a for a in alleles if a),
            })

    fieldnames = ["sample_id", "fasta", "scheme", "st", "mlst_quality", "is_st1_candidate", "allele_profile"]
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with st1_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows([row for row in rows if row["is_st1_candidate"] == "TRUE"])

    print(f"MLST bin: {mlst_bin}")
    print(f"FASTA screened: {len(fastas)}")
    print(f"ST1 candidates: {sum(row['is_st1_candidate'] == 'TRUE' for row in rows)}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {st1_path}")


if __name__ == "__main__":
    main()
