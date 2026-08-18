#!/usr/bin/env python3
"""Rapidly screen Jahnen et al. CC1 FASTQ files for exact diagnostic seh k-mers.

Positive files are streamed only until a diagnostic k-mer is observed. Files with
no hit are read completely, making this suitable for locating rare seh-negative
isolates without retaining the full public read collection.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import time
import urllib.request
from pathlib import Path


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    ).upper()


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def diagnostic_kmers(seq: str, length: int = 31, count: int = 12) -> set[str]:
    last = len(seq) - length
    starts = sorted({round(i * last / (count - 1)) for i in range(count)})
    kmers = {seq[start : start + length] for start in starts}
    return kmers | {reverse_complement(kmer) for kmer in kmers}


def screen_fastq(url: str, kmers: set[str]) -> tuple[bool, int, int]:
    request = urllib.request.Request(url, headers={"User-Agent": "GridION-seh-screen/1.0"})
    compressed_bytes = 0
    reads = 0
    with urllib.request.urlopen(request, timeout=120) as response:
        class CountingReader:
            def read(self, size=-1):
                nonlocal compressed_bytes
                chunk = response.read(size)
                compressed_bytes += len(chunk)
                return chunk

            def close(self):
                return None

        with gzip.GzipFile(fileobj=CountingReader()) as stream:
            while True:
                header = stream.readline()
                if not header:
                    return False, reads, compressed_bytes
                sequence = stream.readline().decode("ascii", "ignore").strip().upper()
                stream.readline()
                stream.readline()
                reads += 1
                if any(kmer in sequence for kmer in kmers):
                    return True, reads, compressed_bytes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accessions", type=Path, required=True)
    parser.add_argument("--ena-files", type=Path, required=True)
    parser.add_argument("--seh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.accessions.open(newline="") as handle:
        accessions = list(csv.DictReader(handle, delimiter="\t"))
    cc1 = [row for row in accessions if row.get("lfd.Nr.", "").rstrip("AB").isdigit()
           and int(row["lfd.Nr."].rstrip("AB")) <= 73]
    with args.ena_files.open(newline="") as handle:
        ena = {row["sample_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}

    kmers = diagnostic_kmers(read_fasta(args.seh))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "paper_id", "short_id", "biosample", "run_accession", "seh_kmer_status",
        "hit_file", "reads_examined", "compressed_bytes_streamed", "note",
    ]
    done = set()
    if args.output.exists():
        with args.output.open(newline="") as handle:
            done = {row["biosample"] for row in csv.DictReader(handle, delimiter="\t")}
    mode = "a" if done else "w"
    with args.output.open(mode, newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        if not done:
            writer.writeheader()
        for index, isolate in enumerate(cc1, start=1):
            biosample = isolate["BioSample"]
            if biosample in done:
                continue
            run = ena.get(biosample)
            result = {
                "paper_id": isolate["ID Paper"], "short_id": isolate["lfd.Nr."],
                "biosample": biosample, "run_accession": run["run_accession"] if run else "",
                "seh_kmer_status": "unresolved", "hit_file": "", "reads_examined": 0,
                "compressed_bytes_streamed": 0, "note": "",
            }
            if not run or not run.get("fastq_ftp"):
                result["note"] = "No public ENA FASTQ URL"
            else:
                total_reads = total_bytes = 0
                for ftp in run["fastq_ftp"].split(";"):
                    url = "https://" + ftp
                    try:
                        hit, reads, nbytes = screen_fastq(url, kmers)
                    except Exception as exc:  # retain resumable audit row
                        result["note"] = f"{type(exc).__name__}: {exc}"
                        break
                    total_reads += reads
                    total_bytes += nbytes
                    if hit:
                        result["seh_kmer_status"] = "positive"
                        result["hit_file"] = ftp.rsplit("/", 1)[-1]
                        break
                else:
                    result["seh_kmer_status"] = "negative"
                result["reads_examined"] = total_reads
                result["compressed_bytes_streamed"] = total_bytes
            writer.writerow(result)
            out_handle.flush()
            print(f"[{index}/{len(cc1)}] {biosample}: {result['seh_kmer_status']}", flush=True)
            time.sleep(0.15)


if __name__ == "__main__":
    main()
