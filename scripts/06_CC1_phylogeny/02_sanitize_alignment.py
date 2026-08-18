#!/usr/bin/env python3
"""Shorten SKA path-based FASTA headers to stable sample identifiers."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    seen = set()
    taxa = 0
    with args.input.open() as source, args.output.open("w") as target:
        for line in source:
            if line.startswith(">"):
                sample = Path(line[1:].strip().split()[0]).stem
                if sample in seen:
                    raise ValueError(f"Duplicate sample after header shortening: {sample}")
                seen.add(sample)
                taxa += 1
                target.write(f">{sample}\n")
            else:
                target.write(line.upper())
    print(f"Wrote {args.output} with {taxa} unique taxa")


if __name__ == "__main__":
    main()
