#!/usr/bin/env python3
"""Minimal predictor stub that outputs an all-unpaired structure."""

import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: echo_predictor.py <fasta> <outdir>")
    fasta = Path(sys.argv[1]).read_text().strip().splitlines()
    if len(fasta) < 2:
        raise SystemExit("Invalid FASTA")
    seq = fasta[1].strip()
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "predictions.db"
    out.write_text(seq + "\n" + ("." * len(seq)) + "\n")


if __name__ == "__main__":
    main()
