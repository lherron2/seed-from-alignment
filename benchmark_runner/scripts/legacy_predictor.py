#!/usr/bin/env python3
"""Wrapper for legacy CaCoFold sampler in this repo."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def load_mapping(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            target_id = row.get("target_id")
            if not target_id:
                continue
            mapping[target_id] = row
    return mapping


def parse_fasta_id(fasta_path: Path) -> str:
    line = fasta_path.read_text().splitlines()[0]
    return line[1:].strip() if line.startswith(">") else line.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description="Legacy CaCoFold predictor wrapper")
    parser.add_argument("--mapping", required=True, help="CSV with target_id, sto_path, cov_path, seq_name")
    parser.add_argument("--n-samples", type=int, default=50)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-loop-sep", type=int, default=3)
    parser.add_argument("fasta", help="Input FASTA path (header contains target_id)")
    parser.add_argument("outdir", help="Output directory")
    args = parser.parse_args()

    mapping = load_mapping(Path(args.mapping))
    target_id = parse_fasta_id(Path(args.fasta))

    if target_id not in mapping:
        raise SystemExit(f"Target {target_id} not found in mapping")

    row = mapping[target_id]
    sto = row.get("sto_path")
    cov = row.get("cov_path")
    seq_name = row.get("seq_name")

    if not sto or not seq_name:
        raise SystemExit("Mapping row must include sto_path and seq_name")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_db = outdir / "predictions.db"

    cmd = [
        sys.executable,
        "-m",
        "src.lib.sample_cacofold_structures",
        sto,
        "--seq",
        seq_name,
        "--n-samples",
        str(args.n_samples),
        "--beta",
        str(args.beta),
        "--seed",
        str(args.seed),
        "--min-loop-sep",
        str(args.min_loop_sep),
        "--out-db",
        str(out_db),
    ]

    if cov:
        cmd.extend([
            "--cov-file",
            cov,
            "--cov-mode",
            "logE_power",
            "--cov-alpha",
            "3.0",
            "--cov-min-power",
            "0.1",
        ])

    subprocess.run(cmd, check=True, cwd=str(Path(__file__).resolve().parents[1]))


if __name__ == "__main__":
    main()
