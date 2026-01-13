#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def bucket_group(bucket: str) -> str:
    bucket = str(bucket)
    if bucket in {"30-79", "80-129"}:
        return "short"
    if bucket in {"130-179", "180-229"}:
        return "medium"
    if bucket in {"230-279", "280-300"}:
        return "long"
    return "unknown"


def main() -> None:
    p = argparse.ArgumentParser(description="Select a stratified subset for hyperparameter ablations.")
    p.add_argument("--in-manifest", required=True, help="Input manifest CSV (must have target_id,bucket).")
    p.add_argument("--out-manifest", required=True, help="Output manifest CSV.")
    p.add_argument("--n-per-group", type=int, default=10, help="Targets per length group (default: 10).")
    args = p.parse_args()

    in_path = Path(args.in_manifest)
    out_path = Path(args.out_manifest)

    df = pd.read_csv(in_path)
    if "target_id" not in df.columns:
        raise SystemExit("manifest missing target_id column")
    if "bucket" not in df.columns:
        raise SystemExit("manifest missing bucket column (expected 50-nt buckets)")

    df = df.copy()
    df["len_group"] = df["bucket"].astype(str).map(bucket_group)
    df = df[df["len_group"].isin(["short", "medium", "long"])].copy()
    if df.empty:
        raise SystemExit("No rows in short/medium/long groups after filtering")

    n = int(args.n_per_group)
    if n <= 0:
        raise SystemExit("--n-per-group must be > 0")

    picked: list[pd.DataFrame] = []
    for g in ["short", "medium", "long"]:
        sub = df[df["len_group"] == g].sort_values(["bucket", "target_id"])
        if len(sub) < n:
            raise SystemExit(f"Not enough targets in group {g}: have {len(sub)}, need {n}")
        picked.append(sub.head(n))

    out = pd.concat(picked, ignore_index=True).sort_values(["len_group", "bucket", "target_id"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} targets to {out_path}")
    print(out["len_group"].value_counts().to_string())


if __name__ == "__main__":
    main()

