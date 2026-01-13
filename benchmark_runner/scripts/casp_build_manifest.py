#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an ssbench manifest for CASP RNA targets")
    parser.add_argument("--targets", required=True, help="CSV with at least: target_id,pdb_id,chain_id")
    parser.add_argument("--out", required=True, help="Write ssbench manifest CSV")
    args = parser.parse_args()

    in_path = Path(args.targets)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    with in_path.open(newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            target_id = str(row.get("target_id", "")).strip()
            pdb_id = str(row.get("pdb_id", "")).strip()
            chain_id = str(row.get("chain_id", "")).strip()
            if not target_id or not pdb_id or not chain_id:
                continue
            rows.append(
                {
                    "target_id": target_id,
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "split": row.get("split", "casp"),
                    "bucket": row.get("bucket", "casp"),
                    "source_header": row.get("source_header", ""),
                    "match_ratio": row.get("match_ratio", "1.0"),
                }
            )

    if not rows:
        raise SystemExit("No rows found in targets CSV (need target_id,pdb_id,chain_id).")

    # Keep columns compatible with other manifests (extra cols are ignored by ssbench).
    fieldnames = [
        "target_id",
        "pdb_id",
        "chain_id",
        "split",
        "bucket",
        "source_header",
        "match_ratio",
    ]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"Wrote {len(rows)} manifest rows to {out_path}")


if __name__ == "__main__":
    main()

