#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter a manifest by truth sequence length")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--truth-dir", required=True)
    parser.add_argument("--max-nt", type=int, default=300)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    truth_dir = Path(args.truth_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    dropped = 0
    rows: list[dict[str, str]] = []
    with manifest_path.open(newline="") as f:
        r = csv.DictReader(f)
        fieldnames = list(r.fieldnames or [])
        for row in r:
            target_id = str(row.get("target_id", "")).strip()
            if not target_id:
                continue
            truth_path = truth_dir / f"{target_id}.json"
            if not truth_path.exists():
                dropped += 1
                continue
            data = json.loads(truth_path.read_text())
            seq = str(data.get("sequence", ""))
            if not seq or len(seq) > int(args.max_nt):
                dropped += 1
                continue
            rows.append({k: str(row.get(k, "")) for k in fieldnames})
            kept += 1

    if not rows:
        raise SystemExit("No rows kept. Did truth build run? Is --max-nt too small?")

    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote filtered manifest to {out_path}: kept={kept} dropped={dropped}")


if __name__ == "__main__":
    main()

