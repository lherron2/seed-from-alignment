#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import pandas as pd


def length_bucket_50(L: int) -> str | None:
    if L < 30 or L > 300:
        return None
    for lo in range(30, 280, 50):
        hi = lo + 49
        if lo <= L <= hi:
            return f"{lo}-{hi}"
    if 280 <= L <= 300:
        return "280-300"
    return None


def main() -> None:
    p = argparse.ArgumentParser(description="Select 50 representative targets (<=300nt) from under400 manifest.")
    p.add_argument("--manifest", required=True, help="Input manifest CSV (under400, no rRNA).")
    p.add_argument("--truth-dir", required=True, help="Truth JSON directory named {target_id}.json.")
    p.add_argument("--out-manifest", required=True, help="Output manifest CSV path.")
    p.add_argument("--seed", type=int, default=0, help="Deterministic RNG seed (default: 0).")
    args = p.parse_args()

    manifest = Path(args.manifest)
    truth_dir = Path(args.truth_dir)
    out_manifest = Path(args.out_manifest)

    df = pd.read_csv(manifest)
    if "target_id" not in df.columns:
        raise SystemExit("manifest missing target_id column")

    records: list[dict] = []
    for _, row in df.iterrows():
        tid = str(row["target_id"])
        truth_path = truth_dir / f"{tid}.json"
        if not truth_path.exists():
            continue
        data = json.loads(truth_path.read_text())
        seq = str(data.get("sequence", ""))
        L = len(seq)
        b = length_bucket_50(L)
        if b is None:
            continue
        rec = row.to_dict()
        rec["truth_len"] = L
        rec["bucket"] = b
        records.append(rec)

    if not records:
        raise SystemExit("No eligible records found (need truth_len in [30,300]).")

    df2 = pd.DataFrame(records)
    rng = random.Random(int(args.seed))

    # Target counts per 50-nt bucket; biased to include all long RNAs.
    desired = {
        "30-79": 12,
        "80-129": 10,
        "130-179": 10,
        "180-229": 9,
        "230-279": 7,
        "280-300": 2,
    }
    selected_rows: list[dict] = []

    for bucket, n in desired.items():
        sub = df2[df2["bucket"].astype(str) == bucket].copy()
        if len(sub) == 0:
            continue
        idxs = list(sub.index)
        rng.shuffle(idxs)
        take = min(int(n), len(idxs))
        chosen = sub.loc[idxs[:take]].to_dict(orient="records")
        selected_rows.extend(chosen)

    # If we couldn't hit 50 exactly (e.g., insufficient long RNAs), fill from the remaining pool.
    selected_ids = {str(r["target_id"]) for r in selected_rows}
    if len(selected_rows) < 50:
        rest = df2[~df2["target_id"].astype(str).isin(selected_ids)].copy()
        idxs = list(rest.index)
        rng.shuffle(idxs)
        need = 50 - len(selected_rows)
        selected_rows.extend(rest.loc[idxs[:need]].to_dict(orient="records"))

    if len(selected_rows) != 50:
        raise SystemExit(f"Expected 50 rows, got {len(selected_rows)} (check bucket availability).")

    out_df = pd.DataFrame(selected_rows)
    out_df["bucket"] = out_df["bucket"].astype(str)
    out_df["len_bucket_50"] = out_df["bucket"].astype(str)
    out_df = out_df.sort_values(["len_bucket_50", "truth_len", "target_id"])
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_manifest, index=False)
    print(f"Wrote {len(out_df)} targets to {out_manifest}")
    print(out_df["len_bucket_50"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
