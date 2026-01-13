#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ssbench.metrics.pair_metrics import compute_pair_metrics
from ssbench.predict.parse_dotbracket import pairs_from_dotbracket, parse_dotbracket


def _read_db_structs(db_path: Path) -> tuple[str, list[str]]:
    if not db_path.exists():
        return "", []
    return parse_dotbracket(db_path.read_text())


def _best_f1(dotbrackets: list[str], ref_pairs: set[tuple[int, int]], length: int) -> float:
    best = 0.0
    for s in dotbrackets:
        pairs = set(pairs_from_dotbracket(s))
        m = compute_pair_metrics(pairs, ref_pairs, length)
        if m.f1 > best:
            best = m.f1
    return best


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute oracle best-F1 over existing candidate DBs (refined/sample)"
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--truth-dir", required=True)
    parser.add_argument("--candidates-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest)
    truth_dir = Path(args.truth_dir)
    candidates_root = Path(args.candidates_dir)

    rows: list[dict[str, object]] = []
    for _, row in manifest.iterrows():
        target_id = str(row["target_id"])
        safe_id = target_id.replace("|", "_")
        truth_path = truth_dir / f"{target_id}.json"
        cand_dir = candidates_root / safe_id
        refined_path = cand_dir / "01_refined.db"
        sampled_path = cand_dir / "02_sampled.db"

        if not truth_path.exists():
            continue

        truth = json.loads(truth_path.read_text())
        seq = str(truth["sequence"])
        ref_pairs = {tuple(p) for p in truth["canonical_pairs"]}

        refined_seq, refined_structs = _read_db_structs(refined_path)
        sampled_seq, sampled_structs = _read_db_structs(sampled_path)

        if refined_seq and len(refined_seq) != len(seq):
            raise SystemExit(
                f"{target_id}: refined seq length mismatch ({len(refined_seq)} vs {len(seq)})"
            )
        if sampled_seq and len(sampled_seq) != len(seq):
            raise SystemExit(
                f"{target_id}: sampled seq length mismatch ({len(sampled_seq)} vs {len(seq)})"
            )

        oracle_refined = _best_f1(refined_structs, ref_pairs, len(seq)) if refined_structs else 0.0
        oracle_sampled = _best_f1(sampled_structs, ref_pairs, len(seq)) if sampled_structs else 0.0
        oracle_union = _best_f1(refined_structs + sampled_structs, ref_pairs, len(seq))

        rows.append(
            {
                "target_id": target_id,
                "split": row.get("split", "train"),
                "bucket": row.get("bucket", "unknown"),
                "n_refined": len(refined_structs),
                "n_sampled": len(sampled_structs),
                "oracle_refined_best_f1": oracle_refined,
                "oracle_sampled_best_f1": oracle_sampled,
                "oracle_union_best_f1": oracle_union,
            }
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise SystemExit("No rows produced. Check --manifest/--truth-dir/--candidates-dir.")

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    mean_union = float(df["oracle_union_best_f1"].mean())
    mean_refined = float(df["oracle_refined_best_f1"].mean())
    mean_sampled = float(df["oracle_sampled_best_f1"].mean())
    print(f"Wrote oracle candidate metrics to {out_path} ({len(df)} rows)")
    print(
        f"Mean oracle F1: union={mean_union:.3f} refined={mean_refined:.3f} sampled={mean_sampled:.3f}"
    )


if __name__ == "__main__":
    main()

