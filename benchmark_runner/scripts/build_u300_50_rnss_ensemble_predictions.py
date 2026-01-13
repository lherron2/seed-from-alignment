#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


def parse_dotbracket(text: str) -> tuple[str, list[str]]:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return "", []
    if lines[0].startswith(">"):
        if len(lines) < 3:
            return "", []
        return lines[1], lines[2:]
    if len(lines) >= 2:
        return lines[0], lines[1:]
    return "", []


def pairs_from_dotbracket(struct: str) -> frozenset[tuple[int, int]]:
    open_to_close = {"(": ")", "[": "]", "{": "}", "<": ">"}
    for i in range(26):
        open_to_close[chr(ord("a") + i)] = chr(ord("A") + i)
    close_to_open = {v: k for k, v in open_to_close.items()}
    stacks: dict[str, list[int]] = {op: [] for op in open_to_close}
    pairs: list[tuple[int, int]] = []
    for idx, ch in enumerate(struct):
        if ch in open_to_close:
            stacks[ch].append(idx)
            continue
        op = close_to_open.get(ch)
        if op is None:
            continue
        if not stacks[op]:
            continue
        i = stacks[op].pop()
        pairs.append((i, idx) if i < idx else (idx, i))
    return frozenset(pairs)


def jaccard_distance(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return 1.0 - (inter / union if union else 0.0)


@dataclass(frozen=True)
class StructItem:
    struct: str
    pairs: frozenset[tuple[int, int]]
    order: int


def unique_structs(structs: list[str], length: int) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for s in structs:
        if len(s) != length:
            continue
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def select_diverse(
    structs: list[str],
    k: int,
    *,
    prefer_early: bool,
    diversity_weight: float,
) -> list[str]:
    if k <= 0 or not structs:
        return []
    if k >= len(structs):
        return list(structs)

    items = [StructItem(s, pairs_from_dotbracket(s), i) for i, s in enumerate(structs)]
    n = len(items)

    selected: list[StructItem] = [items[0]]
    selected_structs: set[str] = {items[0].struct}
    selected_pairs = [items[0].pairs]

    def base_score(order: int) -> float:
        if not prefer_early:
            return 1.0
        if n <= 1:
            return 1.0
        return 1.0 - (order / (n - 1))

    while len(selected) < k:
        best: StructItem | None = None
        best_u: float | None = None
        best_tie: tuple[int, str] | None = None

        for it in items:
            if it.struct in selected_structs:
                continue
            base = base_score(it.order)
            min_dist = min(jaccard_distance(it.pairs, sp) for sp in selected_pairs) if selected_pairs else 1.0
            u = base + diversity_weight * min_dist
            tie = (it.order, it.struct)
            if best_u is None or u > best_u or (u == best_u and (best_tie is None or tie < best_tie)):
                best = it
                best_u = u
                best_tie = tie

        if best is None:
            break
        selected.append(best)
        selected_structs.add(best.struct)
        selected_pairs.append(best.pairs)

    return [it.struct for it in selected]


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def quotas_for_length(length: int, k: int) -> tuple[int, int, int]:
    """Return (q_rs, q_lf, q_ef) summing to k.

    Heuristic for representative50 (<=300nt):
    - Keep a large share of RNss(EF) structures to avoid dropping EF's best-of-100 winner.
    - Always reserve a non-trivial slice for RNss(LF) and RNss(RS) to capture complementary wins.
    - For longer RNAs, upweight RNss(EF) further (LF/EF scale better vs RS).
    """
    s = clamp01((float(length) - 150.0) / 150.0)
    q_ef = int(round(float(k) * (0.70 + 0.20 * s)))  # 70 -> 90
    q_lf = int(round(float(k) * (0.15 - 0.05 * s)))  # 15 -> 10
    q_rs = int(k) - q_ef - q_lf

    # Guardrails: keep some representation for RS/LF even at long lengths.
    q_lf = max(8, q_lf)
    q_rs = max(8, q_rs)
    # Rebalance (take from EF first).
    total = q_rs + q_lf + q_ef
    if total > int(k):
        drop = total - int(k)
        take = min(drop, q_ef - 1)
        q_ef -= take
        drop -= take
        take = min(drop, q_lf - 8)
        q_lf -= take
        drop -= take
        q_rs = int(k) - q_ef - q_lf
    if q_rs + q_lf + q_ef != int(k):
        q_ef = int(k) - q_rs - q_lf
    if q_ef < 1:
        q_ef = 1
        q_rs = max(0, int(k) - q_ef - q_lf)
    return q_rs, q_lf, q_ef


def load_manifest_target_ids(path: Path) -> list[str]:
    with path.open(newline="") as f:
        r = csv.DictReader(f)
        if not r.fieldnames or "target_id" not in r.fieldnames:
            raise SystemExit(f"manifest missing target_id: {path}")
        out: list[str] = []
        for row in r:
            tid = (row.get("target_id") or "").strip()
            if tid:
                out.append(tid)
        return out


def read_predictions_db(path: Path) -> tuple[str, list[str]]:
    seq, structs = parse_dotbracket(path.read_text())
    return seq, structs


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build an ensemble predictions directory by merging per-target predictions from "
            "RNAnneal-ss scaffold variants (RNss(RS)/RNss(LF)/RNss(EF)) with length-adaptive quotas "
            "that preserve most RNss(EF) candidates."
        )
    )
    parser.add_argument("--manifest", required=True, help="CSV with target_id column")
    parser.add_argument("--out-dir", required=True, help="Output predictions directory")
    parser.add_argument("--rnss-rs-preds", required=True, help="RNss(RS) predictions directory")
    parser.add_argument("--rnss-lf-preds", required=True, help="RNss(LF) predictions directory")
    parser.add_argument("--rnss-ef-preds", required=True, help="RNss(EF) predictions directory")
    parser.add_argument("--k", type=int, default=100, help="Number of structures per target (default: 100)")
    args = parser.parse_args()

    manifest = Path(args.manifest)
    out_root = Path(args.out_dir)
    rs_root = Path(args.rnss_rs_preds)
    lf_root = Path(args.rnss_lf_preds)
    ef_root = Path(args.rnss_ef_preds)
    k = int(args.k)
    if k <= 0:
        raise SystemExit("--k must be > 0")

    out_root.mkdir(parents=True, exist_ok=True)
    target_ids = load_manifest_target_ids(manifest)
    missing = 0
    wrote = 0

    for target_id in target_ids:
        safe_id = target_id.replace("|", "_")
        rs_db = rs_root / safe_id / "predictions.db"
        lf_db = lf_root / safe_id / "predictions.db"
        ef_db = ef_root / safe_id / "predictions.db"
        if not (rs_db.exists() and lf_db.exists() and ef_db.exists()):
            missing += 1
            continue

        seq_r, rs_structs = read_predictions_db(rs_db)
        seq_l, lf_structs = read_predictions_db(lf_db)
        seq_e, ef_structs = read_predictions_db(ef_db)
        seq = seq_e or seq_r or seq_l
        if not seq:
            missing += 1
            continue

        L = len(seq)
        rs_structs = unique_structs(rs_structs, L)
        lf_structs = unique_structs(lf_structs, L)
        ef_structs = unique_structs(ef_structs, L)

        q_rs, q_lf, q_ef = quotas_for_length(L, k)

        rs_sel = select_diverse(rs_structs, q_rs, prefer_early=True, diversity_weight=0.6)
        lf_sel = select_diverse(lf_structs, q_lf, prefer_early=True, diversity_weight=0.6)
        ef_sel = select_diverse(ef_structs, q_ef, prefer_early=False, diversity_weight=1.0)

        selected: list[str] = []
        seen: set[str] = set()

        def add_many(structs: list[str]) -> None:
            for s in structs:
                if s in seen:
                    continue
                seen.add(s)
                selected.append(s)

        # Preserve RNss(EF) top-1 as the ensemble top-1.
        if ef_structs:
            add_many([ef_structs[0]])
        add_many(ef_sel)
        add_many(lf_sel)
        add_many(rs_sel)

        # Fill remaining slots via global diversity over the merged pool.
        if len(selected) < k:
            pool: list[str] = []
            pool.extend(ef_structs)
            pool.extend(lf_structs)
            pool.extend(rs_structs)
            pool = unique_structs(pool, L)

            pool_items = [StructItem(s, pairs_from_dotbracket(s), i) for i, s in enumerate(pool)]
            selected_pairs = [pairs_from_dotbracket(s) for s in selected]

            while len(selected) < k and len(seen) < len(pool_items):
                best: StructItem | None = None
                best_u: float | None = None
                best_tie: tuple[int, str] | None = None
                for it in pool_items:
                    if it.struct in seen:
                        continue
                    base = 1.0 if len(pool_items) <= 1 else 1.0 - (it.order / (len(pool_items) - 1))
                    min_dist = min(jaccard_distance(it.pairs, sp) for sp in selected_pairs) if selected_pairs else 1.0
                    u = base + 0.8 * min_dist
                    tie = (it.order, it.struct)
                    if best_u is None or u > best_u or (u == best_u and (best_tie is None or tie < best_tie)):
                        best = it
                        best_u = u
                        best_tie = tie
                if best is None:
                    break
                seen.add(best.struct)
                selected.append(best.struct)
                selected_pairs.append(best.pairs)

        if not selected:
            selected = ["." * L]

        if len(selected) < k:
            pad = selected[0]
            selected.extend([pad] * (k - len(selected)))
        else:
            selected = selected[:k]

        (out_root / f"{safe_id}.fa").write_text(f">{target_id}\n{seq}\n")
        out_dir = out_root / safe_id
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "predictions.db").write_text("\n".join([seq] + selected) + "\n")
        wrote += 1

    print(f"Wrote {wrote} targets to {out_root} (missing={missing})")


if __name__ == "__main__":
    main()

