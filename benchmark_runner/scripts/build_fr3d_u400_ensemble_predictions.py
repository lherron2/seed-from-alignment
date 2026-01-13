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
        if i < idx:
            pairs.append((i, idx))
        else:
            pairs.append((idx, i))
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
    # Maintain each candidate's current min-distance to the selected set.
    # This makes selection O(k*n) jaccard computations instead of O(k*n*|selected|).
    min_dist = [1.0] * n
    for i in range(1, n):
        min_dist[i] = jaccard_distance(items[i].pairs, items[0].pairs)

    def base_score(order: int) -> float:
        if not prefer_early:
            return 1.0
        if n <= 1:
            return 1.0
        return 1.0 - (order / (n - 1))

    while len(selected) < k:
        best_idx: int | None = None
        best_u: float | None = None
        best_tie: tuple[int, str] | None = None
        for idx, it in enumerate(items):
            if it.struct in selected_structs:
                continue
            base = base_score(it.order)
            u = base + diversity_weight * float(min_dist[idx])
            tie = (it.order, it.struct)
            if best_u is None or u > best_u or (u == best_u and (best_tie is None or tie < best_tie)):
                best_idx = idx
                best_u = u
                best_tie = tie

        if best_idx is None:
            break
        best = items[best_idx]
        selected.append(best)
        selected_structs.add(best.struct)
        selected_pairs.append(best.pairs)
        # Update min-distance cache after adding the new structure.
        for i, it in enumerate(items):
            if it.struct in selected_structs:
                continue
            d = jaccard_distance(it.pairs, best.pairs)
            if d < min_dist[i]:
                min_dist[i] = d

    return [it.struct for it in selected]


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def quotas_for_length(length: int, k: int) -> tuple[int, int, int]:
    """Return (q_pipeline, q_linearfold, q_eternafold) summing to k."""
    # Smooth length-adaptive schedule:
    # - short RNAs: keep substantial RNAnneal-ss share, but reserve room for baselines
    # - long RNAs: heavily upweight LF/EF (they scale better and improve @100)
    s = clamp01((float(length) - 150.0) / 150.0)
    q_lf = int(round(float(k) * (0.20 + 0.20 * s)))  # 20 -> 40
    q_ef = int(round(float(k) * (0.25 + 0.20 * s)))  # 25 -> 45
    q_pipe = int(k) - q_lf - q_ef
    if q_pipe < 1:
        # Keep at least one RNAnneal-ss candidate (top-1).
        deficit = 1 - q_pipe
        q_pipe = 1
        # Reduce EF first, then LF.
        take = min(deficit, q_ef)
        q_ef -= take
        deficit -= take
        take = min(deficit, q_lf)
        q_lf -= take
        deficit -= take
    # Guardrails.
    q_pipe = max(0, q_pipe)
    q_lf = max(0, q_lf)
    q_ef = max(0, q_ef)
    # Rebalance to sum exactly.
    total = q_pipe + q_lf + q_ef
    if total != int(k):
        q_pipe += int(k) - total
    return q_pipe, q_lf, q_ef


def quotas_for_length_with_rnastructure(length: int, k: int) -> tuple[int, int, int, int]:
    """Return (q_pipeline, q_linearfold, q_eternafold, q_rnastructure) summing to k.

    Rationale: keep a small but consistent RNAstructure share for diversity/energy-model variety,
    while upweighting LF/EF for longer targets (they typically scale better at best-of-K).
    """
    s = clamp01((float(length) - 150.0) / 150.0)
    q_rs = int(round(float(k) * 0.15))  # constant 15% share
    q_lf = int(round(float(k) * (0.20 + 0.20 * s)))  # 20 -> 40
    q_ef = int(round(float(k) * (0.25 + 0.20 * s)))  # 25 -> 45
    q_pipe = int(k) - q_rs - q_lf - q_ef
    if q_pipe < 1:
        deficit = 1 - q_pipe
        q_pipe = 1
        # Reduce EF first, then LF, then RNAstructure.
        take = min(deficit, q_ef)
        q_ef -= take
        deficit -= take
        take = min(deficit, q_lf)
        q_lf -= take
        deficit -= take
        take = min(deficit, q_rs)
        q_rs -= take
        deficit -= take

    q_pipe = max(0, q_pipe)
    q_lf = max(0, q_lf)
    q_ef = max(0, q_ef)
    q_rs = max(0, q_rs)
    total = q_pipe + q_lf + q_ef + q_rs
    if total != int(k):
        q_pipe += int(k) - total
    return q_pipe, q_lf, q_ef, q_rs


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
            "Build an ensemble predictions directory by merging existing per-target predictions "
            "from RNAnneal-ss, LinearFold-V, and EternaFold, with length-adaptive quotas."
        )
    )
    parser.add_argument("--manifest", required=True, help="CSV with target_id column")
    parser.add_argument("--out-dir", required=True, help="Output predictions directory")
    parser.add_argument("--pipeline-preds", required=True, help="RNAnneal-ss predictions directory")
    parser.add_argument("--linearfold-preds", required=True, help="LinearFold-V predictions directory")
    parser.add_argument("--eternafold-preds", required=True, help="EternaFold predictions directory")
    parser.add_argument(
        "--rnastructure-preds",
        default=None,
        help="Optional RNAstructure AllSub predictions directory (if provided, included in the ensemble).",
    )
    parser.add_argument("--k", type=int, default=100, help="Number of structures per target (default: 100)")
    args = parser.parse_args()

    manifest = Path(args.manifest)
    out_root = Path(args.out_dir)
    pipe_root = Path(args.pipeline_preds)
    lf_root = Path(args.linearfold_preds)
    ef_root = Path(args.eternafold_preds)
    rs_root = Path(args.rnastructure_preds) if args.rnastructure_preds else None
    k = int(args.k)
    if k <= 0:
        raise SystemExit("--k must be > 0")

    out_root.mkdir(parents=True, exist_ok=True)

    target_ids = load_manifest_target_ids(manifest)
    missing = 0
    wrote = 0

    for target_id in target_ids:
        safe_id = target_id.replace("|", "_")
        pipe_db = pipe_root / safe_id / "predictions.db"
        lf_db = lf_root / safe_id / "predictions.db"
        ef_db = ef_root / safe_id / "predictions.db"
        rs_db = rs_root / safe_id / "predictions.db" if rs_root else None
        if not (pipe_db.exists() and lf_db.exists() and ef_db.exists()):
            missing += 1
            continue
        if rs_db is not None and not rs_db.exists():
            missing += 1
            continue

        seq_p, pipe_structs = read_predictions_db(pipe_db)
        seq_l, lf_structs = read_predictions_db(lf_db)
        seq_e, ef_structs = read_predictions_db(ef_db)
        seq_r, rs_structs = ("", [])
        if rs_db is not None:
            seq_r, rs_structs = read_predictions_db(rs_db)
        seq = seq_p or seq_l or seq_e
        if not seq:
            missing += 1
            continue
        if seq_l and seq_l != seq:
            # Prefer the pipeline's sequence line when present.
            pass
        if seq_e and seq_e != seq:
            pass
        if seq_r and seq_r != seq:
            pass

        L = len(seq)
        pipe_structs = unique_structs(pipe_structs, L)
        lf_structs = unique_structs(lf_structs, L)
        ef_structs = unique_structs(ef_structs, L)
        rs_structs = unique_structs(rs_structs, L) if rs_db is not None else []

        if rs_db is None:
            q_pipe, q_lf, q_ef = quotas_for_length(L, k)
            q_rs = 0
        else:
            q_pipe, q_lf, q_ef, q_rs = quotas_for_length_with_rnastructure(L, k)

        # Within-source selection: keep quality by preferring early items for RNAnneal-ss and LF,
        # but strongly prioritize diversity for EternaFold samples (best-of-100 is often late).
        pipe_sel = select_diverse(pipe_structs, q_pipe, prefer_early=True, diversity_weight=0.6)
        lf_sel = select_diverse(lf_structs, q_lf, prefer_early=True, diversity_weight=0.6)
        ef_sel = select_diverse(ef_structs, q_ef, prefer_early=False, diversity_weight=1.0)
        rs_sel = select_diverse(rs_structs, q_rs, prefer_early=True, diversity_weight=0.8) if q_rs else []

        selected: list[str] = []
        seen: set[str] = set()

        def add_many(structs: list[str]) -> None:
            for s in structs:
                if s in seen:
                    continue
                seen.add(s)
                selected.append(s)

        # Preserve RNAnneal-ss top-1 as the ensemble top-1 when available.
        if pipe_structs:
            add_many([pipe_structs[0]])
        add_many(pipe_sel)
        add_many(lf_sel)
        add_many(ef_sel)
        add_many(rs_sel)

        # Fill remaining slots (e.g. after cross-source dedup) by global diversity.
        if len(selected) < k:
            pool: list[str] = []
            pool.extend(pipe_structs)
            pool.extend(lf_structs)
            pool.extend(ef_structs)
            pool.extend(rs_structs)
            pool = unique_structs(pool, L)

            pool_items = [StructItem(s, pairs_from_dotbracket(s), i) for i, s in enumerate(pool)]
            selected_pairs = [pairs_from_dotbracket(s) for s in selected]
            selected_structs_set = set(selected)
            # Cache each pool candidate's min-distance to the current selected set.
            pool_min_dist = [1.0] * len(pool_items)
            if selected_pairs:
                for i, it in enumerate(pool_items):
                    if it.struct in selected_structs_set:
                        pool_min_dist[i] = 0.0
                        continue
                    md = 1.0
                    for sp in selected_pairs:
                        d = jaccard_distance(it.pairs, sp)
                        if d < md:
                            md = d
                            if md <= 0.0:
                                break
                    pool_min_dist[i] = md

            while len(selected) < k and len(seen) < len(pool_items):
                best_idx: int | None = None
                best_u: float | None = None
                best_tie: tuple[int, str] | None = None
                for idx, it in enumerate(pool_items):
                    if it.struct in seen:
                        continue
                    base = 1.0 if len(pool_items) <= 1 else 1.0 - (it.order / (len(pool_items) - 1))
                    u = base + 0.8 * float(pool_min_dist[idx])
                    tie = (it.order, it.struct)
                    if best_u is None or u > best_u or (u == best_u and (best_tie is None or tie < best_tie)):
                        best_idx = idx
                        best_u = u
                        best_tie = tie
                if best_idx is None:
                    break
                best = pool_items[best_idx]
                seen.add(best.struct)
                selected.append(best.struct)
                selected_pairs.append(best.pairs)
                # Update min-distance cache after adding the new structure.
                for i, it in enumerate(pool_items):
                    if it.struct in seen:
                        continue
                    d = jaccard_distance(it.pairs, best.pairs)
                    if d < pool_min_dist[i]:
                        pool_min_dist[i] = d

        if not selected:
            selected = ["." * L]

        # Pad to exactly k (duplicates are fine for @100 metrics).
        if len(selected) < k:
            pad = selected[0]
            selected.extend([pad] * (k - len(selected)))
        else:
            selected = selected[:k]

        # Write per-target output in ssbench format.
        out_fa = out_root / f"{safe_id}.fa"
        out_fa.write_text(f">{target_id}\n{seq}\n")
        out_dir = out_root / safe_id
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "predictions.db").write_text("\n".join([seq] + selected) + "\n")
        wrote += 1

    print(f"Wrote {wrote} targets to {out_root} (missing={missing})")


if __name__ == "__main__":
    main()
