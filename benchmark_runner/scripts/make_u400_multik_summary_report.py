#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


def pairs_from_dotbracket(db: str) -> set[tuple[int, int]]:
    stack: list[int] = []
    out: set[tuple[int, int]] = set()
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                continue
            j = stack.pop()
            if j < i:
                out.add((j, i))
            else:
                out.add((i, j))
    return out


def read_predictions_db(path: Path) -> tuple[str, list[str]]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        return "", []
    seq = lines[0]
    structs = [ln for ln in lines[1:] if len(ln) == len(seq)]
    return seq, structs


def fmt_float(x: float) -> str:
    return f"{x:.3f}"


def render_table_tex(
    *,
    out_path: Path,
    caption: str,
    label: str,
    headers: list[str],
    rows: list[list[str]],
    col_spec: str,
) -> None:
    lines: list[str] = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\small")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}")
    lines.append(r"\toprule")
    lines.append(" & ".join(headers) + r" \\")
    lines.append(r"\midrule")
    for row in rows:
        lines.append(" & ".join(row) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    out_path.write_text("\n".join(lines) + "\n")

def write_wide_summary_csv(
    *,
    out_path: Path,
    summary: pd.DataFrame,
    methods: list[MethodSpec],
    ks: list[int],
) -> None:
    rows: list[dict[str, object]] = []
    for k in ks:
        row: dict[str, object] = {"k": int(k)}
        for m in methods:
            sub = summary[(summary["method"] == m.key) & (summary["k"] == int(k))]
            if len(sub) != 1:
                row[f"{m.key}_mean"] = float("nan")
                row[f"{m.key}_median"] = float("nan")
                continue
            row[f"{m.key}_mean"] = float(sub["mean"].iloc[0])
            row[f"{m.key}_median"] = float(sub["median"].iloc[0])
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, index=False)


@dataclass(frozen=True)
class MethodSpec:
    key: str
    label: str
    pred_dir: Path


def main() -> None:
    p = argparse.ArgumentParser(description="u400 multi-K summary (mean/median F1 and MCC) for multiple methods.")
    p.add_argument("--out-dir", required=True, help="Docs report directory (writes generated/*).")
    p.add_argument("--manifest", required=True, help="Manifest CSV (target_id).")
    p.add_argument("--truth-dir", required=True, help="Truth JSON directory named {target_id}.json.")
    p.add_argument("--ensemble-preds", required=True)
    p.add_argument("--linearfold-preds", required=True)
    p.add_argument("--eternafold-preds", required=True)
    p.add_argument("--rnastructure-preds", default=None, help="Optional RNAstructure predictions directory.")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    gen = out_dir / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(Path(args.manifest))
    truth_dir = Path(args.truth_dir)

    methods = [
        MethodSpec("ens", "Ensemble", Path(args.ensemble_preds)),
        MethodSpec("lf", "LF-V", Path(args.linearfold_preds)),
        MethodSpec("ef", "EFold", Path(args.eternafold_preds)),
    ]
    if args.rnastructure_preds:
        methods.append(MethodSpec("rs", "RNAstr", Path(args.rnastructure_preds)))
    ks = [1, 50, 100, 200, 500]

    from ssbench.metrics.pair_metrics import compute_pair_metrics

    # Determine shared target set (fair comparison).
    all_ids: list[str] = [str(t) for t in manifest["target_id"].astype(str).tolist()]
    shared: list[str] = []
    for tid in all_ids:
        truth_path = truth_dir / f"{tid}.json"
        if not truth_path.exists():
            continue
        safe = tid.replace("|", "_")
        if all((m.pred_dir / safe / "predictions.db").exists() for m in methods):
            shared.append(tid)

    if not shared:
        raise SystemExit("No shared targets across methods. Check predictions/truth paths.")

    # Preload truth pairs.
    truth_pairs: dict[str, set[tuple[int, int]]] = {}
    truth_len: dict[str, int] = {}
    for tid in shared:
        data = json.loads((truth_dir / f"{tid}.json").read_text())
        seq = str(data.get("sequence", ""))
        L = len(seq)
        truth_len[tid] = L
        pairs = {tuple(p) for p in data.get("canonical_pairs", [])}
        truth_pairs[tid] = {(int(i), int(j)) for i, j in pairs if 0 <= int(i) < int(j) < L}

    rows: list[dict[str, object]] = []
    for m in methods:
        for tid in shared:
            safe = tid.replace("|", "_")
            pred_path = m.pred_dir / safe / "predictions.db"
            _seq, structs = read_predictions_db(pred_path)
            L = truth_len[tid]
            structs = [s for s in structs if len(s) == L]
            pred_pairs_list = [pairs_from_dotbracket(s) for s in structs]
            if not pred_pairs_list:
                pred_pairs_list = [set()]

            ref = truth_pairs[tid]
            for k in ks:
                n = min(int(k), len(pred_pairs_list))
                best_idx = 0
                best_f1 = -1.0
                best_metrics = None
                for i in range(n):
                    met = compute_pair_metrics(pred_pairs_list[i], ref, L)
                    if met.f1 > best_f1:
                        best_f1 = met.f1
                        best_idx = i
                        best_metrics = met
                if best_metrics is None:
                    best_metrics = compute_pair_metrics(set(), ref, L)
                rows.append(
                    {
                        "target_id": tid,
                        "method": m.key,
                        "k": int(k),
                        "effective_k": int(n),
                        "f1": float(best_metrics.f1),
                        "mcc": float(best_metrics.mcc),
                    }
                )

    df = pd.DataFrame(rows)
    df.to_csv(gen / "multik_metrics.csv", index=False)

    # Summary tables.
    def summarize(metric: str) -> pd.DataFrame:
        out = (
            df.groupby(["method", "k"], as_index=False)[metric]
            .agg(mean="mean", median="median")
            .sort_values(["k", "method"])
        )
        return out

    f1_sum = summarize("f1")
    mcc_sum = summarize("mcc")

    key_to_label = {m.key: m.label for m in methods}

    def table_rows(summary: pd.DataFrame) -> list[list[str]]:
        rows_out: list[list[str]] = []
        for k in ks:
            sub = summary[summary["k"] == k]
            # start with K, then for each method mean/median
            row = [f"@{k}"]
            for m in methods:
                r = sub[sub["method"] == m.key]
                if len(r) != 1:
                    row.extend(["-", "-"])
                    continue
                row.append(fmt_float(float(r["mean"].iloc[0])))
                row.append(fmt_float(float(r["median"].iloc[0])))
            rows_out.append(row)
        return rows_out

    headers = ["K"]
    for m in methods:
        headers.extend([f"{m.label} mean", f"{m.label} med"])

    render_table_tex(
        out_path=gen / "f1_table.tex",
        caption=f"Best-of-K F1 on u400 non-rRNA benchmark (N={len(shared)} targets; oracle within top-K prefix).",
        label="tab:f1",
        headers=headers,
        rows=table_rows(f1_sum),
        col_spec="l" + "rr" * len(methods),
    )

    render_table_tex(
        out_path=gen / "mcc_table.tex",
        caption=f"Best-of-K MCC on u400 non-rRNA benchmark (N={len(shared)} targets; MCC for the structure achieving best F1 within top-K prefix).",
        label="tab:mcc",
        headers=headers,
        rows=table_rows(mcc_sum),
        col_spec="l" + "rr" * len(methods),
    )

    # Curve sources for LaTeX (no external plotting deps required).
    write_wide_summary_csv(out_path=gen / "f1_wide.csv", summary=f1_sum, methods=methods, ks=ks)
    write_wide_summary_csv(out_path=gen / "mcc_wide.csv", summary=mcc_sum, methods=methods, ks=ks)


if __name__ == "__main__":
    main()
