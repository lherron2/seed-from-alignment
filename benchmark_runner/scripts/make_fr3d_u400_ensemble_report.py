#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class MethodRow:
    key: str
    name: str
    metrics_csv: Path


def fmt_float(x: float) -> str:
    return f"{x:.3f}"


def frac_fmt(x: float) -> str:
    return f"{100.0 * x:.1f}\\%"


def compute_summary(df: pd.DataFrame) -> dict[str, float]:
    out: dict[str, float] = {}
    out["n"] = float(len(df))
    out["mean_top1_f1"] = float(df["f1"].mean())
    out["median_top1_f1"] = float(df["f1"].median())
    out["mean_best100_f1"] = float(df["best_of_k_f1"].mean())
    out["median_best100_f1"] = float(df["best_of_k_f1"].median())
    out["fail_top1"] = float((df["f1"] <= 0.0).mean())
    out["fail_best100"] = float((df["best_of_k_f1"] <= 0.0).mean())
    out["delta_f1_mean"] = float((df["best_of_k_f1"] - df["f1"]).mean())
    return out


def assign_length_bucket(length: int) -> str:
    if 30 <= length <= 80:
        return "30-80"
    if 81 <= length <= 150:
        return "81-150"
    if 151 <= length <= 300:
        return "151-300"
    if 301 <= length <= 400:
        return "301-400"
    return "other"


def render_table_tex(
    *,
    out_path: Path,
    caption: str,
    label: str,
    headers: list[str],
    rows: list[list[str]],
) -> None:
    lines: list[str] = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(r"\begin{tabular}{lrrrrrr}")
    lines.append(r"\toprule")
    lines.append(" & ".join(headers) + r" \\")
    lines.append(r"\midrule")
    for row in rows:
        lines.append(" & ".join(row) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    out_path.write_text("\n".join(lines) + "\n")

def write_cdf_dat(*, out_path: Path, values: pd.Series) -> None:
    vals = values.dropna().astype(float).sort_values().tolist()
    n = len(vals)
    lines = ["f1 cdf", "0 0"]
    if n == 0:
        out_path.write_text("\n".join(lines) + "\n")
        return
    for idx, v in enumerate(vals, start=1):
        lines.append(f"{v:.6f} {idx / n:.6f}")
    out_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate LaTeX tables for FR3D under400 ensemble comparison.")
    parser.add_argument("--out-dir", required=True, help="Output directory (will create generated/*.tex)")
    parser.add_argument("--truth-dir", required=True, help="Truth JSON directory named {target_id}.json")
    parser.add_argument("--pipeline", required=True, help="CSV from ssbench score for RNAnneal-ss pipeline")
    parser.add_argument("--ensemble", required=True, help="CSV from ssbench score for the ensemble method")
    parser.add_argument("--rnastructure", required=True, help="CSV from ssbench score for RNAstructure Fold baseline")
    parser.add_argument("--linearfold", required=True, help="CSV from ssbench score for LinearFold-V baseline")
    parser.add_argument("--eternafold", required=True, help="CSV from ssbench score for EternaFold baseline")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    truth_dir = Path(args.truth_dir)
    gen = out_dir / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    methods = [
        MethodRow("pipeline", "RNAnneal-ss", Path(args.pipeline)),
        MethodRow("ensemble", "Ensemble", Path(args.ensemble)),
        MethodRow("rnastructure", "RNAstructure", Path(args.rnastructure)),
        MethodRow("linearfold", "LinearFold-V", Path(args.linearfold)),
        MethodRow("eternafold", "EternaFold", Path(args.eternafold)),
    ]

    dfs: dict[str, pd.DataFrame] = {}
    for m in methods:
        df = pd.read_csv(m.metrics_csv)
        if "target_id" not in df.columns:
            raise SystemExit(f"missing target_id column: {m.metrics_csv}")
        dfs[m.key] = df

    shared = set.intersection(*(set(df["target_id"].astype(str)) for df in dfs.values()))
    if not shared:
        raise SystemExit("No overlapping target_ids across metrics files.")

    lengths: dict[str, int] = {}
    for target_id in sorted(shared):
        truth_path = truth_dir / f"{target_id}.json"
        if not truth_path.exists():
            continue
        data = json.loads(truth_path.read_text())
        seq = str(data.get("sequence", ""))
        lengths[target_id] = len(seq)

    aligned: dict[str, pd.DataFrame] = {}
    for k, df in dfs.items():
        sub = df[df["target_id"].astype(str).isin(shared)].copy()
        sub = sub.sort_values("target_id")
        sub["truth_len"] = sub["target_id"].astype(str).map(lambda t: lengths.get(t, 0))
        sub["bucket"] = sub["truth_len"].map(assign_length_bucket)
        aligned[k] = sub

    rows = []
    for m in methods:
        s = compute_summary(aligned[m.key])
        rows.append(
            [
                m.name,
                fmt_float(s["mean_top1_f1"]),
                fmt_float(s["median_top1_f1"]),
                fmt_float(s["mean_best100_f1"]),
                fmt_float(s["median_best100_f1"]),
                fmt_float(s["delta_f1_mean"]),
                frac_fmt(s["fail_best100"]),
            ]
        )

    render_table_tex(
        out_path=gen / "overall_table.tex",
        caption=f"Overall performance on FR3D/BGSU under400 non-rRNA benchmark (N={len(shared)}; metrics @1 and best-of-100).",
        label="tab:overall",
        headers=["Method", "Mean F1@1", "Med F1@1", "Mean F1@100", "Med F1@100", r"$\Delta$F1", "Fail@100"],
        rows=rows,
    )

    bucket_order = ["30-80", "81-150", "151-300", "301-400", "other"]
    lines: list[str] = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(rf"\caption{{Mean F1@1 / F1@100 by truth-length bucket (N={len(shared)}).}}")
    lines.append(r"\label{tab:length}")
    lines.append(r"\begin{tabular}{llrrr}")
    lines.append(r"\toprule")
    lines.append(r"Bucket & Method & N & Mean F1@1 & Mean F1@100 \\")
    lines.append(r"\midrule")
    for b in bucket_order:
        n_bucket = int(sum(1 for t in shared if assign_length_bucket(lengths.get(t, 0)) == b))
        if n_bucket == 0:
            continue
        for m in methods:
            df = aligned[m.key]
            sub = df[df["bucket"].astype(str) == b]
            if len(sub) == 0:
                continue
            lines.append(
                " & ".join(
                    [
                        b,
                        m.name,
                        str(int(len(sub))),
                        fmt_float(float(sub["f1"].mean())),
                        fmt_float(float(sub["best_of_k_f1"].mean())),
                    ]
                )
                + r" \\"
            )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    (gen / "length_table.tex").write_text("\n".join(lines) + "\n")

    # CDF data files (top-1 and best-of-100).
    for m in methods:
        df = aligned[m.key]
        write_cdf_dat(out_path=gen / f"cdf_top1_{m.key}.dat", values=df["f1"])
        write_cdf_dat(out_path=gen / f"cdf_best100_{m.key}.dat", values=df["best_of_k_f1"])


if __name__ == "__main__":
    main()
