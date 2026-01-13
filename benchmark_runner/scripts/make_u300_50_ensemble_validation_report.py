#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class MethodSpec:
    key: str
    label: str
    metrics_csv: Path
    predictions_dir: Path
    color: str


def fmt_float(x: float) -> str:
    return f"{x:.3f}"


def frac_fmt(x: float) -> str:
    return f"{100.0 * x:.1f}\\%"


def f1_at_cdf(values: list[float], cdf: float) -> float | None:
    if not values:
        return None
    if cdf <= 0:
        return min(values)
    if cdf >= 1:
        return max(values)
    vals = sorted(float(v) for v in values)
    n = len(vals)
    import math

    k = int(math.ceil(cdf * n))  # 1..n
    k = max(1, min(n, k))
    return vals[k - 1]


def read_predictions_db(path: Path) -> tuple[str, list[str]]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        return "", []
    seq = lines[0]
    structs = [ln for ln in lines[1:] if len(ln) == len(seq)]
    return seq, structs


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
            out.add((j, i) if j < i else (i, j))
    return out


def pair_probs_from_topk(structs: list[str], length: int, k: int = 100) -> dict[tuple[int, int], float]:
    n = min(int(k), len(structs))
    if n <= 0:
        return {}
    counts: dict[tuple[int, int], int] = {}
    for db in structs[:n]:
        for i, j in pairs_from_dotbracket(db):
            if 0 <= i < j < length:
                counts[(i, j)] = counts.get((i, j), 0) + 1
    denom = float(n)
    return {p: c / denom for p, c in counts.items()}


def pr_curve(
    *,
    probs_by_target: dict[str, dict[tuple[int, int], float]],
    truth_pairs_by_target: dict[str, set[tuple[int, int]]],
    thresholds: list[float],
) -> tuple[list[float], list[float]]:
    recalls: list[float] = []
    precisions: list[float] = []
    targets = sorted(probs_by_target.keys())
    for t in thresholds:
        ps: list[float] = []
        rs: list[float] = []
        for tid in targets:
            probs = probs_by_target[tid]
            truth = truth_pairs_by_target.get(tid, set())
            pred = {p for p, v in probs.items() if v >= t}
            tp = len(pred & truth)
            fp = len(pred - truth)
            fn = len(truth - pred)
            pval = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            rval = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            ps.append(pval)
            rs.append(rval)
        precisions.append(float(sum(ps) / len(ps)) if ps else 0.0)
        recalls.append(float(sum(rs) / len(rs)) if rs else 0.0)
    return recalls, precisions


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


def write_xy_dat(*, out_path: Path, xs: list[float], ys: list[float], x_name: str, y_name: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{x_name} {y_name}"]
    for x, y in zip(xs, ys, strict=True):
        lines.append(f"{x:.6f} {y:.6f}")
    out_path.write_text("\n".join(lines) + "\n")


def write_survival_dat(*, out_path: Path, values: list[float]) -> None:
    vals = sorted(float(v) for v in values)
    n = len(vals)
    if n == 0:
        write_xy_dat(out_path=out_path, xs=[0.0], ys=[1.0], x_name="f1", y_name="survival")
        return
    xs = [0.0] + vals
    ys = [1.0] + [1.0 - (i / n) for i in range(1, n + 1)]
    write_xy_dat(out_path=out_path, xs=xs, ys=ys, x_name="f1", y_name="survival")


def latex_color(hex_color: str) -> str:
    c = hex_color.strip().lstrip("#")
    if len(c) != 6:
        raise ValueError(f"Expected #RRGGBB hex color, got {hex_color!r}")
    return c.upper()


def main() -> None:
    p = argparse.ArgumentParser(description="Build report assets for u300 representative50 ensemble validation.")
    p.add_argument("--out-dir", required=True, help="Docs report directory (writes generated/*).")
    p.add_argument("--manifest", required=True, help="Manifest CSV (must include bucket).")
    p.add_argument("--truth-dir", required=True, help="Truth JSON directory named {target_id}.json.")
    p.add_argument("--v3-metrics-dir", required=True, help="u300_50_scaffold_validation_v3 metrics dir.")
    p.add_argument("--v3-pred-root", required=True, help="u300_50_scaffold_validation_v3 predictions root.")
    p.add_argument("--ens-metrics-dir", required=True, help="u300_50_ensemble_validation_v1 metrics dir.")
    p.add_argument("--ens-pred-root", required=True, help="u300_50_ensemble_validation_v1 predictions root.")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    gen = out_dir / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    truth_dir = Path(args.truth_dir)
    v3_pred_root = Path(args.v3_pred_root)
    v3_metrics_dir = Path(args.v3_metrics_dir)
    ens_pred_root = Path(args.ens_pred_root)
    ens_metrics_dir = Path(args.ens_metrics_dir)

    methods = [
        MethodSpec(
            "ens_tuned",
            "Ens(RNss)",
            ens_metrics_dir / "ensemble_merge_rnss_tuned.csv",
            ens_pred_root / "ensemble_merge_rnss_tuned",
            "#000000",
        ),
        MethodSpec(
            "rn_ef",
            "RNss(EF)",
            v3_metrics_dir / "rnanneal_ss_ef_scaff.csv",
            v3_pred_root / "rnanneal_ss_ef_scaff",
            "#D62728",
        ),
        MethodSpec(
            "rn_lf",
            "RNss(LF)",
            v3_metrics_dir / "rnanneal_ss_lf_scaff.csv",
            v3_pred_root / "rnanneal_ss_lf_scaff",
            "#9467BD",
        ),
        MethodSpec(
            "rn_rs",
            "RNss(RS)",
            v3_metrics_dir / "rnanneal_ss_rs_scaff.csv",
            v3_pred_root / "rnanneal_ss_rs_scaff",
            "#1F77B4",
        ),
        MethodSpec(
            "ens3_scaff",
            "RNss(Ens3)",
            ens_metrics_dir / "rnanneal_ss_ens_scaff.csv",
            ens_pred_root / "rnanneal_ss_ens_scaff",
            "#7F7F7F",
        ),
        MethodSpec("lf", "LF-V", v3_metrics_dir / "linearfold.csv", v3_pred_root / "linearfold", "#FF7F0E"),
        MethodSpec("ef", "EFold", v3_metrics_dir / "eternafold.csv", v3_pred_root / "eternafold", "#2CA02C"),
        MethodSpec("rs", "RNAstr", v3_metrics_dir / "rnastructure.csv", v3_pred_root / "rnastructure", "#BCBD22"),
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

    manifest_df = pd.read_csv(Path(args.manifest))
    bucket_map = {str(row["target_id"]): str(row.get("bucket", "unknown")) for _, row in manifest_df.iterrows()}

    aligned: dict[str, pd.DataFrame] = {}
    for k, df in dfs.items():
        sub = df[df["target_id"].astype(str).isin(shared)].copy()
        sub["bucket"] = sub["target_id"].astype(str).map(lambda t: bucket_map.get(t, "unknown"))
        aligned[k] = sub.sort_values("target_id")

    def summary(df: pd.DataFrame) -> dict[str, float]:
        out: dict[str, float] = {}
        out["n"] = float(len(df))
        out["mean_f1_1"] = float(df["f1"].mean())
        out["median_f1_1"] = float(df["f1"].median())
        out["mean_f1_100"] = float(df["best_of_k_f1"].mean())
        out["median_f1_100"] = float(df["best_of_k_f1"].median())
        out["fail_100"] = float((df["best_of_k_f1"] <= 0.0).mean())
        out["delta_mean"] = float((df["best_of_k_f1"] - df["f1"]).mean())
        return out

    overall_rows: list[list[str]] = []
    for m in methods:
        s = summary(aligned[m.key])
        overall_rows.append(
            [
                m.label,
                fmt_float(s["mean_f1_1"]),
                fmt_float(s["median_f1_1"]),
                fmt_float(s["mean_f1_100"]),
                fmt_float(s["median_f1_100"]),
                fmt_float(s["delta_mean"]),
                frac_fmt(s["fail_100"]),
            ]
        )

    render_table_tex(
        out_path=gen / "overall_table.tex",
        caption=f"Overall performance on representative50 (N={len(shared)}; metrics @1 and best-of-100).",
        label="tab:overall",
        headers=["Method", "Mean F1@1", "Med F1@1", "Mean F1@100", "Med F1@100", r"$\Delta$F1", "Fail@100"],
        rows=overall_rows,
        col_spec="lrrrrrr",
    )

    bucket_order = ["30-79", "80-129", "130-179", "180-229", "230-279", "280-300", "unknown"]
    bucket_rows: list[list[str]] = []
    for bucket in bucket_order:
        for m in methods:
            sub = aligned[m.key][aligned[m.key]["bucket"].astype(str) == bucket]
            if len(sub) == 0:
                continue
            s = summary(sub)
            bucket_rows.append(
                [
                    bucket,
                    m.label,
                    str(int(len(sub))),
                    fmt_float(s["mean_f1_1"]),
                    fmt_float(s["mean_f1_100"]),
                    frac_fmt(s["fail_100"]),
                ]
            )

    render_table_tex(
        out_path=gen / "bucket_table.tex",
        caption="Performance by truth-length bucket (50-nt buckets from the manifest).",
        label="tab:buckets",
        headers=["Bucket", "Method", "N", "Mean F1@1", "Mean F1@100", "Fail@100"],
        rows=bucket_rows,
        col_spec="llrrrr",
    )

    # ---- Plot data + TikZ/pgfplots snippets ----
    series_overall: dict[str, dict[str, list[float]]] = {}
    for m in methods:
        df = aligned[m.key]
        series_overall[m.key] = {
            "f1_1": df["f1"].astype(float).tolist(),
            "f1_100": df["best_of_k_f1"].astype(float).tolist(),
        }

    for m in methods:
        write_survival_dat(out_path=gen / f"survival_top1_{m.key}.dat", values=series_overall[m.key]["f1_1"])
        write_survival_dat(out_path=gen / f"survival_best100_{m.key}.dat", values=series_overall[m.key]["f1_100"])

    def cname(m: MethodSpec) -> str:
        return f"c_{m.key}"

    def define_colors() -> list[str]:
        return [rf"\definecolor{{{cname(m)}}}{{HTML}}{{{latex_color(m.color)}}}" for m in methods]

    def write_survival_tex(*, out_tex: Path, title: str, dat_prefix: str) -> None:
        rel = "generated"
        lines: list[str] = []
        lines.append(r"\begin{tikzpicture}")
        lines.extend(define_colors())
        lines.append(r"\begin{axis}[")
        lines.append(r"  width=\linewidth, height=0.62\linewidth,")
        lines.append(r"  xmin=0, xmax=1, ymin=0, ymax=1,")
        lines.append(r"  xlabel={F1}, ylabel={1 - CDF},")
        lines.append(rf"  title={{{title}}},")
        lines.append(r"  grid=both, grid style={black!10},")
        lines.append(r"  legend pos=north east, legend cell align=left,")
        lines.append(r"]")
        for m in methods:
            vals = series_overall[m.key]["f1_1" if dat_prefix == "survival_top1" else "f1_100"]
            med = f1_at_cdf(vals, 0.50)
            label = m.label if med is None else f"{m.label} (p50={med:.3f})"
            lines.append(
                rf"\addplot+[const plot, very thick, draw={cname(m)}] table[x=f1,y=survival] {{{rel}/{dat_prefix}_{m.key}.dat}};"
            )
            lines.append(rf"\addlegendentry{{{label}}}")
        lines.append(r"\end{axis}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    write_survival_tex(out_tex=gen / "survival_top1.tex", title="Top-1 (F1@1)", dat_prefix="survival_top1")
    write_survival_tex(out_tex=gen / "survival_best100.tex", title="Best-of-100 (F1@100)", dat_prefix="survival_best100")

    # Precision/recall curves from top-100 pair frequencies (python-computed).
    truth_pairs_by_target: dict[str, set[tuple[int, int]]] = {}
    for tid in sorted(shared):
        truth_path = truth_dir / f"{tid}.json"
        data = json.loads(truth_path.read_text())
        pairs = {tuple(p) for p in data.get("canonical_pairs", [])}
        seq = str(data.get("sequence", ""))
        L = len(seq)
        truth_pairs_by_target[tid] = {(int(i), int(j)) for i, j in pairs if 0 <= int(i) < int(j) < L}

    def build_probs(method: MethodSpec, target_ids: list[str]) -> dict[str, dict[tuple[int, int], float]]:
        out: dict[str, dict[tuple[int, int], float]] = {}
        for tid in target_ids:
            safe = tid.replace("|", "_")
            db_path = method.predictions_dir / safe / "predictions.db"
            seq, structs = read_predictions_db(db_path)
            out[tid] = pair_probs_from_topk(structs, len(seq), k=100)
        return out

    thresholds = [i / 100 for i in range(0, 101)]
    probs_overall: dict[str, dict[str, dict[tuple[int, int], float]]] = {}
    for m in methods:
        probs_overall[m.key] = build_probs(m, sorted(shared))
        r, p_ = pr_curve(
            probs_by_target=probs_overall[m.key],
            truth_pairs_by_target=truth_pairs_by_target,
            thresholds=thresholds,
        )
        write_xy_dat(out_path=gen / f"pr_overall_{m.key}.dat", xs=r, ys=p_, x_name="recall", y_name="precision")

    buckets = sorted({bucket_map.get(t, "unknown") for t in shared})
    for bucket in buckets:
        bucket_ids = [t for t in sorted(shared) if bucket_map.get(t, "unknown") == bucket]
        for m in methods:
            probs = build_probs(m, bucket_ids)
            r, p_ = pr_curve(
                probs_by_target=probs,
                truth_pairs_by_target=truth_pairs_by_target,
                thresholds=thresholds,
            )
            safe_bucket = bucket.replace("-", "_")
            write_xy_dat(
                out_path=gen / f"pr_bucket_{safe_bucket}_{m.key}.dat",
                xs=r,
                ys=p_,
                x_name="recall",
                y_name="precision",
            )

    def write_pr_tex_overall(out_tex: Path) -> None:
        rel = "generated"
        lines: list[str] = []
        lines.append(r"\begin{tikzpicture}")
        lines.extend(define_colors())
        lines.append(r"\begin{axis}[")
        lines.append(r"  width=0.95\linewidth, height=0.58\linewidth,")
        lines.append(r"  xmin=0, xmax=1, ymin=0, ymax=1,")
        lines.append(r"  xlabel={Recall}, ylabel={Precision},")
        lines.append(r"  title={Pair PR curves from top-100 pair frequencies},")
        lines.append(r"  grid=both, grid style={black!10},")
        lines.append(r"  legend pos=south west, legend cell align=left,")
        lines.append(r"]")
        for m in methods:
            lines.append(
                rf"\addplot+[very thick, draw={cname(m)}] table[x=recall,y=precision] {{{rel}/pr_overall_{m.key}.dat}};"
            )
            lines.append(rf"\addlegendentry{{{m.label}}}")
        lines.append(r"\end{axis}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    def write_pr_tex_by_bucket(out_tex: Path) -> None:
        rel = "generated"
        safe_buckets = [b.replace("-", "_") for b in buckets]
        lines: list[str] = []
        lines.append(r"\begin{tikzpicture}")
        lines.extend(define_colors())
        lines.append(r"\begin{groupplot}[")
        lines.append(r"  group style={group size=3 by 2, horizontal sep=1.2em, vertical sep=1.2em},")
        lines.append(r"  width=0.33\linewidth, height=0.33\linewidth,")
        lines.append(r"  xmin=0, xmax=1, ymin=0, ymax=1,")
        lines.append(r"  xlabel={Recall}, ylabel={Precision},")
        lines.append(r"  grid=both, grid style={black!10},")
        lines.append(r"]")
        for bucket, safe_bucket in zip(buckets, safe_buckets, strict=True):
            lines.append(rf"\nextgroupplot[title={{{bucket}}}]")
            for m in methods:
                lines.append(
                    rf"\addplot+[very thick, draw={cname(m)}] table[x=recall,y=precision] {{{rel}/pr_bucket_{safe_bucket}_{m.key}.dat}};"
                )
        lines.append(r"\end{groupplot}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    write_pr_tex_overall(gen / "pr_overall.tex")
    write_pr_tex_by_bucket(gen / "pr_by_bucket.tex")

    def five_number(values: list[float]) -> tuple[float, float, float, float, float]:
        if not values:
            return 0.0, 0.0, 0.0, 0.0, 0.0
        vs = sorted(float(v) for v in values)
        import numpy as np

        return (
            float(vs[0]),
            float(np.quantile(vs, 0.25)),
            float(np.quantile(vs, 0.50)),
            float(np.quantile(vs, 0.75)),
            float(vs[-1]),
        )

    def write_boxplot_group_tex(*, out_tex: Path, metric: str) -> None:
        lines: list[str] = []
        lines.append(r"\begin{tikzpicture}")
        lines.extend(define_colors())
        lines.append(r"\begin{groupplot}[")
        lines.append(r"  group style={group size=3 by 2, horizontal sep=1.2em, vertical sep=1.2em},")
        lines.append(r"  width=0.33\linewidth, height=0.33\linewidth,")
        lines.append(r"  xmin=0.5, xmax=8.5,")
        lines.append(r"  ymin=0, ymax=1,")
        lines.append(r"  ymajorgrids=true, grid style={black!10},")
        lines.append(r"  xtick={1,2,3,4,5,6,7,8},")
        lines.append(r"  xticklabels={Ens,RNss(EF),RNss(LF),RNss(RS),RNss(Ens3),LF-V,EFold,RNAstr},")
        lines.append(r"  x tick label style={rotate=30,anchor=east,font=\scriptsize},")
        lines.append(r"  title style={font=\small},")
        lines.append(r"]")
        for bucket in buckets:
            lines.append(rf"\nextgroupplot[title={{{bucket}}}]")
            for idx, m in enumerate(methods, start=1):
                sub = aligned[m.key][aligned[m.key]["bucket"].astype(str) == bucket]
                vals = (
                    sub["f1"].astype(float).tolist()
                    if metric == "f1_1"
                    else sub["best_of_k_f1"].astype(float).tolist()
                )
                mn, q1, med, q3, mx = five_number(vals)
                lines.append(
                    rf"\addplot+ [boxplot prepared={{draw position={idx}, lower whisker={mn:.4f}, lower quartile={q1:.4f}, median={med:.4f}, upper quartile={q3:.4f}, upper whisker={mx:.4f}}},"
                    rf" draw={cname(m)}, fill={cname(m)}, fill opacity=0.20, very thick] coordinates {{}};"
                )
        lines.append(r"\end{groupplot}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    write_boxplot_group_tex(out_tex=gen / "dist_f1_top1_by_bucket.tex", metric="f1_1")
    write_boxplot_group_tex(out_tex=gen / "dist_f1_best100_by_bucket.tex", metric="f1_100")


if __name__ == "__main__":
    main()
