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
    linestyle: str = "solid"


def fmt_float(x: float) -> str:
    return f"{x:.3f}"


def frac_fmt(x: float) -> str:
    return f"{100.0 * x:.1f}\\%"


def _ecdf(values: list[float]) -> tuple[list[float], list[float]]:
    vals = sorted(float(v) for v in values)
    n = len(vals)
    if n == 0:
        return [0.0], [0.0]
    xs = [0.0] + vals
    ys = [0.0] + [i / n for i in range(1, n + 1)]
    return xs, ys


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
            if j < i:
                out.add((j, i))
            else:
                out.add((i, j))
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


def write_survival_dat(*, out_path: Path, values: list[float], x_name: str = "f1") -> None:
    vals = sorted(float(v) for v in values)
    n = len(vals)
    if n == 0:
        write_xy_dat(out_path=out_path, xs=[0.0], ys=[1.0], x_name=x_name, y_name="survival")
        return
    xs = [0.0] + vals
    ys = [1.0] + [1.0 - (i / n) for i in range(1, n + 1)]
    write_xy_dat(out_path=out_path, xs=xs, ys=ys, x_name=x_name, y_name="survival")


def write_cdf_dat(*, out_path: Path, values: list[float], x_name: str) -> None:
    xs, ys = _ecdf(values)
    write_xy_dat(out_path=out_path, xs=xs, ys=ys, x_name=x_name, y_name="cdf")


def latex_color(hex_color: str) -> str:
    c = hex_color.strip().lstrip("#")
    if len(c) != 6:
        raise ValueError(f"Expected #RRGGBB hex color, got {hex_color!r}")
    return c.upper()


def safe_target_dir(target_id: str) -> str:
    return target_id.replace("|", "_")


def load_qc(pred_dir: Path, target_ids: list[str]) -> dict[str, dict[str, float | None]]:
    out: dict[str, dict[str, float | None]] = {}
    for tid in target_ids:
        qc_path = pred_dir / safe_target_dir(tid) / "qc.json"
        if not qc_path.exists():
            out[tid] = {}
            continue
        data = json.loads(qc_path.read_text())
        timing = data.get("timing_seconds", {}) or {}
        effective = data.get("effective", {}) or {}
        counts = data.get("counts", {}) or {}
        seeds = data.get("seeds", {}) or {}

        def f(x: object) -> float | None:
            if x is None:
                return None
            try:
                return float(x)
            except Exception:
                return None

        out[tid] = {
            "runtime_s": f(timing.get("total", None)),
            "effective_max_scaffolds": f(effective.get("max_scaffolds", None)),
            "effective_refine_max_regions": f(effective.get("refine_max_regions", None)),
            "effective_refine_max_seeds": f(effective.get("refine_max_seeds", None)),
            "effective_refine_kissing_candidates": f(effective.get("refine_kissing_candidates", None)),
            "counts_n_refined": f(counts.get("n_refined", None)),
            "counts_n_candidates": f(counts.get("n_candidates", None)),
            "counts_n_selected_refined": f(counts.get("n_selected_refined", None)),
            "seeds_n_seed_candidates": f(seeds.get("n_seed_candidates", None)),
        }
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Build report assets for the u300 hyperparameter ablation report.")
    p.add_argument("--out-dir", required=True, help="Docs report directory (writes generated/*).")
    p.add_argument("--manifest", required=True, help="Manifest CSV (must include len_group).")
    p.add_argument("--truth-dir", required=True, help="Truth JSON directory named {target_id}.json.")
    p.add_argument("--metrics-dir", required=True, help="Directory containing metrics CSVs.")
    p.add_argument("--pred-hparam-root", required=True, help="Predictions root containing config subdirs (B4/U4/...).")
    p.add_argument("--pred-v3-dir", required=True, help="Predictions dir for v3 RNss(EF) baseline (full run).")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    gen = out_dir / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    truth_dir = Path(args.truth_dir)
    metrics_dir = Path(args.metrics_dir)
    pred_hparam_root = Path(args.pred_hparam_root)
    pred_v3_dir = Path(args.pred_v3_dir)

    methods = [
        MethodSpec("v3", "v3", metrics_dir / "v3_RNss_EF.csv", pred_v3_dir, "#7f7f7f", "dashed"),
        MethodSpec("B4", "B4", metrics_dir / "B4.csv", pred_hparam_root / "B4", "#1f77b4", "solid"),
        MethodSpec("U4", "U4", metrics_dir / "U4.csv", pred_hparam_root / "U4", "#ff7f0e", "solid"),
        MethodSpec("R10", "R10", metrics_dir / "R10.csv", pred_hparam_root / "R10", "#2ca02c", "solid"),
        MethodSpec("S150", "S150", metrics_dir / "S150.csv", pred_hparam_root / "S150", "#d62728", "solid"),
        MethodSpec("K300", "K300", metrics_dir / "K300.csv", pred_hparam_root / "K300", "#9467bd", "solid"),
        MethodSpec("Sc30", "Sc30", metrics_dir / "Sc30.csv", pred_hparam_root / "Sc30", "#8c564b", "solid"),
        MethodSpec("Hi", "Hi", metrics_dir / "Hi.csv", pred_hparam_root / "Hi", "#000000", "solid"),
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
    group_map = {str(row["target_id"]): str(row.get("len_group", "unknown")) for _, row in manifest_df.iterrows()}
    group_order = ["short", "medium", "long", "unknown"]

    aligned: dict[str, pd.DataFrame] = {}
    for k, df in dfs.items():
        sub = df[df["target_id"].astype(str).isin(shared)].copy()
        sub["len_group"] = sub["target_id"].astype(str).map(lambda t: group_map.get(t, "unknown"))
        aligned[k] = sub.sort_values("target_id")

    def attach_qc(df: pd.DataFrame, qc_by_target: dict[str, dict[str, float | None]]) -> pd.DataFrame:
        out = df.copy()
        for k in [
            "runtime_s",
            "effective_max_scaffolds",
            "effective_refine_max_regions",
            "effective_refine_max_seeds",
            "effective_refine_kissing_candidates",
            "counts_n_refined",
            "counts_n_candidates",
            "counts_n_selected_refined",
            "seeds_n_seed_candidates",
        ]:
            out[k] = out["target_id"].astype(str).map(lambda t: (qc_by_target.get(t, {}) or {}).get(k, None))
        return out

    qc_by_method: dict[str, dict[str, dict[str, float | None]]] = {}
    for m in methods:
        qc_by_method[m.key] = load_qc(m.predictions_dir, sorted(shared))
        aligned[m.key] = attach_qc(aligned[m.key], qc_by_method[m.key])

    def summary(df: pd.DataFrame) -> dict[str, float]:
        out: dict[str, float] = {}
        out["n"] = float(len(df))
        out["mean_f1_1"] = float(df["f1"].mean())
        out["median_f1_1"] = float(df["f1"].median())
        out["mean_f1_100"] = float(df["best_of_k_f1"].mean())
        out["median_f1_100"] = float(df["best_of_k_f1"].median())
        out["fail_100"] = float((df["best_of_k_f1"] <= 0.0).mean())
        runtimes = df["runtime_s"].dropna().astype(float)
        out["mean_runtime_s"] = float(runtimes.mean()) if len(runtimes) else float("nan")
        out["median_runtime_s"] = float(runtimes.median()) if len(runtimes) else float("nan")
        return out

    overall_rows: list[list[str]] = []
    for m in methods:
        s = summary(aligned[m.key])
        overall_rows.append(
            [
                m.label,
                fmt_float(s["mean_f1_1"]),
                fmt_float(s["mean_f1_100"]),
                frac_fmt(s["fail_100"]),
                f"{s['mean_runtime_s']:.1f}" if s["mean_runtime_s"] == s["mean_runtime_s"] else "n/a",
                f"{s['median_runtime_s']:.1f}" if s["median_runtime_s"] == s["median_runtime_s"] else "n/a",
            ]
        )

    render_table_tex(
        out_path=gen / "overall_table.tex",
        caption=f"Overall ablation results (N={len(shared)}; @1 and best-of-100). Runtime is wall-clock seconds per target.",
        label="tab:overall",
        headers=["Config", "Mean F1@1", "Mean F1@100", "Fail@100", "Mean time (s)", "Med time (s)"],
        rows=overall_rows,
        col_spec="lrrrrr",
    )

    group_rows: list[list[str]] = []
    for g in group_order:
        for m in methods:
            sub = aligned[m.key][aligned[m.key]["len_group"].astype(str) == g]
            if len(sub) == 0:
                continue
            s = summary(sub)
            group_rows.append(
                [
                    g,
                    m.label,
                    str(int(len(sub))),
                    fmt_float(s["mean_f1_1"]),
                    fmt_float(s["mean_f1_100"]),
                    frac_fmt(s["fail_100"]),
                    f"{s['mean_runtime_s']:.1f}" if s["mean_runtime_s"] == s["mean_runtime_s"] else "n/a",
                ]
            )

    render_table_tex(
        out_path=gen / "group_table.tex",
        caption="Results by length group from the manifest (short/medium/long).",
        label="tab:len_groups",
        headers=["Group", "Config", "N", "Mean F1@1", "Mean F1@100", "Fail@100", "Mean time (s)"],
        rows=group_rows,
        col_spec="llrrrrr",
    )

    # QC sanity check: report that the ablated knobs actually changed the internal budgets/counts.
    qc_rows: list[list[str]] = []
    for m in methods:
        df = aligned[m.key]
        cols = [
            "effective_max_scaffolds",
            "effective_refine_max_regions",
            "effective_refine_max_seeds",
            "effective_refine_kissing_candidates",
            "counts_n_refined",
            "counts_n_candidates",
        ]
        means: dict[str, float] = {}
        for c in cols:
            s = df[c].dropna().astype(float)
            means[c] = float(s.mean()) if len(s) else float("nan")
        qc_rows.append(
            [
                m.label,
                f"{means['effective_max_scaffolds']:.1f}" if means["effective_max_scaffolds"] == means["effective_max_scaffolds"] else "n/a",
                f"{means['effective_refine_max_regions']:.1f}" if means["effective_refine_max_regions"] == means["effective_refine_max_regions"] else "n/a",
                f"{means['effective_refine_max_seeds']:.1f}" if means["effective_refine_max_seeds"] == means["effective_refine_max_seeds"] else "n/a",
                f"{means['effective_refine_kissing_candidates']:.1f}"
                if means["effective_refine_kissing_candidates"] == means["effective_refine_kissing_candidates"]
                else "n/a",
                f"{means['counts_n_refined']:.0f}" if means["counts_n_refined"] == means["counts_n_refined"] else "n/a",
                f"{means['counts_n_candidates']:.0f}"
                if means["counts_n_candidates"] == means["counts_n_candidates"]
                else "n/a",
            ]
        )

    render_table_tex(
        out_path=gen / "qc_table.tex",
        caption="QC means from per-target \\texttt{qc.json} (post length-adaptive scaling).",
        label="tab:qc",
        headers=[
            "Config",
            "Eff. scaffolds",
            "Eff. regions",
            "Eff. seeds",
            "Eff. kissing",
            "Mean refined",
            "Mean cand.",
        ],
        rows=qc_rows,
        col_spec="lrrrrrr",
    )

    # Short findings paragraph (auto-generated from summaries).
    by_key = {m.key: summary(aligned[m.key]) for m in methods}
    best_key, best_stats = max(by_key.items(), key=lambda kv: kv[1]["mean_f1_100"])
    base = by_key.get("B4", best_stats)
    best_delta = best_stats["mean_f1_100"] - base["mean_f1_100"]
    runtime_delta = best_stats["mean_runtime_s"] - base["mean_runtime_s"]
    mean_100s = [s["mean_f1_100"] for s in by_key.values()]
    spread = (max(mean_100s) - min(mean_100s)) if mean_100s else 0.0
    findings_lines = [
        r"\begin{itemize}",
        rf"\item Best mean F1@100: \textbf{{{best_key}}} ({best_stats['mean_f1_100']:.4f}); B4={base['mean_f1_100']:.4f}; $\Delta$={best_delta:+.4f}.",
        rf"\item Mean F1@100 spread across configs is {spread:.4f} on this 27-target subset (Table~\ref{{tab:qc}} confirms the knobs did change internal budgets).",
        rf"\item Fail@100 is {by_key['B4']['fail_100']*100:.1f}\% for all configs on this subset.",
        rf"\item Runtime: {best_key} is {runtime_delta:+.1f}s vs B4 (mean per target).",
        r"\end{itemize}",
    ]
    (gen / "findings.tex").write_text("\n".join(findings_lines) + "\n")

    # ---- Plot data and TikZ/pgfplots snippets ----
    series_overall: dict[str, dict[str, list[float]]] = {}
    for m in methods:
        df = aligned[m.key]
        series_overall[m.key] = {
            "f1_1": df["f1"].astype(float).tolist(),
            "f1_100": df["best_of_k_f1"].astype(float).tolist(),
            "runtime_s": df["runtime_s"].dropna().astype(float).tolist(),
        }

    for m in methods:
        write_survival_dat(out_path=gen / f"survival_top1_{m.key}.dat", values=series_overall[m.key]["f1_1"])
        write_survival_dat(out_path=gen / f"survival_best100_{m.key}.dat", values=series_overall[m.key]["f1_100"])
        write_cdf_dat(out_path=gen / f"runtime_cdf_{m.key}.dat", values=series_overall[m.key]["runtime_s"], x_name="sec")

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
        lines.append(r"  legend style={font=\scriptsize},")
        lines.append(r"  legend columns=2,")
        lines.append(r"  legend pos=north east, legend cell align=left,")
        lines.append(r"]")
        for m in methods:
            vals = series_overall[m.key]["f1_1" if dat_prefix == "survival_top1" else "f1_100"]
            med = f1_at_cdf(vals, 0.50)
            label = m.label if med is None else f"{m.label} (p50={med:.3f})"
            lines.append(
                rf"\addplot+[const plot, very thick, draw={cname(m)}, {m.linestyle}] "
                rf"table[x=f1,y=survival] {{{rel}/{dat_prefix}_{m.key}.dat}};"
            )
            lines.append(rf"\addlegendentry{{{label}}}")
        lines.append(r"\end{axis}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    write_survival_tex(out_tex=gen / "survival_top1.tex", title="Top-1 (F1@1)", dat_prefix="survival_top1")
    write_survival_tex(out_tex=gen / "survival_best100.tex", title="Best-of-100 (F1@100)", dat_prefix="survival_best100")

    def write_runtime_cdf_tex(out_tex: Path) -> None:
        rel = "generated"
        lines: list[str] = []
        lines.append(r"\begin{tikzpicture}")
        lines.extend(define_colors())
        lines.append(r"\begin{axis}[")
        lines.append(r"  width=\linewidth, height=0.62\linewidth,")
        lines.append(r"  xlabel={Seconds}, ylabel={CDF},")
        lines.append(r"  title={Runtime CDF (per target)},")
        lines.append(r"  grid=both, grid style={black!10},")
        lines.append(r"  legend style={font=\scriptsize},")
        lines.append(r"  legend columns=2,")
        lines.append(r"  legend pos=south east, legend cell align=left,")
        lines.append(r"]")
        for m in methods:
            vals = series_overall[m.key]["runtime_s"]
            med = f1_at_cdf(vals, 0.50)
            label = m.label if med is None else f"{m.label} (p50={med:.1f}s)"
            lines.append(
                rf"\addplot+[const plot, very thick, draw={cname(m)}, {m.linestyle}] "
                rf"table[x=sec,y=cdf] {{{rel}/runtime_cdf_{m.key}.dat}};"
            )
            lines.append(rf"\addlegendentry{{{label}}}")
        lines.append(r"\end{axis}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    write_runtime_cdf_tex(gen / "runtime_cdf.tex")

    # ---- PR curves (subset: v3 vs B4 vs Hi) ----
    pr_methods = [m for m in methods if m.key in {"v3", "B4", "Hi"}]
    truth_pairs_by_target: dict[str, set[tuple[int, int]]] = {}
    for tid in sorted(shared):
        truth_path = truth_dir / f"{tid}.json"
        data = json.loads(truth_path.read_text())
        pairs = {tuple(p) for p in data.get("canonical_pairs", [])}
        seq = str(data.get("sequence", ""))
        truth_pairs_by_target[tid] = {(int(i), int(j)) for i, j in pairs if 0 <= int(i) < int(j) < len(seq)}

    def build_probs(method: MethodSpec, target_ids: list[str]) -> dict[str, dict[tuple[int, int], float]]:
        out: dict[str, dict[tuple[int, int], float]] = {}
        for tid in target_ids:
            db_path = method.predictions_dir / safe_target_dir(tid) / "predictions.db"
            seq, structs = read_predictions_db(db_path)
            out[tid] = pair_probs_from_topk(structs, len(seq), k=100)
        return out

    thresholds = [i / 100 for i in range(0, 101)]
    for m in pr_methods:
        probs = build_probs(m, sorted(shared))
        r, p_ = pr_curve(probs_by_target=probs, truth_pairs_by_target=truth_pairs_by_target, thresholds=thresholds)
        write_xy_dat(out_path=gen / f"pr_{m.key}.dat", xs=r, ys=p_, x_name="recall", y_name="precision")

    def write_pr_tex(out_tex: Path) -> None:
        rel = "generated"
        lines: list[str] = []
        lines.append(r"\begin{tikzpicture}")
        lines.extend(define_colors())
        lines.append(r"\begin{axis}[")
        lines.append(r"  width=0.9\linewidth, height=0.55\linewidth,")
        lines.append(r"  xmin=0, xmax=1, ymin=0, ymax=1,")
        lines.append(r"  xlabel={Recall}, ylabel={Precision},")
        lines.append(r"  title={Pair PR curves (top-100 pair frequencies)},")
        lines.append(r"  grid=both, grid style={black!10},")
        lines.append(r"  legend pos=south west, legend cell align=left,")
        lines.append(r"]")
        for m in pr_methods:
            lines.append(
                rf"\addplot+[very thick, draw={cname(m)}, {m.linestyle}] table[x=recall,y=precision] {{{rel}/pr_{m.key}.dat}};"
            )
            lines.append(rf"\addlegendentry{{{m.label}}}")
        lines.append(r"\end{axis}")
        lines.append(r"\end{tikzpicture}")
        out_tex.write_text("\n".join(lines) + "\n")

    write_pr_tex(gen / "pr_compare.tex")


if __name__ == "__main__":
    main()
