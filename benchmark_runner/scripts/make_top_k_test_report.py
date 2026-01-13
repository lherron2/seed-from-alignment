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


def frac_fmt(x: float) -> str:
    return f"{100.0 * x:.1f}\\%"


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


def _ecdf(values: list[float]) -> tuple[list[float], list[float]]:
    vals = sorted(float(v) for v in values)
    n = len(vals)
    if n == 0:
        return [0.0], [0.0]
    xs = [0.0] + vals
    ys = [0.0] + [i / n for i in range(1, n + 1)]
    return xs, ys


def write_cdf_dat(*, out_path: Path, values: list[float]) -> None:
    xs, ys = _ecdf(values)
    write_xy_dat(out_path=out_path, xs=xs, ys=ys, x_name="f1", y_name="cdf")


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


def latex_color(hex_color: str) -> str:
    c = hex_color.strip().lstrip("#")
    if len(c) != 6:
        raise ValueError(f"Expected #RRGGBB hex color, got {hex_color!r}")
    return c.upper()


@dataclass(frozen=True)
class SeriesSpec:
    k: int
    label: str
    color: str


def main() -> None:
    p = argparse.ArgumentParser(description="Build report assets for the top-K output-size test.")
    p.add_argument("--out-dir", required=True, help="Docs report directory (writes generated/*).")
    p.add_argument("--manifest", required=True, help="Manifest CSV (must include truth_len and len_bucket_50).")
    p.add_argument("--truth-dir", required=True, help="Truth JSON directory named {target_id}.json.")
    p.add_argument("--pred-dir", required=True, help="Predictions directory containing {safe_id}/predictions.db.")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    gen = out_dir / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(Path(args.manifest))
    truth_dir = Path(args.truth_dir)
    pred_dir = Path(args.pred_dir)

    # Keep in ascending order so curves make sense.
    ks = [1, 50, 100, 200, 500]

    # Load metrics helpers from ssbench (already in the venv).
    from ssbench.metrics.pair_metrics import compute_pair_metrics

    rows: list[dict[str, object]] = []
    for _, row in manifest.iterrows():
        target_id = str(row["target_id"])
        safe_id = target_id.replace("|", "_")
        truth_path = truth_dir / f"{target_id}.json"
        pred_path = pred_dir / safe_id / "predictions.db"
        if not truth_path.exists() or not pred_path.exists():
            continue

        truth = json.loads(truth_path.read_text())
        seq = str(truth.get("sequence", ""))
        L = len(seq)
        ref_pairs = {tuple(p) for p in truth.get("canonical_pairs", [])}
        ref_pairs = {(int(i), int(j)) for i, j in ref_pairs if 0 <= int(i) < int(j) < L}

        _pred_seq, structs = read_predictions_db(pred_path)
        structs = [s for s in structs if len(s) == L]
        pred_pairs_list = [pairs_from_dotbracket(s) for s in structs]
        if not pred_pairs_list:
            pred_pairs_list = [set()]

        for k in ks:
            n = min(int(k), len(pred_pairs_list))
            best_idx = 0
            best_f1 = -1.0
            for i in range(n):
                f1 = compute_pair_metrics(pred_pairs_list[i], ref_pairs, L).f1
                if f1 > best_f1:
                    best_f1 = f1
                    best_idx = i
            best = compute_pair_metrics(pred_pairs_list[best_idx], ref_pairs, L)
            rows.append(
                {
                    "target_id": target_id,
                    "len_bucket_50": str(row.get("len_bucket_50", "unknown")),
                    "truth_len": int(row.get("truth_len", L)),
                    "k": int(k),
                    "effective_k": int(n),
                    "precision": float(best.precision),
                    "recall": float(best.recall),
                    "f1": float(best.f1),
                }
            )

    if not rows:
        raise SystemExit("No rows scored. Check manifest/truth/predictions.")

    df = pd.DataFrame(rows)
    out_csv = gen / "topk_metrics.csv"
    df.to_csv(out_csv, index=False)

    # Two-way length stratification for "do longer RNAs benefit more?"
    def len_group(L: int) -> str:
        return "long(151-300)" if int(L) >= 151 else "short(<=150)"

    df["len_group"] = df["truth_len"].astype(int).map(len_group)

    # Summary table.
    overall_rows: list[list[str]] = []
    for k in ks:
        sub = df[df["k"] == int(k)]
        if len(sub) == 0:
            continue
        overall_rows.append(
            [
                f"@{k}",
                fmt_float(float(sub["f1"].mean())),
                fmt_float(float(sub["f1"].median())),
                frac_fmt(float((sub["f1"] <= 0.0).mean())),
                fmt_float(float(sub["precision"].mean())),
                fmt_float(float(sub["recall"].mean())),
                fmt_float(float(sub["effective_k"].mean())),
            ]
        )

    render_table_tex(
        out_path=gen / "overall_table.tex",
        caption=f"Top-K oracle metrics vs output size (N={df['target_id'].nunique()}; higher is better).",
        label="tab:overall",
        headers=["K", "Mean F1", "Med F1", "Fail", "Mean P", "Mean R", "Mean eff.K"],
        rows=overall_rows,
        col_spec="lrrrrrr",
    )

    # Benefit vs K by length group (delta vs @1).
    base = df[df["k"] == 1][["target_id", "f1"]].rename(columns={"f1": "f1_at1"})
    merged = df.merge(base, on="target_id", how="left")
    merged["delta_f1"] = merged["f1"].astype(float) - merged["f1_at1"].astype(float)

    delta_rows: list[list[str]] = []
    for g in ["short(<=150)", "long(151-300)"]:
        for k in [50, 100, 200, 500]:
            sub = merged[(merged["len_group"] == g) & (merged["k"] == k)]
            if len(sub) == 0:
                continue
            delta_rows.append(
                [
                    g,
                    f"@{k}",
                    str(int(len(sub))),
                    fmt_float(float(sub["delta_f1"].mean())),
                    fmt_float(float(sub["delta_f1"].median())),
                ]
            )

    render_table_tex(
        out_path=gen / "delta_by_len_table.tex",
        caption=r"Mean/median $\Delta$F1 (best-of-$K$ minus best-of-1), stratified by RNA length.",
        label="tab:delta_len",
        headers=["Group", "K", "N", "Mean $\\Delta$F1", "Med $\\Delta$F1"],
        rows=delta_rows,
        col_spec="llrrr",
    )

    # Brief findings (written as LaTeX) to answer: do longer RNAs benefit more?
    lines: list[str] = [r"\begin{itemize}"]
    for k in [50, 100, 200, 500]:
        s = float(merged[(merged["len_group"] == "short(<=150)") & (merged["k"] == k)]["delta_f1"].mean())
        l = float(merged[(merged["len_group"] == "long(151-300)") & (merged["k"] == k)]["delta_f1"].mean())
        lines.append(rf"\item Mean $\Delta$F1 @{k}: short={s:.3f}, long={l:.3f} (long-short={l-s:+.3f}).")

    d500 = merged[merged["k"] == 500].copy()
    rho = float(d500["truth_len"].corr(d500["delta_f1"], method="spearman"))
    lines.append(rf"\item Spearman $\rho$ between truth length and $\Delta$F1 @500: {rho:.3f}.")
    lines.append(r"\end{itemize}")
    (gen / "len_findings.tex").write_text("\n".join(lines) + "\n")

    # Delta-vs-K curves.
    delta_series = []
    for g, color in [("short(<=150)", "#1f77b4"), ("long(151-300)", "#d62728")]:
        ys = [float(merged[(merged["len_group"] == g) & (merged["k"] == k)]["delta_f1"].mean()) for k in ks]
        write_xy_dat(out_path=gen / f"delta_vs_k_{g.split('(')[0].strip()}.dat", xs=[float(k) for k in ks], ys=ys, x_name="k", y_name="delta")
        delta_series.append((g, color, f"delta_vs_k_{g.split('(')[0].strip()}.dat"))

    delta_lines: list[str] = []
    delta_lines.append(r"\begin{tikzpicture}")
    delta_lines.append(r"\definecolor{c_short}{HTML}{%s}" % latex_color("#1f77b4"))
    delta_lines.append(r"\definecolor{c_long}{HTML}{%s}" % latex_color("#d62728"))
    delta_lines.append(r"\begin{axis}[")
    delta_lines.append(r"  width=0.9\linewidth, height=0.55\linewidth,")
    delta_lines.append(r"  xmin=0, xlabel={K (structures kept)}, ylabel={Mean $\Delta$F1 vs @1},")
    delta_lines.append(r"  title={Do longer RNAs gain more from larger K?},")
    delta_lines.append(r"  grid=both, grid style={black!10},")
    delta_lines.append(r"  legend pos=north west, legend cell align=left,")
    delta_lines.append(r"  xtick={1,50,100,200,500},")
    delta_lines.append(r"]")
    delta_lines.append(r"\addplot+[very thick, mark=none, draw=c_short] table[x=k,y=delta] {generated/delta_vs_k_short.dat};")
    delta_lines.append(r"\addlegendentry{short (<=150)}")
    delta_lines.append(r"\addplot+[very thick, mark=none, draw=c_long] table[x=k,y=delta] {generated/delta_vs_k_long.dat};")
    delta_lines.append(r"\addlegendentry{long (151-300)}")
    delta_lines.append(r"\end{axis}")
    delta_lines.append(r"\end{tikzpicture}")
    (gen / "delta_vs_k_by_len.tex").write_text("\n".join(delta_lines) + "\n")

    # Mean curves.
    k_xs = [float(k) for k in ks]
    mean_f1 = [float(df[df["k"] == k]["f1"].mean()) for k in ks]
    fail = [float((df[df["k"] == k]["f1"] <= 0.0).mean()) for k in ks]
    write_xy_dat(out_path=gen / "mean_f1_vs_k.dat", xs=k_xs, ys=mean_f1, x_name="k", y_name="f1")
    write_xy_dat(out_path=gen / "fail_vs_k.dat", xs=k_xs, ys=fail, x_name="k", y_name="fail")

    # CDF curves for each K.
    series = [
        SeriesSpec(1, "@1", "#1f77b4"),
        SeriesSpec(50, "@50", "#ff7f0e"),
        SeriesSpec(100, "@100", "#2ca02c"),
        SeriesSpec(200, "@200", "#d62728"),
        SeriesSpec(500, "@500", "#9467bd"),
    ]
    for s in series:
        vals = df[df["k"] == s.k]["f1"].astype(float).tolist()
        write_cdf_dat(out_path=gen / f"cdf_{s.k}.dat", values=vals)

    def cname(k: int) -> str:
        return f"c_k{k}"

    def define_colors() -> list[str]:
        return [rf"\definecolor{{{cname(s.k)}}}{{HTML}}{{{latex_color(s.color)}}}" for s in series]

    # Plots (pgfplots snippets).
    mean_lines: list[str] = []
    mean_lines.append(r"\begin{tikzpicture}")
    mean_lines.append(r"\begin{axis}[")
    mean_lines.append(r"  width=0.9\linewidth, height=0.55\linewidth,")
    mean_lines.append(r"  xmin=0, xlabel={K (structures kept)}, ylabel={Mean F1},")
    mean_lines.append(r"  title={Mean best-of-K F1 vs output size},")
    mean_lines.append(r"  grid=both, grid style={black!10},")
    mean_lines.append(r"  xtick={1,50,100,200,500},")
    mean_lines.append(r"]")
    mean_lines.append(r"\addplot+[very thick, mark=none] table[x=k,y=f1] {generated/mean_f1_vs_k.dat};")
    mean_lines.append(r"\end{axis}")
    mean_lines.append(r"\end{tikzpicture}")
    (gen / "mean_f1_vs_k.tex").write_text("\n".join(mean_lines) + "\n")

    fail_lines: list[str] = []
    fail_lines.append(r"\begin{tikzpicture}")
    fail_lines.append(r"\begin{axis}[")
    fail_lines.append(r"  width=0.9\linewidth, height=0.55\linewidth,")
    fail_lines.append(r"  xmin=0, ymin=0, ymax=1,")
    fail_lines.append(r"  xlabel={K (structures kept)}, ylabel={Fail rate (F1=0)},")
    fail_lines.append(r"  title={Fail rate vs output size},")
    fail_lines.append(r"  grid=both, grid style={black!10},")
    fail_lines.append(r"  xtick={1,50,100,200,500},")
    fail_lines.append(r"]")
    fail_lines.append(r"\addplot+[very thick, mark=none] table[x=k,y=fail] {generated/fail_vs_k.dat};")
    fail_lines.append(r"\end{axis}")
    fail_lines.append(r"\end{tikzpicture}")
    (gen / "fail_vs_k.tex").write_text("\n".join(fail_lines) + "\n")

    cdf_lines: list[str] = []
    cdf_lines.append(r"\begin{tikzpicture}")
    cdf_lines.extend(define_colors())
    cdf_lines.append(r"\begin{axis}[")
    cdf_lines.append(r"  width=\linewidth, height=0.62\linewidth,")
    cdf_lines.append(r"  xmin=0, xmax=1, ymin=0, ymax=1,")
    cdf_lines.append(r"  xlabel={F1}, ylabel={CDF},")
    cdf_lines.append(r"  title={CDF of best-of-K F1},")
    cdf_lines.append(r"  grid=both, grid style={black!10},")
    cdf_lines.append(r"  legend pos=north west, legend cell align=left,")
    cdf_lines.append(r"  legend style={font=\scriptsize},")
    cdf_lines.append(r"]")
    for s in series:
        vals = df[df["k"] == s.k]["f1"].astype(float).tolist()
        med = f1_at_cdf(vals, 0.50)
        label = s.label if med is None else f"{s.label} (p50={med:.3f})"
        cdf_lines.append(
            rf"\addplot+[const plot, very thick, draw={cname(s.k)}] table[x=f1,y=cdf] {{generated/cdf_{s.k}.dat}};"
        )
        cdf_lines.append(rf"\addlegendentry{{{label}}}")
    cdf_lines.append(r"\end{axis}")
    cdf_lines.append(r"\end{tikzpicture}")
    (gen / "cdf.tex").write_text("\n".join(cdf_lines) + "\n")


if __name__ == "__main__":
    main()
