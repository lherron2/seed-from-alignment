#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ssbench.metrics.pair_metrics import compute_pair_metrics
from ssbench.predict.parse_dotbracket import pairs_from_dotbracket, parse_dotbracket


@dataclass(frozen=True)
class TargetMetrics:
    target_id: str
    pdb_id: str
    chain_id: str
    truth_len: int
    ref_pairs: int
    ref_pk_pairs: int
    ref_pair_density: float
    length_bucket: str
    rfam_id: str
    n_predictions: int
    top1_precision: float
    top1_recall: float
    top1_f1: float
    top1_mcc: float
    top100_rank_f1: int
    top100_pred_pairs: int
    top100_precision: float
    top100_recall: float
    top100_f1: float
    top100_mcc: float

    @property
    def delta_f1(self) -> float:
        return float(self.top100_f1 - self.top1_f1)

    @property
    def delta_mcc(self) -> float:
        return float(self.top100_mcc - self.top1_mcc)


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


def _first_rfam_hit_from_tblout(tblout_path: Path) -> str | None:
    if not tblout_path.exists():
        return None
    for line in tblout_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        accession = parts[1].strip()
        if accession and accession != "-":
            return accession
        return None
    return None


def _top_hit_from_tblout(tblout_path: Path) -> tuple[str, str, str] | None:
    if not tblout_path.exists():
        return None
    for line in tblout_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        target_name = parts[0].strip()
        accession = parts[1].strip()
        desc = " ".join(parts[18:]).strip() if len(parts) > 18 else ""
        return target_name, accession, desc
    return None


def is_ribosomal_rna_hit(target_name: str, accession: str, description: str) -> bool:
    _ = accession
    import re

    if re.search(r"(^|_)rRNA(_|$)", target_name):
        return True
    if "ribosomal RNA" in description:
        return True
    return False


def infer_rfam_hit(pred_dir: Path) -> tuple[str, str, str]:
    hit = _top_hit_from_tblout(pred_dir / "cmscan.tblout")
    if hit is not None:
        target_name, accession, desc = hit
        if accession and accession != "-":
            return accession, target_name, desc
    # Fallback: infer from cmfetch output name.
    cms = sorted(pred_dir.glob("RF*.cm"))
    if cms:
        return cms[0].stem, cms[0].stem, ""
    return "no_rfam_hit", "", ""


def best_by_f1(metrics_rows: list[tuple[int, float, float, float, float, int]]) -> tuple[int, float, float, float, float]:
    # rows: (rank, precision, recall, f1, mcc, n_pred_pairs)
    # Prefer higher F1; break ties by MCC then recall then fewer predicted pairs then rank.
    best = max(
        metrics_rows,
        key=lambda r: (r[3], r[4], r[2], -r[5], -r[0]),
    )
    rank, precision, recall, f1, mcc, _n_pred_pairs = best
    return rank, precision, recall, f1, mcc


def render_table_tex(
    out_path: Path,
    *,
    caption: str,
    label: str,
    headers: list[str],
    rows: list[list[str]],
    align: str,
) -> None:
    lines: list[str] = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(rf"\begin{{tabular}}{{{align}}}")
    lines.append(r"\toprule")
    lines.append(" & ".join(headers) + r" \\")
    lines.append(r"\midrule")
    for row in rows:
        lines.append(" & ".join(row) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    out_path.write_text("\n".join(lines) + "\n")


def detok(text: str) -> str:
    return r"\texttt{\detokenize{" + text + "}}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate data + LaTeX snippets for the FR3D/BGSU representative benchmark report.",
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--truth-dir", required=True)
    parser.add_argument("--pred-dir", required=True)
    parser.add_argument("--k", type=int, default=100)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--exclude-ribosomal",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude targets whose cmscan top hit appears to be a ribosomal RNA (rRNA).",
    )
    args = parser.parse_args()

    manifest = Path(args.manifest)
    truth_dir = Path(args.truth_dir)
    pred_dir = Path(args.pred_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_manifest = pd.read_csv(manifest)
    required_cols = {"target_id", "pdb_id", "chain_id"}
    missing_cols = required_cols - set(df_manifest.columns)
    if missing_cols:
        raise SystemExit(f"manifest missing required columns: {sorted(missing_cols)}")

    records: list[TargetMetrics] = []
    for _, row in df_manifest.iterrows():
        target_id = str(row["target_id"])
        pdb_id = str(row["pdb_id"]).lower()
        chain_id = str(row["chain_id"])
        truth_path = truth_dir / f"{target_id}.json"
        if not truth_path.exists():
            continue
        pred_subdir = pred_dir / target_id.replace("|", "_")
        pred_db = pred_subdir / "predictions.db"
        if not pred_db.exists():
            continue

        truth = json.loads(truth_path.read_text())
        seq = str(truth.get("sequence", ""))
        L = len(seq)
        ref_pairs = {tuple(p) for p in truth.get("canonical_pairs", [])}
        ref_pk_pairs = {tuple(p) for p in truth.get("pk_pairs", [])}

        seq2, dotbrs = parse_dotbracket(pred_db.read_text())
        if len(seq2) != L:
            raise SystemExit(
                f"Sequence length mismatch for {target_id}: truth={L} pred={len(seq2)} ({pred_db})"
            )
        if not dotbrs:
            dotbrs = ["." * L]

        k_eff = min(int(args.k), len(dotbrs))
        dotbrs_k = dotbrs[:k_eff]

        metrics_rows = []
        for rank, s in enumerate(dotbrs_k, start=1):
            pred_pairs = set(pairs_from_dotbracket(s))
            m = compute_pair_metrics(pred_pairs, ref_pairs, L)
            metrics_rows.append((rank, m.precision, m.recall, m.f1, m.mcc, len(pred_pairs)))

        top1_rank, top1_precision, top1_recall, top1_f1, top1_mcc, top1_npairs = metrics_rows[0]
        best_rank, best_precision, best_recall, best_f1, best_mcc = best_by_f1(metrics_rows)
        best_npairs = metrics_rows[best_rank - 1][5]

        bucket = assign_length_bucket(L)
        rfam_id, rfam_name, rfam_desc = infer_rfam_hit(pred_subdir)
        if bool(args.exclude_ribosomal) and is_ribosomal_rna_hit(rfam_name, rfam_id, rfam_desc):
            continue
        ref_density = (len(ref_pairs) / L) if L > 0 else 0.0

        records.append(
            TargetMetrics(
                target_id=target_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                truth_len=int(L),
                ref_pairs=int(len(ref_pairs)),
                ref_pk_pairs=int(len(ref_pk_pairs)),
                ref_pair_density=float(ref_density),
                length_bucket=bucket,
                rfam_id=rfam_id,
                n_predictions=int(len(dotbrs_k)),
                top1_precision=float(top1_precision),
                top1_recall=float(top1_recall),
                top1_f1=float(top1_f1),
                top1_mcc=float(top1_mcc),
                top100_rank_f1=int(best_rank),
                top100_pred_pairs=int(best_npairs),
                top100_precision=float(best_precision),
                top100_recall=float(best_recall),
                top100_f1=float(best_f1),
                top100_mcc=float(best_mcc),
            )
        )

    if not records:
        raise SystemExit("No targets produced metrics. Check manifest/truth/pred paths.")

    df = pd.DataFrame([r.__dict__ | {"delta_f1": r.delta_f1, "delta_mcc": r.delta_mcc} for r in records])
    df = df.sort_values(["target_id"]).reset_index(drop=True)

    per_target_csv = out_dir / "per_target_metrics.csv"
    df.to_csv(per_target_csv, index=False)

    # Add derived columns used for analysis.
    df["rfam_hit"] = df["rfam_id"] != "no_rfam_hit"
    df["pk_frac"] = df["ref_pk_pairs"] / df["ref_pairs"].clip(lower=1)

    # Overall summary.
    def summarize(sub: pd.DataFrame) -> dict[str, float | int]:
        return {
            "n": int(len(sub)),
            "top1_precision_mean": float(sub["top1_precision"].mean()),
            "top1_recall_mean": float(sub["top1_recall"].mean()),
            "top1_f1_mean": float(sub["top1_f1"].mean()),
            "top100_precision_mean": float(sub["top100_precision"].mean()),
            "top100_recall_mean": float(sub["top100_recall"].mean()),
            "top100_f1_mean": float(sub["top100_f1"].mean()),
            "top1_mcc_mean": float(sub["top1_mcc"].mean()),
            "top100_mcc_mean": float(sub["top100_mcc"].mean()),
            "top1_precision_median": float(sub["top1_precision"].median()),
            "top1_recall_median": float(sub["top1_recall"].median()),
            "top1_f1_median": float(sub["top1_f1"].median()),
            "top100_precision_median": float(sub["top100_precision"].median()),
            "top100_recall_median": float(sub["top100_recall"].median()),
            "top100_f1_median": float(sub["top100_f1"].median()),
            "delta_f1_mean": float(sub["delta_f1"].mean()),
            "delta_f1_median": float(sub["delta_f1"].median()),
            "p_best100_ge_0p7": float((sub["top100_f1"] >= 0.7).mean()),
            "p_best100_ge_0p9": float((sub["top100_f1"] >= 0.9).mean()),
            "min_n_predictions": int(sub["n_predictions"].min()),
            "p_best_rank_le_10": float((sub["top100_rank_f1"] <= 10).mean()),
            "best_rank_median": float(sub["top100_rank_f1"].median()),
            "best_rank_p90": float(sub["top100_rank_f1"].quantile(0.9)),
            "best_rank_max": int(sub["top100_rank_f1"].max()),
            "p_delta_f1_ge_0p3": float((sub["delta_f1"] >= 0.3).mean()),
        }

    overall = summarize(df)
    (out_dir / "overall_summary.json").write_text(json.dumps(overall, indent=2) + "\n")

    # Group: length buckets.
    by_len = (
        df.groupby("length_bucket", dropna=False)
        .apply(lambda g: pd.Series(summarize(g)))
        .reset_index()
    )
    # Keep a stable and human-friendly ordering.
    bucket_order = {"30-80": 0, "81-150": 1, "151-300": 2, "301-400": 3, "other": 4}
    by_len["_bucket_order"] = by_len["length_bucket"].map(bucket_order).fillna(999).astype(int)
    by_len = by_len.sort_values(["_bucket_order", "length_bucket"]).drop(columns=["_bucket_order"])
    by_len.to_csv(out_dir / "by_length_bucket.csv", index=False)

    # Group: rfam id.
    by_rfam = (
        df.groupby("rfam_id", dropna=False)
        .apply(lambda g: pd.Series(summarize(g)))
        .reset_index()
        .sort_values(["n", "top100_f1_mean"], ascending=[False, True])
    )
    by_rfam.to_csv(out_dir / "by_rfam.csv", index=False)

    # Group: Rfam hit vs fallback.
    by_rfam_hit = (
        df.groupby("rfam_hit", dropna=False)
        .apply(lambda g: pd.Series(summarize(g)))
        .reset_index()
        .sort_values("rfam_hit", ascending=False)
    )
    by_rfam_hit.to_csv(out_dir / "by_rfam_hit.csv", index=False)

    # Group: pseudoknot fraction buckets.
    df["pk_bucket"] = pd.cut(
        df["pk_frac"],
        bins=[-1e-9, 0.0, 0.2, 0.4, 1.0],
        labels=["0", "(0,0.2)", "[0.2,0.4)", "[0.4,1]"],
        include_lowest=True,
    )
    by_pk = (
        df.groupby("pk_bucket", dropna=False)
        .apply(lambda g: pd.Series(summarize(g)))
        .reset_index()
    )
    by_pk.to_csv(out_dir / "by_pk_bucket.csv", index=False)

    # Correlations (Pearson).
    corr_rows = []
    for feature in ["truth_len", "ref_pairs", "ref_pair_density", "ref_pk_pairs", "pk_frac"]:
        x = df[feature].astype(float).to_numpy()
        corr_top1 = float(np.corrcoef(x, df["top1_f1"].astype(float).to_numpy())[0, 1])
        corr_top100 = float(np.corrcoef(x, df["top100_f1"].astype(float).to_numpy())[0, 1])
        corr_rows.append({"feature": feature, "corr_top1_f1": corr_top1, "corr_top100_f1": corr_top100})
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(out_dir / "correlations.csv", index=False)

    # CDF data (rank fraction vs metrics).
    n = len(df)
    top1_sorted = sorted(df["top1_f1"].astype(float).tolist())
    top100_sorted = sorted(df["top100_f1"].astype(float).tolist())
    delta_sorted = sorted(df["delta_f1"].astype(float).tolist())
    cdf = pd.DataFrame(
        {
            "frac": [(i + 1) / n for i in range(n)],
            "top1_f1": top1_sorted,
            "top100_f1": top100_sorted,
            "delta_f1": delta_sorted,
        }
    )
    cdf.to_csv(out_dir / "cdf_f1.csv", index=False)

    # Scatter data.
    df_scatter = df[["truth_len", "ref_pairs", "ref_pk_pairs", "ref_pair_density", "top1_f1", "top100_f1", "delta_f1"]]
    df_scatter.to_csv(out_dir / "scatter_len.csv", index=False)

    # LaTeX snippets.
    tex_dir = out_dir / "tables"
    tex_dir.mkdir(parents=True, exist_ok=True)

    n_total = int(overall["n"])
    k_val = int(args.k)
    render_table_tex(
        tex_dir / "overall.tex",
        caption=f"Overall performance on the FR3D/BGSU representative benchmark (N={n_total}, best-of-{k_val}).",
        label="tab:overall",
        headers=[
            "Metric",
            "Top-1",
            "Best-of-100",
        ],
        rows=[
            [
                "Mean Precision",
                f"{overall['top1_precision_mean']:.3f}",
                f"{overall['top100_precision_mean']:.3f}",
            ],
            [
                "Mean Recall",
                f"{overall['top1_recall_mean']:.3f}",
                f"{overall['top100_recall_mean']:.3f}",
            ],
            ["Mean F1", f"{overall['top1_f1_mean']:.3f}", f"{overall['top100_f1_mean']:.3f}"],
            ["Median F1", f"{overall['top1_f1_median']:.3f}", f"{overall['top100_f1_median']:.3f}"],
            ["Mean MCC", f"{overall['top1_mcc_mean']:.3f}", f"{overall['top100_mcc_mean']:.3f}"],
            ["Mean $\\Delta$F1", "\\multicolumn{2}{c}{" + f"{overall['delta_f1_mean']:.3f}" + "}"],
            ["Min \\#preds", "\\multicolumn{2}{c}{" + f"{overall['min_n_predictions']}" + "}"],
            ["Median rank(best@100)", "\\multicolumn{2}{c}{" + f"{overall['best_rank_median']:.1f}" + "}"],
            [
                "Frac(rank(best@100)$\\le 10$)",
                "\\multicolumn{2}{c}{"
                + f"{float(overall['p_best_rank_le_10'])*100:.2f}\\%"
                + "}",
            ],
            [
                "Frac($\\Delta$F1$\\ge 0.3$)",
                "\\multicolumn{2}{c}{"
                + f"{float(overall['p_delta_f1_ge_0p3'])*100:.2f}\\%"
                + "}",
            ],
            [
                "Frac(best-of-100 F1$\\ge0.7$)",
                "\\multicolumn{2}{c}{"
                + f"{float(overall['p_best100_ge_0p7'])*100:.2f}\\%"
                + "}",
            ],
            [
                "Frac(best-of-100 F1$\\ge0.9$)",
                "\\multicolumn{2}{c}{"
                + f"{float(overall['p_best100_ge_0p9'])*100:.2f}\\%"
                + "}",
            ],
        ],
        align="lcc",
    )

    # Length bucket table (means).
    rows = []
    for _, r in by_len.iterrows():
        rows.append(
            [
                str(r["length_bucket"]),
                str(int(r["n"])),
                f"{float(r['top1_precision_mean']):.3f}",
                f"{float(r['top1_recall_mean']):.3f}",
                f"{float(r['top1_f1_mean']):.3f}",
                f"{float(r['top1_mcc_mean']):.3f}",
                f"{float(r['top100_precision_mean']):.3f}",
                f"{float(r['top100_recall_mean']):.3f}",
                f"{float(r['top100_f1_mean']):.3f}",
                f"{float(r['top100_mcc_mean']):.3f}",
                f"{float(r['delta_f1_mean']):.3f}",
            ]
        )
    render_table_tex(
        tex_dir / "by_length_bucket.tex",
        caption="Performance by truth length bucket (means).",
        label="tab:by_length_bucket",
        headers=[
            "Bucket",
            "N",
            "P@1",
            "R@1",
            "F1@1",
            "MCC@1",
            "P@100",
            "R@100",
            "F1@100",
            "MCC@100",
            "$\\Delta$F1",
        ],
        rows=rows,
        align="lrrrrrrrrrr",
    )

    # Rfam hit vs fallback table.
    rows = []
    for _, r in by_rfam_hit.iterrows():
        label = "rfam\\_hit" if bool(r["rfam_hit"]) else "no\\_rfam\\_hit"
        rows.append(
            [
                label,
                str(int(r["n"])),
                f"{float(r['top1_precision_mean']):.3f}",
                f"{float(r['top1_recall_mean']):.3f}",
                f"{float(r['top1_f1_mean']):.3f}",
                f"{float(r['top1_mcc_mean']):.3f}",
                f"{float(r['top100_precision_mean']):.3f}",
                f"{float(r['top100_recall_mean']):.3f}",
                f"{float(r['top100_f1_mean']):.3f}",
                f"{float(r['top100_mcc_mean']):.3f}",
                f"{float(r['delta_f1_mean']):.3f}",
            ]
        )
    render_table_tex(
        tex_dir / "by_rfam_hit.tex",
        caption="Performance split by whether Infernal found an Rfam CM hit (means).",
        label="tab:by_rfam_hit",
        headers=[
            "Group",
            "N",
            "P@1",
            "R@1",
            "F1@1",
            "MCC@1",
            "P@100",
            "R@100",
            "F1@100",
            "MCC@100",
            "$\\Delta$F1",
        ],
        rows=rows,
        align="lrrrrrrrrrr",
    )

    # Pseudoknot bucket table.
    rows = []
    for _, r in by_pk.iterrows():
        rows.append(
            [
                detok(str(r["pk_bucket"])),
                str(int(r["n"])),
                f"{float(r['top1_f1_mean']):.3f}",
                f"{float(r['top100_f1_mean']):.3f}",
                f"{float(r['delta_f1_mean']):.3f}",
            ]
        )
    render_table_tex(
        tex_dir / "by_pk_bucket.tex",
        caption="Performance vs. pseudoknot fraction in the truth (means).",
        label="tab:by_pk_bucket",
        headers=["PK fraction bucket", "N", "F1@1", "F1@100", "$\\Delta$F1"],
        rows=rows,
        align="lrrrr",
    )

    # Correlation table.
    rows = []
    for _, r in corr_df.iterrows():
        rows.append(
            [
                detok(str(r["feature"])),
                f"{float(r['corr_top1_f1']):.3f}",
                f"{float(r['corr_top100_f1']):.3f}",
            ]
        )
    render_table_tex(
        tex_dir / "correlations.tex",
        caption="Pearson correlation of features with per-target F1 (higher magnitude indicates stronger relationship).",
        label="tab:correlations",
        headers=["Feature", "Corr(F1@1)", "Corr(F1@100)"],
        rows=rows,
        align="lrr",
    )

    # Top Rfam families by count.
    top_rfam = by_rfam.sort_values(["n", "top100_f1_mean"], ascending=[False, True]).head(15)
    rows = []
    for _, r in top_rfam.iterrows():
        rows.append(
            [
                detok(str(r["rfam_id"])),
                str(int(r["n"])),
                f"{float(r['top1_f1_mean']):.3f}",
                f"{float(r['top100_f1_mean']):.3f}",
                f"{float(r['delta_f1_mean']):.3f}",
            ]
        )
    render_table_tex(
        tex_dir / "top_rfam.tex",
        caption="Top Rfam IDs by frequency in this benchmark (means).",
        label="tab:top_rfam",
        headers=["Rfam", "N", "F1@1", "F1@100", "$\\Delta$F1"],
        rows=rows,
        align="lrrrr",
    )

    # Worst targets (best-of-100 F1).
    worst = df.sort_values("top100_f1").head(12)
    rows = []
    for _, r in worst.iterrows():
        rows.append(
            [
                detok(str(r["target_id"])),
                detok(str(r["rfam_id"])),
                f"{int(r['truth_len'])}",
                f"{int(r['ref_pairs'])}",
                f"{float(r['top1_f1']):.3f}",
                f"{float(r['top100_f1']):.3f}",
                f"{int(r['top100_rank_f1'])}",
            ]
        )
    render_table_tex(
        tex_dir / "worst_targets.tex",
        caption="Worst targets by best-of-100 F1 (lower is worse).",
        label="tab:worst_targets",
        headers=["Target", "Rfam", "L", "$|P_{ref}|$", "F1@1", "F1@100", "Rank(best)"],
        rows=rows,
        align="llrrrrr",
    )

    # Biggest improvements (best-of-100 - top-1).
    improv = df.sort_values("delta_f1", ascending=False).head(12)
    rows = []
    for _, r in improv.iterrows():
        rows.append(
            [
                detok(str(r["target_id"])),
                detok(str(r["rfam_id"])),
                f"{int(r['truth_len'])}",
                f"{float(r['top1_f1']):.3f}",
                f"{float(r['top100_f1']):.3f}",
                f"{float(r['delta_f1']):.3f}",
                f"{int(r['top100_rank_f1'])}",
            ]
        )
    render_table_tex(
        tex_dir / "biggest_improvements.tex",
        caption="Targets with the largest improvement from rank-1 to best-of-100 (ranking/sampling sensitivity).",
        label="tab:biggest_improvements",
        headers=["Target", "Rfam", "L", "F1@1", "F1@100", "$\\Delta$F1", "Rank(best)"],
        rows=rows,
        align="llrrrrr",
    )

    # Full Rfam table in LaTeX (for appendix).
    rfam_top1_tex = tex_dir / "by_rfam_top1_longtable.tex"
    rfam_top100_tex = tex_dir / "by_rfam_top100_longtable.tex"

    def write_rfam_longtable(path: Path, *, caption: str, label: str, cols: list[str], rows: list[str]) -> None:
        lines = []
        lines.append(rf"\begin{{longtable}}{{{cols[0]}}}")
        lines.append(rf"\caption{{{caption}}}\label{{{label}}}\\")
        lines.append(r"\toprule")
        lines.append(cols[1] + r" \\")
        lines.append(r"\midrule")
        lines.append(r"\endfirsthead")
        lines.append(r"\toprule")
        lines.append(cols[1] + r" \\")
        lines.append(r"\midrule")
        lines.append(r"\endhead")
        lines.extend(rows)
        lines.append(r"\bottomrule")
        lines.append(r"\end{longtable}")
        path.write_text("\n".join(lines) + "\n")

    by_rfam_sorted = by_rfam.sort_values(["n", "top100_f1_mean"], ascending=[False, True])
    top1_rows = []
    top100_rows = []
    for _, r in by_rfam_sorted.iterrows():
        rfam_label = detok(str(r["rfam_id"]))
        top1_rows.append(
            f"{rfam_label} & {int(r['n'])} & {float(r['top1_precision_mean']):.3f} & {float(r['top1_recall_mean']):.3f} & {float(r['top1_f1_mean']):.3f} & {float(r['top1_mcc_mean']):.3f} \\\\"
        )
        top100_rows.append(
            f"{rfam_label} & {int(r['n'])} & {float(r['top100_precision_mean']):.3f} & {float(r['top100_recall_mean']):.3f} & {float(r['top100_f1_mean']):.3f} & {float(r['top100_mcc_mean']):.3f} \\\\"
        )

    write_rfam_longtable(
        rfam_top1_tex,
        caption="Performance by Rfam ID using the rank-1 prediction (means).",
        label="tab:by_rfam_top1",
        cols=["lrrrrr", "Rfam & N & P@1 & R@1 & F1@1 & MCC@1"],
        rows=top1_rows,
    )

    write_rfam_longtable(
        rfam_top100_tex,
        caption="Performance by Rfam ID using best-of-100 (means).",
        label="tab:by_rfam_top100",
        cols=["lrrrrr", "Rfam & N & P@100 & R@100 & F1@100 & MCC@100"],
        rows=top100_rows,
    )

    print(f"Wrote per-target metrics: {per_target_csv}")
    print(f"Wrote report assets to: {out_dir}")


if __name__ == "__main__":
    main()
