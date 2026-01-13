"""Summarize metrics results."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def summarize_metrics(metrics_csv: str | Path, out_path: str | Path) -> None:
    df = pd.read_csv(metrics_csv)
    summary = []

    for split in sorted(df["split"].unique()):
        sub = df[df["split"] == split]
        summary.append({
            "split": split,
            "mean_f1": sub["f1"].mean(),
            "mean_precision": sub["precision"].mean(),
            "mean_recall": sub["recall"].mean(),
        })

    out = Path(out_path)
    lines = ["# ssbench summary", ""]
    for row in summary:
        lines.append(
            f"- {row['split']}: F1={row['mean_f1']:.3f}, P={row['mean_precision']:.3f}, R={row['mean_recall']:.3f}"
        )
    out.write_text("\n".join(lines) + "\n")
