#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MethodSpec:
    key: str
    label: str
    metrics_csv: Path
    color: str
    linestyle: str = "-"


def read_csv_column(path: Path, col: str) -> list[float]:
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return []
        try:
            idx = header.index(col)
        except ValueError as exc:
            raise SystemExit(f"missing column {col!r} in {path}") from exc

        out: list[float] = []
        for row in reader:
            if not row or idx >= len(row):
                continue
            cell = row[idx].strip()
            if not cell:
                continue
            try:
                out.append(float(cell))
            except ValueError:
                continue
        return out


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
    # ECDF is a right-continuous step function. The "x at CDF==cdf" is the smallest
    # x such that ECDF(x) >= cdf.
    import math

    k = int(math.ceil(cdf * n))  # 1..n
    k = max(1, min(n, k))
    return vals[k - 1]


def _allclose(a: list[float], b: list[float], tol: float = 1e-12) -> bool:
    if len(a) != len(b):
        return False
    return all(abs(x - y) <= tol for x, y in zip(a, b, strict=True))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot ECDFs (top-1 and best-of-100) for the FR3D under400 ensemble_merge report."
    )
    parser.add_argument("--out-dir", required=True, help="Directory to write plots into")
    parser.add_argument("--pipeline", required=True, help="RNAnneal-ss metrics CSV")
    parser.add_argument("--ensemble", required=True, help="Ensemble metrics CSV")
    parser.add_argument("--rnastructure", required=True, help="RNAstructure metrics CSV")
    parser.add_argument("--linearfold", required=True, help="LinearFold-V metrics CSV")
    parser.add_argument("--eternafold", required=True, help="EternaFold metrics CSV")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Matplotlib is an optional dependency for this repo; import lazily so other tooling works.
    # In sandboxed environments, $HOME/.cache may be unwritable; redirect MPL cache to out_dir.
    mpl_cfg = out_dir / ".mplconfig"
    mpl_cfg.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cfg))
    import matplotlib.pyplot as plt  # type: ignore[import-not-found]

    # Keep a stable color palette and use line styles so overlapping curves remain visible.
    # (In particular, Ensemble top-1 can coincide with RNAnneal-ss top-1 by construction.)
    methods = [
        MethodSpec("rnastructure", "RNAstructure", Path(args.rnastructure), "#7f7f7f", "-"),  # gray
        MethodSpec("linearfold", "LinearFold-V", Path(args.linearfold), "#ff7f0e", "-"),  # orange
        MethodSpec("eternafold", "EternaFold", Path(args.eternafold), "#2ca02c", "-"),  # green
        MethodSpec("pipeline", "RNAnneal-ss", Path(args.pipeline), "#1f77b4", "--"),  # blue dashed
        MethodSpec("ensemble", "Ensemble", Path(args.ensemble), "#9467bd", "-"),  # purple
    ]

    series: dict[str, dict[str, list[float]]] = {}
    for m in methods:
        series[m.key] = {
            "f1": read_csv_column(m.metrics_csv, "f1"),
            "best_of_k_f1": read_csv_column(m.metrics_csv, "best_of_k_f1"),
        }

    def plot(which: str, col: str, out_name: str) -> None:
        fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=200)
        plot_methods = methods
        if col == "best_of_k_f1":
            # If a method never produces more than one effective candidate (e.g., RNAstructure MFE
            # padded to 100), its best-of-100 is identical to top-1 and isn't meaningful as a top-K
            # curve. Drop it from the @100 CDF.
            filtered: list[MethodSpec] = []
            for m in methods:
                if _allclose(series[m.key]["best_of_k_f1"], series[m.key]["f1"]):
                    continue
                filtered.append(m)
            plot_methods = filtered

        for m in plot_methods:
            vals = series[m.key][col]
            x, y = _ecdf(vals)
            med = f1_at_cdf(vals, 0.50)
            label = m.label if med is None else f"{m.label} (p50={med:.3f})"
            ax.step(
                x,
                y,
                where="post",
                linewidth=2.0,
                color=m.color,
                linestyle=m.linestyle,
                label=label,
            )
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel("F1")
        ax.set_ylabel("CDF")
        ax.set_title(which)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(loc="lower right", frameon=True, fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / f"{out_name}.pdf")
        fig.savefig(out_dir / f"{out_name}.png")
        plt.close(fig)

    plot("Top-1 (F1@1)", "f1", "cdf_top1")
    plot("Best-of-100 (F1@100)", "best_of_k_f1", "cdf_best100")


if __name__ == "__main__":
    main()
