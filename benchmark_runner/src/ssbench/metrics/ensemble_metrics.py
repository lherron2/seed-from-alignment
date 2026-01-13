"""Ensemble metrics for top-K predictions."""

from __future__ import annotations

from dataclasses import dataclass

from .pair_metrics import compute_pair_metrics


@dataclass
class EnsembleMetrics:
    best_f1: float
    top1_f1: float
    coverage: dict[int, float]


def compute_ensemble_metrics(
    preds: list[set[tuple[int, int]]],
    ref: set[tuple[int, int]],
    length: int,
    coverage_points: list[int] | None = None,
) -> EnsembleMetrics:
    """Compute best-of-K, top-1, and coverage curve metrics."""
    if coverage_points is None:
        coverage_points = [1, 2, 5, 10, 20, 50, 100]

    f1s = [compute_pair_metrics(p, ref, length).f1 for p in preds]
    best_f1 = max(f1s) if f1s else 0.0
    top1_f1 = f1s[0] if f1s else 0.0

    coverage = {}
    for k in coverage_points:
        n = min(k, len(preds))
        if n == 0:
            coverage[k] = 0.0
        else:
            coverage[k] = max(f1s[:n])

    return EnsembleMetrics(best_f1=best_f1, top1_f1=top1_f1, coverage=coverage)
