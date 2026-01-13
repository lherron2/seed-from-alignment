"""Pair-level metrics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class PairMetrics:
    precision: float
    recall: float
    f1: float
    mcc: float


def compute_pair_metrics(pred: set[tuple[int, int]], ref: set[tuple[int, int]], length: int) -> PairMetrics:
    """Compute precision, recall, F1, and MCC for pairs."""
    tp = len(pred & ref)
    fp = len(pred - ref)
    fn = len(ref - pred)
    total_pairs = length * (length - 1) // 2
    tn = total_pairs - tp - fp - fn

    precision = tp / (tp + fp) if tp + fp > 0 else 0.0
    recall = tp / (tp + fn) if tp + fn > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0

    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if denom == 0:
        mcc = 0.0
    else:
        mcc = (tp * tn - fp * fn) / (denom ** 0.5)

    return PairMetrics(precision=precision, recall=recall, f1=f1, mcc=mcc)
