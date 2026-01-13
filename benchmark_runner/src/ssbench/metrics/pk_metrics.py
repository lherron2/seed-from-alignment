"""Pseudoknot-aware metrics."""

from __future__ import annotations

from dataclasses import dataclass

from .pair_metrics import PairMetrics, compute_pair_metrics


@dataclass
class PKMetrics:
    nested: PairMetrics
    pk: PairMetrics


def compute_pk_metrics(
    pred_nested: set[tuple[int, int]],
    pred_pk: set[tuple[int, int]],
    ref_nested: set[tuple[int, int]],
    ref_pk: set[tuple[int, int]],
    length: int,
) -> PKMetrics:
    """Compute metrics for nested vs pseudoknotted pairs separately."""
    nested_metrics = compute_pair_metrics(pred_nested, ref_nested, length)
    pk_metrics = compute_pair_metrics(pred_pk, ref_pk, length)
    return PKMetrics(nested=nested_metrics, pk=pk_metrics)
