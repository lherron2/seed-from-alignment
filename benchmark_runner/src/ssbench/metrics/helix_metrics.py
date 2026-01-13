"""Helix-level metrics."""

from __future__ import annotations

from dataclasses import dataclass

from ..truth.helices import helix_overlap, identify_helices


@dataclass
class HelixMetrics:
    precision: float
    recall: float
    lonely_pair_rate: float


def compute_helix_metrics(
    pred_pairs: list[tuple[int, int]],
    ref_pairs: list[tuple[int, int]],
    overlap_threshold: float = 0.5,
) -> HelixMetrics:
    """Compute helix precision/recall and lonely pair rate."""
    pred_helices = identify_helices(pred_pairs)
    ref_helices = identify_helices(ref_pairs)

    matched_pred = 0
    matched_ref = 0

    for ph in pred_helices:
        for rh in ref_helices:
            if helix_overlap(ph, rh) >= overlap_threshold:
                matched_pred += 1
                break

    for rh in ref_helices:
        for ph in pred_helices:
            if helix_overlap(ph, rh) >= overlap_threshold:
                matched_ref += 1
                break

    precision = matched_pred / len(pred_helices) if pred_helices else 0.0
    recall = matched_ref / len(ref_helices) if ref_helices else 0.0

    # Lonely pair rate: fraction of predicted pairs not in any helix of length >= 2
    lonely = 0
    for helix in pred_helices:
        if len(helix) == 1:
            lonely += 1
    lonely_pairs = lonely
    total_pairs = len(pred_pairs)
    lonely_pair_rate = lonely_pairs / total_pairs if total_pairs else 0.0

    return HelixMetrics(precision=precision, recall=recall, lonely_pair_rate=lonely_pair_rate)
