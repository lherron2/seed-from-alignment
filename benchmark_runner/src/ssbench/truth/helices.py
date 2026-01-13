"""Helix detection utilities."""

from __future__ import annotations


def identify_helices(pairs: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """Identify helices as runs of adjacent stacked pairs."""
    if not pairs:
        return []
    pairs_sorted = sorted(pairs)
    helices: list[list[tuple[int, int]]] = []
    current = [pairs_sorted[0]]
    for i in range(1, len(pairs_sorted)):
        prev_i, prev_j = pairs_sorted[i - 1]
        cur_i, cur_j = pairs_sorted[i]
        if cur_i == prev_i + 1 and cur_j == prev_j - 1:
            current.append((cur_i, cur_j))
        else:
            helices.append(current)
            current = [(cur_i, cur_j)]
    helices.append(current)
    return helices


def helix_overlap(h1: list[tuple[int, int]], h2: list[tuple[int, int]]) -> float:
    """Compute overlap fraction between two helices by pair intersection."""
    set1 = set(h1)
    set2 = set(h2)
    if not set1 or not set2:
        return 0.0
    inter = len(set1 & set2)
    return inter / max(len(set1), len(set2))
