"""
Stem-level statistics for a set of base pairs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class StemStats:
    n_stems: int
    n_len1: int
    n_len2: int
    sum_log1p_len: float
    len1_pairs: set[tuple[int, int]]


def stem_stats(pairs: set[tuple[int, int]]) -> StemStats:
    if not pairs:
        return StemStats(0, 0, 0, 0.0, set())

    normalized = {(min(i, j), max(i, j)) for i, j in pairs}
    len1_pairs: set[tuple[int, int]] = set()
    stem_lengths: list[int] = []

    for i, j in normalized:
        if (i - 1, j + 1) in normalized:
            continue
        length = 1
        while (i + length, j - length) in normalized:
            length += 1
        stem_lengths.append(length)
        if length == 1:
            len1_pairs.add((i, j))

    n_stems = len(stem_lengths)
    n_len1 = sum(1 for l in stem_lengths if l == 1)
    n_len2 = sum(1 for l in stem_lengths if l == 2)
    sum_log1p_len = sum(math.log1p(l) for l in stem_lengths)
    return StemStats(n_stems, n_len1, n_len2, sum_log1p_len, len1_pairs)
