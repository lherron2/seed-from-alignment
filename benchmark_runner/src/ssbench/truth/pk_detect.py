"""Pseudoknot detection utilities."""

from __future__ import annotations


def crosses(p1: tuple[int, int], p2: tuple[int, int]) -> bool:
    """Return True if p1 and p2 cross (i<k<j<l or k<i<l<j)."""
    i, j = p1
    k, l = p2
    if i > j:
        i, j = j, i
    if k > l:
        k, l = l, k
    return (i < k < j < l) or (k < i < l < j)


def split_nested_pk(pairs: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Split pairs into nested and pseudoknot sets based on crossings."""
    nested: list[tuple[int, int]] = []
    pk: list[tuple[int, int]] = []
    for i, p in enumerate(pairs):
        is_pk = False
        for q in pairs[i + 1 :]:
            if crosses(p, q):
                is_pk = True
                break
        if is_pk:
            pk.append(p)
        else:
            nested.append(p)
    return nested, pk
