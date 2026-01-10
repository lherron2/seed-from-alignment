"""
Pseudoknot computation module for MCMC RNA structure sampling.

This module provides efficient pseudoknot detection and counting,
with incremental updates for delta energy computation.

Key functions:
- crosses(p1, p2): Check if two pairs form a pseudoknot
- crossing_partners(e, M): Find all pairs in M that cross e
- PKCache: Cache for efficient incremental updates

See spec/03_PK_MODULE.md for detailed specification.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "crosses",
    "crossing_partners",
    "count_crossings",
    "PKCache",
    "build_pk_cache",
    "update_pk_cache_on_add",
    "update_pk_cache_on_remove",
    "pairs_to_layers",
]


def crosses(p1: tuple[int, int], p2: tuple[int, int]) -> bool:
    """Check if two base pairs cross (form a pseudoknot).

    Two pairs (i, j) and (k, l) cross if and only if:
        i < k < j < l  OR  k < i < l < j

    This is the standard arc-diagram crossing test.

    Args:
        p1: First pair (i, j)
        p2: Second pair (k, l)

    Returns:
        True if the pairs cross
    """
    i, j = p1
    k, l = p2

    # Normalize so i < j and k < l
    if i > j:
        i, j = j, i
    if k > l:
        k, l = l, k

    return (i < k < j < l) or (k < i < l < j)


def crossing_partners(
    edge: tuple[int, int],
    matching: set[tuple[int, int]],
) -> set[tuple[int, int]]:
    """Find all pairs in matching that cross the given edge.

    Args:
        edge: The (i, j) pair to check
        matching: Current set of base pairs

    Returns:
        Set of pairs that cross the edge
    """
    return {p for p in matching if p != edge and crosses(edge, p)}


def count_crossings(matching: set[tuple[int, int]]) -> int:
    """Count total number of crossings in a matching.

    Each crossing (pair of crossing pairs) is counted once.

    Args:
        matching: Set of base pairs

    Returns:
        Number of unique crossings
    """
    pairs = list(matching)
    n = len(pairs)
    count = 0

    for i in range(n):
        for j in range(i + 1, n):
            if crosses(pairs[i], pairs[j]):
                count += 1

    return count


@dataclass
class PKCache:
    """Cache for efficient pseudoknot computations.

    Attributes:
        pk_pairs: Set of pairs involved in at least one crossing
        cross_counts: Number of crossings per pair
        total_crossings: Total crossing count
        max_crossings: Maximum crossings for any pair
    """

    pk_pairs: set[tuple[int, int]] = field(default_factory=set)
    cross_counts: dict[tuple[int, int], int] = field(default_factory=dict)
    total_crossings: int = 0
    max_crossings: int = 0


def build_pk_cache(matching: set[tuple[int, int]]) -> PKCache:
    """Build a PK cache from scratch for a matching.

    Args:
        matching: Set of base pairs

    Returns:
        Initialized PKCache
    """
    cache = PKCache()
    pairs = list(matching)
    n = len(pairs)

    # Initialize cross counts
    for p in pairs:
        cache.cross_counts[p] = 0

    # Count crossings
    for i in range(n):
        for j in range(i + 1, n):
            if crosses(pairs[i], pairs[j]):
                cache.cross_counts[pairs[i]] += 1
                cache.cross_counts[pairs[j]] += 1
                cache.total_crossings += 1

    # Identify PK pairs and find max
    for p, count in cache.cross_counts.items():
        if count > 0:
            cache.pk_pairs.add(p)
            cache.max_crossings = max(cache.max_crossings, count)

    return cache


def update_pk_cache_on_add(
    edge: tuple[int, int],
    matching: set[tuple[int, int]],
    cache: PKCache,
) -> tuple[int, int, PKCache]:
    """Update PK cache after adding an edge.

    Args:
        edge: The (i, j) pair being added
        matching: Matching BEFORE the addition
        cache: Current cache

    Returns:
        Tuple of (delta_pk_pairs, delta_max, new_cache)
    """
    # Placeholder implementation - will be filled in Phase 1
    raise NotImplementedError("update_pk_cache_on_add not yet implemented")


def update_pk_cache_on_remove(
    edge: tuple[int, int],
    matching: set[tuple[int, int]],
    cache: PKCache,
) -> tuple[int, int, PKCache]:
    """Update PK cache after removing an edge.

    Args:
        edge: The (i, j) pair being removed
        matching: Matching BEFORE the removal (includes edge)
        cache: Current cache

    Returns:
        Tuple of (delta_pk_pairs, delta_max, new_cache)
    """
    # Placeholder implementation - will be filled in Phase 1
    raise NotImplementedError("update_pk_cache_on_remove not yet implemented")


def pairs_to_layers(
    pairs: list[tuple[int, int]] | set[tuple[int, int]],
) -> list[set[tuple[int, int]]]:
    """Partition pairs into non-crossing layers.

    Each layer contains pairs that don't cross each other.
    This is used for converting to Rosetta notation.

    Args:
        pairs: List or set of base pairs

    Returns:
        List of layers, each containing non-crossing pairs
    """
    if isinstance(pairs, set):
        pairs = sorted(pairs)
    else:
        pairs = sorted(pairs)

    layers: list[set[tuple[int, int]]] = []

    for p in pairs:
        placed = False
        for layer in layers:
            if not any(crosses(p, q) for q in layer):
                layer.add(p)
                placed = True
                break
        if not placed:
            layers.append({p})

    return layers
