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
    # Normalize edge
    i, j = edge
    if i > j:
        edge = (j, i)

    old_pk_pairs_count = len(cache.pk_pairs)
    old_max = cache.max_crossings

    # Create new cache with copied data
    new_cross_counts = dict(cache.cross_counts)
    new_pk_pairs = set(cache.pk_pairs)
    new_total = cache.total_crossings
    new_max = cache.max_crossings

    # Find all pairs that cross the new edge
    crossing_pairs = crossing_partners(edge, matching)
    new_edge_crossings = len(crossing_pairs)

    # Update cross counts for existing pairs
    for p in crossing_pairs:
        old_count = new_cross_counts.get(p, 0)
        new_count = old_count + 1
        new_cross_counts[p] = new_count

        # If this pair just became a PK pair, add it
        if old_count == 0:
            new_pk_pairs.add(p)

        # Update max
        if new_count > new_max:
            new_max = new_count

    # Add the new edge to cross counts
    new_cross_counts[edge] = new_edge_crossings
    if new_edge_crossings > 0:
        new_pk_pairs.add(edge)
        if new_edge_crossings > new_max:
            new_max = new_edge_crossings

    # Update total crossings
    new_total += new_edge_crossings

    new_cache = PKCache(
        pk_pairs=new_pk_pairs,
        cross_counts=new_cross_counts,
        total_crossings=new_total,
        max_crossings=new_max,
    )

    delta_pk_pairs = len(new_pk_pairs) - old_pk_pairs_count
    delta_max = new_max - old_max

    return delta_pk_pairs, delta_max, new_cache


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
    # Normalize edge
    i, j = edge
    if i > j:
        edge = (j, i)

    old_pk_pairs_count = len(cache.pk_pairs)
    old_max = cache.max_crossings

    # Create new cache with copied data
    new_cross_counts = dict(cache.cross_counts)
    new_pk_pairs = set(cache.pk_pairs)
    new_total = cache.total_crossings

    # Get the crossing count of the edge being removed
    edge_crossings = new_cross_counts.get(edge, 0)

    # Find all pairs that crossed the removed edge (in the remaining matching)
    # These are pairs in matching - {edge} that cross edge
    remaining_matching = matching - {edge}
    crossing_pairs = crossing_partners(edge, remaining_matching)

    # Update cross counts for the pairs that crossed the removed edge
    for p in crossing_pairs:
        old_count = new_cross_counts.get(p, 0)
        new_count = max(0, old_count - 1)
        new_cross_counts[p] = new_count

        # If this pair is no longer a PK pair, remove it
        if new_count == 0:
            new_pk_pairs.discard(p)

    # Remove the edge from cross counts
    if edge in new_cross_counts:
        del new_cross_counts[edge]
    new_pk_pairs.discard(edge)

    # Update total crossings
    new_total = max(0, new_total - edge_crossings)

    # Recompute max (we need to do this since we might have decreased the max)
    new_max = 0
    for count in new_cross_counts.values():
        if count > new_max:
            new_max = count

    new_cache = PKCache(
        pk_pairs=new_pk_pairs,
        cross_counts=new_cross_counts,
        total_crossings=new_total,
        max_crossings=new_max,
    )

    delta_pk_pairs = len(new_pk_pairs) - old_pk_pairs_count
    delta_max = new_max - old_max

    return delta_pk_pairs, delta_max, new_cache


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
