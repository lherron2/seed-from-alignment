"""
Energy module for MCMC RNA structure sampling.

This module provides the energy function and delta energy computation for
the Metropolis-Hastings sampler. The target distribution is:

    π(M) ∝ exp(-β · E(M))

where the energy function is:

    E(M) = -Σ w_ij + α·PK_pairs(M) + γ·PK_max(M) + λ·lonely_pairs(M) + ...

Key requirements:
- Delta computation must match full recompute (within numerical tolerance)
- delta(add) = -delta(remove) for reversibility
- O(|M|) time complexity for delta computation

See spec/02_ENERGY_MODULE.md for detailed specification.
"""

from __future__ import annotations

from dataclasses import dataclass

from .pk import PKCache, build_pk_cache, update_pk_cache_on_add, update_pk_cache_on_remove

__all__ = [
    "EnergyParams",
    "EnergyCache",
    "compute_energy",
    "compute_delta_energy",
    "normalize_weights",
]


@dataclass
class EnergyParams:
    """Parameters for the energy function.

    Attributes:
        beta: Inverse temperature (default 1.0)
        pk_alpha: Per-crossing penalty coefficient
        pk_gamma: Max-crossing penalty coefficient
        lonely_penalty: Penalty for isolated base pairs
        hairpin_min: Minimum hairpin loop length
        hairpin_penalty: Penalty per hairpin violation (0 = hard constraint)
    """

    beta: float = 1.0
    pk_alpha: float = 2.0
    pk_gamma: float = 0.0
    lonely_penalty: float = 0.0
    hairpin_min: int = 3
    hairpin_penalty: float = 0.0  # 0 = hard constraint (reject), >0 = soft penalty


@dataclass
class EnergyCache:
    """Cached values for efficient delta energy computation.

    Attributes:
        total_energy: Current total energy
        weight_sum: Sum of pair weights in matching
        pk_cache: Cached pseudoknot information
        lonely_count: Number of isolated pairs
        hairpin_violations: Number of hairpin violations
    """

    total_energy: float = 0.0
    weight_sum: float = 0.0
    pk_cache: PKCache | None = None
    lonely_count: int = 0
    hairpin_violations: int = 0

    @property
    def pk_pairs_count(self) -> int:
        """Number of pairs involved in crossings."""
        return len(self.pk_cache.pk_pairs) if self.pk_cache else 0

    @property
    def pk_max(self) -> int:
        """Maximum crossings for any single pair."""
        return self.pk_cache.max_crossings if self.pk_cache else 0


def _count_lonely_pairs(matching: set[tuple[int, int]]) -> int:
    """Count isolated (lonely) pairs in a matching.

    A pair is lonely if it has no stacking partner on either side.

    Args:
        matching: Set of base pairs

    Returns:
        Number of lonely pairs
    """
    count = 0
    for i, j in matching:
        has_stack_inside = (i + 1, j - 1) in matching
        has_stack_outside = (i - 1, j + 1) in matching
        if not has_stack_inside and not has_stack_outside:
            count += 1
    return count


def _is_lonely_after_add(
    edge: tuple[int, int],
    matching: set[tuple[int, int]],
) -> tuple[bool, int]:
    """Check if adding an edge creates/removes lonely pairs.

    Args:
        edge: Edge being added (not yet in matching)
        matching: Current matching (before add)

    Returns:
        Tuple of (is_edge_lonely, delta_lonely_count)
    """
    i, j = edge

    # Check if the new edge is lonely
    has_stack_inside = (i + 1, j - 1) in matching
    has_stack_outside = (i - 1, j + 1) in matching
    edge_is_lonely = not has_stack_inside and not has_stack_outside

    delta = 1 if edge_is_lonely else 0

    # Check if adding this edge fixes any adjacent lonely pairs
    # Inside neighbor at (i+1, j-1) - adding edge creates its outer stack
    inside = (i + 1, j - 1)
    if inside in matching:
        # Was inside lonely before? It had no outer stack (edge not there yet)
        # So it was lonely iff it had no inner stack either
        has_inner_of_inside = (i + 2, j - 2) in matching
        was_lonely = not has_inner_of_inside
        # After adding edge, inside has outer stack, so not lonely anymore
        if was_lonely:
            delta -= 1

    # Outside neighbor at (i-1, j+1) - adding edge creates its inner stack
    outside = (i - 1, j + 1)
    if outside in matching:
        # Was outside lonely before? It had no inner stack (edge not there yet)
        # So it was lonely iff it had no outer stack either
        has_outer_of_outside = (i - 2, j + 2) in matching
        was_lonely = not has_outer_of_outside
        # After adding edge, outside has inner stack, so not lonely anymore
        if was_lonely:
            delta -= 1

    return edge_is_lonely, delta


def _is_lonely_after_remove(
    edge: tuple[int, int],
    matching: set[tuple[int, int]],
) -> tuple[bool, int]:
    """Check if removing an edge creates/removes lonely pairs.

    Args:
        edge: Edge being removed (currently in matching)
        matching: Current matching (before remove, includes edge)

    Returns:
        Tuple of (was_edge_lonely, delta_lonely_count)
    """
    i, j = edge

    # Check if the edge being removed was lonely
    has_stack_inside = (i + 1, j - 1) in matching
    has_stack_outside = (i - 1, j + 1) in matching
    edge_was_lonely = not has_stack_inside and not has_stack_outside

    delta = -1 if edge_was_lonely else 0

    # Check if removing this edge makes any adjacent pairs lonely
    # Inside neighbor at (i+1, j-1) - removing edge removes its outer stack
    inside = (i + 1, j - 1)
    if inside in matching:
        # After removal, inside loses outer stack
        # It becomes lonely if it also has no inner stack
        has_inner_of_inside = (i + 2, j - 2) in matching
        becomes_lonely = not has_inner_of_inside
        if becomes_lonely:
            delta += 1

    # Outside neighbor at (i-1, j+1) - removing edge removes its inner stack
    outside = (i - 1, j + 1)
    if outside in matching:
        # After removal, outside loses inner stack
        # It becomes lonely if it also has no outer stack
        has_outer_of_outside = (i - 2, j + 2) in matching
        becomes_lonely = not has_outer_of_outside
        if becomes_lonely:
            delta += 1

    return edge_was_lonely, delta


def compute_energy(
    matching: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    params: EnergyParams,
    partners: list[int] | None = None,
) -> tuple[float, EnergyCache]:
    """Compute total energy for a matching from scratch.

    Energy formula:
        E(M) = -Σ w_ij + α·|PK_pairs| + γ·PK_max + λ·lonely_count

    Args:
        matching: Set of base pairs (i, j) with i < j
        weights: Weight for each candidate pair
        params: Energy function parameters
        partners: Optional partners array for efficiency

    Returns:
        Tuple of (energy, cache) for delta computation
    """
    # Compute weight sum
    weight_sum = sum(weights.get(p, 0.0) for p in matching)

    # Build PK cache
    pk_cache = build_pk_cache(matching)

    # Count lonely pairs
    lonely_count = _count_lonely_pairs(matching) if params.lonely_penalty > 0 else 0

    # Compute total energy
    energy = -weight_sum
    energy += params.pk_alpha * len(pk_cache.pk_pairs)
    energy += params.pk_gamma * pk_cache.max_crossings
    energy += params.lonely_penalty * lonely_count

    cache = EnergyCache(
        total_energy=energy,
        weight_sum=weight_sum,
        pk_cache=pk_cache,
        lonely_count=lonely_count,
        hairpin_violations=0,  # Hairpins are handled as hard constraints
    )

    return energy, cache


def compute_delta_energy(
    edge: tuple[int, int],
    matching: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    params: EnergyParams,
    cache: EnergyCache,
    is_add: bool,
    partners: list[int] | None = None,
) -> tuple[float, EnergyCache]:
    """Compute delta energy for adding or removing an edge.

    This is the critical function for MCMC efficiency. Must satisfy:
    - delta(add) = -delta(remove) for the same edge (reversibility)
    - Result matches full recompute within tolerance (1e-10)
    - O(|M|) time complexity

    Args:
        edge: The (i, j) pair to add or remove
        matching: Current matching (before the move)
        weights: Weight for each candidate pair
        params: Energy function parameters
        cache: Current energy cache
        is_add: True for addition, False for removal
        partners: Optional partners array for efficiency

    Returns:
        Tuple of (delta_energy, new_cache)
    """
    # Normalize edge
    i, j = edge
    if i > j:
        edge = (j, i)

    # Get edge weight
    edge_weight = weights.get(edge, 0.0)

    if is_add:
        # Adding edge: delta_weight = -w (energy decreases with higher weight)
        delta_weight = -edge_weight
        new_weight_sum = cache.weight_sum + edge_weight

        # Update PK cache
        if cache.pk_cache is None:
            cache.pk_cache = build_pk_cache(matching)
        delta_pk_pairs, delta_pk_max, new_pk_cache = update_pk_cache_on_add(
            edge, matching, cache.pk_cache
        )

        # Update lonely count
        if params.lonely_penalty > 0:
            _, delta_lonely = _is_lonely_after_add(edge, matching)
            new_lonely = cache.lonely_count + delta_lonely
        else:
            delta_lonely = 0
            new_lonely = 0

    else:
        # Removing edge: delta_weight = +w (energy increases)
        delta_weight = edge_weight
        new_weight_sum = cache.weight_sum - edge_weight

        # Update PK cache
        if cache.pk_cache is None:
            cache.pk_cache = build_pk_cache(matching)
        delta_pk_pairs, delta_pk_max, new_pk_cache = update_pk_cache_on_remove(
            edge, matching, cache.pk_cache
        )

        # Update lonely count
        if params.lonely_penalty > 0:
            _, delta_lonely = _is_lonely_after_remove(edge, matching)
            new_lonely = cache.lonely_count + delta_lonely
        else:
            delta_lonely = 0
            new_lonely = 0

    # Compute delta energy
    delta_energy = delta_weight
    delta_energy += params.pk_alpha * delta_pk_pairs
    delta_energy += params.pk_gamma * delta_pk_max
    delta_energy += params.lonely_penalty * delta_lonely

    new_cache = EnergyCache(
        total_energy=cache.total_energy + delta_energy,
        weight_sum=new_weight_sum,
        pk_cache=new_pk_cache,
        lonely_count=new_lonely,
        hairpin_violations=cache.hairpin_violations,
    )

    return delta_energy, new_cache


def normalize_weights(
    weights: dict[tuple[int, int], float],
    method: str = "max",
) -> tuple[dict[tuple[int, int], float], float]:
    """Normalize weights to a standard scale.

    Args:
        weights: Original weights
        method: Normalization method ("max", "mean", "sum", or "none")

    Returns:
        Tuple of (normalized_weights, scale_factor)
    """
    if not weights or method == "none":
        return dict(weights), 1.0

    values = list(weights.values())
    if method == "max":
        scale = max(values) if values else 1.0
    elif method == "mean":
        scale = sum(values) / len(values) if values else 1.0
    elif method == "sum":
        scale = sum(values) if values else 1.0
    else:
        raise ValueError(f"Unknown normalization method: {method}")

    if scale == 0:
        scale = 1.0

    normalized = {k: v / scale for k, v in weights.items()}
    return normalized, scale
