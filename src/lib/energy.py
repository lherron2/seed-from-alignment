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
        pk_pairs_count: Number of pairs involved in crossings
        pk_max: Maximum crossings for any single pair
        lonely_count: Number of isolated pairs
        hairpin_violations: Number of hairpin violations
    """

    total_energy: float = 0.0
    weight_sum: float = 0.0
    pk_pairs_count: int = 0
    pk_max: int = 0
    lonely_count: int = 0
    hairpin_violations: int = 0


def compute_energy(
    matching: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    params: EnergyParams,
    partners: list[int] | None = None,
) -> tuple[float, EnergyCache]:
    """Compute total energy for a matching from scratch.

    Args:
        matching: Set of base pairs (i, j) with i < j
        weights: Weight for each candidate pair
        params: Energy function parameters
        partners: Optional partners array for efficiency

    Returns:
        Tuple of (energy, cache) for delta computation
    """
    # Placeholder implementation - will be filled in Phase 1
    raise NotImplementedError("compute_energy not yet implemented")


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
    # Placeholder implementation - will be filled in Phase 1
    raise NotImplementedError("compute_delta_energy not yet implemented")


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
