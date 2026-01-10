"""
Sampler state management for MCMC RNA structure sampling.

This module provides state management with caches for efficient
delta computation. The state must maintain consistency between:
- pair_set: Set of current pairs
- partners: Array mapping positions to partners
- energy_cache: Cached energy components
- pk_cache: Cached pseudoknot information

See spec/06_SAMPLER_MODULE.md for detailed specification.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field

from .energy import EnergyCache, EnergyParams
from .pk import PKCache

__all__ = [
    "SamplerState",
    "create_initial_state",
]


@dataclass
class SamplerState:
    """Complete state of the MCMC sampler.

    Invariants that must always hold:
    - partners[i] == j iff partners[j] == i iff (min(i,j), max(i,j)) in pair_set
    - partners[i] == -1 iff position i is unpaired
    - energy_cache reflects current pair_set and pk_cache
    - pk_cache reflects current pair_set

    Attributes:
        length: Sequence length
        pair_set: Current set of base pairs
        partners: Array mapping position to partner (-1 if unpaired)
        energy_cache: Cached energy components
        pk_cache: Cached pseudoknot information
        fixed_pairs: Pairs that cannot be modified (scaffold)
    """

    length: int
    pair_set: set[tuple[int, int]]
    partners: list[int]
    energy_cache: EnergyCache
    pk_cache: PKCache
    fixed_pairs: set[tuple[int, int]] = field(default_factory=set)

    def copy(self) -> SamplerState:
        """Create a deep copy of the state.

        Returns:
            New SamplerState with copied data
        """
        return SamplerState(
            length=self.length,
            pair_set=set(self.pair_set),
            partners=list(self.partners),
            energy_cache=deepcopy(self.energy_cache),
            pk_cache=deepcopy(self.pk_cache),
            fixed_pairs=set(self.fixed_pairs),
        )

    def add_pair(self, i: int, j: int) -> None:
        """Add a pair to the state.

        Updates pair_set and partners array. Does not update caches.

        Args:
            i: First position
            j: Second position
        """
        if i > j:
            i, j = j, i

        self.pair_set.add((i, j))
        self.partners[i] = j
        self.partners[j] = i

    def remove_pair(self, i: int, j: int) -> None:
        """Remove a pair from the state.

        Updates pair_set and partners array. Does not update caches.

        Args:
            i: First position
            j: Second position
        """
        if i > j:
            i, j = j, i

        self.pair_set.discard((i, j))
        self.partners[i] = -1
        self.partners[j] = -1

    def is_position_free(self, pos: int) -> bool:
        """Check if a position is unpaired.

        Args:
            pos: Position to check

        Returns:
            True if unpaired
        """
        return self.partners[pos] == -1

    def validate_consistency(self) -> list[str]:
        """Check state consistency invariants.

        Returns:
            List of error messages (empty if consistent)
        """
        errors = []

        # Check partners array matches pair_set
        for i, j in self.pair_set:
            if i > j:
                errors.append(f"Pair ({i}, {j}) not normalized (i > j)")
            if self.partners[i] != j:
                errors.append(f"partners[{i}] = {self.partners[i]}, expected {j}")
            if self.partners[j] != i:
                errors.append(f"partners[{j}] = {self.partners[j]}, expected {i}")

        # Check that paired positions have valid partners
        for pos in range(self.length):
            partner = self.partners[pos]
            if partner != -1:
                expected_pair = (min(pos, partner), max(pos, partner))
                if expected_pair not in self.pair_set:
                    errors.append(
                        f"partners[{pos}] = {partner} but {expected_pair} not in pair_set"
                    )

        return errors


def create_initial_state(
    length: int,
    fixed_pairs: set[tuple[int, int]] | None = None,
    weights: dict[tuple[int, int], float] | None = None,
    params: EnergyParams | None = None,
) -> SamplerState:
    """Create initial sampler state.

    Args:
        length: Sequence length
        fixed_pairs: Pairs that must be in every sample (scaffold)
        weights: Weight for each candidate pair
        params: Energy function parameters

    Returns:
        Initialized SamplerState
    """
    if fixed_pairs is None:
        fixed_pairs = set()
    if params is None:
        params = EnergyParams()

    # Initialize partners array
    partners = [-1] * length

    # Normalize fixed pairs and set partners
    normalized_fixed = set()
    for i, j in fixed_pairs:
        if i > j:
            i, j = j, i
        normalized_fixed.add((i, j))
        partners[i] = j
        partners[j] = i

    # Build initial caches
    from .pk import build_pk_cache

    pk_cache = build_pk_cache(normalized_fixed)

    # Compute initial energy cache
    # For now, use default cache - will be updated when energy module is complete
    energy_cache = EnergyCache()

    return SamplerState(
        length=length,
        pair_set=set(normalized_fixed),  # Copy to avoid aliasing
        partners=partners,
        energy_cache=energy_cache,
        pk_cache=pk_cache,
        fixed_pairs=normalized_fixed,
    )
