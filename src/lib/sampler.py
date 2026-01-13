"""
MCMC sampling loop and diagnostics for RNA structure sampling.

This module provides the main sampling interface and diagnostic tools
for monitoring convergence and mixing.

Key components:
- SamplerConfig: Configuration for the sampler
- SamplerDiagnostics: Tracking acceptance rates, energy traces, ESS
- sample_matchings: Main sampling function

See spec/06_SAMPLER_MODULE.md for detailed specification.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

from .candidate_index import CandidateEdgeIndex
from .energy import EnergyParams
from .moves import MoveType
from .segments import SegmentIndex

__all__ = [
    "SamplerConfig",
    "SamplerDiagnostics",
    "sample_matchings_new",
    "enumerate_all_matchings",
    "compute_exact_distribution",
]


@dataclass
class SamplerConfig:
    """Configuration for the MCMC sampler.

    Attributes:
        n_samples: Number of samples to collect
        burn_in: Number of burn-in steps
        thin: Thinning interval (keep 1 in every N)
        seed: Random seed for reproducibility
        energy_params: Parameters for energy function
        move_probs: Probability for each move type
        min_hairpin: Minimum hairpin loop length
        validate_structures: Whether to validate each sample
    """

    n_samples: int = 10
    burn_in: int = 2000
    thin: int = 20
    seed: int | None = None
    energy_params: EnergyParams = field(default_factory=EnergyParams)
    move_probs: dict[MoveType, float] | None = None
    min_hairpin: int = 3
    validate_structures: bool = False
    delayed_accept: bool = True


@dataclass
class SamplerDiagnostics:
    """Diagnostics for monitoring MCMC convergence and mixing.

    Attributes:
        total_proposals: Total number of proposals by move type
        accepted: Number of accepted proposals by move type
        energy_trace: List of energies at each sample
        sample_sizes: Number of pairs in each sample
    """

    total_proposals: dict[MoveType, int] = field(default_factory=dict)
    accepted: dict[MoveType, int] = field(default_factory=dict)
    energy_trace: list[float] = field(default_factory=list)
    sample_sizes: list[int] = field(default_factory=list)
    pk_pairs_count: list[int] = field(default_factory=list)
    pk_max_crossings: list[int] = field(default_factory=list)
    wall_seconds: float = 0.0
    total_steps: int = 0

    def record_move(
        self,
        move_type: MoveType,
        accepted: bool,
    ) -> None:
        """Record a move proposal result.

        Args:
            move_type: Type of move proposed
            accepted: Whether the move was accepted
        """
        self.total_proposals[move_type] = self.total_proposals.get(move_type, 0) + 1
        if accepted:
            self.accepted[move_type] = self.accepted.get(move_type, 0) + 1

    def record_state(
        self,
        energy: float,
        matching_size: int,
        pk_pairs_count: int = 0,
        pk_max_crossings: int = 0,
    ) -> None:
        """Record state at sample time.

        Args:
            energy: Current energy
            matching_size: Number of pairs in matching
        """
        self.energy_trace.append(energy)
        self.sample_sizes.append(matching_size)
        self.pk_pairs_count.append(pk_pairs_count)
        self.pk_max_crossings.append(pk_max_crossings)

    def acceptance_rate(self, move_type: MoveType) -> float:
        """Get acceptance rate for a move type.

        Args:
            move_type: Type of move

        Returns:
            Acceptance rate (0-1)
        """
        total = self.total_proposals.get(move_type, 0)
        if total == 0:
            return 0.0
        return self.accepted.get(move_type, 0) / total

    def estimate_ess(self) -> float:
        """Estimate effective sample size from energy trace.

        Uses the autocorrelation-based ESS estimator.

        Returns:
            Estimated effective sample size
        """
        if len(self.energy_trace) < 2:
            return float(len(self.energy_trace))

        # Autocorrelation-based ESS estimator on the energy trace.
        # For short chains this is noisy; it's mainly a regression signal.
        x = self.energy_trace
        n = len(x)
        if n < 3:
            return float(n)

        mean = sum(x) / n
        var = sum((v - mean) ** 2 for v in x) / (n - 1)
        if var <= 0:
            return float(n)

        def autocov(lag: int) -> float:
            return sum((x[t] - mean) * (x[t + lag] - mean) for t in range(n - lag)) / (n - 1)

        tau = 1.0
        # Stop when correlation becomes negative (initial positive sequence).
        for lag in range(1, min(n // 2, 2000)):
            rho = autocov(lag) / var
            if rho <= 0:
                break
            tau += 2.0 * rho

        ess = n / max(1.0, tau)
        # Guard against numerical weirdness.
        return float(max(1.0, min(n, ess)))

    def summary(self) -> dict[str, Any]:
        """Generate summary statistics.

        Returns:
            Dictionary with summary stats
        """
        return {
            "n_samples": len(self.energy_trace),
            "total_steps": self.total_steps,
            "wall_seconds": self.wall_seconds,
            "acceptance_rates": {
                mt.name: self.acceptance_rate(mt)
                for mt in MoveType
                if self.total_proposals.get(mt, 0) > 0
            },
            "ess": self.estimate_ess(),
            "mean_energy": sum(self.energy_trace) / len(self.energy_trace)
            if self.energy_trace
            else 0.0,
            "mean_size": sum(self.sample_sizes) / len(self.sample_sizes)
            if self.sample_sizes
            else 0.0,
            "mean_pk_pairs": sum(self.pk_pairs_count) / len(self.pk_pairs_count)
            if self.pk_pairs_count
            else 0.0,
            "mean_pk_max_crossings": sum(self.pk_max_crossings) / len(self.pk_max_crossings)
            if self.pk_max_crossings
            else 0.0,
        }


def sample_matchings_new(
    length: int,
    candidate_pairs: list[tuple[int, int, float]],
    config: SamplerConfig,
    scaffold_pairs: set[tuple[int, int]] | None = None,
    initial_pairs: set[tuple[int, int]] | None = None,
    segment_index: SegmentIndex | None = None,
    loop_folds: list[set[tuple[int, int]]] | None = None,
    weights: dict[tuple[int, int], float] | None = None,
) -> tuple[list[set[tuple[int, int]]], SamplerDiagnostics]:
    """Sample matchings using MCMC with correct detailed balance.

    This is the new implementation with proper energy delta computation
    and multiple move types.

    Args:
        length: Sequence length
        candidate_pairs: List of (i, j, weight) candidates
        config: Sampler configuration
        scaffold_pairs: Fixed pairs that must be in every sample
        segment_index: Precomputed segment index (optional)
        loop_folds: Precomputed loop foldings (optional)

    Returns:
        Tuple of (samples, diagnostics)
    """
    import math
    import random

    from .energy import compute_delta_energy, compute_delta_energy_no_pk, compute_energy
    from .moves import select_and_propose_move
    from .sampler_state import create_initial_state

    if scaffold_pairs is None:
        scaffold_pairs = set()

    # Build weight dictionary (can be passed in to avoid repeated work across chains)
    if weights is None:
        weights = {}
        for i, j, w in candidate_pairs:
            if i > j:
                i, j = j, i
            weights[(i, j)] = w

    # Initialize RNG
    rng = random.Random(config.seed)

    # Create initial state
    state = create_initial_state(
        length=length,
        fixed_pairs=scaffold_pairs,
        initial_pairs=initial_pairs,
        weights=weights,
        params=config.energy_params,
    )

    # Compute initial energy
    energy, cache = compute_energy(state.pair_set, weights, config.energy_params, state.partners)
    state.energy_cache = cache
    assert cache.pk_cache is not None  # compute_energy always sets pk_cache
    state.pk_cache = cache.pk_cache

    # Diagnostics
    diagnostics = SamplerDiagnostics()
    candidate_index = CandidateEdgeIndex(
        length=length,
        candidate_pairs=candidate_pairs,
        partners=state.partners,
        fixed_pairs=state.fixed_pairs,
        min_hairpin=config.min_hairpin,
    )

    # Calculate total steps
    total_steps = config.burn_in + config.n_samples * config.thin
    samples: list[set[tuple[int, int]]] = []
    diagnostics.total_steps = total_steps
    t0 = time.perf_counter()

    for step in range(total_steps):
        # Use select_and_propose_move to select from available move types
        move = select_and_propose_move(
            matching=state.pair_set,
            partners=state.partners,
            candidate_pairs=candidate_pairs,
            segment_index=segment_index,
            loop_folds=loop_folds,
            fixed_pairs=state.fixed_pairs,
            weights=weights,
            rng=rng,
            candidate_index=candidate_index,
            move_probs=config.move_probs,
            min_hairpin=config.min_hairpin,
        )

        accepted = False

        if move.is_valid:
            # Optional structural validation on the proposed final matching.
            if config.validate_structures:
                from .validate_structure import pairs_to_structure, validate_structure

                dummy_seq = "A" * length
                temp_matching_final = set(state.pair_set)
                if move.pairs_removed:
                    for i, j in move.pairs_removed:
                        if i > j:
                            i, j = j, i
                        temp_matching_final.discard((i, j))
                if move.pairs_added:
                    for i, j in move.pairs_added:
                        if i > j:
                            i, j = j, i
                        temp_matching_final.add((i, j))

                struct = pairs_to_structure(sorted(temp_matching_final), length)
                is_valid, _ = validate_structure(
                    dummy_seq,
                    struct,
                    allow_pk=True,
                    min_hairpin=config.min_hairpin,
                    allow_lonely_pairs=True,
                    canonical_only=False,
                )
                if not is_valid:
                    move.is_valid = False

            if move.is_valid:
                if config.delayed_accept:
                    # Stage 1: cheap energy (no PK terms)
                    delta_fast = 0.0
                    fast_cache = state.energy_cache
                    temp_matching = set(state.pair_set)

                    if move.pairs_removed:
                        for i, j in move.pairs_removed:
                            if i > j:
                                i, j = j, i
                            edge_delta_fast, fast_cache = compute_delta_energy_no_pk(
                                (i, j),
                                temp_matching,
                                weights,
                                config.energy_params,
                                fast_cache,
                                is_add=False,
                            )
                            delta_fast += edge_delta_fast
                            temp_matching.discard((i, j))

                    if move.pairs_added:
                        for i, j in move.pairs_added:
                            if i > j:
                                i, j = j, i
                            edge_delta_fast, fast_cache = compute_delta_energy_no_pk(
                                (i, j),
                                temp_matching,
                                weights,
                                config.energy_params,
                                fast_cache,
                                is_add=True,
                            )
                            delta_fast += edge_delta_fast
                            temp_matching.add((i, j))

                    log_accept1 = -config.energy_params.beta * delta_fast
                    log_accept1 += (
                        math.log(move.hastings_ratio) if move.hastings_ratio > 0 else float("-inf")
                    )

                    if log_accept1 >= 0 or rng.random() < math.exp(log_accept1):
                        # Stage 2: compute full delta only for stage-1-accepted proposals.
                        delta_e = 0.0
                        current_cache = state.energy_cache
                        temp_matching_full = set(state.pair_set)

                        if move.pairs_removed:
                            for i, j in move.pairs_removed:
                                if i > j:
                                    i, j = j, i
                                edge_delta, current_cache = compute_delta_energy(
                                    (i, j),
                                    temp_matching_full,
                                    weights,
                                    config.energy_params,
                                    current_cache,
                                    is_add=False,
                                )
                                delta_e += edge_delta
                                temp_matching_full.discard((i, j))

                        if move.pairs_added:
                            for i, j in move.pairs_added:
                                if i > j:
                                    i, j = j, i
                                edge_delta, current_cache = compute_delta_energy(
                                    (i, j),
                                    temp_matching_full,
                                    weights,
                                    config.energy_params,
                                    current_cache,
                                    is_add=True,
                                )
                                delta_e += edge_delta
                                temp_matching_full.add((i, j))

                        correction = delta_e - delta_fast
                        log_accept2 = -config.energy_params.beta * correction
                        if log_accept2 >= 0 or rng.random() < math.exp(log_accept2):
                            accepted = True
                else:
                    # Standard Metropolis-Hastings acceptance
                    delta_e = 0.0
                    current_cache = state.energy_cache
                    temp_matching = set(state.pair_set)

                    if move.pairs_removed:
                        for i, j in move.pairs_removed:
                            if i > j:
                                i, j = j, i
                            edge_delta, current_cache = compute_delta_energy(
                                (i, j),
                                temp_matching,
                                weights,
                                config.energy_params,
                                current_cache,
                                is_add=False,
                            )
                            delta_e += edge_delta
                            temp_matching.discard((i, j))

                    if move.pairs_added:
                        for i, j in move.pairs_added:
                            if i > j:
                                i, j = j, i
                            edge_delta, current_cache = compute_delta_energy(
                                (i, j),
                                temp_matching,
                                weights,
                                config.energy_params,
                                current_cache,
                                is_add=True,
                            )
                            delta_e += edge_delta
                            temp_matching.add((i, j))

                    log_accept = -config.energy_params.beta * delta_e
                    log_accept += (
                        math.log(move.hastings_ratio) if move.hastings_ratio > 0 else float("-inf")
                    )
                    if log_accept >= 0 or rng.random() < math.exp(log_accept):
                        accepted = True

                if accepted:
                    # Apply the move
                    if move.pairs_removed:
                        for i, j in move.pairs_removed:
                            state.remove_pair(i, j)
                            candidate_index.apply_remove_pair(i, j, state.partners)
                    if move.pairs_added:
                        for i, j in move.pairs_added:
                            state.add_pair(i, j)
                            candidate_index.apply_add_pair(i, j, state.partners)
                    state.energy_cache = current_cache
                    assert current_cache.pk_cache is not None
                    state.pk_cache = current_cache.pk_cache
                    energy = current_cache.total_energy

        diagnostics.record_move(move.move_type, accepted)

        # Collect sample after burn-in at thinning interval
        if step >= config.burn_in and (step - config.burn_in) % config.thin == 0:
            diagnostics.record_state(
                energy,
                len(state.pair_set),
                pk_pairs_count=len(state.pk_cache.pk_pairs),
                pk_max_crossings=state.pk_cache.max_crossings,
            )
            samples.append(set(state.pair_set))

    diagnostics.wall_seconds = time.perf_counter() - t0
    return samples, diagnostics


def enumerate_all_matchings(
    length: int,
    candidate_pairs: list[tuple[int, int]],
    fixed_pairs: set[tuple[int, int]] | None = None,
    min_hairpin: int = 3,
) -> list[frozenset[tuple[int, int]]]:
    """Enumerate all valid matchings for exact testing.

    Only practical for small systems (L < 15, few candidates).

    Args:
        length: Sequence length
        candidate_pairs: List of candidate (i, j) pairs
        fixed_pairs: Pairs that must be in every matching
        min_hairpin: Minimum hairpin loop length

    Returns:
        List of all valid matchings as frozen sets
    """
    if fixed_pairs is None:
        fixed_pairs = set()

    # Normalize pairs
    normalized_candidates = []
    for i, j in candidate_pairs:
        if i > j:
            i, j = j, i
        # Check hairpin constraint
        if j - i - 1 >= min_hairpin:
            normalized_candidates.append((i, j))

    normalized_fixed = set()
    for i, j in fixed_pairs:
        if i > j:
            i, j = j, i
        normalized_fixed.add((i, j))

    # Remove fixed pairs from candidates (they're always included)
    optional_candidates = [p for p in normalized_candidates if p not in normalized_fixed]

    # Check fixed pairs don't conflict
    fixed_positions: set[int] = set()
    for i, j in normalized_fixed:
        if i in fixed_positions or j in fixed_positions:
            # Fixed pairs conflict - no valid matchings
            return []
        fixed_positions.add(i)
        fixed_positions.add(j)

    # Enumerate all subsets of optional candidates
    results: list[frozenset[tuple[int, int]]] = []

    def is_valid_matching(pairs: set[tuple[int, int]]) -> bool:
        """Check if pairs form a valid matching."""
        positions: set[int] = set(fixed_positions)
        for i, j in pairs:
            if i in positions or j in positions:
                return False
            positions.add(i)
            positions.add(j)
        return True

    n = len(optional_candidates)
    for mask in range(1 << n):
        subset = {optional_candidates[k] for k in range(n) if mask & (1 << k)}
        if is_valid_matching(subset):
            full_matching = normalized_fixed | subset
            results.append(frozenset(full_matching))

    return results


def compute_exact_distribution(
    matchings: list[frozenset[tuple[int, int]]],
    weights: dict[tuple[int, int], float],
    params: EnergyParams,
) -> dict[frozenset[tuple[int, int]], float]:
    """Compute exact Boltzmann distribution for testing.

    Args:
        matchings: List of all valid matchings
        weights: Weight for each candidate pair
        params: Energy function parameters

    Returns:
        Dictionary mapping matching -> probability
    """
    import math

    from .energy import compute_energy

    if not matchings:
        return {}

    # Compute energy for each matching
    energies: dict[frozenset[tuple[int, int]], float] = {}
    for matching in matchings:
        energy, _ = compute_energy(set(matching), weights, params)
        energies[matching] = energy

    # Compute Boltzmann weights: exp(-beta * E)
    # Use log-sum-exp trick for numerical stability
    min_energy = min(energies.values())
    boltzmann_weights: dict[frozenset[tuple[int, int]], float] = {}
    for matching, energy in energies.items():
        boltzmann_weights[matching] = math.exp(-params.beta * (energy - min_energy))

    # Normalize
    total = sum(boltzmann_weights.values())
    return {m: w / total for m, w in boltzmann_weights.items()}
