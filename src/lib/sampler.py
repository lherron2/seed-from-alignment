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

from dataclasses import dataclass, field
from typing import Any

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
    ) -> None:
        """Record state at sample time.

        Args:
            energy: Current energy
            matching_size: Number of pairs in matching
        """
        self.energy_trace.append(energy)
        self.sample_sizes.append(matching_size)

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

        # Placeholder - will implement autocorrelation-based ESS
        return float(len(self.energy_trace))

    def summary(self) -> dict[str, Any]:
        """Generate summary statistics.

        Returns:
            Dictionary with summary stats
        """
        return {
            "n_samples": len(self.energy_trace),
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
        }


def sample_matchings_new(
    length: int,
    candidate_pairs: list[tuple[int, int, float]],
    config: SamplerConfig,
    scaffold_pairs: set[tuple[int, int]] | None = None,
    segment_index: SegmentIndex | None = None,
    loop_folds: list[set[tuple[int, int]]] | None = None,
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
    # Placeholder implementation - will be filled in Phase 2
    raise NotImplementedError("sample_matchings_new not yet implemented")


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
    # Placeholder implementation - will be filled in Phase 2
    raise NotImplementedError("enumerate_all_matchings not yet implemented")


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
    # Placeholder implementation - will be filled in Phase 2
    raise NotImplementedError("compute_exact_distribution not yet implemented")
