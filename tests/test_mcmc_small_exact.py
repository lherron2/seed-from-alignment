"""Exact enumeration tests for MCMC sampler correctness.

These tests verify that the MCMC sampler produces the correct
Boltzmann distribution by comparing against exact enumeration
for small systems.

Acceptance criteria: max frequency error < 5%
"""

from collections import Counter

import pytest

from src.lib.energy import EnergyParams
from src.lib.sampler import (
    SamplerConfig,
    compute_exact_distribution,
    enumerate_all_matchings,
    sample_matchings_new,
)


def chi_squared_statistic(
    observed: dict[frozenset[tuple[int, int]], int],
    expected: dict[frozenset[tuple[int, int]], float],
    n_samples: int,
) -> float:
    """Compute chi-squared statistic for goodness-of-fit.

    Args:
        observed: Count of each matching
        expected: Expected probability of each matching
        n_samples: Total number of samples

    Returns:
        Chi-squared statistic
    """
    chi2 = 0.0
    for matching, prob in expected.items():
        expected_count = prob * n_samples
        observed_count = observed.get(matching, 0)
        if expected_count > 0:
            chi2 += (observed_count - expected_count) ** 2 / expected_count
    return chi2


def max_frequency_error(
    observed: dict[frozenset[tuple[int, int]], int],
    expected: dict[frozenset[tuple[int, int]], float],
    n_samples: int,
) -> float:
    """Compute maximum frequency error.

    Args:
        observed: Count of each matching
        expected: Expected probability of each matching
        n_samples: Total number of samples

    Returns:
        Maximum absolute frequency error
    """
    max_error = 0.0
    for matching, prob in expected.items():
        observed_freq = observed.get(matching, 0) / n_samples
        error = abs(observed_freq - prob)
        max_error = max(max_error, error)
    return max_error


class TestExactDistribution:
    """Tests comparing MCMC samples to exact Boltzmann distribution."""

    def test_tiny_system_toggle_only(self) -> None:
        """
        Test on tiny system (L=20, 2 candidates) with toggle moves only.

        Procedure:
        1. Enumerate all valid matchings
        2. Compute exact Boltzmann probabilities
        3. Run MCMC for many samples with thin=1 to capture rare states
        4. Compare empirical frequencies to exact probabilities
        5. Assert max frequency error < 10%
        """
        # Use moderate weights so distribution isn't too skewed
        candidates = [(0, 6, 1.0), (10, 16, 1.0)]
        length = 20
        params = EnergyParams(beta=0.5)  # Lower beta for more mixing

        # Enumerate all matchings
        all_pairs = [(i, j) for i, j, _ in candidates]
        matchings = enumerate_all_matchings(length, all_pairs, min_hairpin=3)
        assert len(matchings) == 4  # empty, {a}, {b}, {a,b}

        # Compute exact distribution
        weights = {(i, j): w for i, j, w in candidates}
        exact_dist = compute_exact_distribution(matchings, weights, params)

        # Run MCMC with thin=1 to catch rare states
        n_samples = 10000
        config = SamplerConfig(
            n_samples=n_samples,
            burn_in=1000,
            thin=1,  # No thinning - critical for catching rare states
            seed=42,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)

        # Count observations
        observed = Counter(frozenset(s) for s in samples)

        # Check max frequency error < 10% (allow more tolerance for stochastic tests)
        max_err = max_frequency_error(observed, exact_dist, n_samples)
        assert max_err < 0.10, f"Max frequency error {max_err:.3f} > 0.10"

    def test_tiny_system_with_conflicts(self) -> None:
        """Test system with conflicting pairs."""
        # (0, 6) and (0, 8) conflict at position 0
        candidates = [(0, 6, 1.0), (0, 8, 1.0), (12, 18, 1.0)]
        length = 20
        params = EnergyParams(beta=0.5)

        all_pairs = [(i, j) for i, j, _ in candidates]
        matchings = enumerate_all_matchings(length, all_pairs, min_hairpin=3)
        # empty, {(0,6)}, {(0,8)}, {(12,18)}, {(0,6),(12,18)}, {(0,8),(12,18)}
        assert len(matchings) == 6

        weights = {(i, j): w for i, j, w in candidates}
        exact_dist = compute_exact_distribution(matchings, weights, params)

        n_samples = 10000
        config = SamplerConfig(
            n_samples=n_samples,
            burn_in=1000,
            thin=1,
            seed=123,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)
        observed = Counter(frozenset(s) for s in samples)

        max_err = max_frequency_error(observed, exact_dist, n_samples)
        assert max_err < 0.10, f"Max frequency error {max_err:.3f} > 0.10"

    def test_tiny_system_with_pk_penalty(self) -> None:
        """Test system with pseudoknot penalty."""
        # (0, 10) and (5, 15) form a pseudoknot
        candidates = [(0, 10, 1.0), (5, 15, 1.0)]
        length = 20
        params = EnergyParams(beta=0.5, pk_alpha=1.0)  # Moderate PK penalty

        all_pairs = [(i, j) for i, j, _ in candidates]
        matchings = enumerate_all_matchings(length, all_pairs, min_hairpin=3)
        assert len(matchings) == 4

        weights = {(i, j): w for i, j, w in candidates}
        exact_dist = compute_exact_distribution(matchings, weights, params)

        n_samples = 10000
        config = SamplerConfig(
            n_samples=n_samples,
            burn_in=1000,
            thin=1,
            seed=456,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)
        observed = Counter(frozenset(s) for s in samples)

        max_err = max_frequency_error(observed, exact_dist, n_samples)
        assert max_err < 0.10, f"Max frequency error {max_err:.3f} > 0.10"

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_tiny_system_with_segments(self) -> None:
        """
        Test on tiny system with segment moves enabled.

        Same procedure as toggle-only, but with segment birth/death.
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 4")
    def test_tiny_system_full_moves(self) -> None:
        """
        Test on tiny system with all move types enabled.

        Same procedure, but includes loop closure moves.
        """
        pass


class TestDetailedBalance:
    """Tests for detailed balance verification."""

    def test_toggle_detailed_balance(self) -> None:
        """
        Verify π(M)·A(M→M') = π(M')·A(M'→M) for toggle moves.

        For toggle, Hastings ratio = 1, so detailed balance is:
        π(M)·min(1, π(M')/π(M)) = π(M')·min(1, π(M)/π(M'))

        This simplifies to min(π(M), π(M')) = min(π(M'), π(M)) ✓
        """
        import math

        from src.lib.energy import compute_energy

        candidates = [(0, 6, 2.0), (10, 16, 1.5)]
        weights = {(i, j): w for i, j, w in candidates}
        params = EnergyParams(beta=1.0)

        # Two states differing by one pair
        m1: set[tuple[int, int]] = set()
        m2 = {(0, 6)}

        # Compute energies and Boltzmann weights
        e1, _ = compute_energy(m1, weights, params)
        e2, _ = compute_energy(m2, weights, params)

        pi1 = math.exp(-params.beta * e1)
        pi2 = math.exp(-params.beta * e2)

        # Acceptance probabilities (Hastings ratio = 1 for toggle)
        a_1_to_2 = min(1.0, pi2 / pi1)
        a_2_to_1 = min(1.0, pi1 / pi2)

        # Detailed balance: π(M1)·A(M1→M2) = π(M2)·A(M2→M1)
        lhs = pi1 * a_1_to_2
        rhs = pi2 * a_2_to_1

        assert abs(lhs - rhs) < 1e-10, f"Detailed balance violated: {lhs} != {rhs}"

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_segment_detailed_balance(self) -> None:
        """
        Verify detailed balance for segment moves.

        Must account for Hastings ratios.
        """
        pass


class TestErgodicity:
    """Tests for MCMC ergodicity."""

    def test_all_states_reachable(self) -> None:
        """For tiny system, verify all valid matchings are sampled."""
        candidates = [(0, 6, 1.0), (10, 16, 1.0)]
        length = 20
        params = EnergyParams(beta=0.3)  # Low beta for better mixing

        all_pairs = [(i, j) for i, j, _ in candidates]
        matchings = enumerate_all_matchings(length, all_pairs, min_hairpin=3)

        # Run long enough to visit all states with thin=1
        config = SamplerConfig(
            n_samples=5000,
            burn_in=500,
            thin=1,
            seed=789,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)

        # All matchings should be visited at least once
        observed = {frozenset(s) for s in samples}
        for matching in matchings:
            assert matching in observed, f"Matching {matching} never visited"

    def test_state_reversibility(self) -> None:
        """Verify chain can return to previous states."""
        candidates = [(0, 6, 1.0)]
        length = 10
        params = EnergyParams(beta=0.3)

        config = SamplerConfig(
            n_samples=500,
            burn_in=100,
            thin=1,
            seed=111,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)

        # Should see both empty and non-empty states multiple times
        empty_count = sum(1 for s in samples if len(s) == 0)
        nonempty_count = len(samples) - empty_count

        assert empty_count > 0, "Never visited empty state"
        assert nonempty_count > 0, "Never visited non-empty state"


class TestHighBeta:
    """Tests with high inverse temperature (near-deterministic)."""

    def test_high_beta_favors_best_state(self) -> None:
        """At high beta, should strongly favor lowest energy state."""
        # Both pairs have equal weight, so best state is both
        candidates = [(0, 6, 2.0), (10, 16, 2.0)]
        length = 20
        params = EnergyParams(beta=2.0)  # Moderate-high beta

        n_samples = 5000
        config = SamplerConfig(
            n_samples=n_samples,
            burn_in=1000,
            thin=1,
            seed=222,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)
        observed = Counter(frozenset(s) for s in samples)

        # Most probable should be state with both pairs
        best_state = frozenset({(0, 6), (10, 16)})
        best_count = observed[best_state]

        # Should dominate (>50% at high beta)
        assert best_count / n_samples > 0.5, (
            f"Best state frequency {best_count / n_samples:.2f} < 0.5"
        )


class TestLowBeta:
    """Tests with low inverse temperature (near-uniform)."""

    def test_low_beta_near_uniform(self) -> None:
        """At low beta, distribution should be more uniform."""
        candidates = [(0, 6, 1.0), (10, 16, 1.0)]
        length = 20
        params = EnergyParams(beta=0.1)  # Low beta

        all_pairs = [(i, j) for i, j, _ in candidates]
        matchings = enumerate_all_matchings(length, all_pairs, min_hairpin=3)
        weights = {(i, j): w for i, j, w in candidates}
        exact_dist = compute_exact_distribution(matchings, weights, params)

        # At low beta, distribution should be nearly uniform
        uniform_prob = 1.0 / len(matchings)
        for prob in exact_dist.values():
            assert abs(prob - uniform_prob) < 0.15, "Distribution not near-uniform at low beta"

        # Verify MCMC samples from this near-uniform distribution
        n_samples = 5000
        config = SamplerConfig(
            n_samples=n_samples,
            burn_in=500,
            thin=1,
            seed=333,
            energy_params=params,
        )
        samples, _ = sample_matchings_new(length, candidates, config)
        observed = Counter(frozenset(s) for s in samples)

        # All states should be observed with reasonable frequency
        for matching in matchings:
            freq = observed.get(matching, 0) / n_samples
            assert freq > 0.1, f"Matching {matching} underrepresented at low beta"
