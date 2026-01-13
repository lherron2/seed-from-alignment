"""Tests for energy computation module."""

import pytest

from src.lib.energy import (
    EnergyCache,
    EnergyParams,
    compute_delta_energy,
    compute_energy,
    normalize_weights,
)


class TestEnergyParams:
    """Tests for EnergyParams dataclass."""

    def test_default_values(self) -> None:
        """Default parameter values are sensible."""
        params = EnergyParams()
        assert params.beta == 1.0
        assert params.pk_alpha >= 0
        assert params.hairpin_min >= 0

    def test_custom_values(self) -> None:
        """Custom parameter values are stored correctly."""
        params = EnergyParams(
            beta=0.5,
            pk_alpha=3.0,
            pk_gamma=1.0,
            lonely_penalty=2.0,
        )
        assert params.beta == 0.5
        assert params.pk_alpha == 3.0
        assert params.pk_gamma == 1.0
        assert params.lonely_penalty == 2.0


class TestEnergyCache:
    """Tests for EnergyCache dataclass."""

    def test_default_values(self) -> None:
        """Default cache values are zero."""
        cache = EnergyCache()
        assert cache.total_energy == 0.0
        assert cache.weight_sum == 0.0
        assert cache.pk_pairs_count == 0

    def test_pk_max_property(self) -> None:
        """pk_max property works correctly."""
        cache = EnergyCache()
        assert cache.pk_max == 0


class TestNormalizeWeights:
    """Tests for normalize_weights() function."""

    def test_empty_weights(self) -> None:
        """Empty weights return empty dict."""
        normalized, scale = normalize_weights({})
        assert normalized == {}
        assert scale == 1.0

    def test_max_normalization(self) -> None:
        """Max normalization divides by max value."""
        weights = {(0, 5): 2.0, (1, 4): 4.0}
        normalized, scale = normalize_weights(weights, method="max")
        assert scale == 4.0
        assert normalized[(0, 5)] == 0.5
        assert normalized[(1, 4)] == 1.0

    def test_mean_normalization(self) -> None:
        """Mean normalization divides by mean value."""
        weights = {(0, 5): 2.0, (1, 4): 4.0}
        normalized, scale = normalize_weights(weights, method="mean")
        assert scale == 3.0

    def test_sum_normalization(self) -> None:
        """Sum normalization divides by total sum."""
        weights = {(0, 5): 2.0, (1, 4): 4.0}
        normalized, scale = normalize_weights(weights, method="sum")
        assert scale == 6.0

    def test_no_normalization(self) -> None:
        """None method returns unchanged weights."""
        weights = {(0, 5): 2.0, (1, 4): 4.0}
        normalized, scale = normalize_weights(weights, method="none")
        assert scale == 1.0
        assert normalized == weights


class TestComputeEnergy:
    """Tests for compute_energy()."""

    def test_empty_matching(self) -> None:
        """Empty matching has zero energy."""
        matching: set[tuple[int, int]] = set()
        weights: dict[tuple[int, int], float] = {}
        params = EnergyParams()

        energy, cache = compute_energy(matching, weights, params)

        assert energy == 0.0
        assert cache.total_energy == 0.0
        assert cache.weight_sum == 0.0
        assert cache.pk_pairs_count == 0

    def test_single_pair_negative_weight(self) -> None:
        """Single pair energy is negative weight (reward)."""
        matching = {(0, 10)}
        weights = {(0, 10): 5.0}
        params = EnergyParams()

        energy, cache = compute_energy(matching, weights, params)

        assert energy == -5.0
        assert cache.weight_sum == 5.0
        assert cache.pk_pairs_count == 0

    def test_multiple_nested_pairs(self) -> None:
        """Multiple nested pairs sum weights negatively."""
        matching = {(0, 10), (2, 8), (3, 7)}
        weights = {(0, 10): 1.0, (2, 8): 2.0, (3, 7): 3.0}
        params = EnergyParams()

        energy, cache = compute_energy(matching, weights, params)

        assert energy == -6.0  # -(1 + 2 + 3)
        assert cache.weight_sum == 6.0
        assert cache.pk_pairs_count == 0  # No crossings

    def test_pk_penalty(self) -> None:
        """Pseudoknot pairs incur penalty."""
        # (0, 5) and (2, 8) cross: 0 < 2 < 5 < 8
        matching = {(0, 5), (2, 8)}
        weights = {(0, 5): 1.0, (2, 8): 1.0}
        params = EnergyParams(pk_alpha=2.0, pk_gamma=0.0)

        energy, cache = compute_energy(matching, weights, params)

        # Energy = -weights + pk_alpha * pk_pairs_count
        # = -2.0 + 2.0 * 2 = 2.0
        assert cache.pk_pairs_count == 2
        assert energy == -2.0 + 2.0 * 2  # = 2.0

    def test_pk_max_penalty(self) -> None:
        """PK max crossing penalty is applied."""
        # (0, 5) and (2, 8) cross once each
        matching = {(0, 5), (2, 8)}
        weights = {(0, 5): 1.0, (2, 8): 1.0}
        params = EnergyParams(pk_alpha=0.0, pk_gamma=3.0)

        energy, cache = compute_energy(matching, weights, params)

        assert cache.pk_max == 1
        assert energy == -2.0 + 3.0 * 1  # = 1.0

    def test_lonely_pair_penalty(self) -> None:
        """Lonely pairs incur penalty."""
        # Single isolated pair (no stacking partners)
        matching = {(0, 10)}
        weights = {(0, 10): 1.0}
        params = EnergyParams(lonely_penalty=5.0)

        energy, cache = compute_energy(matching, weights, params)

        assert cache.lonely_count == 1
        assert energy == -1.0 + 5.0 * 1  # = 4.0

    def test_stacked_pairs_not_lonely(self) -> None:
        """Stacked pairs are not lonely."""
        # (0, 10) and (1, 9) are stacking partners
        matching = {(0, 10), (1, 9)}
        weights = {(0, 10): 1.0, (1, 9): 1.0}
        params = EnergyParams(lonely_penalty=5.0)

        energy, cache = compute_energy(matching, weights, params)

        assert cache.lonely_count == 0
        assert energy == -2.0  # Just the weights


class TestComputeDeltaEnergy:
    """Tests for compute_delta_energy()."""

    TOLERANCE = 1e-10

    def test_add_single_pair(self) -> None:
        """Adding a pair decreases energy by weight."""
        matching: set[tuple[int, int]] = set()
        weights = {(0, 10): 5.0}
        params = EnergyParams()

        _, cache = compute_energy(matching, weights, params)
        delta, new_cache = compute_delta_energy(
            (0, 10), matching, weights, params, cache, is_add=True
        )

        assert delta == -5.0
        assert new_cache.total_energy == -5.0

    def test_remove_single_pair(self) -> None:
        """Removing a pair increases energy by weight."""
        matching = {(0, 10)}
        weights = {(0, 10): 5.0}
        params = EnergyParams()

        _, cache = compute_energy(matching, weights, params)
        delta, new_cache = compute_delta_energy(
            (0, 10), matching, weights, params, cache, is_add=False
        )

        assert delta == 5.0
        assert new_cache.total_energy == 0.0

    def test_delta_matches_recompute_add(self) -> None:
        """Delta add matches full recompute within tolerance."""
        matching = {(0, 10), (2, 8)}
        weights = {(0, 10): 1.0, (2, 8): 2.0, (3, 7): 3.0}
        params = EnergyParams(pk_alpha=2.0, pk_gamma=1.0, lonely_penalty=0.5)

        _, cache = compute_energy(matching, weights, params)

        # Add a new pair
        edge = (3, 7)
        delta, new_cache = compute_delta_energy(edge, matching, weights, params, cache, is_add=True)

        # Compute full energy after add
        new_matching = matching | {edge}
        expected_energy, _ = compute_energy(new_matching, weights, params)

        assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

    def test_delta_matches_recompute_remove(self) -> None:
        """Delta remove matches full recompute within tolerance."""
        matching = {(0, 10), (2, 8), (3, 7)}
        weights = {(0, 10): 1.0, (2, 8): 2.0, (3, 7): 3.0}
        params = EnergyParams(pk_alpha=2.0, pk_gamma=1.0, lonely_penalty=0.5)

        _, cache = compute_energy(matching, weights, params)

        # Remove a pair
        edge = (3, 7)
        delta, new_cache = compute_delta_energy(
            edge, matching, weights, params, cache, is_add=False
        )

        # Compute full energy after remove
        new_matching = matching - {edge}
        expected_energy, _ = compute_energy(new_matching, weights, params)

        assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

    def test_delta_reversibility(self) -> None:
        """delta(add) = -delta(remove) for same edge."""
        matching = {(0, 10), (2, 8)}
        weights = {(0, 10): 1.0, (2, 8): 2.0, (5, 15): 3.0}
        params = EnergyParams(pk_alpha=2.0, pk_gamma=1.0, lonely_penalty=0.5)

        _, cache = compute_energy(matching, weights, params)

        # Compute delta for add
        edge = (5, 15)
        delta_add, cache_after_add = compute_delta_energy(
            edge, matching, weights, params, cache, is_add=True
        )

        # Compute delta for remove from the new state
        new_matching = matching | {edge}
        delta_remove, _ = compute_delta_energy(
            edge, new_matching, weights, params, cache_after_add, is_add=False
        )

        assert abs(delta_add + delta_remove) < self.TOLERANCE

    def test_delta_pk_crossing_add(self) -> None:
        """Adding crossing pair correctly updates PK penalties."""
        # Start with nested pair
        matching = {(0, 10)}
        weights = {(0, 10): 1.0, (5, 15): 2.0}
        params = EnergyParams(pk_alpha=3.0, pk_gamma=0.0)

        _, cache = compute_energy(matching, weights, params)

        # Add crossing pair: (5, 15) crosses (0, 10)
        edge = (5, 15)
        delta, new_cache = compute_delta_energy(edge, matching, weights, params, cache, is_add=True)

        # Verify with full recompute
        new_matching = matching | {edge}
        expected_energy, expected_cache = compute_energy(new_matching, weights, params)

        assert new_cache.pk_pairs_count == expected_cache.pk_pairs_count
        assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

    def test_delta_pk_crossing_remove(self) -> None:
        """Removing crossing pair correctly updates PK penalties."""
        # Start with crossing pairs
        matching = {(0, 5), (2, 8)}  # These cross
        weights = {(0, 5): 1.0, (2, 8): 2.0}
        params = EnergyParams(pk_alpha=3.0, pk_gamma=0.0)

        _, cache = compute_energy(matching, weights, params)
        assert cache.pk_pairs_count == 2

        # Remove one crossing pair
        edge = (2, 8)
        delta, new_cache = compute_delta_energy(
            edge, matching, weights, params, cache, is_add=False
        )

        # After removal, no crossings
        assert new_cache.pk_pairs_count == 0

        # Verify with full recompute
        new_matching = matching - {edge}
        expected_energy, _ = compute_energy(new_matching, weights, params)
        assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

    def test_delta_lonely_pair_add(self) -> None:
        """Adding pair correctly updates lonely pair count."""
        # Start with lonely pair
        matching = {(0, 10)}
        weights = {(0, 10): 1.0, (1, 9): 2.0}
        params = EnergyParams(lonely_penalty=5.0)

        _, cache = compute_energy(matching, weights, params)
        assert cache.lonely_count == 1

        # Add stacking partner
        edge = (1, 9)
        delta, new_cache = compute_delta_energy(edge, matching, weights, params, cache, is_add=True)

        # Now neither is lonely
        assert new_cache.lonely_count == 0

        # Verify with full recompute
        new_matching = matching | {edge}
        expected_energy, _ = compute_energy(new_matching, weights, params)
        assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

    def test_delta_lonely_pair_remove(self) -> None:
        """Removing pair correctly updates lonely pair count."""
        # Start with stacked pairs
        matching = {(0, 10), (1, 9)}
        weights = {(0, 10): 1.0, (1, 9): 2.0}
        params = EnergyParams(lonely_penalty=5.0)

        _, cache = compute_energy(matching, weights, params)
        assert cache.lonely_count == 0

        # Remove one, making other lonely
        edge = (1, 9)
        delta, new_cache = compute_delta_energy(
            edge, matching, weights, params, cache, is_add=False
        )

        # Now the remaining is lonely
        assert new_cache.lonely_count == 1

        # Verify with full recompute
        new_matching = matching - {edge}
        expected_energy, _ = compute_energy(new_matching, weights, params)
        assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

    def test_edge_normalization(self) -> None:
        """Edge (j, i) is normalized to (i, j)."""
        matching = {(0, 10)}
        weights = {(0, 10): 5.0}
        params = EnergyParams()

        _, cache = compute_energy(matching, weights, params)

        # Use reversed edge
        delta, _ = compute_delta_energy((10, 0), matching, weights, params, cache, is_add=False)

        assert delta == 5.0

    def test_complex_scenario_reversibility(self) -> None:
        """Complex scenario maintains reversibility."""
        # Multiple crossing pairs
        matching = {(0, 10), (5, 15), (8, 20)}
        weights = {
            (0, 10): 1.0,
            (5, 15): 2.0,
            (8, 20): 3.0,
            (3, 18): 4.0,
        }
        params = EnergyParams(pk_alpha=2.0, pk_gamma=1.0, lonely_penalty=0.5)

        _, cache = compute_energy(matching, weights, params)

        # Add then remove multiple times
        edge = (3, 18)
        delta_add, cache2 = compute_delta_energy(
            edge, matching, weights, params, cache, is_add=True
        )

        matching2 = matching | {edge}
        delta_remove, cache3 = compute_delta_energy(
            edge, matching2, weights, params, cache2, is_add=False
        )

        # Should return to original
        assert abs(delta_add + delta_remove) < self.TOLERANCE
        assert abs(cache3.total_energy - cache.total_energy) < self.TOLERANCE


class TestComputeDeltaEnergyRandomized:
    """Randomized tests for delta energy computation."""

    TOLERANCE = 1e-10

    @pytest.mark.parametrize("seed", range(5))
    def test_random_add_remove_sequence(self, seed: int) -> None:
        """Random add/remove operations maintain consistency."""
        import random

        rng = random.Random(seed)

        # Generate random weights
        n = 30
        weights: dict[tuple[int, int], float] = {}
        for i in range(n):
            for j in range(i + 5, n):
                if rng.random() < 0.3:
                    weights[(i, j)] = rng.uniform(0.1, 5.0)

        params = EnergyParams(pk_alpha=2.0, pk_gamma=1.0, lonely_penalty=0.5)

        # Start empty
        matching: set[tuple[int, int]] = set()
        _, cache = compute_energy(matching, weights, params)

        # Perform random operations
        available = list(weights.keys())
        for _ in range(20):
            if not available and not matching:
                break

            if not matching or (available and rng.random() < 0.6):
                # Add
                edge = rng.choice(available)
                available.remove(edge)

                delta, new_cache = compute_delta_energy(
                    edge, matching, weights, params, cache, is_add=True
                )
                matching.add(edge)

                # Verify
                expected_energy, _ = compute_energy(matching, weights, params)
                assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

                cache = new_cache
            else:
                # Remove
                edge = rng.choice(list(matching))

                delta, new_cache = compute_delta_energy(
                    edge, matching, weights, params, cache, is_add=False
                )
                matching.remove(edge)
                available.append(edge)

                # Verify
                expected_energy, _ = compute_energy(matching, weights, params)
                assert abs(new_cache.total_energy - expected_energy) < self.TOLERANCE

                cache = new_cache
