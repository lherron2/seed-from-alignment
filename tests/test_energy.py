"""Tests for energy computation module."""

import pytest

from src.lib.energy import (
    EnergyCache,
    EnergyParams,
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


class TestEnergyCache:
    """Tests for EnergyCache dataclass."""

    def test_default_values(self) -> None:
        """Default cache values are zero."""
        cache = EnergyCache()
        assert cache.total_energy == 0.0
        assert cache.weight_sum == 0.0
        assert cache.pk_pairs_count == 0


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


# Placeholder tests for functions to be implemented in Phase 1


class TestComputeEnergy:
    """Tests for compute_energy() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 1")
    def test_empty_matching(self) -> None:
        """Empty matching has zero energy."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 1")
    def test_single_pair(self) -> None:
        """Single pair energy is negative weight."""
        pass


class TestComputeDeltaEnergy:
    """Tests for compute_delta_energy() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 1")
    def test_delta_matches_recompute(self) -> None:
        """Delta must match full recompute within tolerance."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 1")
    def test_delta_reversibility(self) -> None:
        """delta(add) = -delta(remove) for same edge."""
        pass
