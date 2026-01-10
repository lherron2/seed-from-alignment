"""Tests for MCMC sampler."""

import pytest

from src.lib.moves import MoveType
from src.lib.sampler import SamplerConfig, SamplerDiagnostics


class TestSamplerConfig:
    """Tests for SamplerConfig dataclass."""

    def test_default_values(self) -> None:
        """Default configuration values are sensible."""
        config = SamplerConfig()
        assert config.n_samples > 0
        assert config.burn_in >= 0
        assert config.thin >= 1


class TestSamplerDiagnostics:
    """Tests for SamplerDiagnostics class."""

    def test_record_move(self) -> None:
        """Recording moves updates counts."""
        diag = SamplerDiagnostics()
        diag.record_move(MoveType.TOGGLE, accepted=True)
        diag.record_move(MoveType.TOGGLE, accepted=False)
        diag.record_move(MoveType.TOGGLE, accepted=True)

        assert diag.total_proposals[MoveType.TOGGLE] == 3
        assert diag.accepted[MoveType.TOGGLE] == 2

    def test_acceptance_rate(self) -> None:
        """Acceptance rate calculation."""
        diag = SamplerDiagnostics()
        diag.record_move(MoveType.TOGGLE, accepted=True)
        diag.record_move(MoveType.TOGGLE, accepted=False)

        rate = diag.acceptance_rate(MoveType.TOGGLE)
        assert rate == 0.5

    def test_acceptance_rate_no_proposals(self) -> None:
        """Acceptance rate is 0 when no proposals."""
        diag = SamplerDiagnostics()
        rate = diag.acceptance_rate(MoveType.TOGGLE)
        assert rate == 0.0

    def test_record_state(self) -> None:
        """Recording state updates traces."""
        diag = SamplerDiagnostics()
        diag.record_state(energy=-10.5, matching_size=5)
        diag.record_state(energy=-12.0, matching_size=6)

        assert len(diag.energy_trace) == 2
        assert len(diag.sample_sizes) == 2

    def test_summary(self) -> None:
        """Summary generation."""
        diag = SamplerDiagnostics()
        diag.record_move(MoveType.TOGGLE, accepted=True)
        diag.record_state(energy=-10.0, matching_size=5)

        summary = diag.summary()
        assert "n_samples" in summary
        assert "acceptance_rates" in summary
        assert "ess" in summary


# Placeholder tests for sampler functions


class TestSampleMatchingsNew:
    """Tests for sample_matchings_new() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_deterministic_seed(self) -> None:
        """Same seed produces same output."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_correct_sample_count(self) -> None:
        """Returns requested number of samples."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_valid_matchings(self) -> None:
        """All samples are valid matchings."""
        pass


class TestEnumerateAllMatchings:
    """Tests for enumerate_all_matchings() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_tiny_system(self) -> None:
        """Enumerate matchings for tiny system."""
        pass
