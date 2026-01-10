"""Tests for MCMC sampler."""

from src.lib.energy import EnergyParams
from src.lib.moves import MoveType
from src.lib.sampler import (
    SamplerConfig,
    SamplerDiagnostics,
    compute_exact_distribution,
    enumerate_all_matchings,
    sample_matchings_new,
)


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


class TestEnumerateAllMatchings:
    """Tests for enumerate_all_matchings()."""

    def test_empty_candidates(self) -> None:
        """No candidates gives single empty matching."""
        matchings = enumerate_all_matchings(10, [], min_hairpin=3)
        assert len(matchings) == 1
        assert matchings[0] == frozenset()

    def test_single_candidate(self) -> None:
        """Single valid candidate gives two matchings."""
        matchings = enumerate_all_matchings(10, [(0, 6)], min_hairpin=3)
        assert len(matchings) == 2
        assert frozenset() in matchings
        assert frozenset({(0, 6)}) in matchings

    def test_hairpin_filter(self) -> None:
        """Pairs violating hairpin are excluded."""
        matchings = enumerate_all_matchings(10, [(0, 2)], min_hairpin=3)
        # (0, 2) has only 1 unpaired base, so filtered out
        assert len(matchings) == 1
        assert matchings[0] == frozenset()

    def test_conflicting_pairs(self) -> None:
        """Conflicting pairs generate correct matchings."""
        # (0, 6) and (0, 8) conflict at position 0
        matchings = enumerate_all_matchings(10, [(0, 6), (0, 8)], min_hairpin=3)
        # Should have: empty, {(0,6)}, {(0,8)}
        assert len(matchings) == 3
        assert frozenset() in matchings
        assert frozenset({(0, 6)}) in matchings
        assert frozenset({(0, 8)}) in matchings

    def test_independent_pairs(self) -> None:
        """Independent pairs give all combinations."""
        # (0, 6) and (8, 14) don't conflict (positions are distinct)
        matchings = enumerate_all_matchings(20, [(0, 6), (8, 14)], min_hairpin=3)
        # Should have: empty, {(0,6)}, {(8,14)}, {(0,6), (8,14)}
        assert len(matchings) == 4

    def test_fixed_pairs(self) -> None:
        """Fixed pairs appear in all matchings."""
        matchings = enumerate_all_matchings(
            20, [(0, 6), (8, 14)], fixed_pairs={(0, 6)}, min_hairpin=3
        )
        # (0, 6) is fixed, so only {(0,6)} and {(0,6), (8,14)}
        assert len(matchings) == 2
        for m in matchings:
            assert (0, 6) in m

    def test_fixed_conflict(self) -> None:
        """Fixed pairs conflicting with optional gives correct result."""
        matchings = enumerate_all_matchings(
            20, [(0, 6), (0, 8)], fixed_pairs={(0, 6)}, min_hairpin=3
        )
        # (0, 6) fixed, (0, 8) conflicts with it
        assert len(matchings) == 1
        assert matchings[0] == frozenset({(0, 6)})


class TestComputeExactDistribution:
    """Tests for compute_exact_distribution()."""

    def test_empty_matchings(self) -> None:
        """Empty list gives empty distribution."""
        dist = compute_exact_distribution([], {}, EnergyParams())
        assert dist == {}

    def test_single_matching(self) -> None:
        """Single matching has probability 1."""
        matchings = [frozenset()]
        dist = compute_exact_distribution(matchings, {}, EnergyParams())
        assert abs(dist[frozenset()] - 1.0) < 1e-10

    def test_uniform_weights(self) -> None:
        """Equal weights give equal probabilities."""
        matchings = [frozenset(), frozenset({(0, 6)})]
        weights = {(0, 6): 0.0}  # Zero weight means no energy contribution
        dist = compute_exact_distribution(matchings, weights, EnergyParams(beta=1.0))
        assert abs(dist[frozenset()] - 0.5) < 1e-10
        assert abs(dist[frozenset({(0, 6)})] - 0.5) < 1e-10

    def test_high_weight_favored(self) -> None:
        """Higher weight pairs more likely."""
        matchings = [frozenset(), frozenset({(0, 6)})]
        weights = {(0, 6): 5.0}  # High weight = lower energy = more probable
        params = EnergyParams(beta=1.0)
        dist = compute_exact_distribution(matchings, weights, params)
        # Energy of empty = 0, energy of {(0,6)} = -5
        # p(empty) ∝ exp(0) = 1
        # p({(0,6)}) ∝ exp(5)
        assert dist[frozenset({(0, 6)})] > dist[frozenset()]


class TestSampleMatchingsNew:
    """Tests for sample_matchings_new()."""

    def test_deterministic_seed(self) -> None:
        """Same seed produces same output."""
        candidate_pairs = [(0, 6, 1.0), (8, 14, 1.0)]
        config = SamplerConfig(n_samples=5, burn_in=100, thin=10, seed=42)

        samples1, _ = sample_matchings_new(20, candidate_pairs, config)
        samples2, _ = sample_matchings_new(20, candidate_pairs, config)

        assert len(samples1) == len(samples2)
        for s1, s2 in zip(samples1, samples2):
            assert s1 == s2

    def test_correct_sample_count(self) -> None:
        """Returns requested number of samples."""
        candidate_pairs = [(0, 6, 1.0)]
        config = SamplerConfig(n_samples=10, burn_in=50, thin=5, seed=42)

        samples, _ = sample_matchings_new(10, candidate_pairs, config)

        assert len(samples) == 10

    def test_valid_matchings(self) -> None:
        """All samples are valid matchings (no position conflicts)."""
        candidate_pairs = [(0, 6, 1.0), (1, 8, 1.0), (2, 9, 1.0)]
        config = SamplerConfig(n_samples=20, burn_in=100, thin=10, seed=42)

        samples, _ = sample_matchings_new(15, candidate_pairs, config)

        for sample in samples:
            positions_used: set[int] = set()
            for i, j in sample:
                assert i not in positions_used, f"Position {i} used twice"
                assert j not in positions_used, f"Position {j} used twice"
                positions_used.add(i)
                positions_used.add(j)

    def test_respects_scaffold(self) -> None:
        """Scaffold pairs appear in all samples."""
        candidate_pairs = [(0, 6, 1.0), (8, 14, 1.0)]
        scaffold = {(0, 6)}
        config = SamplerConfig(n_samples=10, burn_in=100, thin=10, seed=42)

        samples, _ = sample_matchings_new(20, candidate_pairs, config, scaffold_pairs=scaffold)

        for sample in samples:
            assert (0, 6) in sample

    def test_diagnostics_recorded(self) -> None:
        """Diagnostics are properly recorded."""
        candidate_pairs = [(0, 6, 1.0)]
        config = SamplerConfig(n_samples=5, burn_in=100, thin=10, seed=42)

        samples, diagnostics = sample_matchings_new(10, candidate_pairs, config)

        assert len(diagnostics.energy_trace) == 5
        assert MoveType.TOGGLE in diagnostics.total_proposals
        assert diagnostics.total_proposals[MoveType.TOGGLE] > 0
