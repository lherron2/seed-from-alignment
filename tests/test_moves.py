"""Tests for MCMC move implementations."""

import random

import pytest

from src.lib.moves import MoveResult, MoveType, propose_toggle


class TestMoveResult:
    """Tests for MoveResult dataclass."""

    def test_default_values(self) -> None:
        """Default MoveResult values."""
        result = MoveResult()
        assert result.accepted is False
        assert result.move_type == MoveType.TOGGLE
        assert result.hastings_ratio == 1.0
        assert result.is_valid is True


class TestMoveType:
    """Tests for MoveType enum."""

    def test_move_types_exist(self) -> None:
        """All expected move types exist."""
        assert MoveType.TOGGLE
        assert MoveType.SEGMENT_BIRTH
        assert MoveType.SEGMENT_DEATH
        assert MoveType.SEGMENT_SWAP
        assert MoveType.LOOP_CLOSURE


class TestProposeToggle:
    """Tests for propose_toggle()."""

    def test_toggle_add(self) -> None:
        """Toggle can add a pair when positions are free."""
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 10
        candidate_pairs = [(0, 6, 1.0), (1, 8, 1.0)]
        fixed_pairs: set[tuple[int, int]] = set()
        weights = {(0, 6): 1.0, (1, 8): 1.0}
        rng = random.Random(42)

        result = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )

        assert result.move_type == MoveType.TOGGLE
        assert result.is_valid is True
        assert result.pairs_added is not None
        assert len(result.pairs_added) == 1
        assert result.pairs_removed is None
        assert result.hastings_ratio == 1.0

    def test_toggle_remove(self) -> None:
        """Toggle can remove a pair."""
        matching = {(0, 6)}
        partners = [-1] * 10
        partners[0] = 6
        partners[6] = 0
        candidate_pairs = [(0, 6, 1.0)]
        fixed_pairs: set[tuple[int, int]] = set()
        weights = {(0, 6): 1.0}
        rng = random.Random(42)

        result = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )

        assert result.move_type == MoveType.TOGGLE
        assert result.is_valid is True
        assert result.pairs_removed is not None
        assert len(result.pairs_removed) == 1
        assert result.pairs_added is None
        assert result.hastings_ratio == 1.0

    def test_toggle_respects_fixed(self) -> None:
        """Toggle cannot remove fixed pairs."""
        matching = {(0, 6)}
        partners = [-1] * 10
        partners[0] = 6
        partners[6] = 0
        candidate_pairs = [(0, 6, 1.0)]
        fixed_pairs = {(0, 6)}
        weights = {(0, 6): 1.0}
        rng = random.Random(42)

        result = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )

        assert result.move_type == MoveType.TOGGLE
        assert result.is_valid is False  # Can't remove fixed pair

    def test_toggle_respects_conflicts(self) -> None:
        """Toggle cannot add when position is already paired."""
        matching = {(0, 6)}
        partners = [-1] * 10
        partners[0] = 6
        partners[6] = 0
        candidate_pairs = [(0, 8, 1.0)]  # Position 0 already paired
        fixed_pairs: set[tuple[int, int]] = set()
        weights = {(0, 8): 1.0}
        rng = random.Random(42)

        result = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )

        assert result.is_valid is False

    def test_toggle_respects_hairpin(self) -> None:
        """Toggle cannot add pairs violating hairpin constraint."""
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 10
        candidate_pairs = [(0, 2, 1.0)]  # Only 1 unpaired between, need 3
        fixed_pairs: set[tuple[int, int]] = set()
        weights = {(0, 2): 1.0}
        rng = random.Random(42)

        result = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )

        assert result.is_valid is False

    def test_toggle_empty_candidates(self) -> None:
        """Toggle with no candidates returns invalid."""
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 10
        candidate_pairs: list[tuple[int, int, float]] = []
        fixed_pairs: set[tuple[int, int]] = set()
        weights: dict[tuple[int, int], float] = {}
        rng = random.Random(42)

        result = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )

        assert result.is_valid is False

    def test_toggle_symmetric(self) -> None:
        """Toggle Hastings ratio is always 1."""
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 10
        candidate_pairs = [(0, 6, 1.0), (1, 8, 2.0)]
        fixed_pairs: set[tuple[int, int]] = set()
        weights = {(0, 6): 1.0, (1, 8): 2.0}

        # Run many toggles and verify Hastings ratio
        for seed in range(20):
            rng = random.Random(seed)
            result = propose_toggle(
                matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
            )
            if result.is_valid:
                assert result.hastings_ratio == 1.0


class TestProposeSegmentBirthDeath:
    """Tests for propose_segment_birth_death() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_birth_hastings_ratio(self) -> None:
        """Birth move has correct Hastings ratio."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_death_hastings_ratio(self) -> None:
        """Death move has correct Hastings ratio."""
        pass


class TestMoveReversibility:
    """Tests for move reversibility (detailed balance)."""

    def test_toggle_reversibility(self) -> None:
        """Toggle move satisfies detailed balance (symmetric)."""
        # For toggle, the move is symmetric so detailed balance
        # is trivially satisfied with Hastings ratio = 1
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 10
        candidate_pairs = [(0, 6, 1.0)]
        fixed_pairs: set[tuple[int, int]] = set()
        weights = {(0, 6): 1.0}
        rng = random.Random(42)

        # Propose add
        result_add = propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )
        assert result_add.is_valid
        assert result_add.hastings_ratio == 1.0
        assert result_add.pairs_added == {(0, 6)}

        # After adding, propose remove
        new_matching = {(0, 6)}
        new_partners = [-1] * 10
        new_partners[0] = 6
        new_partners[6] = 0

        result_remove = propose_toggle(
            new_matching, new_partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin=3
        )
        assert result_remove.is_valid
        assert result_remove.hastings_ratio == 1.0
        assert result_remove.pairs_removed == {(0, 6)}

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_segment_birth_death_reversibility(self) -> None:
        """Segment birth/death satisfies detailed balance."""
        pass
