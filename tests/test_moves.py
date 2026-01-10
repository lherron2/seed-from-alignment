"""Tests for MCMC move implementations."""

import pytest

from src.lib.moves import MoveResult, MoveType


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


# Placeholder tests for move functions


class TestProposeToggle:
    """Tests for propose_toggle() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_toggle_add(self) -> None:
        """Toggle can add a pair."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_toggle_remove(self) -> None:
        """Toggle can remove a pair."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_toggle_respects_fixed(self) -> None:
        """Toggle cannot remove fixed pairs."""
        pass


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

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_toggle_reversibility(self) -> None:
        """Toggle move satisfies detailed balance."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_segment_birth_death_reversibility(self) -> None:
        """Segment birth/death satisfies detailed balance."""
        pass
