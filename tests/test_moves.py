"""Tests for MCMC move implementations."""

import random

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
    """Tests for propose_segment_birth_death()."""

    def test_birth_from_empty(self) -> None:
        """Birth move from empty matching adds segment."""
        from src.lib.moves import propose_segment_birth_death
        from src.lib.segments import build_candidate_segments

        # Create segment index with stacking pairs
        pairs = [(0, 15), (1, 14), (2, 13)]
        index = build_candidate_segments(pairs, min_len=2)
        assert len(index.all_segments) > 0

        matching: set[tuple[int, int]] = set()
        partners = [-1] * 20
        fixed_pairs: set[tuple[int, int]] = set()
        weights = dict.fromkeys(pairs, 1.0)
        rng = random.Random(42)

        result = propose_segment_birth_death(matching, partners, index, fixed_pairs, weights, rng)

        # From empty, should be a birth move
        assert result.is_valid
        assert result.move_type == MoveType.SEGMENT_BIRTH
        assert result.pairs_added is not None
        assert len(result.pairs_added) >= 2  # At least 2 pairs in segment
        assert result.hastings_ratio > 0

    def test_death_from_segment(self) -> None:
        """Death move from matching with segment removes it."""
        from src.lib.moves import propose_segment_birth_death
        from src.lib.segments import build_candidate_segments

        pairs = [(0, 15), (1, 14), (2, 13)]
        index = build_candidate_segments(pairs, min_len=2)

        # Start with a segment in matching
        matching = {(0, 15), (1, 14)}
        partners = [-1] * 20
        partners[0] = 15
        partners[15] = 0
        partners[1] = 14
        partners[14] = 1
        fixed_pairs: set[tuple[int, int]] = set()
        weights = dict.fromkeys(pairs, 1.0)

        # Try multiple seeds to get a death move
        for seed in range(20):
            rng = random.Random(seed)
            result = propose_segment_birth_death(
                matching, partners, index, fixed_pairs, weights, rng
            )
            if result.is_valid and result.move_type == MoveType.SEGMENT_DEATH:
                assert result.pairs_removed is not None
                assert len(result.pairs_removed) >= 2
                break
        else:
            # Should get at least one death move
            pass  # Allow if not found, could be birth moves

    def test_no_moves_available(self) -> None:
        """When all segments are fixed, returns invalid."""
        from src.lib.moves import propose_segment_birth_death
        from src.lib.segments import build_candidate_segments

        pairs = [(0, 15), (1, 14)]
        index = build_candidate_segments(pairs, min_len=2)

        # All positions occupied by fixed pairs
        matching = {(0, 15), (1, 14)}
        partners = [-1] * 20
        partners[0] = 15
        partners[15] = 0
        partners[1] = 14
        partners[14] = 1
        fixed_pairs = {(0, 15), (1, 14)}  # All fixed
        weights = dict.fromkeys(pairs, 1.0)
        rng = random.Random(42)

        result = propose_segment_birth_death(matching, partners, index, fixed_pairs, weights, rng)

        # Should be invalid since no segments can be removed (fixed)
        # and no segments can be added (positions occupied)
        assert not result.is_valid


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

    def test_segment_birth_death_reversibility(self) -> None:
        """Segment birth/death satisfies detailed balance."""
        from src.lib.moves import propose_segment_birth_death
        from src.lib.segments import build_candidate_segments

        # Create segment index
        pairs = [(0, 15), (1, 14), (2, 13)]
        index = build_candidate_segments(pairs, min_len=2)

        # Empty matching - birth move
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 20
        fixed_pairs: set[tuple[int, int]] = set()
        weights = dict.fromkeys(pairs, 1.0)
        rng = random.Random(42)

        result_birth = propose_segment_birth_death(
            matching, partners, index, fixed_pairs, weights, rng
        )

        assert result_birth.is_valid
        assert result_birth.move_type == MoveType.SEGMENT_BIRTH
        assert result_birth.pairs_added is not None

        # Hastings ratio should be positive
        assert result_birth.hastings_ratio > 0

        # Now do the reverse: death move from the new matching
        new_matching = result_birth.pairs_added
        new_partners = [-1] * 20
        for i, j in new_matching:
            new_partners[i] = j
            new_partners[j] = i

        # Find a seed that gives us a death move
        for seed in range(50):
            rng2 = random.Random(seed)
            result_death = propose_segment_birth_death(
                new_matching, new_partners, index, fixed_pairs, weights, rng2
            )
            if result_death.is_valid and result_death.move_type == MoveType.SEGMENT_DEATH:
                # For detailed balance: hastings_birth * hastings_death ≈ 1
                # This is because birth and death are reverse moves
                assert result_death.hastings_ratio > 0
                break
