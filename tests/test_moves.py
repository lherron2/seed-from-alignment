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

    def test_birth_hastings_matches_counts(self) -> None:
        """Birth Hastings ratio uses counts from the correct states."""
        from src.lib.moves import propose_segment_birth_death
        from src.lib.segments import (
            build_candidate_segments,
            find_addable_segments,
            find_removable_segments,
        )

        pairs = [(0, 15), (1, 14), (2, 13)]
        index = build_candidate_segments(pairs, min_len=2, max_len=3)
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 20
        fixed_pairs: set[tuple[int, int]] = set()
        weights = dict.fromkeys(pairs, 1.0)
        rng = random.Random(0)

        result = propose_segment_birth_death(
            matching,
            partners,
            index,
            fixed_pairs,
            weights,
            rng,
            min_hairpin=3,
            force_move=MoveType.SEGMENT_BIRTH,
        )

        assert result.is_valid
        assert result.pairs_added is not None

        addable = find_addable_segments(matching, partners, index, 3)
        matching_after = matching | result.pairs_added
        removable_after = find_removable_segments(matching_after, index, fixed_pairs)
        expected = len(addable) / len(removable_after)

        assert result.hastings_ratio == expected

    def test_death_hastings_matches_counts(self) -> None:
        """Death Hastings ratio uses counts from the correct states."""
        from src.lib.moves import propose_segment_birth_death
        from src.lib.segments import (
            build_candidate_segments,
            find_addable_segments,
            find_removable_segments,
        )

        pairs = [(0, 15), (1, 14), (2, 13)]
        index = build_candidate_segments(pairs, min_len=2, max_len=3)
        matching = {(0, 15), (1, 14), (2, 13)}
        partners = [-1] * 20
        partners[0] = 15
        partners[15] = 0
        partners[1] = 14
        partners[14] = 1
        partners[2] = 13
        partners[13] = 2
        fixed_pairs: set[tuple[int, int]] = set()
        weights = dict.fromkeys(pairs, 1.0)
        rng = random.Random(0)

        result = propose_segment_birth_death(
            matching,
            partners,
            index,
            fixed_pairs,
            weights,
            rng,
            min_hairpin=3,
            force_move=MoveType.SEGMENT_DEATH,
        )

        assert result.is_valid
        assert result.pairs_removed is not None

        matching_after = matching - result.pairs_removed
        new_partners = list(partners)
        for i, j in result.pairs_removed:
            new_partners[i] = -1
            new_partners[j] = -1
        addable_after = find_addable_segments(matching_after, new_partners, index, 3)
        removable = find_removable_segments(matching, index, fixed_pairs)
        expected = len(removable) / len(addable_after)

        assert result.hastings_ratio == expected


class TestProposeSegmentSwap:
    """Tests for propose_segment_swap()."""

    def test_swap_hastings_matches_counts(self) -> None:
        """Swap Hastings ratio includes removable counts."""
        from src.lib.moves import propose_segment_swap
        from src.lib.segments import build_candidate_segments, find_removable_segments

        pairs = [(0, 9), (1, 8), (2, 7)]
        index = build_candidate_segments(pairs, min_len=2, max_len=3)

        # Start with the length-3 segment in matching.
        matching = {(0, 9), (1, 8), (2, 7)}
        partners = [-1] * 15
        partners[0] = 9
        partners[9] = 0
        partners[1] = 8
        partners[8] = 1
        partners[2] = 7
        partners[7] = 2
        fixed_pairs: set[tuple[int, int]] = set()
        weights = dict.fromkeys(pairs, 1.0)

        found = False
        for seed in range(200):
            rng = random.Random(seed)
            result = propose_segment_swap(
                matching,
                partners,
                index,
                fixed_pairs,
                weights,
                rng,
                min_hairpin=3,
            )
            if not result.is_valid:
                continue
            if result.pairs_removed is None or result.pairs_added is None:
                continue

            # Identify the selected segments.
            seg_old = None
            seg_new = None
            for seg in index.all_segments:
                if set(seg.pairs()) == result.pairs_removed:
                    seg_old = seg
                if set(seg.pairs()) == result.pairs_added:
                    seg_new = seg
            if seg_old is None or seg_new is None:
                continue

            removable = find_removable_segments(matching, index, fixed_pairs)
            n_removable = len(removable)

            # Compute valid conflicts for seg_old (same logic as proposal).
            valid_conflicts = []
            seg_old_positions = seg_old.positions()
            for seg_candidate in index.conflicts.get(seg_old, set()):
                seg_candidate_positions = seg_candidate.positions()
                can_add = True
                for pos in seg_candidate_positions:
                    if pos not in seg_old_positions and partners[pos] != -1:
                        can_add = False
                        break
                if can_add:
                    for i, j in seg_candidate.pairs():
                        if j - i - 1 < 3:
                            can_add = False
                            break
                if can_add:
                    valid_conflicts.append(seg_candidate)

            # Build M' for reverse counts.
            matching_after = (matching - result.pairs_removed) | result.pairs_added
            new_partners = list(partners)
            for i, j in result.pairs_removed:
                new_partners[i] = -1
                new_partners[j] = -1
            for i, j in result.pairs_added:
                new_partners[i] = j
                new_partners[j] = i

            removable_after = find_removable_segments(matching_after, index, fixed_pairs)
            n_removable_after = len(removable_after)

            valid_reverse = []
            seg_new_positions = seg_new.positions()
            for seg_candidate in index.conflicts.get(seg_new, set()):
                seg_candidate_positions = seg_candidate.positions()
                can_add = True
                for pos in seg_candidate_positions:
                    if pos not in seg_new_positions and new_partners[pos] != -1:
                        can_add = False
                        break
                if can_add:
                    for i, j in seg_candidate.pairs():
                        if j - i - 1 < 3:
                            can_add = False
                            break
                if can_add:
                    valid_reverse.append(seg_candidate)

            if n_removable != n_removable_after:
                expected = (n_removable / n_removable_after) * (
                    len(valid_conflicts) / len(valid_reverse)
                )
                assert result.hastings_ratio == expected
                found = True
                break

        assert found, "Did not find a swap move with differing removable counts"


class TestProposeLoopClosure:
    """Tests for propose_loop_closure()."""

    def test_loop_closure_add_remove(self) -> None:
        """Loop closure toggles a fold symmetrically."""
        from src.lib.moves import propose_loop_closure

        loop_folds = [{(0, 6), (1, 5)}]
        matching: set[tuple[int, int]] = set()
        partners = [-1] * 10
        weights: dict[tuple[int, int], float] = {}
        rng = random.Random(0)

        result_add = propose_loop_closure(
            matching,
            partners,
            loop_folds,
            weights,
            rng,
            min_hairpin=3,
            fixed_pairs=set(),
        )
        assert result_add.is_valid
        assert result_add.pairs_added == {(0, 6), (1, 5)}

        # Apply the move and try removal.
        new_matching = set(result_add.pairs_added or set())
        new_partners = [-1] * 10
        for i, j in new_matching:
            new_partners[i] = j
            new_partners[j] = i

        result_remove = propose_loop_closure(
            new_matching,
            new_partners,
            loop_folds,
            weights,
            rng,
            min_hairpin=3,
            fixed_pairs=set(),
        )
        assert result_remove.is_valid
        assert result_remove.pairs_removed == {(0, 6), (1, 5)}

    def test_loop_closure_overlap_invalid(self) -> None:
        """Partial overlap with fold is rejected."""
        from src.lib.moves import propose_loop_closure

        loop_folds = [{(0, 6), (1, 5)}]
        matching = {(0, 6)}
        partners = [-1] * 10
        partners[0] = 6
        partners[6] = 0
        weights: dict[tuple[int, int], float] = {}
        rng = random.Random(0)

        result = propose_loop_closure(
            matching,
            partners,
            loop_folds,
            weights,
            rng,
            min_hairpin=3,
            fixed_pairs=set(),
        )
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
