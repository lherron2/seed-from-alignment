"""Tests for segment/helix data structures."""

import pytest

from src.lib.segments import (
    Segment,
    build_candidate_segments,
    find_addable_segments,
    find_removable_segments,
    identify_segments_in_matching,
    segments_conflict,
)


class TestSegment:
    """Tests for Segment dataclass."""

    def test_basic_segment(self) -> None:
        """Basic segment creation."""
        seg = Segment(start5=2, start3=15, length=4)
        assert seg.start5 == 2
        assert seg.start3 == 15
        assert seg.length == 4

    def test_pairs(self) -> None:
        """Segment pairs generation."""
        seg = Segment(start5=2, start3=15, length=4)
        pairs = seg.pairs()
        expected = [(2, 15), (3, 14), (4, 13), (5, 12)]
        assert pairs == expected

    def test_positions(self) -> None:
        """Segment positions set."""
        seg = Segment(start5=2, start3=15, length=2)
        positions = seg.positions()
        assert positions == {2, 3, 14, 15}

    def test_invalid_length(self) -> None:
        """Segment with length < 1 raises error."""
        with pytest.raises(ValueError):
            Segment(start5=2, start3=15, length=0)

    def test_invalid_order(self) -> None:
        """Segment with start5 >= start3 raises error."""
        with pytest.raises(ValueError):
            Segment(start5=10, start3=5, length=2)

    def test_segment_is_frozen(self) -> None:
        """Segment should be hashable (frozen)."""
        seg = Segment(start5=2, start3=15, length=4)
        # Should be able to add to set
        s = {seg}
        assert seg in s


class TestSegmentsConflict:
    """Tests for segments_conflict() function."""

    def test_same_segment_conflicts(self) -> None:
        """A segment conflicts with itself."""
        seg = Segment(start5=2, start3=15, length=4)
        assert segments_conflict(seg, seg) is True

    def test_overlapping_segments_conflict(self) -> None:
        """Overlapping segments conflict."""
        seg1 = Segment(start5=2, start3=15, length=4)
        seg2 = Segment(start5=4, start3=13, length=2)
        assert segments_conflict(seg1, seg2) is True

    def test_disjoint_segments_no_conflict(self) -> None:
        """Disjoint segments don't conflict."""
        seg1 = Segment(start5=0, start3=10, length=2)
        seg2 = Segment(start5=20, start3=30, length=2)
        assert segments_conflict(seg1, seg2) is False


class TestBuildCandidateSegments:
    """Tests for build_candidate_segments()."""

    def test_empty_pairs(self) -> None:
        """Empty candidate pairs gives empty index."""
        index = build_candidate_segments([])
        assert len(index.all_segments) == 0

    def test_single_pair_too_short(self) -> None:
        """Single pair doesn't form a segment (min_len=2)."""
        index = build_candidate_segments([(0, 10)], min_len=2)
        assert len(index.all_segments) == 0

    def test_single_pair_min_len_1(self) -> None:
        """Single pair forms segment if min_len=1."""
        index = build_candidate_segments([(0, 10)], min_len=1)
        assert len(index.all_segments) == 1
        seg = list(index.all_segments)[0]
        assert seg.pairs() == [(0, 10)]

    def test_stacking_pairs(self) -> None:
        """Stacking pairs form segments."""
        # (0, 10), (1, 9), (2, 8) are stacking
        pairs = [(0, 10), (1, 9), (2, 8)]
        index = build_candidate_segments(pairs, min_len=2, max_len=3)

        # Should have segments of length 2 and 3
        lengths = {seg.length for seg in index.all_segments}
        assert 2 in lengths
        assert 3 in lengths

    def test_non_stacking_pairs(self) -> None:
        """Non-stacking pairs don't form segments."""
        # (0, 10) and (5, 15) don't stack
        pairs = [(0, 10), (5, 15)]
        index = build_candidate_segments(pairs, min_len=2)
        assert len(index.all_segments) == 0

    def test_conflict_detection(self) -> None:
        """Conflicting segments are detected."""
        # Two overlapping helices: (0, 10), (1, 9) and (0, 12), (1, 11)
        pairs = [(0, 10), (1, 9), (0, 12), (1, 11)]
        index = build_candidate_segments(pairs, min_len=2)

        # Find segments and check conflicts
        for seg in index.all_segments:
            conflicts = index.conflicts.get(seg, set())
            # Segments sharing position 0 or 1 should conflict
            for other in conflicts:
                assert seg.positions() & other.positions()

    def test_segment_by_pair(self) -> None:
        """segment_by_pair mapping is correct."""
        pairs = [(0, 10), (1, 9), (2, 8)]
        index = build_candidate_segments(pairs, min_len=2, max_len=3)

        # Each pair should map to at least one segment
        for pair in [(0, 10), (1, 9), (2, 8)]:
            assert pair in index.segment_by_pair
            assert len(index.segment_by_pair[pair]) > 0


class TestIdentifySegmentsInMatching:
    """Tests for identify_segments_in_matching()."""

    def test_empty_matching(self) -> None:
        """Empty matching has no segments."""
        pairs = [(0, 10), (1, 9)]
        index = build_candidate_segments(pairs, min_len=2)
        result = identify_segments_in_matching(set(), index)
        assert len(result) == 0

    def test_complete_segment_in_matching(self) -> None:
        """Complete segment is identified."""
        pairs = [(0, 10), (1, 9), (2, 8)]
        index = build_candidate_segments(pairs, min_len=2, max_len=3)

        # Matching has the full helix
        matching = {(0, 10), (1, 9), (2, 8)}
        result = identify_segments_in_matching(matching, index)

        # Should find the length-3 segment
        assert any(seg.length == 3 for seg in result)

    def test_partial_segment_not_counted(self) -> None:
        """Partial segment is not identified."""
        pairs = [(0, 10), (1, 9), (2, 8)]
        index = build_candidate_segments(pairs, min_len=3, max_len=3)

        # Matching has only first two pairs
        matching = {(0, 10), (1, 9)}
        result = identify_segments_in_matching(matching, index)

        # Should not find the length-3 segment
        assert len(result) == 0


class TestFindAddableSegments:
    """Tests for find_addable_segments()."""

    def test_empty_matching_all_addable(self) -> None:
        """Empty matching allows adding any segment."""
        pairs = [(0, 10), (1, 9)]
        index = build_candidate_segments(pairs, min_len=2)
        partners = [-1] * 15

        result = find_addable_segments(set(), partners, index)
        assert len(result) > 0

    def test_occupied_positions_not_addable(self) -> None:
        """Segments with occupied positions are not addable."""
        pairs = [(0, 10), (1, 9), (5, 15), (6, 14)]
        index = build_candidate_segments(pairs, min_len=2)

        # Set up matching with first segment
        matching = {(0, 10), (1, 9)}
        partners = [-1] * 20
        partners[0] = 10
        partners[10] = 0
        partners[1] = 9
        partners[9] = 1

        result = find_addable_segments(matching, partners, index)

        # Only segment at (5, 15), (6, 14) should be addable
        for seg in result:
            seg_pairs = set(seg.pairs())
            assert (0, 10) not in seg_pairs
            assert (1, 9) not in seg_pairs


class TestFindRemovableSegments:
    """Tests for find_removable_segments()."""

    def test_empty_matching_none_removable(self) -> None:
        """Empty matching has no removable segments."""
        pairs = [(0, 10), (1, 9)]
        index = build_candidate_segments(pairs, min_len=2)

        result = find_removable_segments(set(), index)
        assert len(result) == 0

    def test_complete_segment_removable(self) -> None:
        """Complete segment in matching is removable."""
        pairs = [(0, 10), (1, 9)]
        index = build_candidate_segments(pairs, min_len=2)

        matching = {(0, 10), (1, 9)}
        result = find_removable_segments(matching, index)
        assert len(result) > 0

    def test_fixed_pairs_not_removable(self) -> None:
        """Segments containing fixed pairs are not removable."""
        pairs = [(0, 10), (1, 9)]
        index = build_candidate_segments(pairs, min_len=2)

        matching = {(0, 10), (1, 9)}
        fixed = {(0, 10)}  # First pair is fixed

        result = find_removable_segments(matching, index, fixed)
        # No segments should be removable because they all contain the fixed pair
        assert len(result) == 0
