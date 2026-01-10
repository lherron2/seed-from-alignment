"""Tests for segment/helix data structures."""

import pytest

from src.lib.segments import (
    Segment,
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


# Placeholder tests for functions to be implemented in Phase 3


class TestBuildCandidateSegments:
    """Tests for build_candidate_segments() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_basic_segments(self) -> None:
        """Build segments from candidate pairs."""
        pass


class TestFindAddableSegments:
    """Tests for find_addable_segments() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_find_addable(self) -> None:
        """Find segments that can be added."""
        pass


class TestFindRemovableSegments:
    """Tests for find_removable_segments() - placeholder."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_find_removable(self) -> None:
        """Find segments that can be removed."""
        pass
