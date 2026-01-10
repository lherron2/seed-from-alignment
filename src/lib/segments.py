"""
Segment/helix data structures for MCMC RNA structure sampling.

This module provides:
- Segment dataclass for representing helical segments
- SegmentIndex for precomputing candidate segments and conflicts
- Functions for finding addable/removable segments

A segment is a contiguous helix (stem) of stacked base pairs.
Segment-level moves improve MCMC mixing by proposing coordinated
changes to multiple pairs at once.

See spec/04_SEGMENTS_MODULE.md for detailed specification.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "Segment",
    "SegmentIndex",
    "build_candidate_segments",
    "segments_conflict",
    "identify_segments_in_matching",
    "find_addable_segments",
    "find_removable_segments",
]


@dataclass(frozen=True)
class Segment:
    """A contiguous helix (stem) of stacked base pairs.

    A segment is defined by its start and end positions on both
    strands, forming a set of consecutive stacked pairs.

    Example: Segment(start5=2, start3=15, length=4) represents:
        pairs: [(2,15), (3,14), (4,13), (5,12)]

    Attributes:
        start5: 5' strand start position
        start3: 3' strand start position (paired with start5)
        length: Number of stacked pairs
    """

    start5: int
    start3: int
    length: int

    def __post_init__(self) -> None:
        if self.length < 1:
            raise ValueError("Segment length must be >= 1")
        if self.start5 >= self.start3:
            raise ValueError("start5 must be less than start3")

    def pairs(self) -> list[tuple[int, int]]:
        """Return list of base pairs in this segment.

        Returns:
            List of (i, j) pairs with i < j
        """
        return [(self.start5 + k, self.start3 - k) for k in range(self.length)]

    def positions(self) -> set[int]:
        """Return set of all positions occupied by this segment.

        Returns:
            Set of position indices
        """
        result = set()
        for i, j in self.pairs():
            result.add(i)
            result.add(j)
        return result

    def as_frozenset(self) -> frozenset[tuple[int, int]]:
        """Return pairs as a frozen set for hashing/comparison.

        Returns:
            Frozen set of pairs
        """
        return frozenset(self.pairs())


@dataclass
class SegmentIndex:
    """Index of candidate segments for efficient lookup.

    Precomputes:
    - All candidate segments from the candidate pairs
    - Conflict relationships between segments

    Attributes:
        all_segments: Set of all candidate segments
        segment_by_pair: Map from pair to segments containing it
        conflicts: Map from segment to conflicting segments
        length: Sequence length
    """

    all_segments: set[Segment] = field(default_factory=set)
    segment_by_pair: dict[tuple[int, int], set[Segment]] = field(default_factory=dict)
    conflicts: dict[Segment, set[Segment]] = field(default_factory=dict)
    length: int = 0


def build_candidate_segments(
    candidate_pairs: list[tuple[int, int]] | set[tuple[int, int]],
    min_len: int = 2,
    max_len: int = 10,
) -> SegmentIndex:
    """Build index of candidate segments from candidate pairs.

    Args:
        candidate_pairs: List or set of candidate (i, j) pairs
        min_len: Minimum segment length
        max_len: Maximum segment length

    Returns:
        SegmentIndex with all candidate segments
    """
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("build_candidate_segments not yet implemented")


def segments_conflict(seg1: Segment, seg2: Segment) -> bool:
    """Check if two segments conflict (share positions).

    Args:
        seg1: First segment
        seg2: Second segment

    Returns:
        True if segments share any positions
    """
    return bool(seg1.positions() & seg2.positions())


def identify_segments_in_matching(
    matching: set[tuple[int, int]],
    index: SegmentIndex,
) -> set[Segment]:
    """Identify complete segments present in a matching.

    A segment is "in" the matching if ALL its pairs are present.

    Args:
        matching: Current set of base pairs
        index: Precomputed segment index

    Returns:
        Set of complete segments in the matching
    """
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("identify_segments_in_matching not yet implemented")


def find_addable_segments(
    matching: set[tuple[int, int]],
    partners: list[int],
    index: SegmentIndex,
) -> list[Segment]:
    """Find segments that can be added to the current matching.

    A segment is addable if:
    - None of its positions are currently paired
    - Adding it would not violate hairpin constraints

    Args:
        matching: Current set of base pairs
        partners: Partners array (partners[i] = j if paired, -1 if unpaired)
        index: Precomputed segment index

    Returns:
        List of addable segments
    """
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("find_addable_segments not yet implemented")


def find_removable_segments(
    matching: set[tuple[int, int]],
    index: SegmentIndex,
    fixed_pairs: set[tuple[int, int]] | None = None,
) -> list[Segment]:
    """Find complete segments that can be removed from matching.

    A segment is removable if:
    - ALL its pairs are in the matching
    - None of its pairs are in fixed_pairs

    Args:
        matching: Current set of base pairs
        index: Precomputed segment index
        fixed_pairs: Pairs that cannot be removed (scaffold)

    Returns:
        List of removable segments
    """
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("find_removable_segments not yet implemented")
