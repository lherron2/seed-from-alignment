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

    A segment is a contiguous helix: [(i, j), (i+1, j-1), ...].
    We find all maximal segments where each consecutive pair stacks.

    Args:
        candidate_pairs: List or set of candidate (i, j) pairs
        min_len: Minimum segment length
        max_len: Maximum segment length

    Returns:
        SegmentIndex with all candidate segments
    """
    # Normalize pairs to (i, j) with i < j
    pairs_set: set[tuple[int, int]] = set()
    for pair in candidate_pairs:
        i, j = pair[0], pair[1]
        if i > j:
            i, j = j, i
        pairs_set.add((i, j))

    if not pairs_set:
        return SegmentIndex()

    # Find connected components under stacking relation
    # Two pairs (i, j) and (i', j') stack if i' = i+1 and j' = j-1
    visited: set[tuple[int, int]] = set()
    segments: set[Segment] = set()
    segment_by_pair: dict[tuple[int, int], set[Segment]] = {}

    # Build stacking adjacency for efficient lookup
    stacking_neighbors: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for i, j in pairs_set:
        neighbors = []
        inner = (i + 1, j - 1)
        outer = (i - 1, j + 1)
        if inner in pairs_set:
            neighbors.append(inner)
        if outer in pairs_set:
            neighbors.append(outer)
        stacking_neighbors[(i, j)] = neighbors

    # Find all maximal contiguous chains (segments)
    for start_pair in pairs_set:
        if start_pair in visited:
            continue

        # BFS/DFS to find all connected pairs under stacking
        chain: list[tuple[int, int]] = []
        stack = [start_pair]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            chain.append(curr)
            for neighbor in stacking_neighbors.get(curr, []):
                if neighbor not in visited:
                    stack.append(neighbor)

        # Sort chain by 5' position to get ordered segment
        chain.sort(key=lambda p: p[0])

        # The chain is a maximal stacking chain
        # Create segments of various lengths from this chain
        chain_len = len(chain)
        if chain_len >= min_len:
            # For each starting position and valid length, create a segment
            for start_idx in range(chain_len):
                for seg_len in range(min_len, min(max_len + 1, chain_len - start_idx + 1)):
                    if start_idx + seg_len > chain_len:
                        break
                    # Verify the pairs are truly consecutive stacking
                    valid = True
                    for k in range(seg_len - 1):
                        p1 = chain[start_idx + k]
                        p2 = chain[start_idx + k + 1]
                        if p2[0] != p1[0] + 1 or p2[1] != p1[1] - 1:
                            valid = False
                            break
                    if valid:
                        first_pair = chain[start_idx]
                        seg = Segment(
                            start5=first_pair[0],
                            start3=first_pair[1],
                            length=seg_len,
                        )
                        segments.add(seg)

    # Build segment_by_pair mapping
    for seg in segments:
        for pair in seg.pairs():
            if pair not in segment_by_pair:
                segment_by_pair[pair] = set()
            segment_by_pair[pair].add(seg)

    # Build conflict graph
    conflicts: dict[Segment, set[Segment]] = {seg: set() for seg in segments}
    segment_list = list(segments)
    for i, seg1 in enumerate(segment_list):
        for seg2 in segment_list[i + 1 :]:
            if segments_conflict(seg1, seg2):
                conflicts[seg1].add(seg2)
                conflicts[seg2].add(seg1)

    # Determine max length for sequence length field
    max_pos = 0
    for i, j in pairs_set:
        max_pos = max(max_pos, i, j)

    return SegmentIndex(
        all_segments=segments,
        segment_by_pair=segment_by_pair,
        conflicts=conflicts,
        length=max_pos + 1,
    )


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
    result: set[Segment] = set()
    for seg in index.all_segments:
        # Check if all pairs of the segment are in the matching
        seg_pairs = seg.pairs()
        if all((i, j) in matching for i, j in seg_pairs):
            result.add(seg)
    return result


def find_addable_segments(
    matching: set[tuple[int, int]],
    partners: list[int],
    index: SegmentIndex,
    min_hairpin: int = 3,
) -> list[Segment]:
    """Find segments that can be added to the current matching.

    A segment is addable if:
    - None of its positions are currently paired
    - Adding it would not violate hairpin constraints

    Args:
        matching: Current set of base pairs
        partners: Partners array (partners[i] = j if paired, -1 if unpaired)
        index: Precomputed segment index
        min_hairpin: Minimum hairpin loop size

    Returns:
        List of addable segments
    """
    result: list[Segment] = []
    for seg in index.all_segments:
        # Check if the segment is already in the matching
        seg_pairs = seg.pairs()
        if any((i, j) in matching for i, j in seg_pairs):
            continue

        # Check if all positions are free
        can_add = True
        for i, j in seg_pairs:
            if partners[i] != -1 or partners[j] != -1:
                can_add = False
                break
            # Check hairpin constraint
            if j - i - 1 < min_hairpin:
                can_add = False
                break

        if can_add:
            result.append(seg)

    return result


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
    if fixed_pairs is None:
        fixed_pairs = set()

    result: list[Segment] = []
    for seg in index.all_segments:
        seg_pairs = seg.pairs()

        # Check if all pairs are in matching
        if not all((i, j) in matching for i, j in seg_pairs):
            continue

        # Check if none of the pairs are fixed
        if any((i, j) in fixed_pairs for i, j in seg_pairs):
            continue

        result.append(seg)

    return result
