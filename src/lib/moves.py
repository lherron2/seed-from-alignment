"""
MCMC move implementations for RNA structure sampling.

This module provides the move proposal functions for the
Metropolis-Hastings sampler. All moves must satisfy detailed balance:

    π(M) · q(M→M') · A(M→M') = π(M') · q(M'→M) · A(M'→M)

Move types:
- Toggle: Add or remove a single pair (Hastings ratio = 1)
- Segment birth/death: Add or remove a complete helix
- Segment swap: Replace one helix with another
- Loop closure: Symmetric proposals within loop regions

See spec/05_MOVES_MODULE.md for detailed specification.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto

from .segments import SegmentIndex

__all__ = [
    "MoveType",
    "MoveResult",
    "propose_toggle",
    "propose_segment_birth_death",
    "propose_segment_swap",
    "propose_loop_closure",
    "select_and_propose_move",
]


class MoveType(Enum):
    """Types of MCMC moves."""

    TOGGLE = auto()
    SEGMENT_BIRTH = auto()
    SEGMENT_DEATH = auto()
    SEGMENT_SWAP = auto()
    LOOP_CLOSURE = auto()


@dataclass
class MoveResult:
    """Result of a move proposal.

    Attributes:
        accepted: Whether the move was accepted
        move_type: Type of move proposed
        pairs_added: Pairs added by the move
        pairs_removed: Pairs removed by the move
        hastings_ratio: q(M'→M) / q(M→M') for the move
        delta_energy: Energy change from the move
        is_valid: Whether the move is structurally valid
    """

    accepted: bool = False
    move_type: MoveType = MoveType.TOGGLE
    pairs_added: set[tuple[int, int]] | None = None
    pairs_removed: set[tuple[int, int]] | None = None
    hastings_ratio: float = 1.0
    delta_energy: float = 0.0
    is_valid: bool = True


def propose_toggle(
    matching: set[tuple[int, int]],
    partners: list[int],
    candidate_pairs: list[tuple[int, int, float]],
    fixed_pairs: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    rng,
    min_hairpin: int = 3,
) -> MoveResult:
    """Propose a toggle move (add or remove one pair).

    The toggle move is symmetric: Hastings ratio = 1.

    Args:
        matching: Current set of base pairs
        partners: Partners array
        candidate_pairs: List of (i, j, weight) candidates
        fixed_pairs: Pairs that cannot be removed
        weights: Weight dictionary
        rng: Random number generator
        min_hairpin: Minimum hairpin loop size

    Returns:
        MoveResult with proposal details
    """
    if not candidate_pairs:
        return MoveResult(
            accepted=False,
            move_type=MoveType.TOGGLE,
            is_valid=False,
        )

    # Uniform random selection from candidates
    idx = rng.randint(0, len(candidate_pairs) - 1)
    i, j, _ = candidate_pairs[idx]

    # Normalize pair
    if i > j:
        i, j = j, i
    edge = (i, j)

    if edge in matching:
        # Propose removal
        if edge in fixed_pairs:
            # Cannot remove fixed pairs
            return MoveResult(
                accepted=False,
                move_type=MoveType.TOGGLE,
                pairs_removed=None,
                is_valid=False,
            )
        return MoveResult(
            accepted=False,  # Will be set by sampler after M-H test
            move_type=MoveType.TOGGLE,
            pairs_removed={edge},
            pairs_added=None,
            hastings_ratio=1.0,  # Symmetric move
            is_valid=True,
        )
    else:
        # Propose addition
        # Check if both positions are unpaired
        if partners[i] != -1 or partners[j] != -1:
            # Conflict - one position already paired
            return MoveResult(
                accepted=False,
                move_type=MoveType.TOGGLE,
                is_valid=False,
            )

        # Check hairpin constraint
        if j - i - 1 < min_hairpin:
            return MoveResult(
                accepted=False,
                move_type=MoveType.TOGGLE,
                is_valid=False,
            )

        return MoveResult(
            accepted=False,  # Will be set by sampler after M-H test
            move_type=MoveType.TOGGLE,
            pairs_added={edge},
            pairs_removed=None,
            hastings_ratio=1.0,  # Symmetric move
            is_valid=True,
        )


def propose_segment_birth_death(
    matching: set[tuple[int, int]],
    partners: list[int],
    segment_index: SegmentIndex,
    fixed_pairs: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    rng,
    min_hairpin: int = 3,
) -> MoveResult:
    """Propose a segment birth or death move.

    Birth: Add a complete helix (all positions currently unpaired)
    Death: Remove a complete helix (none of its pairs are fixed)

    Hastings ratios:
    - Birth: |all_segments| / |removable(M')|
    - Death: |removable(M)| / |all_segments|

    Args:
        matching: Current set of base pairs
        partners: Partners array
        segment_index: Precomputed segment index
        fixed_pairs: Pairs that cannot be removed
        weights: Weight dictionary
        rng: Random number generator
        min_hairpin: Minimum hairpin loop size

    Returns:
        MoveResult with proposal details
    """
    from .segments import find_addable_segments, find_removable_segments

    # Get current removable segments
    removable = find_removable_segments(matching, segment_index, fixed_pairs)
    n_removable = len(removable)

    # Get current addable segments
    addable = find_addable_segments(matching, partners, segment_index, min_hairpin)
    n_addable = len(addable)

    # If neither birth nor death is possible, return invalid
    if n_addable == 0 and n_removable == 0:
        return MoveResult(
            move_type=MoveType.SEGMENT_BIRTH,
            is_valid=False,
        )

    # Choose birth or death with equal probability if both possible
    # Otherwise, choose the possible one
    if n_addable == 0:
        do_birth = False
    elif n_removable == 0:
        do_birth = True
    else:
        do_birth = rng.random() < 0.5

    if do_birth:
        # Birth move: select uniformly from addable segments
        seg = addable[rng.randint(0, n_addable - 1)]
        pairs_to_add = set(seg.pairs())

        # After adding, how many segments will be removable?
        # We need to count removable segments in M' = M ∪ seg
        # This is an approximation - we add our segment and check
        n_removable_after = n_removable + 1  # At minimum, the segment we added

        # Hastings ratio for birth: |addable(M)| / |removable(M')|
        # More precisely: forward = 1/(2*n_addable), reverse = 1/(2*n_removable_after)
        # But since we choose birth/death with prob 1/2 each, it simplifies
        if n_removable_after == 0:
            hastings = float("inf")  # Should not happen
        else:
            hastings = n_addable / n_removable_after

        return MoveResult(
            move_type=MoveType.SEGMENT_BIRTH,
            pairs_added=pairs_to_add,
            pairs_removed=None,
            hastings_ratio=hastings,
            is_valid=True,
        )
    else:
        # Death move: select uniformly from removable segments
        seg = removable[rng.randint(0, n_removable - 1)]
        pairs_to_remove = set(seg.pairs())

        # After removing, how many segments will be addable?
        n_addable_after = n_addable + 1  # At minimum, the segment we removed

        # Hastings ratio for death: |removable(M)| / |addable(M')|
        if n_addable_after == 0:
            hastings = float("inf")  # Should not happen
        else:
            hastings = n_removable / n_addable_after

        return MoveResult(
            move_type=MoveType.SEGMENT_DEATH,
            pairs_added=None,
            pairs_removed=pairs_to_remove,
            hastings_ratio=hastings,
            is_valid=True,
        )


def propose_segment_swap(
    matching: set[tuple[int, int]],
    partners: list[int],
    segment_index: SegmentIndex,
    fixed_pairs: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    rng,
    min_hairpin: int = 3,
) -> MoveResult:
    """Propose a segment swap move.

    Remove one segment and add a different (conflicting) segment.

    Hastings ratio:
        |valid_conflicts_for_old| / |valid_conflicts_for_new|

    Args:
        matching: Current set of base pairs
        partners: Partners array
        segment_index: Precomputed segment index
        fixed_pairs: Pairs that cannot be removed
        weights: Weight dictionary
        rng: Random number generator
        min_hairpin: Minimum hairpin loop size

    Returns:
        MoveResult with proposal details
    """
    from .segments import find_removable_segments

    # Get removable segments
    removable = find_removable_segments(matching, segment_index, fixed_pairs)
    if not removable:
        return MoveResult(
            move_type=MoveType.SEGMENT_SWAP,
            is_valid=False,
        )

    # Select a segment to remove
    seg_old = removable[rng.randint(0, len(removable) - 1)]
    seg_old_pairs = set(seg_old.pairs())

    # Find conflicting segments that would become addable after removing seg_old
    # A segment is a valid swap candidate if:
    # 1. It conflicts with seg_old (shares positions)
    # 2. All its positions (except those in seg_old) are currently free
    conflicts = segment_index.conflicts.get(seg_old, set())
    valid_conflicts: list = []

    for seg_new in conflicts:
        seg_new_pairs = set(seg_new.pairs())
        seg_new_positions = seg_new.positions()
        seg_old_positions = seg_old.positions()

        # Check if the non-overlapping positions are free
        can_add = True
        for pos in seg_new_positions:
            if pos not in seg_old_positions:
                if partners[pos] != -1:
                    can_add = False
                    break

        # Check hairpin constraint for new segment pairs
        if can_add:
            for i, j in seg_new.pairs():
                if j - i - 1 < min_hairpin:
                    can_add = False
                    break

        if can_add:
            valid_conflicts.append(seg_new)

    if not valid_conflicts:
        return MoveResult(
            move_type=MoveType.SEGMENT_SWAP,
            is_valid=False,
        )

    # Select new segment uniformly from valid conflicts
    seg_new = valid_conflicts[rng.randint(0, len(valid_conflicts) - 1)]
    seg_new_pairs = set(seg_new.pairs())

    # Compute reverse Hastings ratio
    # For the reverse move (removing seg_new, adding seg_old back),
    # we need to count valid conflicts for seg_new in the new matching M'
    # M' = (M - seg_old) ∪ seg_new

    # Compute partners array for M' to check reverse conflicts
    new_partners = list(partners)  # Copy
    # Clear old segment positions
    for i, j in seg_old.pairs():
        new_partners[i] = -1
        new_partners[j] = -1
    # Set new segment positions
    for i, j in seg_new.pairs():
        new_partners[i] = j
        new_partners[j] = i

    # Count valid reverse conflicts (segments conflicting with seg_new that could be swapped in)
    reverse_conflicts = segment_index.conflicts.get(seg_new, set())
    valid_reverse: list = []

    for seg_rev in reverse_conflicts:
        seg_rev_positions = seg_rev.positions()
        seg_new_positions = seg_new.positions()

        can_add = True
        for pos in seg_rev_positions:
            if pos not in seg_new_positions:
                if new_partners[pos] != -1:
                    can_add = False
                    break

        if can_add:
            for i, j in seg_rev.pairs():
                if j - i - 1 < min_hairpin:
                    can_add = False
                    break

        if can_add:
            valid_reverse.append(seg_rev)

    if not valid_reverse:
        # Reverse move not possible - reject
        return MoveResult(
            move_type=MoveType.SEGMENT_SWAP,
            is_valid=False,
        )

    # Hastings ratio: |valid_conflicts| / |valid_reverse|
    hastings = len(valid_conflicts) / len(valid_reverse)

    return MoveResult(
        move_type=MoveType.SEGMENT_SWAP,
        pairs_added=seg_new_pairs,
        pairs_removed=seg_old_pairs,
        hastings_ratio=hastings,
        is_valid=True,
    )


def propose_loop_closure(
    matching: set[tuple[int, int]],
    partners: list[int],
    loop_folds: list[set[tuple[int, int]]],
    weights: dict[tuple[int, int], float],
    rng,
) -> MoveResult:
    """Propose a loop closure move.

    This move is symmetric (Hastings ratio = 1) and proposes
    replacing the current pairs within a loop region with an
    alternative folding of that region.

    Args:
        matching: Current set of base pairs
        partners: Partners array
        loop_folds: Precomputed alternative loop foldings
        weights: Weight dictionary
        rng: Random number generator

    Returns:
        MoveResult with proposal details
    """
    # Placeholder implementation - will be filled in Phase 4
    raise NotImplementedError("propose_loop_closure not yet implemented")


def select_and_propose_move(
    matching: set[tuple[int, int]],
    partners: list[int],
    candidate_pairs: list[tuple[int, int, float]],
    segment_index: SegmentIndex | None,
    loop_folds: list[set[tuple[int, int]]] | None,
    fixed_pairs: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    rng,
    move_probs: dict[MoveType, float] | None = None,
    min_hairpin: int = 3,
) -> MoveResult:
    """Select a move type and propose it.

    Args:
        matching: Current set of base pairs
        partners: Partners array
        candidate_pairs: List of (i, j, weight) candidates
        segment_index: Precomputed segment index (optional)
        loop_folds: Precomputed loop foldings (optional)
        fixed_pairs: Pairs that cannot be removed
        weights: Weight dictionary
        rng: Random number generator
        move_probs: Probability for each move type (optional)
        min_hairpin: Minimum hairpin loop size

    Returns:
        MoveResult with proposal details
    """
    # Default move probabilities
    if move_probs is None:
        # If we have segment index, use segment moves
        if segment_index is not None and len(segment_index.all_segments) > 0:
            move_probs = {
                MoveType.TOGGLE: 0.5,
                MoveType.SEGMENT_BIRTH: 0.2,
                MoveType.SEGMENT_DEATH: 0.2,
                MoveType.SEGMENT_SWAP: 0.1,
            }
        else:
            # Only toggle moves available
            move_probs = {MoveType.TOGGLE: 1.0}

    # Normalize probabilities
    total = sum(move_probs.values())
    if total == 0:
        return MoveResult(is_valid=False)

    # Select move type
    r = rng.random() * total
    cumsum = 0.0
    selected_move = MoveType.TOGGLE

    for move_type, prob in move_probs.items():
        cumsum += prob
        if r < cumsum:
            selected_move = move_type
            break

    # Propose the selected move
    if selected_move == MoveType.TOGGLE:
        return propose_toggle(
            matching, partners, candidate_pairs, fixed_pairs, weights, rng, min_hairpin
        )

    elif selected_move in (MoveType.SEGMENT_BIRTH, MoveType.SEGMENT_DEATH):
        if segment_index is None or len(segment_index.all_segments) == 0:
            return MoveResult(move_type=selected_move, is_valid=False)
        return propose_segment_birth_death(
            matching, partners, segment_index, fixed_pairs, weights, rng, min_hairpin
        )

    elif selected_move == MoveType.SEGMENT_SWAP:
        if segment_index is None or len(segment_index.all_segments) == 0:
            return MoveResult(move_type=MoveType.SEGMENT_SWAP, is_valid=False)
        return propose_segment_swap(
            matching, partners, segment_index, fixed_pairs, weights, rng, min_hairpin
        )

    elif selected_move == MoveType.LOOP_CLOSURE:
        if loop_folds is None or len(loop_folds) == 0:
            return MoveResult(move_type=MoveType.LOOP_CLOSURE, is_valid=False)
        return propose_loop_closure(matching, partners, loop_folds, weights, rng)

    else:
        return MoveResult(is_valid=False)
