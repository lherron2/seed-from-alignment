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

from .candidate_index import CandidateEdgeIndex
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
    candidate_index: CandidateEdgeIndex | None = None,
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
    if candidate_index is not None:
        counts_before = candidate_index.counts()
        p_add, p_rem = candidate_index.toggle_choice_probs(
            counts_before.addable, counts_before.removable
        )
        if p_add == 0.0 and p_rem == 0.0:
            return MoveResult(
                accepted=False,
                move_type=MoveType.TOGGLE,
                is_valid=False,
            )

        do_add = p_add > 0 and (p_rem == 0 or rng.random() < (p_add / (p_add + p_rem)))
        if do_add:
            if counts_before.addable <= 0:
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)
            edge_idx = candidate_index.addable.random(rng)
            i = candidate_index.edges_i[edge_idx]
            j = candidate_index.edges_j[edge_idx]
            if partners[i] != -1 or partners[j] != -1:
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)
            if j - i - 1 < min_hairpin:
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)

            counts_after = candidate_index.estimate_counts_after_add(edge_idx, partners)
            p_add_after, p_rem_after = candidate_index.toggle_choice_probs(
                counts_after.addable, counts_after.removable
            )
            if counts_after.removable <= 0 or p_rem_after <= 0 or p_add <= 0:
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)

            q_fwd = p_add / float(counts_before.addable)
            q_rev = p_rem_after / float(counts_after.removable)
            hastings = q_rev / q_fwd

            return MoveResult(
                accepted=False,
                move_type=MoveType.TOGGLE,
                pairs_added={(i, j)},
                pairs_removed=None,
                hastings_ratio=hastings,
                is_valid=True,
            )
        else:
            if counts_before.removable <= 0:
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)
            edge_idx = candidate_index.removable.random(rng)
            i = candidate_index.edges_i[edge_idx]
            j = candidate_index.edges_j[edge_idx]
            edge = (i, j)
            if edge in fixed_pairs:
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)
            if not (partners[i] == j and partners[j] == i):
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)

            counts_after = candidate_index.estimate_counts_after_remove(edge_idx, partners)
            p_add_after, p_rem_after = candidate_index.toggle_choice_probs(
                counts_after.addable, counts_after.removable
            )
            if (
                counts_before.removable <= 0
                or p_rem <= 0
                or counts_after.addable <= 0
                or p_add_after <= 0
            ):
                return MoveResult(move_type=MoveType.TOGGLE, is_valid=False)

            q_fwd = p_rem / float(counts_before.removable)
            q_rev = p_add_after / float(counts_after.addable)
            hastings = q_rev / q_fwd

            return MoveResult(
                accepted=False,
                move_type=MoveType.TOGGLE,
                pairs_removed={edge},
                pairs_added=None,
                hastings_ratio=hastings,
                is_valid=True,
            )

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
    force_move: MoveType | None = None,
    fast_uniform_pool: bool = False,
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
        force_move: If set, force birth or death move type

    Returns:
        MoveResult with proposal details
    """
    from .segments import find_addable_segments, find_removable_segments

    segments = segment_index.all_segments_list
    if not segments:
        return MoveResult(move_type=MoveType.SEGMENT_BIRTH, is_valid=False)

    if fast_uniform_pool:
        # Uniform-from-global-pool proposal:
        # q selects a segment uniformly from all segments, then proposes add/remove if valid.
        # This avoids O(#segments) scans and keeps detailed balance (Hastings within type = 1).
        seg = segments[rng.randrange(len(segments))]
        seg_pairs = seg.pairs()
        pairs_set = set(seg_pairs)

        want_birth = force_move == MoveType.SEGMENT_BIRTH
        want_death = force_move == MoveType.SEGMENT_DEATH
        if force_move is None:
            # Default: decide based on whether the segment is currently present.
            # If present and removable -> death; else attempt birth.
            want_death = all(p in matching for p in seg_pairs)
            want_birth = not want_death

        if want_birth:
            # Add segment if all positions are free and hairpins OK.
            for i, j in seg_pairs:
                if partners[i] != -1 or partners[j] != -1:
                    return MoveResult(move_type=MoveType.SEGMENT_BIRTH, is_valid=False)
                if j - i - 1 < min_hairpin:
                    return MoveResult(move_type=MoveType.SEGMENT_BIRTH, is_valid=False)
            return MoveResult(
                move_type=MoveType.SEGMENT_BIRTH,
                pairs_added=pairs_set,
                pairs_removed=None,
                hastings_ratio=1.0,
                is_valid=True,
            )

        if want_death:
            # Remove segment if fully present and none fixed.
            if not all(p in matching for p in seg_pairs):
                return MoveResult(move_type=MoveType.SEGMENT_DEATH, is_valid=False)
            if any(p in fixed_pairs for p in seg_pairs):
                return MoveResult(move_type=MoveType.SEGMENT_DEATH, is_valid=False)
            return MoveResult(
                move_type=MoveType.SEGMENT_DEATH,
                pairs_added=None,
                pairs_removed=pairs_set,
                hastings_ratio=1.0,
                is_valid=True,
            )

        return MoveResult(move_type=MoveType.SEGMENT_BIRTH, is_valid=False)

    # Slow exact mode (legacy): compute counts of addable/removable segments.
    # Get current removable segments
    removable = find_removable_segments(matching, segment_index, fixed_pairs)
    n_removable = len(removable)

    # Get current addable segments
    addable = find_addable_segments(matching, partners, segment_index, min_hairpin)
    n_addable = len(addable)

    # If neither birth nor death is possible, return invalid
    if n_addable == 0 and n_removable == 0:
        return MoveResult(move_type=MoveType.SEGMENT_BIRTH, is_valid=False)

    # Choose birth or death, optionally forcing a specific type
    if force_move == MoveType.SEGMENT_BIRTH:
        do_birth = True
    elif force_move == MoveType.SEGMENT_DEATH:
        do_birth = False
    else:
        # Choose birth or death with equal probability if both possible
        # Otherwise, choose the possible one
        if n_addable == 0:
            do_birth = False
        elif n_removable == 0:
            do_birth = True
        else:
            do_birth = rng.random() < 0.5

    if do_birth:
        if n_addable == 0:
            return MoveResult(
                move_type=MoveType.SEGMENT_BIRTH,
                is_valid=False,
            )
        # Birth move: select uniformly from addable segments
        seg = addable[rng.randint(0, n_addable - 1)]
        pairs_to_add = set(seg.pairs())

        # After adding, count removable segments in M' = M ∪ seg
        matching_after = matching | pairs_to_add
        removable_after = find_removable_segments(
            matching_after,
            segment_index,
            fixed_pairs,
        )
        n_removable_after = len(removable_after)

        # Hastings ratio for birth: |addable(M)| / |removable(M')|
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
        if n_removable == 0:
            return MoveResult(
                move_type=MoveType.SEGMENT_DEATH,
                is_valid=False,
            )
        # Death move: select uniformly from removable segments
        seg = removable[rng.randint(0, n_removable - 1)]
        pairs_to_remove = set(seg.pairs())

        # After removing, count addable segments in M' = M \ seg
        matching_after = matching - pairs_to_remove
        new_partners = list(partners)
        for i, j in pairs_to_remove:
            new_partners[i] = -1
            new_partners[j] = -1
        addable_after = find_addable_segments(
            matching_after,
            new_partners,
            segment_index,
            min_hairpin,
        )
        n_addable_after = len(addable_after)

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
    n_removable = len(removable)
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

    # Count removable segments in M' for the reverse move
    matching_after = (matching - seg_old_pairs) | seg_new_pairs
    removable_after = find_removable_segments(matching_after, segment_index, fixed_pairs)
    n_removable_after = len(removable_after)
    if n_removable_after == 0:
        return MoveResult(
            move_type=MoveType.SEGMENT_SWAP,
            is_valid=False,
        )

    # Hastings ratio: (|removable(M)| / |removable(M')|)
    #                   * (|valid_conflicts| / |valid_reverse|)
    hastings = (n_removable / n_removable_after) * (len(valid_conflicts) / len(valid_reverse))

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
    min_hairpin: int = 3,
    fixed_pairs: set[tuple[int, int]] | None = None,
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
        min_hairpin: Minimum hairpin loop size
        fixed_pairs: Pairs that cannot be removed

    Returns:
        MoveResult with proposal details
    """
    if not loop_folds:
        return MoveResult(
            move_type=MoveType.LOOP_CLOSURE,
            is_valid=False,
        )

    if fixed_pairs is None:
        fixed_pairs = set()

    # Select uniformly from all precomputed loop folds.
    fold_pairs = set(loop_folds[rng.randint(0, len(loop_folds) - 1)])

    if fold_pairs.issubset(matching):
        # Propose removal of the fold.
        if any(pair in fixed_pairs for pair in fold_pairs):
            return MoveResult(
                move_type=MoveType.LOOP_CLOSURE,
                is_valid=False,
            )
        return MoveResult(
            move_type=MoveType.LOOP_CLOSURE,
            pairs_removed=fold_pairs,
            pairs_added=None,
            hastings_ratio=1.0,
            is_valid=True,
        )

    if fold_pairs & matching:
        # Overlap would break reversibility (partial fold).
        return MoveResult(
            move_type=MoveType.LOOP_CLOSURE,
            is_valid=False,
        )

    # Addition: ensure all positions are free and hairpin constraints are met.
    for i, j in fold_pairs:
        if partners[i] != -1 or partners[j] != -1:
            return MoveResult(
                move_type=MoveType.LOOP_CLOSURE,
                is_valid=False,
            )
        if j - i - 1 < min_hairpin:
            return MoveResult(
                move_type=MoveType.LOOP_CLOSURE,
                is_valid=False,
            )

    return MoveResult(
        move_type=MoveType.LOOP_CLOSURE,
        pairs_added=fold_pairs,
        pairs_removed=None,
        hastings_ratio=1.0,
        is_valid=True,
    )


def select_and_propose_move(
    matching: set[tuple[int, int]],
    partners: list[int],
    candidate_pairs: list[tuple[int, int, float]],
    segment_index: SegmentIndex | None,
    loop_folds: list[set[tuple[int, int]]] | None,
    fixed_pairs: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    rng,
    candidate_index: CandidateEdgeIndex | None = None,
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
            matching,
            partners,
            candidate_pairs,
            fixed_pairs,
            weights,
            rng,
            min_hairpin,
            candidate_index=candidate_index,
        )

    elif selected_move in (MoveType.SEGMENT_BIRTH, MoveType.SEGMENT_DEATH):
        if segment_index is None or len(segment_index.all_segments) == 0:
            return MoveResult(move_type=selected_move, is_valid=False)
        move = propose_segment_birth_death(
            matching,
            partners,
            segment_index,
            fixed_pairs,
            weights,
            rng,
            min_hairpin,
            force_move=selected_move,
            fast_uniform_pool=True,
        )
        # Include move selection probability ratio if birth/death probs differ.
        if move.is_valid and move_probs is not None:
            p_birth = move_probs.get(MoveType.SEGMENT_BIRTH, 0.0)
            p_death = move_probs.get(MoveType.SEGMENT_DEATH, 0.0)
            if move.move_type == MoveType.SEGMENT_BIRTH:
                if p_birth <= 0 or p_death <= 0:
                    move.is_valid = False
                else:
                    move.hastings_ratio *= p_death / p_birth
            elif move.move_type == MoveType.SEGMENT_DEATH:
                if p_birth <= 0 or p_death <= 0:
                    move.is_valid = False
                else:
                    move.hastings_ratio *= p_birth / p_death
        return move

    elif selected_move == MoveType.SEGMENT_SWAP:
        if segment_index is None or len(segment_index.all_segments) == 0:
            return MoveResult(move_type=MoveType.SEGMENT_SWAP, is_valid=False)
        return propose_segment_swap(
            matching, partners, segment_index, fixed_pairs, weights, rng, min_hairpin
        )

    elif selected_move == MoveType.LOOP_CLOSURE:
        if loop_folds is None or len(loop_folds) == 0:
            return MoveResult(move_type=MoveType.LOOP_CLOSURE, is_valid=False)
        return propose_loop_closure(
            matching,
            partners,
            loop_folds,
            weights,
            rng,
            min_hairpin,
            fixed_pairs=fixed_pairs,
        )

    else:
        return MoveResult(is_valid=False)
