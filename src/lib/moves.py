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

    Returns:
        MoveResult with proposal details
    """
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("propose_segment_birth_death not yet implemented")


def propose_segment_swap(
    matching: set[tuple[int, int]],
    partners: list[int],
    segment_index: SegmentIndex,
    fixed_pairs: set[tuple[int, int]],
    weights: dict[tuple[int, int], float],
    rng,
) -> MoveResult:
    """Propose a segment swap move.

    Remove one segment and add a different (conflicting) segment.

    Hastings ratio:
        (|removable(M)| · |addable_after_A|) / (|removable(M')| · |addable_after_B|)

    Args:
        matching: Current set of base pairs
        partners: Partners array
        segment_index: Precomputed segment index
        fixed_pairs: Pairs that cannot be removed
        weights: Weight dictionary
        rng: Random number generator

    Returns:
        MoveResult with proposal details
    """
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("propose_segment_swap not yet implemented")


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
    # Placeholder implementation - will be filled in Phase 3
    raise NotImplementedError("select_and_propose_move not yet implemented")
