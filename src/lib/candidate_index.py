"""
Candidate edge indexing for fast MCMC proposals.

Provides:
- O(1) uniform sampling from addable/removable edge sets
- O(deg(i)+deg(j)) updates after accepting a move that changes partners[i/j]
"""

from __future__ import annotations

import random
from dataclasses import dataclass


class IndexedSet:
    """Set-like container with O(1) add/remove and O(1) uniform random sampling."""

    __slots__ = ("items", "pos")

    def __init__(self) -> None:
        self.items: list[int] = []
        self.pos: dict[int, int] = {}

    def __len__(self) -> int:  # pragma: no cover
        return len(self.items)

    def __contains__(self, x: int) -> bool:  # pragma: no cover
        return x in self.pos

    def add(self, x: int) -> None:
        if x in self.pos:
            return
        self.pos[x] = len(self.items)
        self.items.append(x)

    def discard(self, x: int) -> None:
        i = self.pos.get(x)
        if i is None:
            return
        last = self.items[-1]
        self.items[i] = last
        self.pos[last] = i
        self.items.pop()
        del self.pos[x]

    def random(self, rng: random.Random) -> int:
        return self.items[rng.randrange(len(self.items))]


@dataclass
class ToggleCounts:
    addable: int
    removable: int


class CandidateEdgeIndex:
    """
    Maintains dynamic addable/removable candidate edges for smart toggle proposals.

    Definitions:
    - addable: edge (i,j) not currently present where partners[i]=partners[j]=-1 and hairpin ok
    - removable: edge (i,j) currently present and not fixed
    """

    def __init__(
        self,
        length: int,
        candidate_pairs: list[tuple[int, int, float]],
        partners: list[int],
        fixed_pairs: set[tuple[int, int]] | None,
        min_hairpin: int,
    ) -> None:
        self.length = length
        self.min_hairpin = min_hairpin
        self.edges_i: list[int] = []
        self.edges_j: list[int] = []
        self.edge_to_idx: dict[tuple[int, int], int] = {}
        self.adj: list[list[int]] = [[] for _ in range(length)]
        self.fixed_pairs = fixed_pairs or set()
        self.addable = IndexedSet()
        self.removable = IndexedSet()
        self._mark: list[int] = []
        self._epoch = 1

        for idx, (i, j, _w) in enumerate(candidate_pairs):
            if i > j:
                i, j = j, i
            self.edges_i.append(i)
            self.edges_j.append(j)
            self.edge_to_idx[(i, j)] = idx
            if 0 <= i < length:
                self.adj[i].append(idx)
            if 0 <= j < length:
                self.adj[j].append(idx)

        self._mark = [0] * len(self.edges_i)
        self._initialize_sets(partners)

    def _initialize_sets(self, partners: list[int]) -> None:
        for idx in range(len(self.edges_i)):
            i = self.edges_i[idx]
            j = self.edges_j[idx]
            in_matching = partners[i] == j and partners[j] == i
            if in_matching:
                if (i, j) not in self.fixed_pairs:
                    self.removable.add(idx)
                continue
            if self._is_addable_idx(idx, partners):
                self.addable.add(idx)

    def _is_addable_idx(self, idx: int, partners: list[int]) -> bool:
        i = self.edges_i[idx]
        j = self.edges_j[idx]
        if partners[i] != -1 or partners[j] != -1:
            return False
        if (j - i - 1) < self.min_hairpin:
            return False
        return True

    def counts(self) -> ToggleCounts:
        return ToggleCounts(addable=len(self.addable.items), removable=len(self.removable.items))

    def _unique_affected_edge_indices(self, i: int, j: int) -> list[int]:
        self._epoch += 1
        if self._epoch >= 2_000_000_000:
            self._epoch = 1
            self._mark = [0] * len(self._mark)

        affected: list[int] = []
        for pos in (i, j):
            if not (0 <= pos < self.length):
                continue
            for eidx in self.adj[pos]:
                if self._mark[eidx] == self._epoch:
                    continue
                self._mark[eidx] = self._epoch
                affected.append(eidx)
        return affected

    @staticmethod
    def toggle_choice_probs(addable: int, removable: int) -> tuple[float, float]:
        if addable > 0 and removable > 0:
            return 0.5, 0.5
        if addable > 0:
            return 1.0, 0.0
        if removable > 0:
            return 0.0, 1.0
        return 0.0, 0.0

    def estimate_counts_after_add(self, edge_idx: int, partners: list[int]) -> ToggleCounts:
        i = self.edges_i[edge_idx]
        j = self.edges_j[edge_idx]
        addable_before = len(self.addable.items)
        removable_before = len(self.removable.items)

        affected = self._unique_affected_edge_indices(i, j)
        delta_addable = 0

        # Hypothetical partners after add
        # (only i and j change, so recompute locally)
        for aidx in affected:
            was = aidx in self.addable.pos
            # After adding, i and j become paired, so any edge touching i/j is not addable.
            now = False
            if aidx != edge_idx:
                ai = self.edges_i[aidx]
                aj = self.edges_j[aidx]
                # If the edge doesn't touch i/j, its addability won't change; but affected edges always touch.
                if ai != i and ai != j and aj != i and aj != j:
                    now = was
            if now != was:
                delta_addable += 1 if now else -1

        addable_after = max(0, addable_before + delta_addable)
        removable_after = removable_before + (0 if (i, j) in self.fixed_pairs else 1)
        return ToggleCounts(addable=addable_after, removable=removable_after)

    def estimate_counts_after_remove(self, edge_idx: int, partners: list[int]) -> ToggleCounts:
        i = self.edges_i[edge_idx]
        j = self.edges_j[edge_idx]
        addable_before = len(self.addable.items)
        removable_before = len(self.removable.items)

        affected = self._unique_affected_edge_indices(i, j)
        delta_addable = 0

        for aidx in affected:
            was = aidx in self.addable.pos
            ai = self.edges_i[aidx]
            aj = self.edges_j[aidx]

            # Hypothetical partners after remove: i and j become free.
            # Edge is addable if both endpoints are free under that hypothetical.
            free_ai = partners[ai] == -1 or ai in (i, j)
            free_aj = partners[aj] == -1 or aj in (i, j)
            in_matching_after = False
            if (ai, aj) == (i, j):
                in_matching_after = False
            else:
                in_matching_after = partners[ai] == aj and partners[aj] == ai

            now = (
                (not in_matching_after)
                and free_ai
                and free_aj
                and (aj - ai - 1) >= self.min_hairpin
            )
            if now != was:
                delta_addable += 1 if now else -1

        addable_after = max(0, addable_before + delta_addable)
        removable_after = max(0, removable_before - 1)
        return ToggleCounts(addable=addable_after, removable=removable_after)

    def apply_add_pair(self, i: int, j: int, partners: list[int]) -> None:
        if i > j:
            i, j = j, i
        idx = self.edge_to_idx.get((i, j))
        if idx is None:
            return

        # Edge becomes removable if not fixed.
        self.addable.discard(idx)
        if (i, j) not in self.fixed_pairs:
            self.removable.add(idx)

        # Any edge incident to i or j cannot be addable now.
        affected = self._unique_affected_edge_indices(i, j)
        for eidx in affected:
            if eidx in self.addable.pos:
                self.addable.discard(eidx)

    def apply_remove_pair(self, i: int, j: int, partners: list[int]) -> None:
        if i > j:
            i, j = j, i
        idx = self.edge_to_idx.get((i, j))
        if idx is None:
            return

        self.removable.discard(idx)

        # Edges incident to i or j may become addable now.
        affected = self._unique_affected_edge_indices(i, j)
        for eidx in affected:
            if eidx == idx:
                # The removed edge itself can become addable now.
                if self._is_addable_idx(eidx, partners):
                    self.addable.add(eidx)
                else:
                    self.addable.discard(eidx)
                continue
            if self._is_addable_idx(eidx, partners):
                self.addable.add(eidx)
            else:
                self.addable.discard(eidx)
