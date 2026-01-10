"""Tests for pseudoknot computation module."""

from src.lib.pk import (
    build_pk_cache,
    count_crossings,
    crosses,
    crossing_partners,
    pairs_to_layers,
)


class TestCrosses:
    """Tests for the crosses() function."""

    def test_nested_pairs_dont_cross(self) -> None:
        """Nested pairs should not cross."""
        assert crosses((0, 10), (2, 8)) is False
        assert crosses((2, 8), (0, 10)) is False

    def test_disjoint_pairs_dont_cross(self) -> None:
        """Disjoint pairs should not cross."""
        assert crosses((0, 5), (6, 10)) is False
        assert crosses((6, 10), (0, 5)) is False

    def test_crossing_pairs_detected(self) -> None:
        """Crossing pairs (pseudoknots) should be detected."""
        # i < k < j < l pattern
        assert crosses((0, 5), (2, 8)) is True
        assert crosses((2, 8), (0, 5)) is True

    def test_adjacent_pairs_dont_cross(self) -> None:
        """Adjacent pairs should not cross."""
        assert crosses((0, 5), (5, 10)) is False

    def test_same_pair_doesnt_cross_itself(self) -> None:
        """A pair should not cross itself."""
        assert crosses((0, 5), (0, 5)) is False


class TestCrossingPartners:
    """Tests for crossing_partners() function."""

    def test_empty_matching(self) -> None:
        """Empty matching has no crossing partners."""
        result = crossing_partners((0, 5), set())
        assert result == set()

    def test_no_crossings(self) -> None:
        """Matching with no crossings."""
        matching = {(0, 10), (2, 8), (12, 15)}
        result = crossing_partners((0, 10), matching)
        assert result == set()

    def test_single_crossing(self) -> None:
        """Single crossing pair."""
        matching = {(0, 5), (2, 8)}
        result = crossing_partners((0, 5), matching)
        assert result == {(2, 8)}


class TestCountCrossings:
    """Tests for count_crossings() function."""

    def test_empty_matching(self) -> None:
        """Empty matching has zero crossings."""
        assert count_crossings(set()) == 0

    def test_single_pair(self) -> None:
        """Single pair has zero crossings."""
        assert count_crossings({(0, 5)}) == 0

    def test_nested_pairs(self) -> None:
        """Nested pairs have zero crossings."""
        assert count_crossings({(0, 10), (2, 8)}) == 0

    def test_one_crossing(self) -> None:
        """One crossing pair."""
        assert count_crossings({(0, 5), (2, 8)}) == 1


class TestPKCache:
    """Tests for PKCache and build_pk_cache()."""

    def test_empty_matching(self) -> None:
        """Empty matching produces empty cache."""
        cache = build_pk_cache(set())
        assert cache.pk_pairs == set()
        assert cache.total_crossings == 0
        assert cache.max_crossings == 0

    def test_non_crossing_matching(self) -> None:
        """Non-crossing matching has empty pk_pairs."""
        matching = {(0, 10), (2, 8), (12, 15)}
        cache = build_pk_cache(matching)
        assert cache.pk_pairs == set()
        assert cache.total_crossings == 0

    def test_crossing_matching(self) -> None:
        """Crossing matching identifies PK pairs."""
        matching = {(0, 5), (2, 8)}
        cache = build_pk_cache(matching)
        assert cache.pk_pairs == matching
        assert cache.total_crossings == 1
        assert cache.max_crossings == 1


class TestPairsToLayers:
    """Tests for pairs_to_layers() function."""

    def test_empty(self) -> None:
        """Empty pairs give empty layers."""
        assert pairs_to_layers([]) == []

    def test_single_pair(self) -> None:
        """Single pair is one layer."""
        layers = pairs_to_layers([(0, 5)])
        assert len(layers) == 1
        assert (0, 5) in layers[0]

    def test_nested_pairs_same_layer(self) -> None:
        """Nested pairs go in same layer."""
        layers = pairs_to_layers([(0, 10), (2, 8)])
        assert len(layers) == 1

    def test_crossing_pairs_different_layers(self) -> None:
        """Crossing pairs go in different layers."""
        layers = pairs_to_layers([(0, 5), (2, 8)])
        assert len(layers) == 2
