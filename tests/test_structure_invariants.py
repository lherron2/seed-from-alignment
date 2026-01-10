"""Tests for RNA structure validation."""

from src.lib.validate_structure import (
    pairs_to_structure,
    parse_structure_to_pairs,
    validate_structure,
)


class TestParseStructureToPairs:
    """Tests for parse_structure_to_pairs()."""

    def test_empty_structure(self) -> None:
        """Empty structure has no pairs."""
        assert parse_structure_to_pairs("....") == []

    def test_simple_hairpin(self) -> None:
        """Simple hairpin structure."""
        pairs = parse_structure_to_pairs("((...))")
        assert pairs == [(0, 6), (1, 5)]

    def test_pseudoknot(self) -> None:
        """Pseudoknot with brackets and square brackets."""
        pairs = parse_structure_to_pairs("((..[..))..]..")
        # (0,8), (1,7), [4,11]
        assert (0, 8) in pairs
        assert (1, 7) in pairs
        assert (4, 11) in pairs


class TestPairsToStructure:
    """Tests for pairs_to_structure()."""

    def test_empty_pairs(self) -> None:
        """No pairs gives dots."""
        struct = pairs_to_structure([], 5)
        assert struct == "....."

    def test_simple_pairs(self) -> None:
        """Simple nested pairs."""
        struct = pairs_to_structure([(0, 6), (1, 5)], 7)
        assert struct == "((...))"

    def test_pseudoknot_layering(self) -> None:
        """Crossing pairs use different bracket types."""
        # (0,5) and (2,8) cross
        struct = pairs_to_structure([(0, 5), (2, 8)], 10)
        # Should use () for one layer and [] for another
        assert "(" in struct
        assert "[" in struct or "{" in struct


class TestValidateStructure:
    """Tests for validate_structure()."""

    def test_valid_simple_structure(self) -> None:
        """Valid simple structure passes."""
        seq = "GGGGAAAACCCC"
        struct = "((((....))))"
        is_valid, errors = validate_structure(seq, struct)
        assert is_valid is True
        assert errors == []

    def test_length_mismatch(self) -> None:
        """Length mismatch is detected."""
        is_valid, errors = validate_structure("AAAA", "((...))")
        assert is_valid is False
        assert any("Length" in e for e in errors)

    def test_unbalanced_brackets(self) -> None:
        """Unbalanced brackets are detected."""
        is_valid, errors = validate_structure("AAAAA", "((..)")
        assert is_valid is False
        assert any("Unmatched" in e for e in errors)

    def test_hairpin_too_small(self) -> None:
        """Small hairpin is detected."""
        seq = "GGCC"
        struct = "(())"  # Only 0 unpaired bases in loop
        is_valid, errors = validate_structure(seq, struct, min_hairpin=3)
        assert is_valid is False
        assert any("Hairpin" in e or "small" in e.lower() for e in errors)

    def test_hairpin_allowed_when_large_enough(self) -> None:
        """Hairpin passes when large enough."""
        seq = "GGAAACC"
        struct = "((...)))"
        # This has 3 unpaired in loop, but structure is invalid
        # Let me use a valid structure
        seq = "GGAAACC"
        struct = "((...))"
        is_valid, errors = validate_structure(seq, struct, min_hairpin=3)
        assert is_valid is True

    def test_lonely_pair_detected(self) -> None:
        """Lonely (isolated) pairs are detected."""
        seq = "G...C"  # 5 characters
        struct = "(...)"  # No stacking partner - hairpin with 3 unpaired
        is_valid, errors = validate_structure(seq, struct, allow_lonely_pairs=False, min_hairpin=3)
        assert is_valid is False
        assert any("Lonely" in e for e in errors)

    def test_lonely_pair_allowed(self) -> None:
        """Lonely pairs pass when allowed."""
        seq = "G...C"  # 5 characters
        struct = "(...)"
        is_valid, errors = validate_structure(seq, struct, allow_lonely_pairs=True, min_hairpin=3)
        assert is_valid is True

    def test_canonical_only(self) -> None:
        """Non-canonical pairs detected when canonical_only=True."""
        seq = "A...A"  # 5 characters, A-A not canonical
        struct = "(...)"
        is_valid, errors = validate_structure(
            seq, struct, canonical_only=True, allow_lonely_pairs=True, min_hairpin=3
        )
        assert is_valid is False
        assert any("Non-canonical" in e for e in errors)

    def test_canonical_pair_passes(self) -> None:
        """Canonical pairs pass."""
        seq = "A...U"  # 5 characters, A-U is canonical
        struct = "(...)"
        is_valid, errors = validate_structure(
            seq, struct, canonical_only=True, allow_lonely_pairs=True, min_hairpin=3
        )
        assert is_valid is True

    def test_pseudoknot_disallowed(self) -> None:
        """Pseudoknots detected when disallowed."""
        # (0,4) and [2,6] cross: 0 < 2 < 4 < 6 is a pseudoknot
        seq = "AAAAAAA"
        struct = "(.[(]).."  # (0,4) [2,6] - this is wrong, let me fix
        # Actually need: (.[.).]
        struct = "(.[.).]"
        is_valid, errors = validate_structure(
            seq, struct, allow_pk=False, allow_lonely_pairs=True, min_hairpin=0
        )
        assert is_valid is False
        assert any("Pseudoknot" in e for e in errors)

    def test_pseudoknot_allowed(self) -> None:
        """Pseudoknots pass when allowed."""
        # Same structure, but PK allowed
        seq = "AAAAAAA"
        struct = "(.[.).]"
        is_valid, errors = validate_structure(
            seq, struct, allow_pk=True, allow_lonely_pairs=True, min_hairpin=0
        )
        assert is_valid is True

    def test_lowercase_forced_unpaired(self) -> None:
        """Lowercase bases must be unpaired when flag is set."""
        seq = "GGaaaCC"  # 7 chars, lowercase in unpaired region
        struct = "((...))"  # 7 chars
        is_valid, errors = validate_structure(
            seq, struct, lowercase_forced_unpaired=True, min_hairpin=3
        )
        assert is_valid is True

        # Now pair a lowercase base
        seq = "a...A"  # 5 chars, lowercase 'a' at position 0 will be paired with A at 4
        struct = "(...)"  # 5 chars, (0,4) pair
        is_valid, errors = validate_structure(
            seq, struct, lowercase_forced_unpaired=True, allow_lonely_pairs=True, min_hairpin=3
        )
        assert is_valid is False
        assert any("lowercase" in e.lower() for e in errors)
