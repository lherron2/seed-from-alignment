"""Regression tests to ensure output format is preserved.

These tests verify that the refactored sampler produces output
in the same format as the original implementation.
"""

import pytest


class TestOutputFormat:
    """Tests for output format preservation."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 4")
    def test_db_file_format(self) -> None:
        """
        Output .db file format matches expected structure:
        - Line 1: ungapped sequence
        - Lines 2+: one PK-annotated dot-bracket string per sample
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 4")
    def test_bracket_notation(self) -> None:
        """
        Bracket notation uses Rosetta-safe hierarchy:
        - Layer 0: ()
        - Layer 1: []
        - Layer 2: {}
        - Layers 3+: aA, bB, cC, ...
        """
        pass


class TestDeterminism:
    """Tests for deterministic output."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 4")
    def test_same_seed_same_output(self) -> None:
        """Same seed produces byte-identical output."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 4")
    def test_different_seeds_different_output(self) -> None:
        """Different seeds produce different output (statistical sanity)."""
        pass


class TestBackwardCompatibility:
    """Tests for backward compatibility with existing code."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 6")
    def test_legacy_interface(self) -> None:
        """Legacy sample_matchings() interface still works."""
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 6")
    def test_cli_arguments(self) -> None:
        """All existing CLI arguments are supported."""
        pass
