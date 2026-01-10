import textwrap
from pathlib import Path

import pytest

from src.lib import refine_unpaired_regions as rup


def test_read_fasta_sequence_first_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(textwrap.dedent(">seq1\nACGUACGU\n>seq2 other stuff\nGGGG\n"))
    seq = rup.read_fasta_sequence(fasta)
    assert seq == "ACGUACGU"


def test_read_fasta_sequence_named(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(textwrap.dedent(">seq1\nACGU\n>seq2 my_favorite\nGGGGCCCC\n"))
    seq = rup.read_fasta_sequence(fasta, seq_name="seq2")
    assert seq == "GGGGCCCC"


def test_read_fasta_sequence_missing_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">seq1\nACGU\n")
    with pytest.raises(ValueError):
        rup.read_fasta_sequence(fasta, seq_name="nope")


def test_read_db_structures_basic(tmp_path: Path) -> None:
    db = tmp_path / "structs.db"
    db.write_text(textwrap.dedent("....((..))...\n\n....((..))...\n"))
    structs = rup.read_db_structures(db)
    assert structs == ["....((..))...", "....((..))..."]


def test_read_db_structures_different_lengths_accepted(tmp_path: Path) -> None:
    # Current implementation doesn't enforce equal structure lengths
    db = tmp_path / "structs_mixed.db"
    db.write_text("....\n.....\n")
    structs = rup.read_db_structures(db)
    # Both lines are structures (dots are structure chars), no length validation
    assert structs == ["....", "....."]


def test_find_unpaired_runs_returns_expected_ranges() -> None:
    struct = "....((..))..."
    runs = rup.find_unpaired_runs(struct, min_len=3)
    # 5' run: [0,4), 3' run: [10,13)
    assert runs == [(0, 4), (10, 13)]


def test_find_unpaired_runs_none_when_below_threshold() -> None:
    struct = "..((..)).."
    runs = rup.find_unpaired_runs(struct, min_len=3)
    assert runs == []


# Tests for functions that exist in the current implementation


def test_pairs_cross_basic() -> None:
    # Non-crossing pairs
    assert rup.pairs_cross((0, 5), (1, 4)) is False  # nested
    assert rup.pairs_cross((0, 3), (4, 7)) is False  # disjoint
    # Crossing pairs (pseudoknot)
    assert rup.pairs_cross((0, 5), (2, 7)) is True  # i < k < j < l


def test_is_canonical_pairs() -> None:
    assert rup.is_canonical("A", "U") is True
    assert rup.is_canonical("G", "C") is True
    assert rup.is_canonical("G", "U") is True  # wobble
    assert rup.is_canonical("A", "A") is False
    assert rup.is_canonical("C", "U") is False
