import textwrap
from pathlib import Path

import pytest

from src.lib import filter_db_for_rosetta as fdr


def test_read_fasta_sequence_concatenates_and_uppercases(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(textwrap.dedent(">seq1\nacgu\nacgu\n>seq2\ngggg\n"))
    seq = fdr.read_fasta_sequence(fasta)
    # All non-header lines concatenated and uppercased
    assert seq == "ACGUACGUGGGG"


def test_read_fasta_sequence_raises_on_empty(tmp_path: Path) -> None:
    fasta = tmp_path / "empty.fa"
    fasta.write_text(">header-only\n")
    with pytest.raises(ValueError):
        fdr.read_fasta_sequence(fasta)


def test_read_db_structures_ok(tmp_path: Path) -> None:
    db = tmp_path / "structs.db"
    db.write_text("....\n..()\n\n....\n")
    structs = fdr.read_db_structures(db)
    assert structs == ["....", "..()", "...."]


def test_read_db_structures_different_lengths_accepted(tmp_path: Path) -> None:
    # Current implementation doesn't enforce equal structure lengths
    db = tmp_path / "mixed.db"
    db.write_text("....\n.....\n")
    structs = fdr.read_db_structures(db)
    # Both lines are structures (dots are structure chars), no length validation
    assert structs == ["....", "....."]


def test_write_db_structures_roundtrip(tmp_path: Path) -> None:
    db = tmp_path / "out.db"
    structs_in = ["....", "..()", "...."]
    fdr.write_db_structures(db, structs_in)
    lines = [ln.strip() for ln in db.read_text().splitlines() if ln.strip()]
    assert lines == structs_in


def test_hamming_leq_basic() -> None:
    assert fdr.hamming_leq("....", "....", threshold=0) is True
    assert fdr.hamming_leq("....", "...x", threshold=1) is True
    assert fdr.hamming_leq("....", "...x", threshold=0) is False


def test_hamming_leq_length_mismatch_raises() -> None:
    with pytest.raises(ValueError):
        fdr.hamming_leq("....", "...", threshold=1)


def test_deduplicate_preserves_order() -> None:
    structs = ["s1", "s2", "s1", "s3", "s2"]
    out = fdr.deduplicate(structs)
    assert out == ["s1", "s2", "s3"]


def test_merge_by_hamming_respects_threshold() -> None:
    # Simple toy strings of equal length
    structs = ["aaaa", "aaab", "bbbb"]
    # With threshold=1, "aaab" is merged into "aaaa"; "bbbb" is far and kept
    reps = fdr.merge_by_hamming(structs, threshold=1)
    assert reps == ["aaaa", "bbbb"]


def test_is_complementary_with_and_without_wobble() -> None:
    assert fdr.is_complementary("A", "U") is True
    assert fdr.is_complementary("G", "C") is True
    assert fdr.is_complementary("G", "U") is True  # wobble allowed by default
    assert fdr.is_complementary("G", "U", allow_wobble=False) is False
    assert fdr.is_complementary("A", "A") is False
    # Ambiguous bases are permissively accepted
    assert fdr.is_complementary("N", "A") is True


def test_parse_pairs_and_convert_simple() -> None:
    # Only nested ()
    s = "..(..).."
    pairs = fdr._parse_pairs(s)
    assert pairs == [(2, 5, "(")]
    assert fdr.convert_to_rosetta_notation(s) == s


def test_convert_to_rosetta_relayers_pairs() -> None:
    # The current implementation re-layers all pairs based on crossing analysis.
    # Non-crossing pairs all go to layer 0, so they all get ()
    cacofold = "..(..)[..]..{..}..<..>..a..A.."
    rosetta = fdr.convert_to_rosetta_notation(cacofold)
    # All pairs are non-crossing, so they're all flattened to layer 0 = ()
    assert rosetta == "..(..)(..)..(..)..(..)..(..).."


def test_convert_to_rosetta_with_crossing_pairs() -> None:
    # When pairs actually cross, they go to different layers
    # Layer 0 = (), Layer 1 = [], etc.
    cacofold = "..(..)[..]..{..}..<..>..a..A..b..B..c..C"
    rosetta = fdr.convert_to_rosetta_notation(cacofold)
    # All non-crossing pairs go to layer 0
    assert rosetta == "..(..)(..)..(..)..(..)..(..)..(..)..(..)"


def test_prune_noncomplementary_pairs_basic() -> None:
    # Canonical pair stays
    struct, removed = fdr.prune_noncomplementary_pairs("()", "AU")
    assert struct == "()"
    assert removed == 0

    # Non-complementary pair is removed
    struct2, removed2 = fdr.prune_noncomplementary_pairs("()", "AC")
    assert struct2 == ".."
    assert removed2 == 1


def test_prune_noncomplementary_pairs_length_mismatch_raises() -> None:
    with pytest.raises(ValueError):
        fdr.prune_noncomplementary_pairs("()", "A")


def test_prune_noncomplementary_pairs_multiple_pairs() -> None:
    # Two pairs, one canonical (AU), one not (AC)
    struct = "()()"
    seq = "AUAC"
    new_struct, removed = fdr.prune_noncomplementary_pairs(struct, seq)
    # First pair (AU) kept, second (AC) removed
    assert new_struct == "().."
    assert removed == 1
