import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from src.lib import refine_unpaired_regions as rup


def test_read_fasta_sequence_first_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(
        textwrap.dedent(
            ">seq1\n"
            "ACGUACGU\n"
            ">seq2 other stuff\n"
            "GGGG\n"
        )
    )
    seq = rup.read_fasta_sequence(fasta)
    assert seq == "ACGUACGU"


def test_read_fasta_sequence_named(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(
        textwrap.dedent(
            ">seq1\n"
            "ACGU\n"
            ">seq2 my_favorite\n"
            "GGGGCCCC\n"
        )
    )
    seq = rup.read_fasta_sequence(fasta, seq_name="seq2")
    assert seq == "GGGGCCCC"


def test_read_fasta_sequence_missing_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">seq1\nACGU\n")
    with pytest.raises(ValueError):
        rup.read_fasta_sequence(fasta, seq_name="nope")


def test_read_db_structures_basic(tmp_path: Path) -> None:
    db = tmp_path / "structs.db"
    db.write_text(
        textwrap.dedent(
            "....((..))...\n"
            "\n"
            "....((..))...\n"
        )
    )
    structs = rup.read_db_structures(db)
    assert structs == ["....((..))...", "....((..))..."]


def test_read_db_structures_length_mismatch_raises(tmp_path: Path) -> None:
    db = tmp_path / "structs_bad.db"
    db.write_text("....\n.....\n")
    with pytest.raises(ValueError):
        rup.read_db_structures(db)


def test_find_unpaired_runs_returns_expected_ranges() -> None:
    struct = "....((..))..."
    runs = rup.find_unpaired_runs(struct, min_len=3)
    # 5' run: [0,4), 3' run: [10,13)
    assert runs == [(0, 4), (10, 13)]


def test_find_unpaired_runs_none_when_below_threshold() -> None:
    struct = "..((..)).."
    runs = rup.find_unpaired_runs(struct, min_len=3)
    assert runs == []


def test_find_terminal_unpaired_ends_basic() -> None:
    struct = "....((..))...."
    len_5, len_3 = rup.find_terminal_unpaired_ends(struct, min_len=3)
    assert (len_5, len_3) == (4, 4)


def test_find_terminal_unpaired_ends_too_short_returns_zero() -> None:
    struct = "..((..))...."
    len_5, len_3 = rup.find_terminal_unpaired_ends(struct, min_len=3)
    assert (len_5, len_3) == (0, 0)


def test_find_terminal_unpaired_ends_all_unpaired_returns_zero() -> None:
    struct = "............"
    len_5, len_3 = rup.find_terminal_unpaired_ends(struct, min_len=3)
    assert (len_5, len_3) == (0, 0)


def test_check_balanced_parentheses() -> None:
    assert rup._check_balanced_parentheses("((..))") is True
    assert rup._check_balanced_parentheses("(()") is False
    assert rup._check_balanced_parentheses("())") is False
    assert rup._check_balanced_parentheses("..()..") is True


def test_refine_structure_calls_fold_on_internal_region() -> None:
    full_seq = "ACGUACGU"
    struct = "........"  # one long unpaired run, treated as internal
    fake_result = "()..().."

    call_args: dict = {}

    def fake_fold(fold_exe, subseq, temperature=None, extra_args=None):
        # Capture args for assertion
        call_args["fold_exe"] = fold_exe
        call_args["subseq"] = subseq
        call_args["temperature"] = temperature
        call_args["extra_args"] = extra_args
        assert len(subseq) == len(full_seq)
        return fake_result

    with patch(
        "src.lib.refine_unpaired_regions.call_rnastructure_fold",
        side_effect=fake_fold,
    ) as mock_fold:
        refined = rup.refine_structure(
            struct=struct,
            full_seq=full_seq,
            fold_exe=Path("/fake/Fold"),
            min_unpaired_len=3,
            temperature=None,
            extra_args=["-d"],
            cache={},
        )

    assert refined == fake_result
    mock_fold.assert_called_once()
    assert call_args["subseq"] == full_seq
    assert call_args["extra_args"] == ["-d"]


def test_refine_structure_respects_cache() -> None:
    full_seq = "ACGUACGU"
    struct = "........"
    cache = {("ACGUACGU", None, tuple()): "()..().."}

    with patch("src.lib.refine_unpaired_regions.call_rnastructure_fold") as mock_fold:
        refined = rup.refine_structure(
            struct=struct,
            full_seq=full_seq,
            fold_exe=Path("/fake/Fold"),
            min_unpaired_len=3,
            temperature=None,
            extra_args=None,
            cache=cache,
        )

    # Should not call Fold because cache already has the entry
    mock_fold.assert_not_called()
    assert refined == "()..().."


def test_refine_terminal_ends_with_duplex_applies_pairs() -> None:
    # Two-dot 5' and 3' ends so each has length 2
    struct = "..(()).."
    full_seq = "ACGUACGU"[: len(struct)]

    def fake_duplex(duplex_exe, seq5, seq3, temperature=None, extra_args=None):
        # Expect the end subsequences
        assert seq5 == full_seq[:2]
        assert seq3 == full_seq[-2:]
        # DuplexFold coordinates: 5' strand first (1..n1), then 3' (n1+1..n1+n2)
        # Pair first 5' base with last 3' base
        return [(1, 4)]

    with patch(
        "src.lib.refine_unpaired_regions.call_rnastructure_duplexfold",
        side_effect=fake_duplex,
    ):
        refined = rup.refine_terminal_ends_with_duplex(
            struct=struct,
            full_seq=full_seq,
            duplex_exe=Path("/fake/DuplexFold"),
            min_terminal_unpaired_len=2,
            temperature=None,
            extra_args=None,
            cache={},
        )

    assert refined[0] == "("
    assert refined[-1] == ")"
    # Internal region should remain untouched
    assert refined[2:6] == "(())"

