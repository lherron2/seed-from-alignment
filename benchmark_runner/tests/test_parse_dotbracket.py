from ssbench.predict.parse_dotbracket import pairs_from_dotbracket, parse_dotbracket


def test_parse_dotbracket_two_line() -> None:
    seq, structs = parse_dotbracket("ACGU\n(..)\n")
    assert seq == "ACGU"
    assert structs == ["(..)"]


def test_pairs_from_dotbracket() -> None:
    pairs = pairs_from_dotbracket("(())")
    assert pairs == [(1, 2), (0, 3)]
