from ssbench.truth.helices import helix_overlap, identify_helices


def test_identify_helices() -> None:
    pairs = [(0, 9), (1, 8), (5, 12)]
    helices = identify_helices(pairs)
    assert len(helices) == 2


def test_helix_overlap() -> None:
    h1 = [(0, 9), (1, 8)]
    h2 = [(1, 8), (2, 7)]
    assert helix_overlap(h1, h2) == 0.5
