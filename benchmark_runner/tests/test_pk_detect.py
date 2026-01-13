from ssbench.truth.pk_detect import crosses, split_nested_pk


def test_crosses() -> None:
    assert crosses((0, 5), (2, 8))
    assert not crosses((0, 5), (6, 10))


def test_split_nested_pk() -> None:
    pairs = [(0, 5), (2, 8)]
    nested, pk = split_nested_pk(pairs)
    assert (0, 5) in pk or (2, 8) in pk
