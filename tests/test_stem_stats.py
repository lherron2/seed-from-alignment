from src.lib import stem_stats as ss


def test_stem_stats_basic() -> None:
    pairs = {(0, 5), (1, 4), (3, 8)}
    stats = ss.stem_stats(pairs)
    assert stats.n_stems == 2
    assert stats.n_len1 == 1
    assert stats.n_len2 == 1
    assert (3, 8) in stats.len1_pairs
