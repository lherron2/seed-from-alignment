from ssbench.metrics.pair_metrics import compute_pair_metrics


def test_pair_metrics() -> None:
    pred = {(0, 5), (2, 7)}
    ref = {(0, 5)}
    metrics = compute_pair_metrics(pred, ref, 10)
    assert metrics.precision == 0.5
    assert metrics.recall == 1.0
