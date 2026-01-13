from src.lib import weight_calibration as wc


def test_normalize_component_robust_z_constant() -> None:
    comp = {(0, 1): 1.0, (1, 2): 1.0}
    norm = wc.normalize_component(comp, wc.NormalizationConfig(method="robust_z", zmax=3.0))
    assert set(norm.values()) == {0.0}


def test_blend_components_none() -> None:
    comps = {"core": {(0, 1): 2.0}, "cov": {(0, 1): 1.0, (1, 2): 3.0}}
    alphas = {"core": 1.0, "cov": 2.0}
    blended = wc.blend_components(comps, alphas, wc.NormalizationConfig(method="none", zmax=3.0))
    assert blended[(0, 1)] == 4.0
    assert blended[(1, 2)] == 6.0


def test_quantile_threshold() -> None:
    values = {(0, 1): 1.0, (1, 2): 2.0, (2, 3): 3.0, (3, 4): 4.0}
    assert wc.quantile_threshold(values, 0.0) == 1.0
    assert wc.quantile_threshold(values, 1.0) == 4.0
