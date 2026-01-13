"""
Utilities for normalizing and blending per-pair evidence components.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass


@dataclass(frozen=True)
class NormalizationConfig:
    method: str = "none"  # none | robust_z
    zmax: float = 3.0


def _robust_center_scale(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values])
    return med, mad


def _minmax_scale(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    vmin = min(values)
    vmax = max(values)
    return vmin, vmax


def normalize_component(
    component: dict[tuple[int, int], float],
    config: NormalizationConfig,
) -> dict[tuple[int, int], float]:
    if config.method == "none" or not component:
        return dict(component)

    values = list(component.values())
    if config.method == "robust_z":
        median, mad = _robust_center_scale(values)
        if mad <= 0:
            vmin, vmax = _minmax_scale(values)
            if vmax <= vmin:
                return dict.fromkeys(component, 0.0)
            span = vmax - vmin
            return {
                k: max(
                    -config.zmax, min(config.zmax, (2.0 * (v - vmin) / span - 1.0) * config.zmax)
                )
                for k, v in component.items()
            }
        inv = 1.0 / mad
        return {
            k: max(-config.zmax, min(config.zmax, (v - median) * inv)) for k, v in component.items()
        }

    raise ValueError(f"Unknown normalization method: {config.method}")


def blend_components(
    components: dict[str, dict[tuple[int, int], float]],
    alphas: dict[str, float],
    config: NormalizationConfig,
) -> dict[tuple[int, int], float]:
    blended: dict[tuple[int, int], float] = {}
    for name, comp in components.items():
        alpha = float(alphas.get(name, 0.0))
        if alpha == 0.0:
            continue
        normed = normalize_component(comp, config)
        for pair, val in normed.items():
            blended[pair] = blended.get(pair, 0.0) + alpha * float(val)
    return blended


def weight_scale(values: dict[tuple[int, int], float]) -> float:
    if not values:
        return 0.0
    try:
        return float(statistics.median(values.values()))
    except statistics.StatisticsError:
        return 0.0


def quantile_threshold(values: dict[tuple[int, int], float], q: float) -> float:
    if not values:
        return 0.0
    q = max(0.0, min(1.0, float(q)))
    sorted_vals = sorted(values.values())
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    idx = int(math.floor(q * (len(sorted_vals) - 1)))
    return float(sorted_vals[idx])
