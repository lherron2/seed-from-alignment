"""Configuration utilities for ssbench."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

_DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[2] / "defaults.yaml"


@dataclass
class BenchmarkConfig:
    """Configuration for dataset selection and scoring."""

    release_id: str = "current"
    resolution: str = "4.0"
    min_length: int | None = None
    max_length: int | None = None
    length_buckets: list[tuple[int, int]] = ((30, 80), (81, 150), (151, 300), (301, 1000))
    split_seed: int = 0
    drop_missing_rfam: bool = True
    fallback_to_ec: bool = False
    helix_overlap_threshold: float = 0.5
    k: int = 1


def load_config(path: str | Path | None) -> BenchmarkConfig:
    """Load YAML config if provided; otherwise return defaults."""

    if path is None:
        if _DEFAULT_CONFIG_PATH.exists():
            path = _DEFAULT_CONFIG_PATH
        else:
            return BenchmarkConfig()

    data: dict[str, Any] = {}
    with open(path) as fh:
        data = yaml.safe_load(fh) or {}

    dataset = data.get("dataset", data)
    score = data.get("score", data)
    predict = data.get("predict", data)

    return BenchmarkConfig(
        release_id=str(dataset.get("release_id", "current")),
        resolution=str(dataset.get("resolution", "4.0")),
        min_length=dataset.get("min_length"),
        max_length=dataset.get("max_length"),
        length_buckets=[
            tuple(b)
            for b in dataset.get("length_buckets", [(30, 80), (81, 150), (151, 300), (301, 1000)])
        ],
        split_seed=int(dataset.get("split_seed", 0)),
        drop_missing_rfam=bool(dataset.get("drop_missing_rfam", True)),
        fallback_to_ec=bool(dataset.get("fallback_to_ec", False)),
        helix_overlap_threshold=float(score.get("helix_overlap_threshold", 0.5)),
        k=int(score.get("k", predict.get("top_k", 1))),
    )
