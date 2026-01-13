"""Predictor interface types."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class PredictionRecord:
    """Predicted structures for a target."""

    target_id: str
    sequence: str
    dotbrackets: list[str]
