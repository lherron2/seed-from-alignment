"""Table generation helpers."""

from __future__ import annotations

import pandas as pd


def aggregate_by_bucket(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["split", "bucket"]).mean(numeric_only=True).reset_index()
