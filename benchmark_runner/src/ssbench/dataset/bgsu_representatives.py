"""Download and filter BGSU RNA 3D Hub representative sets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

BGSU_CURRENT_URL = "http://rna.bgsu.edu/rna3dhub/nrlist/download/current/{resolution}/csv"
BGSU_FULL_URL = "https://rna.bgsu.edu/rna3dhub/nrlist/download/NR/{release}/{resolution}/csv/full"


@dataclass
class DatasetFilters:
    """Selection filters for representative sets."""

    release_id: str = "current"
    resolution: str = "4.0"
    single_chain_only: bool = True
    min_length: int | None = None
    max_length: int | None = None
    limit: int | None = None


def _download_csv(url: str) -> pd.DataFrame:
    import io

    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text))

    # Some BGSU CSVs have no header; the first data row becomes columns.
    if len(df.columns) <= 3 and any(str(c).startswith("NR_") for c in df.columns):
        df = pd.read_csv(io.StringIO(resp.text), header=None)
        df.columns = ["ec_id", "ife_id", "members"][: len(df.columns)]

    return df


def download_representatives(filters: DatasetFilters, full: bool = False) -> pd.DataFrame:
    """Download BGSU representative set CSV.

    Args:
        filters: Dataset filters.
        full: Whether to download the "full" representative set (with rfam column).

    Returns:
        DataFrame with representative set entries.
    """
    if full:
        url = BGSU_FULL_URL.format(release=filters.release_id, resolution=filters.resolution)
    else:
        url = BGSU_CURRENT_URL.format(resolution=filters.resolution)

    df = _download_csv(url)

    if filters.single_chain_only and "ife_id" in df.columns:
        df = df[~df["ife_id"].astype(str).str.contains("+", regex=False)]

    if filters.min_length is not None and "nts_observed" in df.columns:
        df = df[df["nts_observed"] >= filters.min_length]
    if filters.max_length is not None and "nts_observed" in df.columns:
        df = df[df["nts_observed"] <= filters.max_length]

    if filters.limit is not None:
        df = df.head(filters.limit)

    return df


def write_manifest(df: pd.DataFrame, out_path: str | Path) -> None:
    """Write manifest CSV with standardized columns."""
    cols = [
        "target_id",
        "pdb_id",
        "chain_id",
        "ife_id",
        "ec_id",
        "rfam",
        "nts_observed",
        "pdb_resolution",
    ]

    out = pd.DataFrame()
    if "ife_id" in df.columns:
        out["ife_id"] = df["ife_id"].astype(str)
        out["target_id"] = out["ife_id"]
    if "ec_id" in df.columns:
        out["ec_id"] = df["ec_id"]
    if "rfam" in df.columns:
        out["rfam"] = df["rfam"]
    if "nts_observed" in df.columns:
        out["nts_observed"] = df["nts_observed"]
    if "pdb_resolution" in df.columns:
        out["pdb_resolution"] = df["pdb_resolution"]

    # Derive pdb_id and chain_id from ife_id if missing.
    if "ife_id" in out.columns:
        parts = out["ife_id"].str.split("|", expand=True)
        if parts.shape[1] >= 3:
            if "pdb_id" not in out.columns:
                out["pdb_id"] = parts[0].str.lower()
            if "chain_id" not in out.columns:
                out["chain_id"] = parts[2]
        else:
            if "pdb_id" not in out.columns:
                out["pdb_id"] = pd.NA
            if "chain_id" not in out.columns:
                out["chain_id"] = pd.NA

    for col in cols:
        if col not in out.columns:
            out[col] = pd.NA
    out = out.reindex(columns=cols)
    out.to_csv(out_path, index=False)


def stratify_buckets(lengths: Iterable[int], buckets: list[tuple[int, int]]) -> list[str]:
    """Assign length bucket labels based on length ranges."""
    labels: list[str] = []
    for length in lengths:
        label = "unknown"
        for lo, hi in buckets:
            if lo <= length <= hi:
                label = f"{lo}-{hi}"
                break
        labels.append(label)
    return labels


def split_by_rfam(df: pd.DataFrame, seed: int, drop_missing: bool = True, fallback_to_ec: bool = False) -> pd.DataFrame:
    """Split dataset into train/test/valid based on Rfam family.

    Args:
        df: Input manifest dataframe.
        seed: Random seed.
        drop_missing: Drop rows without rfam.
        fallback_to_ec: Use ec_id to group when rfam missing.
    """
    import random

    if "rfam" not in df.columns:
        if fallback_to_ec and "ec_id" in df.columns:
            df = df.copy()
            df["rfam"] = df["ec_id"]
        else:
            df = df.copy()
            df["split"] = "train"
            return df

    df = df.copy()
    if drop_missing:
        rfam_series = df["rfam"]
        if rfam_series.isna().all():
            # No rfam available; keep all and assign a default split.
            df["split"] = "train"
            return df
        df = df[rfam_series.notna()]

    rfams = sorted(df["rfam"].astype(str).unique())
    rng = random.Random(seed)
    rng.shuffle(rfams)

    n = len(rfams)
    n_train = int(0.8 * n)
    n_valid = int(0.1 * n)
    train = set(rfams[:n_train])
    valid = set(rfams[n_train : n_train + n_valid])
    test = set(rfams[n_train + n_valid :])

    def assign(rfam: str) -> str:
        if rfam in train:
            return "train"
        if rfam in valid:
            return "valid"
        if rfam in test:
            return "test"
        return "train"

    df["split"] = df["rfam"].astype(str).map(assign)
    return df
