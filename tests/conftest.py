# tests/conftest.py
"""Shared test fixtures for MCMC sampler tests."""

import sys
from pathlib import Path

import pytest

# Repo root = parent of this file's directory
ROOT = Path(__file__).resolve().parents[1]

# Ensure the repo root is on sys.path so `import src...` works
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture
def fixtures_dir() -> Path:
    """Path to test fixtures directory."""
    return ROOT / "tests" / "fixtures"


@pytest.fixture
def small_candidates() -> list[tuple[int, int, float]]:
    """Small set of candidate pairs for unit tests.

    These form a simple hairpin-like structure with some alternatives.
    """
    return [
        (0, 11, 2.0),  # Outer pair
        (1, 10, 2.0),  # Next inner
        (2, 9, 2.0),  # Inner
        (3, 8, 1.5),  # Innermost
        (0, 7, 1.0),  # Alternative outer (crosses with (2,9) etc)
        (4, 11, 1.0),  # Alternative (crosses with (0,11))
    ]


@pytest.fixture
def small_weights(small_candidates: list[tuple[int, int, float]]) -> dict[tuple[int, int], float]:
    """Weight dictionary from small_candidates."""
    return {(i, j): w for i, j, w in small_candidates}


@pytest.fixture
def default_energy_params():
    """Default energy parameters."""
    from src.lib.energy import EnergyParams

    return EnergyParams()


@pytest.fixture
def tiny_sto_path(fixtures_dir: Path) -> Path:
    """Path to tiny test Stockholm file."""
    return fixtures_dir / "tiny.sto"


@pytest.fixture
def tiny_pk_sto_path(fixtures_dir: Path) -> Path:
    """Path to tiny pseudoknotted Stockholm file."""
    return fixtures_dir / "tiny_pk.sto"


@pytest.fixture
def tiny_db_path(fixtures_dir: Path) -> Path:
    """Path to tiny test .db file."""
    return fixtures_dir / "tiny.db"
