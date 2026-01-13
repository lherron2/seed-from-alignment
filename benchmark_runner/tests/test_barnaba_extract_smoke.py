import pytest


def test_barnaba_import() -> None:
    try:
        import barnaba  # noqa: F401
    except ImportError:
        pytest.skip("barnaba not installed")
