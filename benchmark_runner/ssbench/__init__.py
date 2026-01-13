"""Shim package to support running without installation."""

from pathlib import Path

# Extend package path to include src/ssbench.
_src_pkg = Path(__file__).resolve().parents[1] / "src" / "ssbench"
if _src_pkg.exists():
    __path__.append(str(_src_pkg))
