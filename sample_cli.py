#!/usr/bin/env python3
"""
CLI entry point for sampling pseudoknotted secondary structures
from a CaCoFold/R-scape Stockholm (.sto) file.

This is a thin wrapper around src.lib.sample_cacofold_structures.main(),
so all of the original algorithms and options remain intact.
"""

from src.lib.sample_cacofold_structures import main


if __name__ == "__main__":
    main()

