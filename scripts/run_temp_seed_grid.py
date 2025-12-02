#!/usr/bin/env python3
"""
Run the CaCoFold → refine → PK → Rosetta pipeline on a grid of
(temperature, seed) combinations, then concatenate all Rosetta-ready
ensembles into a single deduplicated .db file.

Usage example (from repo root):

    python scripts/run_temp_seed_grid.py \
        --config example/pipeline_config.yaml \
        --temps 10,20,30 \
        --seeds 0,1,2,3 \
        --mode refine_first \
        --output example/pipeline_output/9E5I_rosetta_all.db
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path
from typing import List

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.lib.pipeline import load_pipeline_config, run_pipeline  # type: ignore

from src.lib.pipeline import load_pipeline_config, run_pipeline


def _parse_float_list(s: str) -> List[float]:
    return [float(x) for x in s.split(",") if x.strip() != ""]


def _parse_int_list(s: str) -> List[int]:
    return [int(x) for x in s.split(",") if x.strip() != ""]


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the RNA pipeline over a grid of temperatures and seeds, "
            "then concatenate + deduplicate the Rosetta-ready outputs."
        )
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="Base pipeline config (JSON or YAML). Typically example/pipeline_config.yaml.",
    )
    parser.add_argument(
        "--temps",
        required=True,
        help="Comma-separated list of temperatures (passed to RNAstructure -T). "
             "Example: 10,20,30",
    )
    parser.add_argument(
        "--seeds",
        required=True,
        help="Comma-separated list of random seeds for the PK sampler. "
             "Example: 0,1,2,3",
    )
    parser.add_argument(
        "--mode",
        choices=["legacy", "refine_first"],
        default="refine_first",
        help="Which pre-wired pipeline to run (default: refine_first).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Path for the final concatenated, deduplicated Rosetta-ready .db. "
            "If not given, defaults to <work_dir>/<basename>_grid_rosetta_all.db"
        ),
    )
    args = parser.parse_args(argv)

    cfg_path = Path(args.config).resolve()
    base_cfg = load_pipeline_config(cfg_path)

    temps = _parse_float_list(args.temps)
    seeds = _parse_int_list(args.seeds)

    mode = args.mode.lower()

    # Base name for per-run rosetta DBs
    base_rosetta_path = base_cfg.io.rosetta_db
    base_rosetta_stem = base_rosetta_path.stem  # e.g. "9E5I_sampled_rosetta"
    work_dir = base_cfg.io.work_dir

    # Determine final output path
    if args.output is not None:
        final_out = Path(args.output).resolve()
        final_out.parent.mkdir(parents=True, exist_ok=True)
    else:
        final_out = (work_dir / f"{base_rosetta_stem}_grid_rosetta_all.db").resolve()
        final_out.parent.mkdir(parents=True, exist_ok=True)

    sys.stderr.write(
        f"[GRID] Base config: {cfg_path}\n"
        f"[GRID] Work dir   : {work_dir}\n"
        f"[GRID] Mode       : {mode}\n"
        f"[GRID] Temps      : {temps}\n"
        f"[GRID] Seeds      : {seeds}\n"
        f"[GRID] Final out  : {final_out}\n"
    )

    rosetta_paths: List[Path] = []

    # ----------------- run the grid ----------------- #
    for T in temps:
        for seed in seeds:
            cfg = copy.deepcopy(base_cfg)

            # Override temperature + seed in the cloned config
            cfg.refine.temperature = float(T)
            cfg.sample.seed = int(seed)

            # Give each run its own rosetta_db filename
            run_tag = f"T{T:g}_seed{seed}"
            run_rosetta_name = f"{base_rosetta_stem}_{run_tag}.db"
            cfg.io.rosetta_db = (work_dir / run_rosetta_name).resolve()

            sys.stderr.write(
                f"[GRID] Running mode={mode}, T={T:g}, seed={seed} → "
                f"{cfg.io.rosetta_db}\n"
            )

            out_path = run_pipeline(cfg, mode=mode)
            rosetta_paths.append(Path(out_path).resolve())

    sys.stderr.write(
        f"[GRID] Completed {len(rosetta_paths)} runs. "
        f"Collecting and deduplicating structures...\n"
    )

    # ----------------- concatenate + deduplicate ----------------- #
    seen = set()
    unique_count = 0

    with final_out.open("w") as out_f:
        for p in rosetta_paths:
            if not p.is_file():
                sys.stderr.write(f"[GRID] WARNING: missing rosetta file {p}, skipping.\n")
                continue

            sys.stderr.write(f"[GRID] Reading {p}\n")
            with p.open() as fh:
                for line in fh:
                    s = line.strip()
                    if not s:
                        continue
                    if s in seen:
                        continue
                    seen.add(s)
                    out_f.write(s + "\n")
                    unique_count += 1

    sys.stderr.write(
        f"[GRID] Wrote {unique_count} unique structures to {final_out}\n"
    )


if __name__ == "__main__":
    main()

