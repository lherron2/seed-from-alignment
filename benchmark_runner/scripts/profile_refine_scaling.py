#!/usr/bin/env python3
"""
Profile runtime sensitivity of refinement hyperparameters.

Runs the predictor on a single target FASTA repeatedly with different refinement
settings and records wall time and whether refinement timed out.

This is designed to find settings that keep a full selected_rnas run under a
tight wall-clock budget by identifying the most expensive knobs.
"""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class Trial:
    name: str
    overrides: dict[str, Any]


def _read_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _cmd_str(parts: list[str]) -> str:
    return " ".join(shlex.quote(p) for p in parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--params",
        default="refine_max_regions,refine_max_seeds,refine_max_solutions,refine_max_end_mask_len,refine_end_mask_step,refine_max_structures,refine_kissing_candidates,refine_min_unpaired,refine_max_helices_sequential,refine_max_helices_pairwise",
        help="Comma-separated refinement params to profile (order matters).",
    )
    parser.add_argument("--timeout-s", type=float, default=300.0)
    parser.add_argument("--refine-max-seconds", type=float, default=30.0)
    args_ns = parser.parse_args()

    br_root = Path(__file__).resolve().parents[1]
    defaults = _read_yaml(br_root / "defaults.yaml")
    paths = defaults.get("paths", {})
    predict = defaults.get("predict", {})
    rn = paths.get("rnastructure", {})

    py = br_root / ".venv" / "bin" / "python"
    python_bin = str(py) if py.exists() else "python3"
    py_path = str(br_root / "src")

    fasta = br_root / "data" / "predictions" / "selected_rnas" / "ESHAPE_AP_001.fa"
    if not fasta.exists():
        raise SystemExit(f"Missing FASTA fixture: {fasta}")

    out_base = br_root / "data" / "predictions" / "refine_scaling"
    out_base.mkdir(parents=True, exist_ok=True)

    base_args: dict[str, Any] = {
        "cm_db": paths.get("rfam_cm", ""),
        "rscape": paths.get("rscape", ""),
        "infernal_bin": paths.get("infernal_bin", ""),
        "allsub": rn.get("allsub", ""),
        "duplex": rn.get("duplex", ""),
        "fold": rn.get("fold", ""),
        "partition": rn.get("partition", ""),
        "probplot": rn.get("probplot", ""),
        "datapath": rn.get("datapath", ""),
        # Keep non-refinement work small for profiling.
        "top_k": 20,
        "n_samples": 30,
        "burn_in": 50,
        "thin": 10,
        "beta": 1.0,
        "seed": 0,
        "min_loop_sep": int(predict.get("min_loop_sep", 0)),
        "pk_alpha": float(predict.get("pk_alpha", 0.5)),
        "cov_mode": str(predict.get("cov_mode", "logE")),
        "cov_alpha": float(predict.get("cov_alpha", 3.0)),
        "cov_min_power": float(predict.get("cov_min_power", 0.1)),
        "weight_calibration_method": str(predict.get("weight_calibration_method", "robust_z")),
        "weight_calibration_zmax": float(predict.get("weight_calibration_zmax", 3.0)),
        "weight_alpha_core": float(predict.get("weight_alpha_core", 1.0)),
        "weight_alpha_alt": float(predict.get("weight_alpha_alt", 1.0)),
        "weight_alpha_cov": float(predict.get("weight_alpha_cov", 1.0)),
        "weight_alpha_thermo": float(predict.get("weight_alpha_thermo", 1.0)),
        "thermo_mode": "off",
        "thermo_weight": 0.0,
        "stem_start_penalty_scale": float(predict.get("stem_start_penalty_scale", 0.05)),
        "stem_len1_penalty_scale": float(predict.get("stem_len1_penalty_scale", 0.10)),
        "stem_len2_penalty_scale": float(predict.get("stem_len2_penalty_scale", 0.02)),
        "stem_log_reward_scale": float(predict.get("stem_log_reward_scale", 0.02)),
        "stem_support_quantile": float(predict.get("stem_support_quantile", 0.5)),
        "pair_penalty_mode": str(predict.get("pair_penalty_mode", "length_aware")),
        "pair_penalty_c0": float(predict.get("pair_penalty_c0", 0.10)),
        "pair_penalty_c1": float(predict.get("pair_penalty_c1", 0.50)),
        "pair_penalty_scale": float(predict.get("pair_penalty_scale", 0.25)),
        "refine_max_seconds": float(args_ns.refine_max_seconds),
        "refine_max_structures": 120,
        "refine_min_unpaired": 6,
        "refine_end_mask_step": 20,
        "refine_max_end_mask_len": 10,
        "refine_max_helices_sequential": 20,
        "refine_max_helices_pairwise": 20,
        "refine_max_regions": 4,
        "refine_max_seeds": 50,
        "refine_max_solutions": 80,
        "refine_kissing_candidates": 50,
    }

    all_sweep_specs: dict[str, list[int]] = {
        "refine_max_structures": [40, 80, 160],
        "refine_min_unpaired": [4, 6, 10],
        "refine_end_mask_step": [10, 20, 30],
        "refine_max_end_mask_len": [10, 20, 30],
        "refine_max_helices_sequential": [8, 16, 24],
        "refine_max_helices_pairwise": [8, 16, 24],
        "refine_max_regions": [2, 4, 8, 12, 24],
        "refine_max_seeds": [5, 10, 25, 50, 100],
        "refine_max_solutions": [20, 40, 80, 120],
        "refine_kissing_candidates": [0, 20, 50, 100],
    }
    selected_params = [p.strip() for p in str(args_ns.params).split(",") if p.strip()]
    sweep_specs = [(p, all_sweep_specs[p]) for p in selected_params if p in all_sweep_specs]

    trials: list[Trial] = [Trial("baseline", {})]
    for param, values in sweep_specs:
        for v in values:
            trials.append(Trial(f"{param}={v}", {param: v}))

    results_path = br_root / "data" / "results" / "refine_scaling_profile.csv"
    results_path.parent.mkdir(parents=True, exist_ok=True)

    with results_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "name",
                "param",
                "value",
                "seconds",
                "exit_code",
            ],
        )
        writer.writeheader()

        for trial in trials:
            args = dict(base_args)
            args.update(trial.overrides)
            outdir = out_base / trial.name.replace("/", "_")
            if outdir.exists():
                subprocess.run(["rm", "-rf", str(outdir)], check=True)
            outdir.mkdir(parents=True, exist_ok=True)

            cmd = [
                "DATAPATH=" + str(args["datapath"]),
                "PYTHONPATH=" + py_path,
                python_bin,
                "-m",
                "ssbench.predict.cacofold_mcmc_pipeline",
                "--cm-db",
                str(args["cm_db"]),
                "--rscape",
                str(args["rscape"]),
                "--infernal-bin",
                str(args["infernal_bin"]),
                "--allsub-exe",
                str(args["allsub"]),
                "--duplex-exe",
                str(args["duplex"]),
                "--fold-exe",
                str(args["fold"]),
                "--top-k",
                str(args["top_k"]),
                "--n-samples",
                str(args["n_samples"]),
                "--burn-in",
                str(args["burn_in"]),
                "--thin",
                str(args["thin"]),
                "--beta",
                str(args["beta"]),
                "--seed",
                str(args["seed"]),
                "--min-loop-sep",
                str(args["min_loop_sep"]),
                "--pk-alpha",
                str(args["pk_alpha"]),
                "--pair-penalty-mode",
                str(args["pair_penalty_mode"]),
                "--pair-penalty-c0",
                str(args["pair_penalty_c0"]),
                "--pair-penalty-c1",
                str(args["pair_penalty_c1"]),
                "--pair-penalty-scale",
                str(args["pair_penalty_scale"]),
                "--cov-mode",
                str(args["cov_mode"]),
                "--cov-alpha",
                str(args["cov_alpha"]),
                "--cov-min-power",
                str(args["cov_min_power"]),
                "--weight-calibration-method",
                str(args["weight_calibration_method"]),
                "--weight-calibration-zmax",
                str(args["weight_calibration_zmax"]),
                "--weight-alpha-core",
                str(args["weight_alpha_core"]),
                "--weight-alpha-alt",
                str(args["weight_alpha_alt"]),
                "--weight-alpha-cov",
                str(args["weight_alpha_cov"]),
                "--weight-alpha-thermo",
                str(args["weight_alpha_thermo"]),
                "--thermo-mode",
                str(args["thermo_mode"]),
                "--thermo-weight",
                str(args["thermo_weight"]),
                "--stem-start-penalty-scale",
                str(args["stem_start_penalty_scale"]),
                "--stem-len1-penalty-scale",
                str(args["stem_len1_penalty_scale"]),
                "--stem-len2-penalty-scale",
                str(args["stem_len2_penalty_scale"]),
                "--stem-log-reward-scale",
                str(args["stem_log_reward_scale"]),
                "--stem-support-quantile",
                str(args["stem_support_quantile"]),
                "--refine-max-seconds",
                str(args["refine_max_seconds"]),
                "--refine-max-structures",
                str(args["refine_max_structures"]),
                "--refine-min-unpaired",
                str(args["refine_min_unpaired"]),
                "--refine-end-mask-step",
                str(args["refine_end_mask_step"]),
                "--refine-max-end-mask-len",
                str(args["refine_max_end_mask_len"]),
                "--refine-max-helices-sequential",
                str(args["refine_max_helices_sequential"]),
                "--refine-max-helices-pairwise",
                str(args["refine_max_helices_pairwise"]),
                "--refine-max-regions",
                str(args["refine_max_regions"]),
                "--refine-max-seeds",
                str(args["refine_max_seeds"]),
                "--refine-max-solutions",
                str(args["refine_max_solutions"]),
                "--refine-kissing-candidates",
                str(args["refine_kissing_candidates"]),
                str(fasta),
                str(outdir),
            ]

            start = time.time()
            try:
                proc = subprocess.run(
                    _cmd_str(cmd),
                    shell=True,
                    text=True,
                    capture_output=True,
                    timeout=float(args_ns.timeout_s),
                )
                code = proc.returncode
            except subprocess.TimeoutExpired:
                code = 124
            elapsed = time.time() - start

            if trial.name == "baseline":
                param = ""
                value = ""
            else:
                param, value = trial.name.split("=", 1)

            writer.writerow(
                {
                    "name": trial.name,
                    "param": param,
                    "value": value,
                    "seconds": f"{elapsed:.3f}",
                    "exit_code": str(code),
                }
            )
            fh.flush()


if __name__ == "__main__":
    main()
