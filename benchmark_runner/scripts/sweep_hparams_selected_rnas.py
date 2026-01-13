#!/usr/bin/env python3
"""
Systematic hyperparameter sweep for the selected_rnas benchmark.

Runs predict+score+report for a fixed manifest/truth and records summary metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class RunConfig:
    name: str
    overrides: dict[str, Any]


def _read_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _mean_metrics(metrics_path: Path) -> dict[str, float]:
    with metrics_path.open() as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if not rows:
        return {"f1": 0.0, "precision": 0.0, "recall": 0.0}
    f1 = sum(float(r["f1"]) for r in rows) / len(rows)
    precision = sum(float(r["precision"]) for r in rows) / len(rows)
    recall = sum(float(r["recall"]) for r in rows) / len(rows)
    return {"f1": f1, "precision": precision, "recall": recall}


def _build_pred_cmd(base: dict[str, Any], overrides: dict[str, Any]) -> list[str]:
    args = dict(base)
    args.update(overrides)

    cmd = [
        f"DATAPATH={args['datapath']}",
        f"RNASTRUCTURE_DOCKER_IMAGE={args.get('rnastructure_docker_image','')}",
        f"PYTHONPATH={args['py_path']}",
        args["python_bin"],
        "-m",
        "ssbench.predict.cacofold_mcmc_pipeline",
        "--cm-db",
        args["rfam_cm"],
        "--rscape",
        args["rscape"],
        "--infernal-bin",
        args["infernal_bin"],
        "--allsub-exe",
        args["allsub_exe"],
        "--duplex-exe",
        args["duplex_exe"],
        "--fold-exe",
        args["fold_exe"],
    ]

    if args.get("partition_exe"):
        cmd += ["--partition-exe", args["partition_exe"]]
    if args.get("probplot_exe"):
        cmd += ["--probplot-exe", args["probplot_exe"]]

    cmd += [
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
        "--thermo-max-structures",
        str(args["thermo_max_structures"]),
        "--thermo-min-count",
        str(args["thermo_min_count"]),
        "--thermo-min-prob",
        str(args["thermo_min_prob"]),
        "--thermo-log-eps",
        str(args["thermo_log_eps"]),
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
        "{fasta}",
        "{outdir}",
    ]

    if args.get("pair_penalty") is not None:
        cmd += ["--pair-penalty", str(args["pair_penalty"])]
    else:
        cmd += ["--pair-penalty-scale", str(args["pair_penalty_scale"])]

    if args.get("pair_penalty_min") is not None:
        cmd += ["--pair-penalty-min", str(args["pair_penalty_min"])]
    if args.get("pair_penalty_max") is not None:
        cmd += ["--pair-penalty-max", str(args["pair_penalty_max"])]
    if args.get("cov_forbid_negative"):
        cmd += ["--cov-forbid-negative"]

    return cmd


def _cmd_str(cmd: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-runs", type=int, default=None)
    parser.add_argument("--top-k", type=int, default=None)
    parser.add_argument("--n-samples", type=int, default=None)
    parser.add_argument("--burn-in", type=int, default=None)
    parser.add_argument("--thin", type=int, default=None)
    parser.add_argument("--refine-max-seconds", type=float, default=None)
    parser.add_argument("--predict-timeout-s", type=float, default=None)
    parser.add_argument("--score-timeout-s", type=float, default=None)
    parser.add_argument("--report-timeout-s", type=float, default=None)
    args = parser.parse_args()

    br_root = Path(__file__).resolve().parents[1]
    defaults = _read_yaml(br_root / "defaults.yaml")
    paths = defaults.get("paths", {})
    predict = defaults.get("predict", {})
    rnastructure = paths.get("rnastructure", {})

    def resolve_from_root(p: str) -> str:
        if not p:
            return p
        if os.path.isabs(p):
            return p
        return str((br_root / p).resolve())

    base = {
        "python_bin": str((br_root / ".venv" / "bin" / "python").resolve()),
        "py_path": str(br_root / "src"),
        "rfam_cm": paths.get("rfam_cm", ""),
        "infernal_bin": paths.get("infernal_bin", ""),
        "rscape": paths.get("rscape", ""),
        "allsub_exe": resolve_from_root(str(rnastructure.get("allsub", ""))),
        "duplex_exe": resolve_from_root(str(rnastructure.get("duplex", ""))),
        "fold_exe": resolve_from_root(str(rnastructure.get("fold", ""))),
        "partition_exe": resolve_from_root(str(rnastructure.get("partition", ""))),
        "probplot_exe": resolve_from_root(str(rnastructure.get("probplot", ""))),
        "datapath": rnastructure.get("datapath", ""),
        "rnastructure_docker_image": str(rnastructure.get("docker_image", "")),
        "top_k": int(args.top_k if args.top_k is not None else predict.get("top_k", 100)),
        "n_samples": int(
            args.n_samples if args.n_samples is not None else predict.get("n_samples", 1000)
        ),
        "burn_in": int(args.burn_in if args.burn_in is not None else predict.get("burn_in", 1000)),
        "thin": int(args.thin if args.thin is not None else predict.get("thin", 100)),
        "beta": float(predict.get("beta", 1.0)),
        "seed": int(predict.get("seed", 0)),
        "min_loop_sep": int(predict.get("min_loop_sep", 0)),
        "pk_alpha": float(predict.get("pk_alpha", 0.5)),
        "pair_penalty": predict.get("pair_penalty"),
        "pair_penalty_scale": float(predict.get("pair_penalty_scale", 0.25)),
        "pair_penalty_mode": str(predict.get("pair_penalty_mode", "legacy")),
        "pair_penalty_c0": float(predict.get("pair_penalty_c0", 0.10)),
        "pair_penalty_c1": float(predict.get("pair_penalty_c1", 0.50)),
        "pair_penalty_min": predict.get("pair_penalty_min"),
        "pair_penalty_max": predict.get("pair_penalty_max"),
        "cov_mode": str(predict.get("cov_mode", "logE")),
        "cov_alpha": float(predict.get("cov_alpha", 3.0)),
        "cov_min_power": float(predict.get("cov_min_power", 0.1)),
        "cov_forbid_negative": bool(predict.get("cov_forbid_negative", False)),
        "weight_calibration_method": str(predict.get("weight_calibration_method", "none")),
        "weight_calibration_zmax": float(predict.get("weight_calibration_zmax", 3.0)),
        "weight_alpha_core": float(predict.get("weight_alpha_core", 1.0)),
        "weight_alpha_alt": float(predict.get("weight_alpha_alt", 1.0)),
        "weight_alpha_cov": float(predict.get("weight_alpha_cov", 1.0)),
        "weight_alpha_thermo": float(predict.get("weight_alpha_thermo", 1.0)),
        "thermo_mode": str(predict.get("thermo_mode", "allsub")),
        "thermo_weight": float(predict.get("thermo_weight", 1.0)),
        "thermo_max_structures": int(predict.get("thermo_max_structures", 50)),
        "thermo_min_count": int(predict.get("thermo_min_count", 2)),
        "thermo_min_prob": float(predict.get("thermo_min_prob", 0.001)),
        "thermo_log_eps": float(predict.get("thermo_log_eps", 1e-6)),
        "stem_start_penalty_scale": float(predict.get("stem_start_penalty_scale", 0.05)),
        "stem_len1_penalty_scale": float(predict.get("stem_len1_penalty_scale", 0.10)),
        "stem_len2_penalty_scale": float(predict.get("stem_len2_penalty_scale", 0.02)),
        "stem_log_reward_scale": float(predict.get("stem_log_reward_scale", 0.02)),
        "stem_support_quantile": float(predict.get("stem_support_quantile", 0.5)),
        "refine_max_seconds": float(
            args.refine_max_seconds
            if args.refine_max_seconds is not None
            else predict.get("refine_max_seconds", 30.0)
        ),
        "refine_max_structures": int(predict.get("refine_max_structures", 200)),
        "refine_min_unpaired": int(predict.get("refine_min_unpaired", 6)),
        "refine_end_mask_step": int(predict.get("refine_end_mask_step", 20)),
        "refine_max_end_mask_len": int(predict.get("refine_max_end_mask_len", 10)),
        "refine_max_helices_sequential": int(predict.get("refine_max_helices_sequential", 20)),
        "refine_max_helices_pairwise": int(predict.get("refine_max_helices_pairwise", 20)),
        "refine_max_regions": int(predict.get("refine_max_regions", 4)),
        "refine_max_seeds": int(predict.get("refine_max_seeds", 100)),
        "refine_max_solutions": int(predict.get("refine_max_solutions", 100)),
        "refine_kissing_candidates": int(predict.get("refine_kissing_candidates", 100)),
    }

    if not Path(base["python_bin"]).exists():
        base["python_bin"] = "python"

    manifest = br_root / "data" / "manifests" / "selected_rnas_truth.csv"
    truth_dir = br_root / "data" / "truth_selected"
    out_root = br_root / "data" / "predictions" / "selected_rnas_sweep"
    results_dir = br_root / "data" / "results"
    out_root.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    runs: list[RunConfig] = [RunConfig("baseline", {})]

    # Refinement sweep: 10 factors, 2 levels each, using an L12 orthogonal array (12 runs).
    # 1 -> low, 2 -> high
    # Factors are ordered to match the user's requested refinement knobs.
    refine_levels: dict[str, tuple[Any, Any]] = {
        "refine_max_structures": (100, 250),
        "refine_min_unpaired": (4, 8),
        "refine_end_mask_step": (10, 25),
        "refine_max_end_mask_len": (20, 30),
        "refine_max_helices_sequential": (12, 24),
        "refine_max_helices_pairwise": (12, 24),
        "refine_max_regions": (12, 24),
        "refine_max_seeds": (50, 120),
        "refine_max_solutions": (60, 120),
        "refine_kissing_candidates": (50, 150),
    }
    refine_factors = list(refine_levels.keys())
    l12 = [
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
        [1, 1, 1, 2, 2, 2, 1, 1, 1, 2],
        [1, 1, 1, 2, 2, 2, 2, 2, 2, 1],
        [1, 2, 2, 1, 1, 2, 1, 1, 2, 1],
        [1, 2, 2, 1, 1, 2, 2, 2, 1, 2],
        [1, 2, 2, 2, 2, 1, 1, 2, 1, 1],
        [1, 2, 2, 2, 2, 1, 2, 1, 2, 2],
        [2, 1, 2, 1, 2, 1, 1, 2, 2, 1],
        [2, 1, 2, 1, 2, 1, 2, 1, 1, 2],
        [2, 2, 1, 2, 1, 1, 1, 2, 2, 2],
        [2, 2, 1, 2, 1, 1, 2, 1, 1, 1],
    ]

    for idx, row in enumerate(l12, start=1):
        overrides: dict[str, Any] = {}
        for factor, level in zip(refine_factors, row, strict=True):
            low, high = refine_levels[factor]
            overrides[factor] = low if level == 1 else high
        runs.append(RunConfig(f"refine_L12_{idx:02d}", overrides))

    summary_path = results_dir / "sweep_selected_rnas_summary.json"
    csv_path = results_dir / "sweep_selected_rnas_summary.csv"
    all_results: list[dict[str, Any]] = []

    start = time.time()
    for run in runs:
        if args.max_runs is not None and len(all_results) >= args.max_runs:
            break
        print(f"[SWEEP] Starting {run.name} with overrides={run.overrides}", flush=True)
        run_id = f"{int(time.time())}_{run.name}"
        pred_dir = out_root / run_id
        metrics_path = results_dir / f"metrics_selected_rnas_{run_id}.csv"
        summary_md = results_dir / f"summary_selected_rnas_{run_id}.md"
        log_path = results_dir / f"sweep_{run_id}.log"

        cmd = _build_pred_cmd(base, run.overrides)
        pred_cmd = _cmd_str(cmd)

        ok = True
        with log_path.open("w") as log:
            try:
                subprocess.run(
                    [
                        base["python_bin"],
                        "-m",
                        "ssbench.cli",
                        "predict",
                        "run",
                        "--manifest",
                        str(manifest),
                        "--truth",
                        str(truth_dir),
                        "--out-dir",
                        str(pred_dir),
                        "--k",
                        str(base["top_k"]),
                        "--predictor-cmd",
                        pred_cmd,
                    ],
                    check=True,
                    timeout=args.predict_timeout_s,
                    env={**os.environ, "PYTHONPATH": base["py_path"]},
                    stdout=log,
                    stderr=log,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                ok = False
                log.write(f"\n[SWEEP][ERROR] predict failed: {e}\n")

            if ok:
                try:
                    subprocess.run(
                        [
                            base["python_bin"],
                            "-m",
                            "ssbench.cli",
                            "score",
                            "--manifest",
                            str(manifest),
                            "--truth",
                            str(truth_dir),
                            "--predictions",
                            str(pred_dir),
                            "--out",
                            str(metrics_path),
                        ],
                        check=True,
                        timeout=args.score_timeout_s,
                        env={**os.environ, "PYTHONPATH": base["py_path"]},
                        stdout=log,
                        stderr=log,
                    )
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                    ok = False
                    log.write(f"\n[SWEEP][ERROR] score failed: {e}\n")

            if ok:
                try:
                    subprocess.run(
                        [
                            base["python_bin"],
                            "-m",
                            "ssbench.cli",
                            "report",
                            "--metrics",
                            str(metrics_path),
                            "--out",
                            str(summary_md),
                        ],
                        check=True,
                        timeout=args.report_timeout_s,
                        env={**os.environ, "PYTHONPATH": base["py_path"]},
                        stdout=log,
                        stderr=log,
                    )
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                    ok = False
                    log.write(f"\n[SWEEP][ERROR] report failed: {e}\n")

        means = _mean_metrics(metrics_path) if ok and metrics_path.exists() else {"f1": 0.0, "precision": 0.0, "recall": 0.0}
        result = {
            "run_id": run_id,
            "name": run.name,
            "metrics_path": str(metrics_path),
            "summary_path": str(summary_md),
            "pred_dir": str(pred_dir),
            "f1": means["f1"],
            "precision": means["precision"],
            "recall": means["recall"],
            "ok": ok,
            "overrides": run.overrides,
        }
        all_results.append(result)

        summary_path.write_text(json.dumps(all_results, indent=2) + "\n")
        with csv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "run_id",
                    "name",
                    "ok",
                    "f1",
                    "precision",
                    "recall",
                    "metrics_path",
                    "summary_path",
                    "pred_dir",
                ],
            )
            writer.writeheader()
            for row in all_results:
                writer.writerow({k: row[k] for k in writer.fieldnames})

        elapsed = time.time() - start
        best = max(all_results, key=lambda r: float(r["f1"])) if all_results else None
        if best is not None:
            print(
                f"[SWEEP] Completed {run.name}: F1={means['f1']:.4f} P={means['precision']:.4f} R={means['recall']:.4f}. "
                f"Best so far: {best['name']} F1={best['f1']:.4f}",
                flush=True,
            )
        if elapsed > 6 * 3600:
            break


if __name__ == "__main__":
    main()
