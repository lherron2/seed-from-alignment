#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml


@dataclass(frozen=True)
class Tools:
    python: Path
    allsub: Path
    partition: Path | None
    probplot: Path | None
    datapath: str | None


def _load_tools(defaults_file: Path, python_bin: Path) -> Tools:
    data = yaml.safe_load(defaults_file.read_text()) or {}
    rnastructure = (data.get("paths", {}) or {}).get("rnastructure", {}) or {}
    allsub = Path(rnastructure.get("allsub"))
    partition = rnastructure.get("partition")
    probplot = rnastructure.get("probplot")
    datapath = rnastructure.get("datapath")
    if not allsub.exists():
        raise SystemExit(f"RNAstructure AllSub not found: {allsub}")
    return Tools(
        python=python_bin,
        allsub=allsub,
        partition=Path(partition) if partition else None,
        probplot=Path(probplot) if probplot else None,
        datapath=str(datapath) if datapath else None,
    )


def _run(cmd: list[str], *, env: dict[str, str]) -> None:
    subprocess.run(cmd, check=True, env=env)


def _format_args(args: list[str]) -> str:
    return " ".join(shlex.quote(a) for a in args)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rerank-only ablations on frozen candidates (selected_rnas)"
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--truth-dir", required=True)
    parser.add_argument("--candidates-dir", required=True, help="Frozen candidate pool root")
    parser.add_argument("--out-root", required=True, help="Predictions output root for experiments")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--defaults", required=True)
    parser.add_argument("--python-bin", required=True)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--limit", type=int, default=None, help="Optional limit on number of targets")
    args = parser.parse_args()

    bench_root = Path(__file__).resolve().parents[1]
    repo_root = bench_root.parents[0]

    python_bin = Path(args.python_bin)
    defaults_file = Path(args.defaults)
    tools = _load_tools(defaults_file, python_bin=python_bin)

    env = os.environ.copy()
    env["PYTHONPATH"] = f"{bench_root / 'src'}:{repo_root}"
    if tools.datapath:
        env["DATAPATH"] = tools.datapath

    manifest = pd.read_csv(args.manifest)
    if args.limit is not None:
        manifest = manifest.head(int(args.limit))

    candidates_root = Path(args.candidates_dir)
    out_root = Path(args.out_root)
    results_dir = Path(args.results_dir)
    truth_dir = Path(args.truth_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Baseline values should mirror `benchmark_runner/defaults.yaml` (predict block) for reranking.
    # Override per-ablation by specifying only the changes.
    experiments: list[tuple[str, dict[str, str]]] = [
        ("baseline", {}),
        (
            "refined_only",
            {
                "--candidate-sources": "refined",
            },
        ),
        (
            "sampled_only",
            {
                "--candidate-sources": "sampled",
            },
        ),
        (
            "no_helix_priors",
            {
                "--stem-start-penalty-scale": "0",
                "--stem-len1-penalty-scale": "0",
                "--stem-len2-penalty-scale": "0",
                "--stem-log-reward-scale": "0",
            },
        ),
        (
            "no_pair_penalty",
            {
                "--pair-penalty": "0",
            },
        ),
        (
            "no_calibration",
            {
                "--weight-calibration-method": "none",
            },
        ),
        (
            "cov_off",
            {
                "--cov-mode": "off",
            },
        ),
        (
            "thermo_off",
            {
                "--thermo-mode": "off",
                "--thermo-weight": "0",
            },
        ),
        (
            "cov_off_thermo_off",
            {
                "--cov-mode": "off",
                "--thermo-mode": "off",
                "--thermo-weight": "0",
            },
        ),
    ]

    # Load baseline predict settings for reranking.
    cfg = yaml.safe_load(defaults_file.read_text()) or {}
    predict = (cfg.get("predict", {}) or {}).copy()

    def baseline_flag(name: str, default: str | None) -> str | None:
        v = predict.get(name)
        if v is None:
            return default
        return str(v)

    baseline_flags: dict[str, str] = {
        "--pk-alpha": baseline_flag("pk_alpha", "0.5") or "0.5",
        "--cov-mode": baseline_flag("cov_mode", "logE") or "logE",
        "--cov-alpha": baseline_flag("cov_alpha", "3.0") or "3.0",
        "--cov-min-power": baseline_flag("cov_min_power", "0.1") or "0.1",
        "--weight-calibration-method": baseline_flag("weight_calibration_method", "robust_z")
        or "robust_z",
        "--weight-calibration-zmax": baseline_flag("weight_calibration_zmax", "3.0") or "3.0",
        "--weight-alpha-core": baseline_flag("weight_alpha_core", "1.0") or "1.0",
        "--weight-alpha-alt": baseline_flag("weight_alpha_alt", "1.0") or "1.0",
        "--weight-alpha-cov": baseline_flag("weight_alpha_cov", "1.0") or "1.0",
        "--weight-alpha-thermo": baseline_flag("weight_alpha_thermo", "1.0") or "1.0",
        "--thermo-mode": baseline_flag("thermo_mode", "off") or "off",
        "--thermo-weight": baseline_flag("thermo_weight", "0.0") or "0.0",
        "--thermo-max-structures": baseline_flag("thermo_max_structures", "50") or "50",
        "--thermo-min-count": baseline_flag("thermo_min_count", "2") or "2",
        "--thermo-min-prob": baseline_flag("thermo_min_prob", "0.001") or "0.001",
        "--thermo-log-eps": baseline_flag("thermo_log_eps", "1e-6") or "1e-6",
        "--stem-start-penalty-scale": baseline_flag("stem_start_penalty_scale", "0.05") or "0.05",
        "--stem-len1-penalty-scale": baseline_flag("stem_len1_penalty_scale", "0.10") or "0.10",
        "--stem-len2-penalty-scale": baseline_flag("stem_len2_penalty_scale", "0.02") or "0.02",
        "--stem-log-reward-scale": baseline_flag("stem_log_reward_scale", "0.02") or "0.02",
        "--stem-support-quantile": baseline_flag("stem_support_quantile", "0.5") or "0.5",
        "--pair-penalty-mode": baseline_flag("pair_penalty_mode", "length_aware") or "length_aware",
        "--pair-penalty-c0": baseline_flag("pair_penalty_c0", "0.10") or "0.10",
        "--pair-penalty-c1": baseline_flag("pair_penalty_c1", "0.50") or "0.50",
        "--pair-penalty-scale": baseline_flag("pair_penalty_scale", "0.25") or "0.25",
    }

    exp_rows: list[dict[str, object]] = []
    for exp_name, overrides in experiments:
        exp_pred_dir = out_root / exp_name
        exp_pred_dir.mkdir(parents=True, exist_ok=True)

        print(f"[RUN] experiment={exp_name} -> {exp_pred_dir}")
        for _, row in manifest.iterrows():
            target_id = str(row["target_id"])
            safe_id = target_id.replace("|", "_")
            cand_dir = candidates_root / safe_id
            if not cand_dir.exists():
                raise SystemExit(f"Missing candidate dir for {target_id}: {cand_dir}")

            out_dir = exp_pred_dir / safe_id
            out_dir.mkdir(parents=True, exist_ok=True)

            cmd = [
                str(tools.python),
                "-m",
                "ssbench.predict.rerank_existing_candidates",
                "--candidates-dir",
                str(cand_dir),
                "--out-dir",
                str(out_dir),
                "--allsub-exe",
                str(tools.allsub),
                "--top-k",
                str(args.top_k),
            ]
            if tools.partition is not None:
                cmd += ["--partition-exe", str(tools.partition)]
            if tools.probplot is not None:
                cmd += ["--probplot-exe", str(tools.probplot)]

            # Baseline rerank settings.
            for k, v in baseline_flags.items():
                cmd += [k, v]
            # Experiment overrides.
            for k, v in overrides.items():
                cmd += [k, v]

            _run(cmd, env=env)

        metrics_path = results_dir / f"metrics_ablation_{exp_name}.csv"
        report_path = results_dir / f"summary_ablation_{exp_name}.md"
        _run(
            [
                str(tools.python),
                "-m",
                "ssbench.cli",
                "score",
                "--manifest",
                args.manifest,
                "--truth",
                str(truth_dir),
                "--predictions",
                str(exp_pred_dir),
                "--out",
                str(metrics_path),
            ],
            env=env,
        )
        _run(
            [
                str(tools.python),
                "-m",
                "ssbench.cli",
                "report",
                "--metrics",
                str(metrics_path),
                "--out",
                str(report_path),
            ],
            env=env,
        )

        df = pd.read_csv(metrics_path)
        exp_rows.append(
            {
                "experiment": exp_name,
                "mean_top1_f1": float(df["top1_f1"].mean()),
                "mean_best_of_k_f1": float(df["best_of_k_f1"].mean()),
                "mean_f1": float(df["f1"].mean()),
                "mean_precision": float(df["precision"].mean()),
                "mean_recall": float(df["recall"].mean()),
                "mean_lonely_pair_rate": float(df["lonely_pair_rate"].mean()),
                "mean_gap": float((df["best_of_k_f1"] - df["top1_f1"]).mean()),
                "top_k": int(args.top_k),
            }
        )

    summary_csv = results_dir / "ablation_rerank_summary.csv"
    df_sum = pd.DataFrame(exp_rows).sort_values("mean_top1_f1", ascending=False)
    df_sum.to_csv(summary_csv, index=False)
    print(f"[DONE] Wrote summary: {summary_csv}")
    print(df_sum.to_string(index=False))
    print(f"[CMD] candidates_dir={candidates_root}")
    print(f"[CMD] out_root={out_root}")
    print(f"[CMD] baseline_flags={_format_args([item for kv in baseline_flags.items() for item in kv])}")


if __name__ == "__main__":
    main()

