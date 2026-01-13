"""Primary RNAnneal-ss CLI (CaCoFold → ensemble scaffolds → top-K structures).

This is a thin wrapper around `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py`,
using the experiment defaults ("v3-high") that emit a single ranked list of structures where
prefixes (K=1,50,100,200,500,...) can be evaluated for best-of-K oracle metrics.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _repo_root() -> Path:
    # src/lib/rnanneal_ss_cli.py -> repo root
    return Path(__file__).resolve().parents[2]


def _default_path(rel: str) -> str:
    return str((_repo_root() / rel).resolve())


def _add_repo_paths_to_syspath() -> None:
    repo_root = _repo_root()
    bench_src = repo_root / "benchmark_runner" / "src"
    # Ensure we can `import ssbench...` and `import src.lib...` regardless of CWD.
    for p in (bench_src, repo_root):
        sp = str(p)
        if sp not in sys.path:
            sys.path.insert(0, sp)


def _v3_high_args(
    *,
    cm_db: str,
    rscape: str,
    infernal_bin: str | None,
    allsub_exe: str,
    fold_exe: str,
    duplex_exe: str,
    top_k: int,
    reuse_cacofold_root: str | None,
) -> list[str]:
    args: list[str] = [
        "--cm-db",
        cm_db,
        "--rscape",
        rscape,
        "--allsub-exe",
        allsub_exe,
        "--duplex-exe",
        duplex_exe,
        "--fold-exe",
        fold_exe,
        "--top-k",
        str(int(top_k)),
        "--n-samples",
        "400",
        "--burn-in",
        "200",
        "--thin",
        "5",
        "--beta",
        "1.0",
        "--seed",
        "0",
        "--min-loop-sep",
        "0",
        "--pk-alpha",
        "0.5",
        "--pair-penalty-mode",
        "length_aware",
        "--pair-penalty-c0",
        "0.10",
        "--pair-penalty-c1",
        "0.50",
        "--cov-mode",
        "logE_power",
        "--cov-alpha",
        "3.0",
        "--cov-min-power",
        "0.1",
        "--weight-calibration-method",
        "robust_z",
        "--weight-calibration-zmax",
        "3.0",
        "--weight-alpha-core",
        "1.0",
        "--weight-alpha-alt",
        "1.0",
        "--weight-alpha-cov",
        "1.0",
        "--weight-alpha-thermo",
        "1.0",
        "--thermo-mode",
        "allsub",
        "--thermo-weight",
        "1.0",
        "--thermo-max-structures",
        "120",
        "--thermo-min-count",
        "2",
        "--thermo-min-prob",
        "0.001",
        "--thermo-log-eps",
        "1e-6",
        "--stem-start-penalty-scale",
        "0.05",
        "--stem-len1-penalty-scale",
        "0.10",
        "--stem-len2-penalty-scale",
        "0.02",
        "--stem-log-reward-scale",
        "0.02",
        "--stem-support-quantile",
        "0.5",
        "--refine-max-seconds",
        "6",
        "--refine-max-structures",
        "200",
        "--refine-min-unpaired",
        "6",
        "--refine-end-mask-step",
        "10",
        "--refine-max-end-mask-len",
        "30",
        "--refine-max-helices-sequential",
        "20",
        "--refine-max-helices-pairwise",
        "20",
        "--refine-max-regions",
        "4",
        "--refine-max-seeds",
        "50",
        "--refine-max-solutions",
        "50",
        "--refine-kissing-candidates",
        "50",
        "--max-scaffolds",
        "15",
        "--max-samples-per-scaffold",
        "120",
        "--length-adaptive",
        "--no-include-unfixed-sampling",
        "--inject-allsub-scaffolds",
        "--inject-allsub-scaffolds-max",
        "25",
        "--inject-allsub-timeout-s",
        "30",
        # Ensures a minimum quota of thermo suboptimals in the final list (empirically helpful).
        "--force-allsub-output",
        "95",
    ]
    if infernal_bin:
        args += ["--infernal-bin", infernal_bin]
    if reuse_cacofold_root:
        args += ["--reuse-cacofold-root", reuse_cacofold_root]
    return args


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="rnanneal-ss",
        description=(
            "RNAnneal-ss (recommended): Infernal+R-scape CaCoFold → ensemble scaffolds "
            "(RNAstructure+LinearFold+EternaFold) → MCMC sampling → ranked top-K dot-brackets."
        ),
    )
    parser.add_argument("fasta", help="Input FASTA")
    parser.add_argument("outdir", help="Output directory (writes predictions.db)")
    parser.add_argument("--top-k", type=int, default=500, help="Number of structures to emit (default: 500)")
    parser.add_argument(
        "--cm-db",
        default=_default_path("benchmark_runner/data/rfam/Rfam.cm"),
        help="Path to Rfam CM database (default: benchmark_runner/data/rfam/Rfam.cm)",
    )
    parser.add_argument(
        "--infernal-bin",
        default=_default_path("infernal-1.1.5/src"),
        help="Path to Infernal bin dir (default: infernal-1.1.5/src)",
    )
    parser.add_argument(
        "--rscape",
        default=_default_path("rscape_v2.6.4/src/R-scape"),
        help="Path to R-scape binary (default: rscape_v2.6.4/src/R-scape)",
    )
    parser.add_argument(
        "--reuse-cacofold-root",
        default=None,
        help="Optional directory with cached per-target CaCoFold outputs to reuse.",
    )
    parser.add_argument(
        "--allsub-exe",
        default=_default_path("benchmark_runner/tools/scaffold_backends/ensemble3/AllSub"),
        help="AllSub backend (default: benchmark_runner/tools/scaffold_backends/ensemble3/AllSub)",
    )
    parser.add_argument(
        "--fold-exe",
        default=_default_path("benchmark_runner/tools/scaffold_backends/ensemble3/Fold"),
        help="Fold backend (default: benchmark_runner/tools/scaffold_backends/ensemble3/Fold)",
    )
    parser.add_argument(
        "--duplex-exe",
        default=_default_path("RNAstructure/exe/DuplexFold"),
        help="Path to RNAstructure DuplexFold (default: RNAstructure/exe/DuplexFold)",
    )
    args = parser.parse_args(argv)

    _add_repo_paths_to_syspath()
    from ssbench.predict import cacofold_mcmc_pipeline  # type: ignore[import-not-found]

    forwarded = _v3_high_args(
        cm_db=str(args.cm_db),
        rscape=str(args.rscape),
        infernal_bin=str(args.infernal_bin) if args.infernal_bin else None,
        allsub_exe=str(args.allsub_exe),
        fold_exe=str(args.fold_exe),
        duplex_exe=str(args.duplex_exe),
        top_k=int(args.top_k),
        reuse_cacofold_root=str(args.reuse_cacofold_root) if args.reuse_cacofold_root else None,
    )
    # Delegate implementation and file formats to the benchmark-proven pipeline.
    sys.argv = ["python", *forwarded, str(args.fasta), str(args.outdir)]
    cacofold_mcmc_pipeline.main()


if __name__ == "__main__":
    main()

