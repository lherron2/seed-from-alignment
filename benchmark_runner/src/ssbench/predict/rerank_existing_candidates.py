from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path


def _default_cov_path(candidate_dir: Path) -> Path | None:
    for name in ("aln_1.cacofold.cov", "aln.cacofold.cov"):
        p = candidate_dir / name
        if p.exists():
            return p
    return None


def _default_sto_path(candidate_dir: Path) -> Path | None:
    for name in ("aln_1.cacofold.sto", "aln.sto", "fallback.sto"):
        p = candidate_dir / name
        if p.exists():
            return p
    return None


def rerank_existing_candidates(
    *,
    sto: Path,
    cov: Path | None,
    seq_name: str | None,
    allsub_exe: Path,
    partition_exe: Path | None,
    probplot_exe: Path | None,
    candidates_dir: Path,
    out_dir: Path,
    top_k: int,
    pk_alpha: float,
    pair_penalty: float | None,
    pair_penalty_scale: float,
    pair_penalty_mode: str,
    pair_penalty_c0: float,
    pair_penalty_c1: float,
    pair_penalty_min: float | None,
    pair_penalty_max: float | None,
    cov_mode: str,
    cov_alpha: float,
    cov_min_power: float,
    cov_forbid_negative: bool,
    weight_calibration_method: str,
    weight_calibration_zmax: float,
    weight_alpha_core: float,
    weight_alpha_alt: float,
    weight_alpha_cov: float,
    weight_alpha_thermo: float,
    thermo_mode: str,
    thermo_weight: float,
    thermo_max_structures: int,
    thermo_min_count: int,
    thermo_min_prob: float,
    thermo_log_eps: float,
    stem_start_penalty_scale: float,
    stem_len1_penalty_scale: float,
    stem_len2_penalty_scale: float,
    stem_log_reward_scale: float,
    stem_support_quantile: float,
    candidate_sources: str,
) -> None:
    repo_root = Path(__file__).resolve().parents[4]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from src.lib import energy as legacy_energy
    from src.lib import refine_unpaired_regions as rup
    from src.lib import sample_cacofold_structures as scs
    from src.lib import stem_stats as ss
    from src.lib import weight_calibration as wc

    seqs, ss_tracks = scs.parse_stockholm_single(str(sto))
    if not seqs:
        raise SystemExit(f"No sequences found in {sto}")

    if seq_name is None:
        if len(seqs) == 1:
            seq_name = next(iter(seqs))
        else:
            raise SystemExit("--seq-name is required when .sto contains multiple sequences")
    if seq_name not in seqs:
        raise SystemExit(f"Sequence {seq_name!r} not found in {sto}")

    aligned_seq = seqs[seq_name]
    aln2seq, _L = scs.aln_to_seq_map(aligned_seq)
    cov_stats = None
    if cov_mode != "off" and cov is not None and cov.exists():
        cov_stats = scs.load_cov_stats(str(cov), aln2seq)

    L, _ungapped_seq, components = scs.extract_candidate_pair_components(
        seq_name,
        seqs,
        ss_tracks,
        cov=cov_stats,
        cov_mode=cov_mode,
        cov_alpha=cov_alpha,
        cov_min_power=cov_min_power,
        cov_forbid_negative=cov_forbid_negative,
    )

    thermo_component: dict[tuple[int, int], float] = {}
    full_seq = rup.read_ungapped_seq_from_sto(sto, seq_name)

    if float(thermo_weight) > 0 and thermo_mode != "off":
        if thermo_mode == "pf" and partition_exe is not None and probplot_exe is not None:
            from src.lib import rnastructure_pf as rpf

            try:
                probs = rpf.call_rnastructure_pf_probs(
                    partition_exe,
                    probplot_exe,
                    full_seq,
                    extra_args=None,
                    env=os.environ.copy(),
                )
            except Exception:
                probs = {}

            for (i, j), p in probs.items():
                if p < float(thermo_min_prob):
                    continue
                if not scs.is_canonical_pair(full_seq[i], full_seq[j]):
                    continue
                thermo_component[(i, j)] = math.log(p + float(thermo_log_eps))

        if not thermo_component and thermo_mode in {"pf", "allsub"}:
            try:
                thermo_structs = rup.call_rnastructure_allsub(allsub_exe, full_seq)
            except Exception:
                thermo_structs = []

            thermo_structs = thermo_structs[: max(0, int(thermo_max_structures))]
            if thermo_structs:
                counts: dict[tuple[int, int], int] = {}
                for s in thermo_structs:
                    for p in rup._pairs_from_struct(s):
                        counts[p] = counts.get(p, 0) + 1
                for p, c in counts.items():
                    if c < int(thermo_min_count):
                        continue
                    i, j = p
                    if not scs.is_canonical_pair(full_seq[i], full_seq[j]):
                        continue
                    thermo_component[p] = float(c)

    components = dict(components)
    if thermo_component:
        components["thermo"] = thermo_component

    if cov_mode != "off" and cov_alpha != 0.0:
        # Back-compat behavior: interpret cov_alpha as the cov coefficient when cov evidence is enabled.
        alphas = {
            "core": 1.0,
            "alt": 1.0,
            "cov": float(cov_alpha),
            "thermo": float(thermo_weight),
        }
    else:
        alphas = {
            "core": float(weight_alpha_core),
            "alt": float(weight_alpha_alt),
            "cov": float(weight_alpha_cov),
            "thermo": float(weight_alpha_thermo),
        }

    norm_cfg = wc.NormalizationConfig(
        method=str(weight_calibration_method),
        zmax=float(weight_calibration_zmax),
    )
    weights = wc.blend_components(components, alphas=alphas, config=norm_cfg)

    w_scale = wc.weight_scale(weights)
    if pair_penalty is None:
        if pair_penalty_mode == "length_aware":
            length = max(1, int(L))
            pair_penalty = w_scale * (
                float(pair_penalty_c0) + float(pair_penalty_c1) / math.sqrt(length)
            )
            if pair_penalty_min is not None:
                pair_penalty = max(float(pair_penalty_min), pair_penalty)
            if pair_penalty_max is not None:
                pair_penalty = min(float(pair_penalty_max), pair_penalty)
        else:
            pair_penalty = float(pair_penalty_scale) * w_scale if w_scale > 0 else 0.0

    params = legacy_energy.EnergyParams(
        pk_alpha=float(pk_alpha),
        pk_gamma=0.0,
        lonely_penalty=0.0,
    )

    stem_start_penalty = float(stem_start_penalty_scale) * w_scale
    stem_len1_penalty = float(stem_len1_penalty_scale) * w_scale
    stem_len2_penalty = float(stem_len2_penalty_scale) * w_scale
    stem_log_reward = float(stem_log_reward_scale) * w_scale
    support_threshold = wc.quantile_threshold(weights, float(stem_support_quantile))

    def score_struct(struct: str) -> float:
        pairs = rup._pairs_from_struct(struct)
        energy, _ = legacy_energy.compute_energy(pairs, weights, params)
        if pair_penalty:
            energy += float(pair_penalty) * len(pairs)
        if any(
            val != 0.0
            for val in (stem_start_penalty, stem_len1_penalty, stem_len2_penalty, stem_log_reward)
        ):
            stats = ss.stem_stats(pairs)
            unsupported_len1 = sum(
                1 for p in stats.len1_pairs if weights.get(p, 0.0) < support_threshold
            )
            energy += stem_start_penalty * stats.n_stems
            energy += stem_len1_penalty * unsupported_len1
            energy += stem_len2_penalty * stats.n_len2
            energy -= stem_log_reward * stats.sum_log1p_len
        return energy

    sources = {s.strip() for s in candidate_sources.split(",") if s.strip()}
    if not sources.issubset({"refined", "sampled"}):
        raise SystemExit("--candidate-sources must be a comma list of refined,sampled")

    structs: list[str] = []
    if "refined" in sources:
        structs.extend(rup.read_db_structures(candidates_dir / "01_refined.db"))
    if "sampled" in sources:
        structs.extend(rup.read_db_structures(candidates_dir / "02_sampled.db"))

    scored = [(s, score_struct(s)) for s in structs]
    scored.sort(key=lambda x: x[1])

    seen = set()
    top: list[str] = []
    for s, _score in scored:
        if s in seen:
            continue
        seen.add(s)
        top.append(s)
        if len(top) >= int(top_k):
            break

    out_dir.mkdir(parents=True, exist_ok=True)
    out_db = out_dir / "predictions.db"
    out_db.write_text("\n".join([full_seq] + top) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rerank an existing candidate pool (01_refined.db/02_sampled.db) into predictions.db"
    )
    parser.add_argument("--candidates-dir", required=True, help="Per-target candidate directory")
    parser.add_argument("--out-dir", required=True, help="Directory to write predictions.db")
    parser.add_argument("--sto", default=None, help="Stockholm alignment (default: from candidates-dir)")
    parser.add_argument("--cov", default=None, help="Covariation file (default: from candidates-dir)")
    parser.add_argument("--seq-name", default=None, help="Sequence name inside .sto (default: if unique)")
    parser.add_argument("--allsub-exe", required=True, help="Path to RNAstructure AllSub")
    parser.add_argument("--partition-exe", default=None, help="Path to RNAstructure Partition")
    parser.add_argument("--probplot-exe", default=None, help="Path to RNAstructure ProbabilityPlot")
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--pk-alpha", type=float, default=0.5)
    parser.add_argument("--pair-penalty", type=float, default=None)
    parser.add_argument("--pair-penalty-scale", type=float, default=0.25)
    parser.add_argument("--pair-penalty-mode", choices=["legacy", "length_aware"], default="legacy")
    parser.add_argument("--pair-penalty-c0", type=float, default=0.10)
    parser.add_argument("--pair-penalty-c1", type=float, default=0.50)
    parser.add_argument("--pair-penalty-min", type=float, default=None)
    parser.add_argument("--pair-penalty-max", type=float, default=None)
    parser.add_argument("--cov-mode", default="logE_power")
    parser.add_argument("--cov-alpha", type=float, default=3.0)
    parser.add_argument("--cov-min-power", type=float, default=0.1)
    parser.add_argument("--cov-forbid-negative", action="store_true")
    parser.add_argument("--weight-calibration-method", choices=["none", "robust_z"], default="none")
    parser.add_argument("--weight-calibration-zmax", type=float, default=3.0)
    parser.add_argument("--weight-alpha-core", type=float, default=1.0)
    parser.add_argument("--weight-alpha-alt", type=float, default=1.0)
    parser.add_argument("--weight-alpha-cov", type=float, default=1.0)
    parser.add_argument("--weight-alpha-thermo", type=float, default=1.0)
    parser.add_argument("--thermo-mode", choices=["allsub", "pf", "off"], default="off")
    parser.add_argument("--thermo-weight", type=float, default=0.0)
    parser.add_argument("--thermo-max-structures", type=int, default=50)
    parser.add_argument("--thermo-min-count", type=int, default=2)
    parser.add_argument("--thermo-min-prob", type=float, default=0.001)
    parser.add_argument("--thermo-log-eps", type=float, default=1e-6)
    parser.add_argument("--stem-start-penalty-scale", type=float, default=0.0)
    parser.add_argument("--stem-len1-penalty-scale", type=float, default=0.0)
    parser.add_argument("--stem-len2-penalty-scale", type=float, default=0.0)
    parser.add_argument("--stem-log-reward-scale", type=float, default=0.0)
    parser.add_argument("--stem-support-quantile", type=float, default=0.5)
    parser.add_argument(
        "--candidate-sources",
        default="refined,sampled",
        help="Comma list: refined,sampled (default: both).",
    )
    args = parser.parse_args()

    candidates_dir = Path(args.candidates_dir)
    sto = Path(args.sto) if args.sto else _default_sto_path(candidates_dir)
    if sto is None:
        raise SystemExit("Missing .sto; pass --sto or place aln_1.cacofold.sto in --candidates-dir")
    cov = Path(args.cov) if args.cov else _default_cov_path(candidates_dir)

    rerank_existing_candidates(
        sto=sto,
        cov=cov,
        seq_name=args.seq_name,
        allsub_exe=Path(args.allsub_exe),
        partition_exe=Path(args.partition_exe) if args.partition_exe else None,
        probplot_exe=Path(args.probplot_exe) if args.probplot_exe else None,
        candidates_dir=candidates_dir,
        out_dir=Path(args.out_dir),
        top_k=args.top_k,
        pk_alpha=args.pk_alpha,
        pair_penalty=args.pair_penalty,
        pair_penalty_scale=args.pair_penalty_scale,
        pair_penalty_mode=str(args.pair_penalty_mode),
        pair_penalty_c0=float(args.pair_penalty_c0),
        pair_penalty_c1=float(args.pair_penalty_c1),
        pair_penalty_min=args.pair_penalty_min,
        pair_penalty_max=args.pair_penalty_max,
        cov_mode=str(args.cov_mode),
        cov_alpha=float(args.cov_alpha),
        cov_min_power=float(args.cov_min_power),
        cov_forbid_negative=bool(args.cov_forbid_negative),
        weight_calibration_method=str(args.weight_calibration_method),
        weight_calibration_zmax=float(args.weight_calibration_zmax),
        weight_alpha_core=float(args.weight_alpha_core),
        weight_alpha_alt=float(args.weight_alpha_alt),
        weight_alpha_cov=float(args.weight_alpha_cov),
        weight_alpha_thermo=float(args.weight_alpha_thermo),
        thermo_mode=str(args.thermo_mode),
        thermo_weight=float(args.thermo_weight),
        thermo_max_structures=int(args.thermo_max_structures),
        thermo_min_count=int(args.thermo_min_count),
        thermo_min_prob=float(args.thermo_min_prob),
        thermo_log_eps=float(args.thermo_log_eps),
        stem_start_penalty_scale=float(args.stem_start_penalty_scale),
        stem_len1_penalty_scale=float(args.stem_len1_penalty_scale),
        stem_len2_penalty_scale=float(args.stem_len2_penalty_scale),
        stem_log_reward_scale=float(args.stem_log_reward_scale),
        stem_support_quantile=float(args.stem_support_quantile),
        candidate_sources=str(args.candidate_sources),
    )


if __name__ == "__main__":
    main()

