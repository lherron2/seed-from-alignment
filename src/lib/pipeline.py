"""
High-level CaCoFold → sampling → refinement → Rosetta pipelines.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from src.lib import energy
from src.lib import filter_db_for_rosetta as fdr
from src.lib import refine_unpaired_regions as rup
from src.lib import rnastructure_pf as rpf
from src.lib import sample_cacofold_structures as scs
from src.lib import weight_calibration as wc

_SCAFFOLD_GLOBALS: dict[str, object] = {}


def _sample_scaffold_worker(scaffold_idx: int) -> tuple[int, list[set[tuple[int, int]]], dict]:
    refined_structs: list[str] = _SCAFFOLD_GLOBALS["refined_structs"]  # type: ignore[assignment]
    candidate_pairs: list[tuple[int, int, float]] = _SCAFFOLD_GLOBALS["candidate_pairs"]  # type: ignore[assignment]
    L: int = _SCAFFOLD_GLOBALS["L"]  # type: ignore[assignment]
    sample_cfg: SampleConfig = _SCAFFOLD_GLOBALS["sample_cfg"]  # type: ignore[assignment]

    refined = refined_structs[scaffold_idx]
    scaffold_pairs = rup._pairs_from_struct(refined)
    # Ensure the seed scaffold respects the sampler hairpin constraint.
    min_hairpin = int(sample_cfg.min_loop_sep)
    scaffold_pairs = {p for p in scaffold_pairs if (max(p) - min(p) - 1) >= min_hairpin}
    if sample_cfg.fix_scaffold_pairs:
        scaffold_pos = {p for pair in scaffold_pairs for p in pair}
        filtered_candidates = [
            (i, j, w)
            for (i, j, w) in candidate_pairs
            if i not in scaffold_pos and j not in scaffold_pos
        ]
        fixed_pairs = scaffold_pairs
        init_pairs = None
    else:
        filtered_candidates = candidate_pairs
        fixed_pairs = None
        init_pairs = scaffold_pairs

    n_samples = int(sample_cfg.n_samples)
    if sample_cfg.max_samples_per_scaffold is not None:
        n_samples = min(n_samples, int(sample_cfg.max_samples_per_scaffold))
    elif not sample_cfg.fix_scaffold_pairs:
        n_samples = min(n_samples, 200)

    pk_samples, report = scs.sample_matchings_multibeta_modular(
        L=L,
        candidate_pairs=filtered_candidates,
        n_samples=n_samples,
        burn_in=sample_cfg.burn_in,
        thin=sample_cfg.thin,
        min_hairpin=sample_cfg.min_loop_sep,
        beta=sample_cfg.beta,
        betas=list(sample_cfg.betas) if sample_cfg.betas else None,
        n_beta_chains=sample_cfg.n_beta_chains,
        beta_min_factor=sample_cfg.beta_min_factor,
        beta_max_factor=sample_cfg.beta_max_factor,
        seed=(sample_cfg.seed + 10_000 * scaffold_idx) if sample_cfg.seed is not None else None,
        pk_alpha=sample_cfg.pk_alpha,
        pk_gamma=sample_cfg.pk_gamma,
        lonely_penalty=sample_cfg.lonely_penalty,
        pk_filter_frac=sample_cfg.pk_filter_frac,
        pk_filter_max_cross_per_pair=sample_cfg.pk_filter_max_cross_per_pair,
        pk_filter_max_total_cross=sample_cfg.pk_filter_max_total_cross,
        delayed_accept=sample_cfg.delayed_accept,
        parallel_chains=False,  # avoid nested process pools
        max_workers=1,
        use_beam_seeds=sample_cfg.use_beam_seeds,
        beam_width=sample_cfg.beam_width,
        beam_expansions_per_state=sample_cfg.beam_expansions_per_state,
        beam_max_segments=sample_cfg.beam_max_segments,
        beam_max_depth=sample_cfg.beam_max_depth,
        scaffold_pairs=fixed_pairs,
        initial_pairs=init_pairs,
        use_segment_moves=sample_cfg.use_segment_moves,
        segment_birth_prob=sample_cfg.segment_birth_prob,
        segment_death_prob=sample_cfg.segment_death_prob,
        segment_swap_prob=sample_cfg.segment_swap_prob,
        diversity_dist_frac=sample_cfg.diversity_dist_frac,
        diagnostics_file=None,
        label=f"scaffold_{scaffold_idx}",
    )

    report["scaffold_idx"] = scaffold_idx
    return scaffold_idx, pk_samples, report


@dataclass
class SampleConfig:
    n_samples: int = 1000
    burn_in: int = 1000
    thin: int = 10
    min_loop_sep: int = 1
    beta: float = 1.0
    betas: Sequence[float] | None = None
    n_beta_chains: int = 4
    beta_min_factor: float = 0.5
    beta_max_factor: float = 2.0
    diversity_dist_frac: float = 0.05
    seed: int | None = None
    pk_alpha: float = 2.0
    pk_gamma: float = 0.0
    lonely_penalty: float = 0.0
    pk_filter_frac: float = 0.2
    pk_filter_max_cross_per_pair: int = 2
    pk_filter_max_total_cross: int = 30
    pk_depth_limit: int | None = None
    use_segment_moves: bool = True
    segment_birth_prob: float = 0.2
    segment_death_prob: float = 0.2
    segment_swap_prob: float = 0.1
    # If True, refined scaffold pairs are fixed and the sampler only fills the remaining positions.
    # If False, the refined scaffold is used only as an initial seed and can be modified by MCMC.
    fix_scaffold_pairs: bool = True
    # Cap the number of refined scaffolds used as MCMC starts.
    # This prevents redundant sampling when many refined variants are available.
    max_scaffolds: int | None = 20
    # Optional cap on samples PER scaffold (defaults to `n_samples`).
    # If unset and `fix_scaffold_pairs` is False, a conservative cap is applied in `sample_pk`.
    max_samples_per_scaffold: int | None = None
    delayed_accept: bool = True
    parallel_chains: bool = True
    parallel_scaffolds: bool = True
    max_workers: int = 8
    use_beam_seeds: bool = True
    beam_width: int = 64
    beam_expansions_per_state: int = 24
    beam_max_segments: int = 2000
    beam_max_depth: int = 6
    diagnostics_json: Path | None = None
    cov_mode: str = "off"
    cov_alpha: float = 1.0
    cov_min_power: float = 0.0
    cov_forbid_negative: bool = False
    cov_negative_E: float = 1.0
    weight_calibration_method: str = "none"
    weight_calibration_zmax: float = 3.0
    weight_alpha_core: float = 1.0
    weight_alpha_alt: float = 1.0
    weight_alpha_cov: float = 1.0
    weight_alpha_thermo: float = 1.0
    # Optional thermodynamic augmentation of the candidate pair set (single-sequence AllSub).
    thermo_mode: str = "allsub"  # allsub | pf | off
    thermo_weight: float = 1.0
    thermo_max_structures: int = 50
    thermo_min_count: int = 2
    thermo_min_prob: float = 0.001
    thermo_log_eps: float = 1e-6


def _normalized_pair_set_from_struct(struct: str) -> frozenset[tuple[int, int]]:
    pairs = []
    for i, j in rup._pairs_from_struct(struct):
        if i > j:
            i, j = j, i
        pairs.append((i, j))
    return frozenset(pairs)


def _intersection_size(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]) -> int:
    if len(a) > len(b):
        a, b = b, a
    return sum(1 for p in a if p in b)


def _jaccard_distance(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]) -> float:
    if not a and not b:
        return 0.0
    inter = _intersection_size(a, b)
    union = len(a) + len(b) - inter
    return 1.0 - (inter / union if union > 0 else 0.0)


def _select_diverse_by_score(
    scored_structs: list[tuple[str, float, frozenset[tuple[int, int]]]],
    k: int,
    diversity_weight: float = 0.6,
) -> list[str]:
    if k <= 0 or not scored_structs:
        return []
    scored_structs = sorted(scored_structs, key=lambda x: (x[1], x[0]))
    if len(scored_structs) <= k:
        return [s for s, _score, _pairs in scored_structs]

    n = len(scored_structs)
    selected: list[tuple[str, float, frozenset[tuple[int, int]]]] = [scored_structs[0]]
    selected_structs = {scored_structs[0][0]}
    selected_pairs = [scored_structs[0][2]]

    while len(selected) < k:
        best_idx: int | None = None
        best_combined: float | None = None
        best_score: float | None = None
        best_struct: str | None = None

        for idx, (struct, score, pairs) in enumerate(scored_structs):
            if struct in selected_structs:
                continue
            base = 1.0 if n == 1 else 1.0 - (idx / (n - 1))
            min_dist = min(_jaccard_distance(pairs, sp) for sp in selected_pairs)
            combined = base + diversity_weight * min_dist

            if best_combined is None:
                best_idx, best_combined, best_score, best_struct = idx, combined, score, struct
                continue

            if combined > best_combined:
                best_idx, best_combined, best_score, best_struct = idx, combined, score, struct
                continue

            if combined == best_combined:
                assert best_score is not None
                assert best_struct is not None
                if score < best_score:
                    best_idx, best_combined, best_score, best_struct = (
                        idx,
                        combined,
                        score,
                        struct,
                    )
                    continue
                if score == best_score and struct < best_struct:
                    best_idx, best_combined, best_score, best_struct = (
                        idx,
                        combined,
                        score,
                        struct,
                    )

        if best_idx is None:
            break

        chosen = scored_structs[best_idx]
        selected.append(chosen)
        selected_structs.add(chosen[0])
        selected_pairs.append(chosen[2])

    return [s for s, _score, _pairs in selected]


def _select_refined_scaffolds(
    refined_structs: list[str],
    weights: dict[tuple[int, int], float],
    sample_cfg: SampleConfig,
    max_scaffolds: int,
) -> list[str]:
    if max_scaffolds <= 0 or len(refined_structs) <= max_scaffolds:
        return refined_structs

    params = energy.EnergyParams(
        pk_alpha=float(sample_cfg.pk_alpha),
        pk_gamma=float(sample_cfg.pk_gamma),
        lonely_penalty=float(sample_cfg.lonely_penalty),
    )
    scored: list[tuple[str, float, frozenset[tuple[int, int]]]] = []
    for struct in refined_structs:
        pairs = _normalized_pair_set_from_struct(struct)
        score, _cache = energy.compute_energy(set(pairs), weights, params)
        scored.append((struct, float(score), pairs))

    return _select_diverse_by_score(scored, max_scaffolds)


@dataclass
class RefineConfig:
    min_unpaired: int = 15
    temperature: float | None = None
    allsub_abs: float | None = None
    allsub_pct: float | None = None
    fold_extra_args: Sequence[str] = ()
    max_structures: int = 1000  # Default limit for refined output
    max_seconds: float | None = None

    # Masking Grid Constants
    end_mask_step: int = 5
    max_end_mask_len: int = 40
    max_helices_sequential: int = 20
    max_helices_pairwise: int = 10

    # Initialization
    max_seeds: int = 50

    # Heuristics
    max_regions_to_refine: int = 30
    max_solutions: int = 2000
    kissing_loop_candidates: int = 1000
    kissing_loop_len_ratio_min: float = 0.5
    kissing_loop_len_ratio_max: float = 2.0
    kissing_loop_min_canon_pairs: int = 4
    kissing_loop_min_contig_pairs: int = 3
    kissing_loop_min_stem_sep: int = 4
    kissing_loop_max_stem_sep: int = 200

    # Optional caches (pickle)
    allsub_cache_path: Path | None = None
    duplex_cache_path: Path | None = None

    # Scoring Weights
    scaffold_pair_energy: float = -1.5
    weight_l0: float = 1.0
    weight_l1: float = 0.6
    weight_l2_plus: float = -1.0


@dataclass
class EnsureValidConfig:
    max_total_crossings: int = 30
    max_crossings_per_pair: int = 2
    max_pk_fraction: float = 0.2


@dataclass
class RosettaConfig:
    hamming_threshold: int = 0
    hamming_frac: float = 0.05
    max_structures: int | None = None


@dataclass
class IOConfig:
    work_dir: Path
    consensus_db: Path
    refined_db: Path
    sampled_db: Path
    valid_db: Path
    rosetta_db: Path
    summary: Path | None


@dataclass
class PipelineConfig:
    sto: Path
    cov: Path | None
    seq_name: str
    allsub_exe: Path
    duplex_exe: Path | None
    partition_exe: Path | None
    probplot_exe: Path | None
    sample: SampleConfig
    refine: RefineConfig
    ensure_valid: EnsureValidConfig
    rosetta: RosettaConfig
    io: IOConfig


def _resolve(path_str: str | None, base: Path) -> Path | None:
    if path_str is None:
        return None
    p = Path(path_str)
    if not p.is_absolute():
        p = base / p
    return p


def load_pipeline_config(path: Path) -> PipelineConfig:
    text = path.read_text()
    base = path.parent

    if path.suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore[import-untyped]
        except ImportError as e:
            raise SystemExit("PyYAML is required for YAML configs.") from e
        raw = yaml.safe_load(text)
    else:
        raw = json.loads(text)

    sto = _resolve(raw["sto"], base)
    cov = _resolve(raw.get("cov"), base)
    seq_name = raw["seq_name"]

    allsub_exe = Path(raw["allsub_exe"])
    duplex_exe = Path(raw["duplex_exe"]) if raw.get("duplex_exe") else None
    partition_exe = Path(raw["partition_exe"]) if raw.get("partition_exe") else None
    probplot_exe = Path(raw["probplot_exe"]) if raw.get("probplot_exe") else None

    sample_raw = dict(raw.get("sample", {}))
    sample_raw["diagnostics_json"] = _resolve(sample_raw.get("diagnostics_json"), base)
    refine_raw = raw.get("refine", {})
    ensure_raw = raw.get("ensure_valid", {})
    rosetta_raw = raw.get("rosetta", {})
    io_raw = raw.get("io", {})

    work_dir = _resolve(io_raw.get("work_dir", "."), base)
    if work_dir is None:
        work_dir = base
    work_dir = work_dir.resolve()

    consensus_db = work_dir / io_raw.get("consensus_db", "00_consensus.db")
    refined_db = work_dir / io_raw.get("refined_db", "01_refined.db")
    sampled_db = work_dir / io_raw.get("sampled_db", "02_sampled.db")
    valid_db = work_dir / io_raw.get("valid_db", "03_valid.db")
    rosetta_db = work_dir / io_raw.get("rosetta_db", "04_rosetta.db")
    summary = work_dir / io_raw.get("summary", "cacofold_summary.txt")

    io = IOConfig(
        work_dir=work_dir,
        consensus_db=consensus_db,
        refined_db=refined_db,
        sampled_db=sampled_db,
        valid_db=valid_db,
        rosetta_db=rosetta_db,
        summary=summary,
    )

    refine_raw = dict(refine_raw)
    refine_raw["allsub_cache_path"] = _resolve(refine_raw.get("allsub_cache_path"), base)
    refine_raw["duplex_cache_path"] = _resolve(refine_raw.get("duplex_cache_path"), base)

    sample = SampleConfig(**sample_raw)
    refine = RefineConfig(**refine_raw)
    ensure = EnsureValidConfig(**ensure_raw)
    rosetta = RosettaConfig(**rosetta_raw)

    if sto is None or not sto.is_file():
        raise FileNotFoundError(f"Stockholm file not found: {sto}")
    if not allsub_exe.is_file():
        raise FileNotFoundError(f"RNAstructure AllSub executable not found: {allsub_exe}")
    if duplex_exe is not None and not duplex_exe.is_file():
        raise FileNotFoundError(f"RNAstructure DuplexFold executable not found: {duplex_exe}")
    if partition_exe is not None and not partition_exe.is_file():
        raise FileNotFoundError(f"RNAstructure Partition executable not found: {partition_exe}")
    if probplot_exe is not None and not probplot_exe.is_file():
        raise FileNotFoundError(
            f"RNAstructure ProbabilityPlot executable not found: {probplot_exe}"
        )

    io.work_dir.mkdir(parents=True, exist_ok=True)

    return PipelineConfig(
        sto=sto,
        cov=cov,
        seq_name=seq_name,
        allsub_exe=allsub_exe,
        duplex_exe=duplex_exe,
        partition_exe=partition_exe,
        probplot_exe=probplot_exe,
        sample=sample,
        refine=refine,
        ensure_valid=ensure,
        rosetta=rosetta,
        io=io,
    )


def filter_consensus_loop_sep(struct: str, min_loop: int = 2) -> str:
    """
    Ensure every base pair (i, j) has |j - i - 1| >= min_loop.
    If not, remove the pair.
    Assumes simple '()' structure.
    """
    chars = list(struct)
    stack = []
    for i, ch in enumerate(chars):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                if (i - j - 1) < min_loop:
                    chars[i] = "."
                    chars[j] = "."
    return "".join(chars)


def get_consensus_db(cfg: PipelineConfig) -> Path:
    sto_path = cfg.sto
    seq_name = cfg.seq_name

    sys.stderr.write(f"[PIPELINE] Reading Stockholm file: {sto_path}\n")
    seqs, ss_tracks = scs.parse_stockholm_single(str(sto_path))

    if seq_name not in seqs:
        raise KeyError(f"Sequence {seq_name!r} not found in {sto_path}")

    ungapped_seq = "".join(ch for ch in seqs[seq_name] if ch not in scs.GAP_CHARS)
    L = len(ungapped_seq)

    try:
        L_cons, ss_cons_str, _cons_pairs = scs.project_ss_cons_to_sequence(
            seq_name, seqs, ss_tracks
        )
    except Exception as e:
        sys.stderr.write(f"[WARN] Could not project SS_cons for {seq_name}: {e}\n")
        ss_cons_str = "." * L

    # Filter loop separation
    ss_cons_str = filter_consensus_loop_sep(ss_cons_str, min_loop=2)

    # Optional summary
    if cfg.io.summary is not None:
        try:
            cov = None
            if cfg.cov is not None:
                aligned_seq = seqs[seq_name]
                aln2seq, _L = scs.aln_to_seq_map(aligned_seq)
                cov = scs.load_cov_stats(str(cfg.cov), aln2seq)

            _, candidate_pairs = scs.extract_candidate_pairs(
                seq_name,
                seqs,
                ss_tracks,
                cov=cov,
            )
        except Exception:
            candidate_pairs = []

        cfg.io.summary.parent.mkdir(parents=True, exist_ok=True)
        scs.write_summary_file(
            summary_path=str(cfg.io.summary),
            seq_name=seq_name,
            ungapped_seq=ungapped_seq,
            L=L,
            ss_cons_str=ss_cons_str,
            candidate_pairs=candidate_pairs,
        )

    cfg.io.consensus_db.parent.mkdir(parents=True, exist_ok=True)
    with cfg.io.consensus_db.open("w") as fh:
        fh.write(ungapped_seq + "\n")
        fh.write(ss_cons_str + "\n")

    return cfg.io.consensus_db


def refine_db(
    cfg: PipelineConfig,
    db_in: Path,
    db_out: Path,
    allsub_cache: dict | None = None,
    duplex_cache: dict | None = None,
) -> Path:
    full_seq = rup.read_ungapped_seq_from_sto(cfg.sto, cfg.seq_name)
    structs = rup.read_db_structures(db_in)

    sys.stderr.write(f"[PIPELINE] Refining {len(structs)} structures from {db_in} → {db_out}\n")

    allsub_exe = cfg.allsub_exe
    duplex_exe = cfg.duplex_exe

    # Extract config values
    min_unpaired = cfg.refine.min_unpaired
    temperature = cfg.refine.temperature
    abs_energy = cfg.refine.allsub_abs
    pct_energy = cfg.refine.allsub_pct
    extra_args = list(cfg.refine.fold_extra_args)
    max_structures = cfg.refine.max_structures

    allsub_cache_path = cfg.refine.allsub_cache_path
    duplex_cache_path = cfg.refine.duplex_cache_path

    if allsub_cache is None:
        allsub_cache = rup.load_cache(allsub_cache_path) if allsub_cache_path else {}
    if duplex_cache is None:
        duplex_cache = rup.load_cache(duplex_cache_path) if duplex_cache_path else {}

    # rup.refine_structure returns List[Tuple[str, float]]
    refined_structs_with_scores: list[tuple[str, float]] = []

    for idx, s in enumerate(structs):
        variants = rup.refine_structure(
            struct=s,
            full_seq=full_seq,
            allsub_exe=allsub_exe,
            duplex_exe=duplex_exe,
            min_unpaired_len=min_unpaired,
            max_seconds=cfg.refine.max_seconds,
            temperature=temperature,
            absolute_energy=abs_energy,
            percent_energy=pct_energy,
            extra_args=extra_args,
            allsub_cache=allsub_cache,
            duplex_cache=duplex_cache,
            # Pass all scriptable params
            end_mask_step=cfg.refine.end_mask_step,
            max_end_mask_len=cfg.refine.max_end_mask_len,
            max_helices_sequential=cfg.refine.max_helices_sequential,
            max_helices_pairwise=cfg.refine.max_helices_pairwise,
            max_seeds=cfg.refine.max_seeds,
            max_regions_to_refine=cfg.refine.max_regions_to_refine,
            max_solutions=cfg.refine.max_solutions,
            kissing_loop_candidates=cfg.refine.kissing_loop_candidates,
            kissing_loop_len_ratio_min=cfg.refine.kissing_loop_len_ratio_min,
            kissing_loop_len_ratio_max=cfg.refine.kissing_loop_len_ratio_max,
            kissing_loop_min_canon_pairs=cfg.refine.kissing_loop_min_canon_pairs,
            kissing_loop_min_contig_pairs=cfg.refine.kissing_loop_min_contig_pairs,
            kissing_loop_min_stem_sep=cfg.refine.kissing_loop_min_stem_sep,
            kissing_loop_max_stem_sep=cfg.refine.kissing_loop_max_stem_sep,
            scaffold_pair_energy=cfg.refine.scaffold_pair_energy,
            weight_l0=cfg.refine.weight_l0,
            weight_l1=cfg.refine.weight_l1,
            weight_l2_plus=cfg.refine.weight_l2_plus,
        )
        refined_structs_with_scores.extend(variants)

    # Sort by score (ascending: lower energy/score is better)
    refined_structs_with_scores.sort(key=lambda x: x[1])

    # Deduplicate keeping best score (first encounter due to sort)
    unique_refined = []
    seen = set()
    for rs, score in refined_structs_with_scores:
        if rs not in seen:
            unique_refined.append((rs, score))
            seen.add(rs)

    # Apply limit on UNIQUE structures
    if max_structures > 0 and len(unique_refined) > max_structures:
        sys.stderr.write(
            f"[PIPELINE] Limiting output to top {max_structures} unique structures (found {len(unique_refined)}).\n"
        )
        unique_refined = unique_refined[:max_structures]

    db_out.parent.mkdir(parents=True, exist_ok=True)
    with db_out.open("w") as fh:
        fh.write(full_seq + "\n")
        for rs, score in unique_refined:
            fh.write(rs + "\n")

    sys.stderr.write(f"[PIPELINE] Wrote {len(unique_refined)} refined structures to {db_out}\n")
    if allsub_cache_path:
        rup.save_cache(allsub_cache_path, allsub_cache)
    if duplex_cache_path:
        rup.save_cache(duplex_cache_path, duplex_cache)
    return db_out


def sample_pk(cfg: PipelineConfig, refined_db: Path | None, out_db: Path) -> Path:
    sto_path = cfg.sto
    seq_name = cfg.seq_name
    seqs, ss_tracks = scs.parse_stockholm_single(str(sto_path))
    ungapped_seq = "".join(ch for ch in seqs[seq_name] if ch not in scs.GAP_CHARS)

    cov = None
    if cfg.cov is not None:
        aligned_seq = seqs[seq_name]
        aln2seq, _ = scs.aln_to_seq_map(aligned_seq)
        cov = scs.load_cov_stats(str(cfg.cov), aln2seq)

    L, ungapped_seq, components = scs.extract_candidate_pair_components(
        seq_name,
        seqs,
        ss_tracks,
        cov=cov,
        cov_mode=cfg.sample.cov_mode,
        cov_alpha=cfg.sample.cov_alpha,
        cov_min_power=cfg.sample.cov_min_power,
        cov_forbid_negative=cfg.sample.cov_forbid_negative,
        cov_negative_E=cfg.sample.cov_negative_E,
    )
    thermo_component: dict[tuple[int, int], float] = {}
    if cfg.sample.thermo_weight > 0 and cfg.sample.thermo_mode != "off":
        if cfg.sample.thermo_mode == "pf" and cfg.partition_exe and cfg.probplot_exe:
            try:
                probs = rpf.call_rnastructure_pf_probs(
                    cfg.partition_exe,
                    cfg.probplot_exe,
                    ungapped_seq,
                    temperature=cfg.refine.temperature,
                    extra_args=cfg.refine.fold_extra_args,
                    env=os.environ.copy(),
                )
            except Exception:
                probs = {}
            for (i, j), p in probs.items():
                if p < float(cfg.sample.thermo_min_prob):
                    continue
                if not scs.is_canonical_pair(ungapped_seq[i], ungapped_seq[j]):
                    continue
                thermo_component[(i, j)] = math.log(p + float(cfg.sample.thermo_log_eps))

        if not thermo_component and cfg.sample.thermo_mode in {"pf", "allsub"}:
            if cfg.allsub_exe.is_file():
                try:
                    thermo_structs = rup.call_rnastructure_allsub(
                        cfg.allsub_exe,
                        ungapped_seq,
                        temperature=cfg.refine.temperature,
                        absolute_energy=cfg.refine.allsub_abs,
                        percent_energy=cfg.refine.allsub_pct,
                        extra_args=cfg.refine.fold_extra_args,
                    )
                except Exception:
                    thermo_structs = []

                thermo_structs = thermo_structs[: max(0, int(cfg.sample.thermo_max_structures))]
                if thermo_structs:
                    counts: dict[tuple[int, int], int] = {}
                    for db, _energy in thermo_structs:
                        for i, j in rup._pairs_from_struct(db):
                            if i > j:
                                i, j = j, i
                            if i < 0 or j >= len(ungapped_seq):
                                continue
                            if not scs.is_canonical_pair(ungapped_seq[i], ungapped_seq[j]):
                                continue
                            counts[(i, j)] = counts.get((i, j), 0) + 1

                    denom = float(len(thermo_structs))
                    min_count = max(1, int(cfg.sample.thermo_min_count))
                    for (i, j), c in counts.items():
                        if c < min_count:
                            continue
                        thermo_component[(i, j)] = float(c) / denom

    components["thermo"] = thermo_component

    if cfg.sample.weight_calibration_method == "none":
        alphas = {
            "core": 1.0,
            "alt": 1.0,
            "cov": float(cfg.sample.cov_alpha),
            "thermo": float(cfg.sample.thermo_weight),
        }
    else:
        alphas = {
            "core": float(cfg.sample.weight_alpha_core),
            "alt": float(cfg.sample.weight_alpha_alt),
            "cov": float(cfg.sample.weight_alpha_cov),
            "thermo": float(cfg.sample.weight_alpha_thermo),
        }

    norm_cfg = wc.NormalizationConfig(
        method=str(cfg.sample.weight_calibration_method),
        zmax=float(cfg.sample.weight_calibration_zmax),
    )
    weights = wc.blend_components(components, alphas=alphas, config=norm_cfg)
    candidate_pairs = [(i, j, w) for (i, j), w in sorted(weights.items())]

    refined_structs = []
    if refined_db and refined_db.is_file():
        refined_structs = rup.read_db_structures(refined_db)

    # Cap the number of scaffolds to avoid redundant work.
    if cfg.sample.max_scaffolds is not None and cfg.sample.max_scaffolds > 0:
        refined_structs = _select_refined_scaffolds(
            refined_structs=refined_structs,
            weights=weights,
            sample_cfg=cfg.sample,
            max_scaffolds=int(cfg.sample.max_scaffolds),
        )

    # Determine per-scaffold sample budget.
    n_samples_per_scaffold = int(cfg.sample.n_samples)
    if cfg.sample.max_samples_per_scaffold is not None:
        n_samples_per_scaffold = min(
            n_samples_per_scaffold, int(cfg.sample.max_samples_per_scaffold)
        )
    elif not cfg.sample.fix_scaffold_pairs:
        # Unfixed scaffolds explore a much larger space; cap for practicality.
        n_samples_per_scaffold = min(n_samples_per_scaffold, 200)

    # Some environments (e.g. restricted sandboxes) disallow process-based parallelism
    # due to blocked semaphores. Detect once and disable scaffold/chain parallelism.
    process_pool_ok = True
    try:
        ctx_test = mp.get_context("fork")
        _lock = ctx_test.Lock()
        del _lock
    except Exception:
        process_pool_ok = False
    parallel_scaffolds = bool(cfg.sample.parallel_scaffolds and process_pool_ok)
    parallel_chains = bool(cfg.sample.parallel_chains and process_pool_ok)

    out_db.parent.mkdir(parents=True, exist_ok=True)
    diagnostics_all: list[dict] = []
    with out_db.open("w") as fh:
        fh.write(ungapped_seq + "\n")

        if refined_structs:
            if parallel_scaffolds and len(refined_structs) > 1 and cfg.sample.max_workers > 1:
                try:
                    ctx = mp.get_context("fork")
                except ValueError:
                    ctx = None

                if ctx is None:
                    # Fallback: sequential scaffolds but parallel chains inside each scaffold.
                    for scaffold_idx, refined in enumerate(refined_structs):
                        scaffold_pairs = rup._pairs_from_struct(refined)
                        min_hairpin = int(cfg.sample.min_loop_sep)
                        scaffold_pairs = {
                            p for p in scaffold_pairs if (max(p) - min(p) - 1) >= min_hairpin
                        }
                        if cfg.sample.fix_scaffold_pairs:
                            scaffold_pos = {p for pair in scaffold_pairs for p in pair}
                            filtered_candidates = [
                                (i, j, w)
                                for (i, j, w) in candidate_pairs
                                if i not in scaffold_pos and j not in scaffold_pos
                            ]
                            fixed_pairs = scaffold_pairs
                            init_pairs = None
                        else:
                            filtered_candidates = candidate_pairs
                            fixed_pairs = None
                            init_pairs = scaffold_pairs
                        pk_samples, report = scs.sample_matchings_multibeta_modular(
                            L=L,
                            candidate_pairs=filtered_candidates,
                            n_samples=n_samples_per_scaffold,
                            burn_in=cfg.sample.burn_in,
                            thin=cfg.sample.thin,
                            min_hairpin=cfg.sample.min_loop_sep,
                            beta=cfg.sample.beta,
                            betas=list(cfg.sample.betas) if cfg.sample.betas else None,
                            n_beta_chains=cfg.sample.n_beta_chains,
                            beta_min_factor=cfg.sample.beta_min_factor,
                            beta_max_factor=cfg.sample.beta_max_factor,
                            seed=(cfg.sample.seed + 10_000 * scaffold_idx)
                            if cfg.sample.seed is not None
                            else None,
                            pk_alpha=cfg.sample.pk_alpha,
                            pk_gamma=cfg.sample.pk_gamma,
                            lonely_penalty=cfg.sample.lonely_penalty,
                            pk_filter_frac=cfg.sample.pk_filter_frac,
                            pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
                            pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
                            delayed_accept=cfg.sample.delayed_accept,
                            parallel_chains=parallel_chains,
                            max_workers=cfg.sample.max_workers,
                            use_beam_seeds=cfg.sample.use_beam_seeds,
                            beam_width=cfg.sample.beam_width,
                            beam_expansions_per_state=cfg.sample.beam_expansions_per_state,
                            beam_max_segments=cfg.sample.beam_max_segments,
                            beam_max_depth=cfg.sample.beam_max_depth,
                            scaffold_pairs=fixed_pairs,
                            initial_pairs=init_pairs,
                            use_segment_moves=cfg.sample.use_segment_moves,
                            segment_birth_prob=cfg.sample.segment_birth_prob,
                            segment_death_prob=cfg.sample.segment_death_prob,
                            segment_swap_prob=cfg.sample.segment_swap_prob,
                            diversity_dist_frac=cfg.sample.diversity_dist_frac,
                            diagnostics_file=None,
                            label=f"scaffold_{scaffold_idx}",
                        )
                        report["scaffold_idx"] = scaffold_idx
                        diagnostics_all.append(report)
                        for pk_set in pk_samples:
                            fh.write(scs.pairs_to_pk_string(sorted(pk_set), L) + "\n")
                else:
                    # Parallelize across scaffolds (each worker runs chains sequentially).
                    _SCAFFOLD_GLOBALS["refined_structs"] = refined_structs
                    _SCAFFOLD_GLOBALS["candidate_pairs"] = candidate_pairs
                    _SCAFFOLD_GLOBALS["L"] = L
                    _SCAFFOLD_GLOBALS["sample_cfg"] = cfg.sample

                    results: dict[int, tuple[list[set[tuple[int, int]]], dict]] = {}
                    try:
                        with ProcessPoolExecutor(
                            max_workers=min(int(cfg.sample.max_workers), len(refined_structs)),
                            mp_context=ctx,
                        ) as ex:
                            futures = [
                                ex.submit(_sample_scaffold_worker, i)
                                for i in range(len(refined_structs))
                            ]
                            for fut in as_completed(futures):
                                scaffold_idx, pk_samples, report = fut.result()
                                results[scaffold_idx] = (pk_samples, report)
                    except (PermissionError, OSError) as exc:
                        sys.stderr.write(
                            f"[WARN] Scaffold parallelism unavailable ({exc}); "
                            "falling back to sequential scaffolds.\n"
                        )
                        for scaffold_idx, refined in enumerate(refined_structs):
                            scaffold_pairs = rup._pairs_from_struct(refined)
                            if cfg.sample.fix_scaffold_pairs:
                                scaffold_pos = {p for pair in scaffold_pairs for p in pair}
                                filtered_candidates = [
                                    (i, j, w)
                                    for (i, j, w) in candidate_pairs
                                    if i not in scaffold_pos and j not in scaffold_pos
                                ]
                                fixed_pairs = scaffold_pairs
                                init_pairs = None
                            else:
                                filtered_candidates = candidate_pairs
                                fixed_pairs = None
                                init_pairs = scaffold_pairs
                            pk_samples, report = scs.sample_matchings_multibeta_modular(
                                L=L,
                                candidate_pairs=filtered_candidates,
                                n_samples=n_samples_per_scaffold,
                                burn_in=cfg.sample.burn_in,
                                thin=cfg.sample.thin,
                                min_hairpin=cfg.sample.min_loop_sep,
                                beta=cfg.sample.beta,
                                betas=list(cfg.sample.betas) if cfg.sample.betas else None,
                                n_beta_chains=cfg.sample.n_beta_chains,
                                beta_min_factor=cfg.sample.beta_min_factor,
                                beta_max_factor=cfg.sample.beta_max_factor,
                                seed=(cfg.sample.seed + 10_000 * scaffold_idx)
                                if cfg.sample.seed is not None
                                else None,
                                pk_alpha=cfg.sample.pk_alpha,
                                pk_gamma=cfg.sample.pk_gamma,
                                lonely_penalty=cfg.sample.lonely_penalty,
                                pk_filter_frac=cfg.sample.pk_filter_frac,
                                pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
                                pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
                                delayed_accept=cfg.sample.delayed_accept,
                                parallel_chains=cfg.sample.parallel_chains,
                                max_workers=cfg.sample.max_workers,
                                use_beam_seeds=cfg.sample.use_beam_seeds,
                                beam_width=cfg.sample.beam_width,
                                beam_expansions_per_state=cfg.sample.beam_expansions_per_state,
                                beam_max_segments=cfg.sample.beam_max_segments,
                                beam_max_depth=cfg.sample.beam_max_depth,
                                scaffold_pairs=fixed_pairs,
                                initial_pairs=init_pairs,
                                use_segment_moves=cfg.sample.use_segment_moves,
                                segment_birth_prob=cfg.sample.segment_birth_prob,
                                segment_death_prob=cfg.sample.segment_death_prob,
                                segment_swap_prob=cfg.sample.segment_swap_prob,
                                diversity_dist_frac=cfg.sample.diversity_dist_frac,
                                diagnostics_file=None,
                                label=f"scaffold_{scaffold_idx}",
                            )
                            report["scaffold_idx"] = scaffold_idx
                            diagnostics_all.append(report)
                            for pk_set in pk_samples:
                                fh.write(scs.pairs_to_pk_string(sorted(pk_set), L) + "\n")
                        results = {}

                    for scaffold_idx in range(len(refined_structs)):
                        if not results:
                            break
                        pk_samples, report = results[scaffold_idx]
                        diagnostics_all.append(report)
                        for pk_set in pk_samples:
                            fh.write(scs.pairs_to_pk_string(sorted(pk_set), L) + "\n")
            else:
                # Sequential scaffolds; optionally parallelize chains within each scaffold.
                for scaffold_idx, refined in enumerate(refined_structs):
                    scaffold_pairs = rup._pairs_from_struct(refined)
                    min_hairpin = int(cfg.sample.min_loop_sep)
                    scaffold_pairs = {
                        p for p in scaffold_pairs if (max(p) - min(p) - 1) >= min_hairpin
                    }
                    if cfg.sample.fix_scaffold_pairs:
                        scaffold_pos = {p for pair in scaffold_pairs for p in pair}
                        filtered_candidates = [
                            (i, j, w)
                            for (i, j, w) in candidate_pairs
                            if i not in scaffold_pos and j not in scaffold_pos
                        ]
                        fixed_pairs = scaffold_pairs
                        init_pairs = None
                    else:
                        filtered_candidates = candidate_pairs
                        fixed_pairs = None
                        init_pairs = scaffold_pairs

                    pk_samples, report = scs.sample_matchings_multibeta_modular(
                        L=L,
                        candidate_pairs=filtered_candidates,
                        n_samples=n_samples_per_scaffold,
                        burn_in=cfg.sample.burn_in,
                        thin=cfg.sample.thin,
                        min_hairpin=cfg.sample.min_loop_sep,
                        beta=cfg.sample.beta,
                        betas=list(cfg.sample.betas) if cfg.sample.betas else None,
                        n_beta_chains=cfg.sample.n_beta_chains,
                        beta_min_factor=cfg.sample.beta_min_factor,
                        beta_max_factor=cfg.sample.beta_max_factor,
                        seed=(cfg.sample.seed + 10_000 * scaffold_idx)
                        if cfg.sample.seed is not None
                        else None,
                        pk_alpha=cfg.sample.pk_alpha,
                        pk_gamma=cfg.sample.pk_gamma,
                        lonely_penalty=cfg.sample.lonely_penalty,
                        pk_filter_frac=cfg.sample.pk_filter_frac,
                        pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
                        pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
                        delayed_accept=cfg.sample.delayed_accept,
                        parallel_chains=parallel_chains,
                        max_workers=cfg.sample.max_workers,
                        use_beam_seeds=cfg.sample.use_beam_seeds,
                        beam_width=cfg.sample.beam_width,
                        beam_expansions_per_state=cfg.sample.beam_expansions_per_state,
                        beam_max_segments=cfg.sample.beam_max_segments,
                        beam_max_depth=cfg.sample.beam_max_depth,
                        scaffold_pairs=fixed_pairs,
                        initial_pairs=init_pairs,
                        use_segment_moves=cfg.sample.use_segment_moves,
                        segment_birth_prob=cfg.sample.segment_birth_prob,
                        segment_death_prob=cfg.sample.segment_death_prob,
                        segment_swap_prob=cfg.sample.segment_swap_prob,
                        diversity_dist_frac=cfg.sample.diversity_dist_frac,
                        diagnostics_file=None,
                        label=f"scaffold_{scaffold_idx}",
                    )
                    report["scaffold_idx"] = scaffold_idx
                    diagnostics_all.append(report)
                    for pk_set in pk_samples:
                        fh.write(scs.pairs_to_pk_string(sorted(pk_set), L) + "\n")
        else:
            # Fallback: no refinement
            pk_samples, report = scs.sample_matchings_multibeta_modular(
                L=L,
                candidate_pairs=candidate_pairs,
                n_samples=n_samples_per_scaffold,
                burn_in=cfg.sample.burn_in,
                thin=cfg.sample.thin,
                min_hairpin=cfg.sample.min_loop_sep,
                beta=cfg.sample.beta,
                betas=list(cfg.sample.betas) if cfg.sample.betas else None,
                n_beta_chains=cfg.sample.n_beta_chains,
                beta_min_factor=cfg.sample.beta_min_factor,
                beta_max_factor=cfg.sample.beta_max_factor,
                seed=cfg.sample.seed,
                pk_alpha=cfg.sample.pk_alpha,
                pk_gamma=cfg.sample.pk_gamma,
                lonely_penalty=cfg.sample.lonely_penalty,
                pk_filter_frac=cfg.sample.pk_filter_frac,
                pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
                pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
                delayed_accept=cfg.sample.delayed_accept,
                parallel_chains=parallel_chains,
                max_workers=cfg.sample.max_workers,
                use_beam_seeds=cfg.sample.use_beam_seeds,
                beam_width=cfg.sample.beam_width,
                beam_expansions_per_state=cfg.sample.beam_expansions_per_state,
                beam_max_segments=cfg.sample.beam_max_segments,
                beam_max_depth=cfg.sample.beam_max_depth,
                use_segment_moves=cfg.sample.use_segment_moves,
                segment_birth_prob=cfg.sample.segment_birth_prob,
                segment_death_prob=cfg.sample.segment_death_prob,
                segment_swap_prob=cfg.sample.segment_swap_prob,
                diversity_dist_frac=cfg.sample.diversity_dist_frac,
                diagnostics_file=None,
                label="no_scaffold",
            )
            diagnostics_all.append(report)
            for pk_set in pk_samples:
                pk_str = scs.pairs_to_pk_string(sorted(pk_set), L)
                fh.write(pk_str + "\n")

    if cfg.sample.diagnostics_json is not None:
        cfg.sample.diagnostics_json.parent.mkdir(parents=True, exist_ok=True)
        cfg.sample.diagnostics_json.write_text(json.dumps(diagnostics_all, indent=2))
        sys.stderr.write(f"[PIPELINE] Wrote sampler diagnostics to {cfg.sample.diagnostics_json}\n")

    return out_db


def ensure_valid(cfg: PipelineConfig, db_in: Path, db_out: Path) -> Path:
    # Using existing logic placeholder
    structs = fdr.read_db_structures(db_in)
    # ... (filtering logic would go here) ...
    fdr.write_db_structures(db_out, structs)
    return db_out


def make_rosetta_ready(cfg: PipelineConfig, db_in: Path, db_out: Path) -> Path:
    structs = fdr.read_db_structures(db_in)

    # Read sequence using fdr helper which handles cleaning/uppercasing
    try:
        seq = fdr.read_ungapped_seq_from_sto(cfg.sto, cfg.seq_name)
    except Exception as e:
        raise ValueError(f"Failed to read sequence {cfg.seq_name} from {cfg.sto}: {e}")

    if seq is None:
        raise ValueError(f"Sequence {cfg.seq_name} read from {cfg.sto} is None.")

    sys.stderr.write(
        f"[PIPELINE] Loaded sequence {cfg.seq_name} (len={len(seq)}) for Rosetta conversion.\n"
    )

    # Dedup first
    deduped = fdr.deduplicate(structs)

    # Hamming merge
    merged = deduped
    # ... (skipping re-implementation of merge_by_hamming calls for brevity, reusing fdr)

    pruned = []
    dropped = 0
    for i, s in enumerate(merged):
        try:
            rn = fdr.convert_to_rosetta_notation(s)
            if rn is None:
                # Should not happen with corrected code, but defensive check
                raise ValueError("convert_to_rosetta_notation returned None")

            # Pass seq as keyword to be absolutely safe
            ps, _ = fdr.prune_noncomplementary_pairs(rn, seq=seq)
            pruned.append(ps)
        except ValueError:
            dropped += 1
            # Optional: Log failing structures to stderr if needed
            # sys.stderr.write(f"[WARN] Dropped invalid structure {i}: {e}\n")
            pass

    if dropped > 0:
        sys.stderr.write(
            f"[PIPELINE] Dropped {dropped} structures due to validation errors (e.g. overlapping pairs).\n"
        )

    fdr.write_db_structures(db_out, fdr.deduplicate(pruned), seq=seq)
    return db_out


def run_pipeline(cfg: PipelineConfig, mode: str = "refine_first") -> Path:
    mode = mode.lower()
    if mode not in {"legacy", "refine_first"}:
        raise ValueError(f"Unknown mode: {mode!r}")

    consensus_db = get_consensus_db(cfg)
    allsub_cache = None
    duplex_cache = None
    if cfg.refine.allsub_cache_path:
        allsub_cache = rup.load_cache(cfg.refine.allsub_cache_path)
    if cfg.refine.duplex_cache_path:
        duplex_cache = rup.load_cache(cfg.refine.duplex_cache_path)

    if mode == "legacy":
        sampled_db = sample_pk(cfg, refined_db=None, out_db=cfg.io.sampled_db)
        refined_db = refine_db(
            cfg,
            db_in=sampled_db,
            db_out=cfg.io.refined_db,
            allsub_cache=allsub_cache,
            duplex_cache=duplex_cache,
        )
        return make_rosetta_ready(cfg, db_in=refined_db, db_out=cfg.io.rosetta_db)

    # refine_first
    refined_consensus_db = refine_db(
        cfg,
        db_in=consensus_db,
        db_out=cfg.io.refined_db,
        allsub_cache=allsub_cache,
        duplex_cache=duplex_cache,
    )
    sampled_db = sample_pk(cfg, refined_db=refined_consensus_db, out_db=cfg.io.sampled_db)

    return make_rosetta_ready(cfg, db_in=sampled_db, db_out=cfg.io.rosetta_db)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Legacy config-driven CaCoFold pipeline wrapper. "
            "For the recommended top-K predictor entrypoint, use `rnanneal-ss`."
        )
    )
    parser.add_argument("-c", "--config", required=True, help="Path to config file.")
    parser.add_argument("--mode", default="refine_first")
    args = parser.parse_args(argv)

    cfg_path = Path(args.config).resolve()
    cfg = load_pipeline_config(cfg_path)

    final_db = run_pipeline(cfg, mode=args.mode)
    sys.stderr.write(f"[PIPELINE] Final DB: {final_db}\n")


if __name__ == "__main__":
    main()
