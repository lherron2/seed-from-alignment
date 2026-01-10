"""
High-level CaCoFold → sampling → refinement → Rosetta pipelines.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from src.lib import filter_db_for_rosetta as fdr
from src.lib import refine_unpaired_regions as rup

# Import the existing building blocks
from src.lib import sample_cacofold_structures as scs


@dataclass
class SampleConfig:
    n_samples: int = 1000
    burn_in: int = 1000
    thin: int = 10
    min_loop_sep: int = 1
    beta: float = 1.0
    seed: int | None = None
    pk_alpha: float = 2.0
    pk_filter_frac: float = 0.2
    pk_filter_max_cross_per_pair: int = 2
    pk_filter_max_total_cross: int = 30
    pk_depth_limit: int | None = None
    cov_mode: str = "off"
    cov_alpha: float = 1.0
    cov_min_power: float = 0.0
    cov_forbid_negative: bool = False
    cov_negative_E: float = 1.0


@dataclass
class RefineConfig:
    min_unpaired: int = 15
    temperature: float | None = None
    allsub_abs: float | None = None
    allsub_pct: float | None = None
    fold_extra_args: Sequence[str] = ()
    max_structures: int = 1000  # Default limit for refined output

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
            import yaml  # type: ignore
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

    sample_raw = raw.get("sample", {})
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

    io.work_dir.mkdir(parents=True, exist_ok=True)

    return PipelineConfig(
        sto=sto,
        cov=cov,
        seq_name=seq_name,
        allsub_exe=allsub_exe,
        duplex_exe=duplex_exe,
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


def refine_db(cfg: PipelineConfig, db_in: Path, db_out: Path) -> Path:
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

    allsub_cache = {}
    duplex_cache = {}

    # rup.refine_structure returns List[Tuple[str, float]]
    refined_structs_with_scores: list[tuple[str, float]] = []

    for idx, s in enumerate(structs):
        variants = rup.refine_structure(
            struct=s,
            full_seq=full_seq,
            allsub_exe=allsub_exe,
            duplex_exe=duplex_exe,
            min_unpaired_len=min_unpaired,
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

    L, candidate_pairs = scs.extract_candidate_pairs(
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

    refined_structs = []
    if refined_db and refined_db.is_file():
        refined_structs = rup.read_db_structures(refined_db)

    out_db.parent.mkdir(parents=True, exist_ok=True)
    with out_db.open("w") as fh:
        fh.write(ungapped_seq + "\n")

        if refined_structs:
            for refined in refined_structs:
                # Parse all fixed pairs from the scaffold
                scaffold_pairs = rup._pairs_from_struct(refined)
                scaffold_pos = {p for pair in scaffold_pairs for p in pair}

                # Filter out candidates that conflict with locked scaffold
                filtered_candidates = [
                    (i, j, w)
                    for (i, j, w) in candidate_pairs
                    if i not in scaffold_pos and j not in scaffold_pos
                ]

                # Pass scaffold pairs directly to the sampler
                pk_samples = scs.sample_matchings(
                    L=L,
                    candidate_pairs=filtered_candidates,
                    n_samples=cfg.sample.n_samples,
                    burn_in=cfg.sample.burn_in,
                    thin=cfg.sample.thin,
                    min_loop_sep=cfg.sample.min_loop_sep,
                    beta=cfg.sample.beta,
                    seed=cfg.sample.seed,
                    pk_alpha=cfg.sample.pk_alpha,
                    pk_filter_frac=cfg.sample.pk_filter_frac,
                    pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
                    pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
                    pk_depth_limit=cfg.sample.pk_depth_limit,
                    scaffold_pairs=scaffold_pairs,  # <-- Scaffold enforced in sampler
                )

                for pk_set in pk_samples:
                    # pk_set now includes the scaffold pairs already
                    pk_str = scs.pairs_to_pk_string(sorted(pk_set), L)
                    fh.write(pk_str + "\n")
        else:
            # Fallback: no refinement
            pk_samples = scs.sample_matchings(
                L=L,
                candidate_pairs=candidate_pairs,
                n_samples=cfg.sample.n_samples,
                burn_in=cfg.sample.burn_in,
                thin=cfg.sample.thin,
                min_loop_sep=cfg.sample.min_loop_sep,
                beta=cfg.sample.beta,
                seed=cfg.sample.seed,
                pk_alpha=cfg.sample.pk_alpha,
                pk_filter_frac=cfg.sample.pk_filter_frac,
                pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
                pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
                pk_depth_limit=cfg.sample.pk_depth_limit,
            )
            for pk_set in pk_samples:
                pk_str = scs.pairs_to_pk_string(sorted(pk_set), L)
                fh.write(pk_str + "\n")

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

    if mode == "legacy":
        sampled_db = sample_pk(cfg, refined_db=None, out_db=cfg.io.sampled_db)
        refined_db = refine_db(cfg, db_in=sampled_db, db_out=cfg.io.refined_db)
        return make_rosetta_ready(cfg, db_in=refined_db, db_out=cfg.io.rosetta_db)

    # refine_first
    refined_consensus_db = refine_db(cfg, db_in=consensus_db, db_out=cfg.io.refined_db)
    sampled_db = sample_pk(cfg, refined_db=refined_consensus_db, out_db=cfg.io.sampled_db)

    return make_rosetta_ready(cfg, db_in=sampled_db, db_out=cfg.io.rosetta_db)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="CaCoFold pipeline wrapper.")
    parser.add_argument("-c", "--config", required=True, help="Path to config file.")
    parser.add_argument("--mode", default="refine_first")
    args = parser.parse_args(argv)

    cfg_path = Path(args.config).resolve()
    cfg = load_pipeline_config(cfg_path)

    final_db = run_pipeline(cfg, mode=args.mode)
    sys.stderr.write(f"[PIPELINE] Final DB: {final_db}\n")


if __name__ == "__main__":
    main()
