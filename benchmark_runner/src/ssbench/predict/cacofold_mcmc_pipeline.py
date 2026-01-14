"""Infernal + CaCoFold + domain masking + MCMC sampling pipeline.

Requires local Rfam CM database, Infernal, R-scape, and RNAstructure AllSub.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path


def run(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, check=True, cwd=cwd)


def parse_fasta_id(fasta_path: Path) -> str:
    line = fasta_path.read_text().splitlines()[0]
    return line[1:].strip().split()[0] if line.startswith(">") else line.strip()


def read_fasta_sequence(fasta_path: Path) -> str:
    lines = [ln.strip() for ln in fasta_path.read_text().splitlines()]
    if not lines:
        return ""
    if lines[0].startswith(">"):
        return "".join(lines[1:]).strip()
    return "".join(lines).strip()


def sanitize_rna_sequence(seq: str) -> tuple[str, bool]:
    """Return (sanitized_seq, changed) with T->U and illegal characters replaced.

    For the under400 benchmark, the only observed non-ACGUT symbols are:
    - 'X': unknown/modified base; tools generally accept it, so we preserve it
    - '&': Barnaba fragment separator; external tools reject it, so map to 'N'
    """
    seq = seq.strip()
    out: list[str] = []
    changed = False
    for ch in seq.upper():
        if ch == "T":
            out.append("U")
            changed = True
            continue
        if ch in {"A", "C", "G", "U", "X", "N"}:
            out.append(ch)
            continue
        # Keep sequence length stable while ensuring a tool-compatible alphabet.
        out.append("N")
        changed = True
    return "".join(out), changed


def write_sanitized_fasta(*, fasta_in: Path, out_dir: Path) -> tuple[Path, bool]:
    seq_id = parse_fasta_id(fasta_in)
    raw_seq = read_fasta_sequence(fasta_in)
    sanitized, changed = sanitize_rna_sequence(raw_seq)
    out_path = out_dir / "input_sanitized.fa"
    out_path.write_text(f">{seq_id}\n{sanitized}\n")
    return out_path, changed


def ensure_rnastructure_datapath(*, exe_hint: Path | None = None) -> None:
    """Best-effort ensure DATAPATH is set so RNAstructure binaries can find thermodynamic tables.

    The benchmark runners usually export DATAPATH already; this is a defensive fallback for
    ad-hoc usage where it's missing.
    """
    # Treat an explicitly exported (possibly empty) DATAPATH as "set".
    # This supports docker-wrapped RNAstructure executables, which provide their own in-container DATAPATH.
    if "DATAPATH" in os.environ:
        return
    candidates: list[Path] = []
    if exe_hint is not None:
        # Common install layout: <root>/exe/{AllSub,Fold,...} and <root>/data_tables.
        candidates.append(exe_hint.parent / "data_tables")
        candidates.append(exe_hint.parent.parent / "data_tables")
        candidates.append(exe_hint.parent.parent.parent / "data_tables")

    repo_root = Path(__file__).resolve().parents[4]
    candidates.append(repo_root / "RNAstructure" / "data_tables")

    for path in candidates:
        if path.is_dir():
            os.environ["DATAPATH"] = str(path)
            return


def parse_tblout_top_hit(tblout_path: Path) -> str | None:
    for line in tblout_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        return parts[1]
    return None


def find_cacofold_sto(out_dir: Path) -> Path | None:
    for path in out_dir.glob("*.cacofold.sto"):
        return path
    for path in out_dir.glob("*.sto"):
        if "cacofold" in path.name:
            return path
    return None


def find_cacofold_cov(out_dir: Path) -> Path | None:
    for path in out_dir.glob("*.cov"):
        return path
    return None


def run_cacofold(
    fasta: Path,
    out_dir: Path,
    cm_db: Path,
    rscape: str,
    infernal_bin: Path | None,
) -> tuple[Path | None, Path | None, str, str]:
    cmscan = "cmscan"
    cmfetch = "cmfetch"
    cmalign = "cmalign"
    if infernal_bin is not None:
        cmscan = str(infernal_bin / "cmscan")
        cmfetch = str(infernal_bin / "cmfetch")
        cmalign = str(infernal_bin / "cmalign")

    seq_id = parse_fasta_id(fasta)
    tblout = out_dir / "cmscan.tblout"
    run([cmscan, "--tblout", str(tblout), str(cm_db), str(fasta)])
    model = parse_tblout_top_hit(tblout)
    if model is None:
        raise SystemExit("No Rfam hit found for input FASTA")

    cm_path = out_dir / f"{model}.cm"
    with cm_path.open("w") as fh:
        subprocess.run([cmfetch, str(cm_db), model], check=True, stdout=fh)

    aln_sto = out_dir / "aln.sto"
    with aln_sto.open("w") as fh:
        subprocess.run([cmalign, str(cm_path), str(fasta)], check=True, stdout=fh)

    try:
        # R-scape can be extremely verbose (ASCII diagrams, etc). Keep the global benchmark logs
        # readable by capturing its output to a per-target file.
        rscape_log = out_dir / "rscape.log"
        with rscape_log.open("w") as fh:
            subprocess.run(
                [rscape, "--cacofold", "--nofigures", "--outdir", str(out_dir), str(aln_sto)],
                check=True,
                stdout=fh,
                stderr=fh,
            )
    except subprocess.CalledProcessError:
        # Fall back to MFE-only pipeline when CaCoFold fails (e.g., no covariation).
        return None, None, seq_id, model

    cacofold_sto = find_cacofold_sto(out_dir)
    cov = find_cacofold_cov(out_dir)
    return cacofold_sto, cov, seq_id, model


def build_fallback_sto(
    fasta: Path,
    out_dir: Path,
    fold_exe: Path,
) -> tuple[Path, str]:
    # Build a minimal Stockholm file with MFE SS_cons when no covariation hit exists.
    from src.lib import refine_unpaired_regions as rup

    seq_id = parse_fasta_id(fasta)
    seq = read_fasta_sequence(fasta)
    ct_path = out_dir / "mfe.ct"
    run([str(fold_exe), str(fasta), str(ct_path)])
    mfe_structs = rup.parse_ct_file(ct_path)
    if not mfe_structs:
        raise SystemExit("Fold produced no structures")
    mfe_struct = mfe_structs[0][0]
    sto_path = out_dir / "fallback.sto"
    sto_path.write_text(
        "# STOCKHOLM 1.0\n"
        f"{seq_id} {seq}\n"
        f"#=GC SS_cons {mfe_struct}\n"
        "//\n"
    )
    return sto_path, seq_id


def build_topk_predictions(
    sto: Path,
    cov: Path | None,
    seq_name: str,
    rfam_id: str,
    fasta_had_nonstandard_bases: bool,
    allsub_exe: Path,
    duplex_exe: Path | None,
    fold_exe: Path | None,
    partition_exe: Path | None,
    probplot_exe: Path | None,
    fasta: Path,
    out_dir: Path,
    top_k: int,
    n_samples: int,
    burn_in: int,
    thin: int,
    beta: float,
    seed: int,
    min_loop_sep: int,
    pk_alpha: float,
    pair_penalty: float | None,
    pair_penalty_scale: float,
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
    pair_penalty_mode: str,
    pair_penalty_c0: float,
    pair_penalty_c1: float,
    pair_penalty_min: float | None,
    pair_penalty_max: float | None,
    refine_max_structures: int,
    refine_min_unpaired: int,
    refine_end_mask_step: int,
    refine_max_end_mask_len: int,
    refine_max_helices_sequential: int,
    refine_max_helices_pairwise: int,
    refine_max_regions: int,
    refine_max_seeds: int,
    refine_max_solutions: int,
    refine_kissing_candidates: int,
    refine_max_seconds: float | None,
    max_scaffolds: int,
    max_samples_per_scaffold: int,
    length_adaptive: bool,
    include_unfixed_sampling: bool,
    inject_allsub_scaffolds: bool,
    inject_allsub_scaffolds_max: int,
    inject_allsub_timeout_s: float | None,
    force_allsub_output: int,
) -> list[str]:
    # Import legacy modules from repo root.
    repo_root = Path(__file__).resolve().parents[4]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from src.lib import energy as legacy_energy
    from src.lib import pipeline as legacy_pipeline
    from src.lib import refine_unpaired_regions as rup
    from src.lib import sample_cacofold_structures as scs

    t0 = time.perf_counter()
    seqs, ss_tracks = scs.parse_stockholm_single(str(sto))
    if seq_name not in seqs:
        if len(seqs) == 1:
            seq_name = next(iter(seqs))
        else:
            raise KeyError(f"Sequence {seq_name!r} not found in {sto}")

    work_dir = out_dir
    io = legacy_pipeline.IOConfig(
        work_dir=work_dir,
        consensus_db=work_dir / "00_consensus.db",
        refined_db=work_dir / "01_refined.db",
        sampled_db=work_dir / "02_sampled.db",
        valid_db=work_dir / "03_valid.db",
        rosetta_db=work_dir / "04_rosetta.db",
        summary=None,
    )

    fasta_seq = read_fasta_sequence(fasta)
    L_fasta = len(fasta_seq)

    def clamp01(x: float) -> float:
        if x <= 0.0:
            return 0.0
        if x >= 1.0:
            return 1.0
        return x

    # Length-adaptive schedule: allocate more compute and reduce reliance on early refined ordering
    # for longer RNAs (151-300), where we empirically see weaker best-of-100 and late ranks.
    len_scale = clamp01((float(L_fasta) - 150.0) / 150.0) if length_adaptive else 0.0

    max_scaffolds_eff = int(max(1, round(float(max_scaffolds) * (1.0 + len_scale))))
    max_samples_per_scaffold_eff = int(
        max(1, round(float(max_samples_per_scaffold) * (1.0 + 2.0 * len_scale)))
    )

    burn_in_eff = int(max(0, round(float(burn_in) * (1.0 + 0.5 * len_scale))))
    n_samples_eff = int(max(0, round(float(n_samples) * (1.0 + 0.5 * len_scale))))
    thin_eff = int(thin)
    if thin_eff <= 0:
        thin_eff = 1

    refine_max_seconds_eff = refine_max_seconds
    if refine_max_seconds_eff is not None:
        refine_max_seconds_eff = float(refine_max_seconds_eff) * (1.0 + 9.0 * len_scale)
        refine_max_seconds_eff = min(float(refine_max_seconds_eff), 300.0)

    refine_max_structures_eff = int(
        max(1, round(float(refine_max_structures) * (1.0 + 1.0 * len_scale)))
    )
    refine_max_regions_eff = int(max(1, round(float(refine_max_regions) * (1.0 + 1.0 * len_scale))))
    refine_max_seeds_eff = int(max(1, round(float(refine_max_seeds) * (1.0 + 1.0 * len_scale))))
    refine_max_solutions_eff = int(
        max(1, round(float(refine_max_solutions) * (1.0 + 1.0 * len_scale)))
    )
    refine_kissing_candidates_eff = int(
        max(1, round(float(refine_kissing_candidates) * (1.0 + 1.0 * len_scale)))
    )

    sample_cfg = legacy_pipeline.SampleConfig(
        n_samples=n_samples_eff,
        burn_in=burn_in_eff,
        thin=thin_eff,
        min_loop_sep=min_loop_sep,
        beta=beta,
        seed=seed,
        pk_alpha=float(pk_alpha),
        cov_mode=cov_mode,
        cov_alpha=cov_alpha,
        cov_min_power=cov_min_power,
        cov_forbid_negative=cov_forbid_negative,
        # Default behavior: keep refined scaffold pairs fixed; sample only in
        # the remaining degrees of freedom (important when covariation is absent).
        fix_scaffold_pairs=True,
        # Keep sampling practical: cap scaffolds and per-scaffold samples (scaled for long RNAs).
        max_scaffolds=max_scaffolds_eff,
        max_samples_per_scaffold=max_samples_per_scaffold_eff,
    )

    refine_cfg = legacy_pipeline.RefineConfig()
    refine_cfg.max_structures = refine_max_structures_eff
    refine_cfg.min_unpaired = refine_min_unpaired
    refine_cfg.end_mask_step = refine_end_mask_step
    refine_cfg.max_end_mask_len = refine_max_end_mask_len
    refine_cfg.max_helices_sequential = refine_max_helices_sequential
    refine_cfg.max_helices_pairwise = refine_max_helices_pairwise
    refine_cfg.max_regions_to_refine = refine_max_regions_eff
    refine_cfg.max_seeds = refine_max_seeds_eff
    refine_cfg.max_solutions = refine_max_solutions_eff
    refine_cfg.kissing_loop_candidates = refine_kissing_candidates_eff
    refine_cfg.max_seconds = refine_max_seconds_eff
    ensure_cfg = legacy_pipeline.EnsureValidConfig()
    rosetta_cfg = legacy_pipeline.RosettaConfig()

    cfg = legacy_pipeline.PipelineConfig(
        sto=sto,
        cov=cov,
        seq_name=seq_name,
        allsub_exe=allsub_exe,
        duplex_exe=duplex_exe,
        partition_exe=partition_exe,
        probplot_exe=probplot_exe,
        sample=sample_cfg,
        refine=refine_cfg,
        ensure_valid=ensure_cfg,
        rosetta=rosetta_cfg,
        io=io,
    )

    # Build weights for scoring and detect covariation availability.
    aligned_seq = seqs[seq_name]
    aln2seq, _L = scs.aln_to_seq_map(aligned_seq)
    cov_stats = None
    if cov is not None and cov.exists():
        cov_stats = scs.load_cov_stats(str(cov), aln2seq)

    # Always start from the CaCoFold/SS_cons scaffold.
    consensus_db = legacy_pipeline.get_consensus_db(cfg)
    consensus_structs = rup.read_db_structures(consensus_db)
    consensus_pairs_n = len(rup._pairs_from_struct(consensus_structs[0])) if consensus_structs else 0
    cov_has_pairs = bool(cov_stats)
    cov_empty = bool(cov is not None and cov.exists() and not cov_has_pairs)
    seed_boost = bool(
        cov is None
        or rfam_id == "no_rfam_hit"
        or fasta_had_nonstandard_bases
        or consensus_pairs_n <= 6
    )

    # Also compute an MFE structure (when possible) but do not send it through the expensive
    # refinement stage; instead, inject it as an extra scaffold for sampling/ranking.
    mfe_struct: str | None = None
    if fold_exe is not None and fold_exe.exists():
        ct_path = work_dir / "mfe.ct"
        run([str(fold_exe), str(fasta), str(ct_path)])
        mfe_structs = rup.parse_ct_file(ct_path)
        if mfe_structs:
            mfe_struct = mfe_structs[0][0]

    t_refine0 = time.perf_counter()
    refined_db = legacy_pipeline.refine_db(cfg, db_in=consensus_db, db_out=cfg.io.refined_db)
    t_refine1 = time.perf_counter()
    # Defensive: refinement can occasionally yield no structures (e.g. fallback.sto + strict timeouts).
    # Downstream stages (injection/sampling/ranking) require at least one scaffold; fall back to the
    # consensus structure (or all-dots) to keep the run valid rather than crashing.
    try:
        refined_existing = rup.read_db_structures(refined_db)
    except Exception:
        refined_existing = []
    if not refined_existing:
        try:
            seq = rup.read_ungapped_seq_from_sto(sto, seq_name)
        except Exception:
            seq = read_fasta_sequence(fasta)
        fallback_struct = consensus_structs[0] if consensus_structs else ("." * len(seq))
        refined_db.write_text(seq + "\n" + fallback_struct + "\n")
        refined_existing = [fallback_struct]

    if mfe_struct is not None:
        try:
            seq = rup.read_ungapped_seq_from_sto(sto, seq_name)
        except Exception:
            seq = read_fasta_sequence(fasta)
        if mfe_struct not in set(refined_existing):
            refined_existing = list(refined_existing) + [mfe_struct]
            refined_db.write_text(seq + "\n" + "\n".join(refined_existing) + "\n")

    # Inject additional backend-generated scaffolds *before* sampling so they can seed MCMC.
    # This is the key change vs. the post-hoc "ensemble merge": these candidates become actual
    # MCMC starting scaffolds in `sample_pk`.
    if inject_allsub_scaffolds and inject_allsub_scaffolds_max > 0:
        try:
            full_seq_for_inject = rup.read_ungapped_seq_from_sto(sto, seq_name)
        except Exception:
            full_seq_for_inject = read_fasta_sequence(fasta)

        def pairs_from_struct(struct: str) -> frozenset[tuple[int, int]]:
            stack: list[int] = []
            out: list[tuple[int, int]] = []
            for i, ch in enumerate(struct):
                if ch == "(":
                    stack.append(i)
                elif ch == ")":
                    if not stack:
                        continue
                    j = stack.pop()
                    out.append((j, i) if j < i else (i, j))
            return frozenset(out)

        def jaccard_dist(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]) -> float:
            if not a and not b:
                return 0.0
            inter = len(a & b)
            union = len(a | b)
            return 1.0 - (inter / union if union else 0.0)

        def select_diverse(structs: list[str], k: int) -> list[str]:
            if k <= 0:
                return []
            if len(structs) <= k:
                return list(structs)
            items = [(s, pairs_from_struct(s), idx) for idx, s in enumerate(structs)]
            chosen = [items[0]]
            chosen_pairs = [items[0][1]]
            chosen_set = {items[0][0]}
            while len(chosen) < k:
                best = None
                best_u = None
                for s, ps, order in items:
                    if s in chosen_set:
                        continue
                    base = 1.0 - (order / max(1, len(items) - 1))
                    min_d = min(jaccard_dist(ps, cp) for cp in chosen_pairs) if chosen_pairs else 1.0
                    u = base + 0.8 * min_d
                    if best_u is None or u > best_u:
                        best = (s, ps, order)
                        best_u = u
                if best is None:
                    break
                chosen.append(best)
                chosen_set.add(best[0])
                chosen_pairs.append(best[1])
            return [s for s, _ps, _ord in chosen]

        try:
            allsub_structs = rup.call_rnastructure_allsub(
                allsub_exe,
                full_seq_for_inject,
                absolute_energy=None,
                percent_energy=None,
                extra_args=None,
                timeout_s=inject_allsub_timeout_s,
            )
        except Exception:
            allsub_structs = []

        pool: list[str] = []
        seen_pool: set[str] = set()
        for db, _e in allsub_structs:
            if len(db) != len(full_seq_for_inject):
                continue
            if db in seen_pool:
                continue
            seen_pool.add(db)
            pool.append(db)
            if len(pool) >= max(200, int(inject_allsub_scaffolds_max) * 5):
                break

        if pool:
            selected = select_diverse(pool, int(inject_allsub_scaffolds_max))
            existing = rup.read_db_structures(refined_db)
            existing_set = set(existing)
            merged = list(existing)
            for struct_db in selected:
                if struct_db in existing_set:
                    continue
                merged.append(struct_db)
                existing_set.add(struct_db)
            refined_db.write_text(full_seq_for_inject + "\n" + "\n".join(merged) + "\n")
    t_sample_fixed0 = time.perf_counter()
    sampled_db = legacy_pipeline.sample_pk(cfg, refined_db=refined_db, out_db=cfg.io.sampled_db)
    t_sample_fixed1 = time.perf_counter()

    # Optional extra sampling pass that allows modifying scaffold pairs.
    # This targets the long-RNA regime where the (single-sequence) scaffold can be a poor global
    # fit to the target's protein-stabilized 3D conformation.
    unfixed_db: Path | None = None
    if include_unfixed_sampling and L_fasta >= 80:
        thermo_pool_n = int(thermo_max_structures)
        thermo_min_count_eff = max(1, int(thermo_min_count))
        if cov_empty:
            thermo_pool_n = max(thermo_pool_n, 200)
            thermo_min_count_eff = 1

        pk_filter_frac = 0.35
        pk_filter_max_cross_per_pair = min(30, max(6, int(round(0.20 * float(L_fasta))))) if L_fasta else 6
        pk_filter_max_total_cross = min(500, max(30, int(round(1.50 * float(L_fasta))))) if L_fasta else 30

        unfixed_sample_cfg = dataclasses.replace(
            sample_cfg,
            fix_scaffold_pairs=False,
            seed=int(seed) + 1,
            max_scaffolds=max(5, int(round(max_scaffolds_eff * 0.25))),
            max_samples_per_scaffold=max(100, int(round(max_samples_per_scaffold_eff * 0.33))),
            pk_alpha=float(sample_cfg.pk_alpha) * 0.5,
            pk_filter_frac=float(pk_filter_frac),
            pk_filter_max_cross_per_pair=int(pk_filter_max_cross_per_pair),
            pk_filter_max_total_cross=int(pk_filter_max_total_cross),
            weight_calibration_method=str(weight_calibration_method),
            weight_calibration_zmax=float(weight_calibration_zmax),
            weight_alpha_core=float(weight_alpha_core),
            weight_alpha_alt=float(weight_alpha_alt),
            weight_alpha_cov=float(weight_alpha_cov),
            weight_alpha_thermo=float(weight_alpha_thermo),
            thermo_mode=str(thermo_mode),
            thermo_weight=float(thermo_weight),
            thermo_max_structures=int(thermo_pool_n),
            thermo_min_count=int(thermo_min_count_eff),
            thermo_min_prob=float(thermo_min_prob),
            thermo_log_eps=float(thermo_log_eps),
        )
        unfixed_io = dataclasses.replace(cfg.io, sampled_db=work_dir / "02_sampled_unfixed.db")
        unfixed_cfg = dataclasses.replace(cfg, sample=unfixed_sample_cfg, io=unfixed_io)
        t_sample_unfixed0 = time.perf_counter()
        try:
            unfixed_db = legacy_pipeline.sample_pk(
                unfixed_cfg, refined_db=refined_db, out_db=unfixed_cfg.io.sampled_db
            )
        except Exception:
            unfixed_db = None
        t_sample_unfixed1 = time.perf_counter()
    else:
        t_sample_unfixed0 = None
        t_sample_unfixed1 = None

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

    from src.lib import rnastructure_pf as rpf
    from src.lib import stem_stats as ss
    from src.lib import weight_calibration as wc

    thermo_component: dict[tuple[int, int], float] = {}
    try:
        full_seq = rup.read_ungapped_seq_from_sto(sto, seq_name)
    except Exception:
        full_seq = read_fasta_sequence(fasta)

    thermo_structs_full: list[tuple[str, float]] = []

    if thermo_weight > 0 and thermo_mode != "off":
        if thermo_mode == "pf" and partition_exe and probplot_exe:
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

        thermo_pool_n = int(thermo_max_structures)
        if seed_boost:
            thermo_pool_n = max(thermo_pool_n, 200)

        if not thermo_component and thermo_mode in {"pf", "allsub"}:
            try:
                thermo_structs_full = rup.call_rnastructure_allsub(allsub_exe, full_seq)
            except Exception:
                thermo_structs_full = []
            thermo_structs = thermo_structs_full[: max(0, thermo_pool_n)]
            if thermo_structs:
                counts: dict[tuple[int, int], int] = {}
                for db, _energy in thermo_structs:
                    for i, j in rup._pairs_from_struct(db):
                        if i > j:
                            i, j = j, i
                        if i < 0 or j >= len(full_seq):
                            continue
                        if not scs.is_canonical_pair(full_seq[i], full_seq[j]):
                            continue
                        counts[(i, j)] = counts.get((i, j), 0) + 1
                denom = float(len(thermo_structs))
                # In weak-evidence regimes we want coverage more than consensus frequency; keep rare pairs too.
                min_count = 1 if seed_boost else max(1, int(thermo_min_count))
                for (i, j), c in counts.items():
                    if c < min_count:
                        continue
                    thermo_component[(i, j)] = float(c) / denom

    components["thermo"] = thermo_component

    if weight_calibration_method == "none":
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
    weights_base = wc.blend_components(components, alphas=alphas, config=norm_cfg)

    w_scale = wc.weight_scale(weights_base)
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
        pk_alpha=sample_cfg.pk_alpha,
        pk_gamma=sample_cfg.pk_gamma,
        lonely_penalty=sample_cfg.lonely_penalty,
    )

    stem_start_penalty = float(stem_start_penalty_scale) * w_scale
    stem_len1_penalty = float(stem_len1_penalty_scale) * w_scale
    stem_len2_penalty = float(stem_len2_penalty_scale) * w_scale
    stem_log_reward = float(stem_log_reward_scale) * w_scale
    support_threshold = wc.quantile_threshold(weights_base, float(stem_support_quantile))

    weights = weights_base

    def normalized_pairs(struct: str) -> frozenset[tuple[int, int]]:
        pairs = []
        for i, j in rup._pairs_from_struct(struct):
            if i > j:
                i, j = j, i
            pairs.append((i, j))
        return frozenset(pairs)

    def compute_pk_stats(pairs: frozenset[tuple[int, int]]) -> tuple[int, int, int]:
        """Return (pk_pairs_count, total_crossings, max_crossings_per_pair)."""
        pairs_list = sorted(pairs)
        cross_counts = {p: 0 for p in pairs_list}
        total_crossings = 0
        for a in range(len(pairs_list)):
            i, j = pairs_list[a]
            for b in range(a + 1, len(pairs_list)):
                k, l = pairs_list[b]
                if (i < k < j < l) or (k < i < l < j):
                    cross_counts[(i, j)] += 1
                    cross_counts[(k, l)] += 1
                    total_crossings += 1
        pk_pairs_count = sum(1 for c in cross_counts.values() if c > 0)
        max_cross = max(cross_counts.values()) if cross_counts else 0
        return pk_pairs_count, total_crossings, max_cross

    def intersection_size(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]) -> int:
        if len(a) > len(b):
            a, b = b, a
        return sum(1 for p in a if p in b)

    def jaccard_distance(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]) -> float:
        if not a and not b:
            return 0.0
        inter = intersection_size(a, b)
        union = len(a) + len(b) - inter
        return 1.0 - (inter / union if union > 0 else 0.0)

    def score_pairs(
        pairs: frozenset[tuple[int, int]],
        pk_alpha_override: float | None,
        pair_penalty_multiplier: float,
    ) -> float:
        effective_params = legacy_energy.EnergyParams(
            pk_alpha=float(params.pk_alpha if pk_alpha_override is None else pk_alpha_override),
            pk_gamma=float(params.pk_gamma),
            lonely_penalty=float(params.lonely_penalty),
        )
        pair_set = set(pairs)
        energy_val, _cache = legacy_energy.compute_energy(pair_set, weights, effective_params)
        if pair_penalty:
            energy_val += float(pair_penalty) * float(pair_penalty_multiplier) * len(pairs)
        if any(
            val != 0.0
            for val in (stem_start_penalty, stem_len1_penalty, stem_len2_penalty, stem_log_reward)
        ):
            stats = ss.stem_stats(set(pairs))
            unsupported_len1 = sum(1 for p in stats.len1_pairs if weights.get(p, 0.0) < support_threshold)
            energy_val += stem_start_penalty * stats.n_stems
            energy_val += stem_len1_penalty * unsupported_len1
            energy_val += stem_len2_penalty * stats.n_len2
            energy_val -= stem_log_reward * stats.sum_log1p_len
        return float(energy_val)

    refined_structs = rup.read_db_structures(refined_db)
    sampled_structs_fixed = rup.read_db_structures(sampled_db)
    n_sampled_fixed = int(len(sampled_structs_fixed))
    n_sampled_unfixed = 0
    sampled_structs = list(sampled_structs_fixed)
    if unfixed_db is not None and unfixed_db.exists():
        sampled_structs_unfixed = rup.read_db_structures(unfixed_db)
        n_sampled_unfixed = int(len(sampled_structs_unfixed))
        sampled_structs += sampled_structs_unfixed

    # --- Candidate pool construction (cheap proxies first) ---
    # Always consider all refined structures (already capped by refine_max_structures).
    refined_candidates: dict[str, dict[str, object]] = {}
    for refined_order, struct in enumerate(refined_structs):
        pairs = normalized_pairs(struct)
        pk_pairs_count, total_crossings, max_cross = compute_pk_stats(pairs)
        support_sum = float(sum(weights.get(p, 0.0) for p in pairs))
        refined_candidates[struct] = {
            "struct": struct,
            "pairs": pairs,
            "is_refined": True,
            "refined_order": refined_order,
            "pk_pairs_count": pk_pairs_count,
            "total_crossings": total_crossings,
            "max_cross": max_cross,
            "support_sum": support_sum,
        }

    # For sampled candidates, keep a pooled subset that covers:
    # - high-support candidates (likely good under nested evidence)
    # - PK-heavy candidates (to avoid dropping PK regimes due to a single scorer)
    sampled_records: list[tuple[str, frozenset[tuple[int, int]], int, int, int, float]] = []
    for struct in sampled_structs:
        pairs = normalized_pairs(struct)
        pk_pairs_count, total_crossings, max_cross = compute_pk_stats(pairs)
        support_sum = float(sum(weights.get(p, 0.0) for p in pairs))
        sampled_records.append((struct, pairs, pk_pairs_count, total_crossings, max_cross, support_sum))

    sampled_records.sort(key=lambda x: x[5], reverse=True)
    base_pool = max(2000, int(top_k) * 20)
    max_sampled_pool = (
        int(round(float(base_pool) * (1.0 + 2.0 * len_scale))) if len_scale > 0 else base_pool
    )
    top_support = sampled_records[: max_sampled_pool // 2]
    top_pk = sorted(sampled_records, key=lambda x: (x[2], x[3], x[5]), reverse=True)[: max_sampled_pool // 2]

    pooled_structs: dict[str, dict[str, object]] = {}
    for struct, pairs, pk_pairs_count, total_crossings, max_cross, support_sum in top_support + top_pk:
        pooled_structs.setdefault(struct, {})
        pooled_structs[struct].update(
            {
                "struct": struct,
                "pairs": pairs,
                "is_refined": False,
                "pk_pairs_count": pk_pairs_count,
                "total_crossings": total_crossings,
                "max_cross": max_cross,
                "support_sum": support_sum,
            }
        )

    candidates: list[dict[str, object]] = list(refined_candidates.values()) + list(pooled_structs.values())

    # --- Seed injection (early diversity without relying on MCMC) ---
    # Add a bounded set of AllSub-derived global folds for weak-evidence cases.  This is most
    # impactful when cmscan has no hit (fallback.sto) and the refined variants are local edits of an
    # incorrect MFE scaffold.
    seed_candidates: list[dict[str, object]] = []
    if thermo_mode != "off" and thermo_weight > 0 and (seed_boost or force_allsub_output > 0):
        if not thermo_structs_full and allsub_exe.exists():
            try:
                thermo_structs_full = rup.call_rnastructure_allsub(allsub_exe, full_seq)
            except Exception:
                thermo_structs_full = []
        # Pull a much larger AllSub prefix when we intend to force AllSub candidates into the
        # final top-K. With ensemble backends (RS/LF/EF), the AllSub stream is interleaved across
        # sources, so the first ~3*K entries roughly correspond to the first ~K entries of each
        # source.
        allsub_cap = max(200, min(3000, int(round(3 * float(top_k)))))
        if force_allsub_output > 0:
            allsub_cap = max(allsub_cap, min(3000, int(round(3 * float(force_allsub_output)))))
        env_max = (os.environ.get("RNANNEAL_SS_SUBOPT_MAX") or "").strip()
        if env_max:
            try:
                allsub_cap = max(allsub_cap, min(3000, int(env_max)))
            except ValueError:
                pass
        allsub_pool = thermo_structs_full[:allsub_cap]
        for seed_order, (struct, energy) in enumerate(allsub_pool, start=1):
            # Normally we skip AllSub structures already present in the refined/sampled pool to
            # avoid ballooning the candidate list. When force_allsub_output>0, we intentionally
            # keep them as explicit seed candidates so the final top-K selection can include them
            # even if refined-order/MMR heuristics would otherwise drop them.
            if struct in refined_candidates:
                if force_allsub_output <= 0:
                    continue
                existing = refined_candidates[struct]
                existing.setdefault("is_seed", True)
                existing.setdefault("seed_kind", "allsub")
                existing.setdefault("seed_order", int(seed_order))
                existing.setdefault("seed_energy", float(energy))
                continue
            if struct in pooled_structs:
                if force_allsub_output <= 0:
                    continue
                existing = pooled_structs[struct]
                existing.setdefault("is_seed", True)
                existing.setdefault("seed_kind", "allsub")
                existing.setdefault("seed_order", int(seed_order))
                existing.setdefault("seed_energy", float(energy))
                continue
            pairs = normalized_pairs(struct)
            pk_pairs_count, total_crossings, max_cross = compute_pk_stats(pairs)
            support_sum = float(sum(weights.get(p, 0.0) for p in pairs))
            seed_candidates.append(
                {
                    "struct": struct,
                    "pairs": pairs,
                    "is_refined": False,
                    "is_seed": True,
                    "seed_kind": "allsub",
                    "seed_order": int(seed_order),
                    "seed_energy": float(energy),
                    "pk_pairs_count": pk_pairs_count,
                    "total_crossings": total_crossings,
                    "max_cross": max_cross,
                    "support_sum": support_sum,
                }
            )

    # If the scaffold is extremely sparse (or AllSub produced very few candidates), add helix-only
    # seeds derived from simple complementarity scanning to cover AU-rich ambiguity.
    allsub_small = bool(thermo_structs_full) and len(thermo_structs_full) <= 10
    add_helix_seeds = bool(consensus_pairs_n <= 6) or allsub_small
    if add_helix_seeds:
        # DP over (i,j): maximum contiguous helix length starting at (i,j).
        seq = full_seq
        Ls = len(seq)
        dp = [[0] * Ls for _ in range(Ls)]
        for i in range(Ls - 1, -1, -1):
            for j in range(i + 1, Ls):
                if scs.is_canonical_pair(seq[i], seq[j]):
                    dp[i][j] = 1 + (dp[i + 1][j - 1] if i + 1 < j else 0)
        helix_structs: list[str] = []
        seen_pairs: set[frozenset[tuple[int, int]]] = set()
        candidates_ij: list[tuple[int, int, int]] = []
        for i in range(Ls):
            for j in range(i + 1, Ls):
                l = dp[i][j]
                if l < 3:
                    continue
                # Only keep outermost stems to avoid duplicates.
                if i > 0 and j + 1 < Ls and dp[i - 1][j + 1] >= l + 1:
                    continue
                candidates_ij.append((i, j, min(l, 12)))
        candidates_ij.sort(key=lambda x: (-x[2], x[0], -x[1]))
        for i, j, l in candidates_ij[:500]:
            ps = frozenset((i + k, j - k) for k in range(l) if (i + k) < (j - k))
            if ps in seen_pairs or len(ps) < 3:
                continue
            seen_pairs.add(ps)
            helix_structs.append(scs.pairs_to_pk_string(sorted(ps), Ls))
            if len(helix_structs) >= 50:
                break

        for struct in helix_structs:
            if struct in refined_candidates or struct in pooled_structs:
                continue
            pairs = normalized_pairs(struct)
            pk_pairs_count, total_crossings, max_cross = compute_pk_stats(pairs)
            support_sum = float(sum(weights.get(p, 0.0) for p in pairs))
            seed_candidates.append(
                {
                    "struct": struct,
                    "pairs": pairs,
                    "is_refined": False,
                    "is_seed": True,
                    "seed_kind": "helix",
                    "pk_pairs_count": pk_pairs_count,
                    "total_crossings": total_crossings,
                    "max_cross": max_cross,
                    "support_sum": support_sum,
                }
            )

    # tRNA-like end-stem register ambiguity: add a few terminal-stem variants.
    if rfam_id == "RF00005" and consensus_structs:
        try:
            base_struct = consensus_structs[0]
            base_pairs = set(normalized_pairs(base_struct))
            Ls = len(full_seq)
            window = 10
            internal = {p for p in base_pairs if not (min(p) < window and max(p) >= Ls - window)}

            # Reuse the dp from helix seeding when available; otherwise compute quickly.
            if not add_helix_seeds:
                seq = full_seq
                dp = [[0] * Ls for _ in range(Ls)]
                for i in range(Ls - 1, -1, -1):
                    for j in range(i + 1, Ls):
                        if scs.is_canonical_pair(seq[i], seq[j]):
                            dp[i][j] = 1 + (dp[i + 1][j - 1] if i + 1 < j else 0)

            term_structs: list[str] = []
            seen = set()
            for i in range(0, window):
                for j in range(Ls - window, Ls):
                    l = dp[i][j] if i < Ls and j < Ls else 0
                    if l < 3:
                        continue
                    l = min(l, window)
                    stem = {(i + k, j - k) for k in range(l) if (i + k) < (j - k)}
                    used = {x for p in stem for x in p}
                    if any(x in {y for p in internal for y in p} for x in used):
                        continue
                    ps = frozenset(internal | stem)
                    if ps in seen:
                        continue
                    seen.add(ps)
                    term_structs.append(scs.pairs_to_pk_string(sorted(ps), Ls))
                    if len(term_structs) >= 20:
                        break
                if len(term_structs) >= 20:
                    break

            for struct in term_structs:
                if struct in refined_candidates or struct in pooled_structs:
                    continue
                pairs = normalized_pairs(struct)
                pk_pairs_count, total_crossings, max_cross = compute_pk_stats(pairs)
                support_sum = float(sum(weights.get(p, 0.0) for p in pairs))
                seed_candidates.append(
                    {
                        "struct": struct,
                        "pairs": pairs,
                        "is_refined": False,
                        "is_seed": True,
                        "seed_kind": "terminal",
                        "pk_pairs_count": pk_pairs_count,
                        "total_crossings": total_crossings,
                        "max_cross": max_cross,
                        "support_sum": support_sum,
                    }
                )
        except Exception:
            pass

    if seed_candidates:
        # Keep the candidate pool bounded.
        seed_seen = set()
        for cand in seed_candidates:
            struct = str(cand["struct"])
            if struct in seed_seen:
                continue
            seed_seen.add(struct)
            candidates.append(cand)
    if not candidates:
        return []

    # --- Scorer ensemble ---
    pk_base = float(params.pk_alpha)
    scorers: dict[str, tuple[float | None, float]] = {
        "base": (None, 1.0),
        "pk0": (0.0, 1.0),
        "pk_low": (0.25 * pk_base, 1.0),
        "no_pair_pen": (None, 0.0),
    }
    if len_scale > 0 and pair_penalty:
        # Encourage inclusion of sparser candidates for long sequences, where over-pairing is a
        # common failure mode and F1@100 is sensitive to FP-heavy structures.
        scorers["pair_pen2"] = (None, 2.0)

    for cand in candidates:
        pairs = cand["pairs"]
        assert isinstance(pairs, frozenset)
        cand_scores: dict[str, float] = {}
        for name, (pk_override, pair_pen_mult) in scorers.items():
            cand_scores[name] = score_pairs(
                pairs=pairs,
                pk_alpha_override=pk_override,
                pair_penalty_multiplier=pair_pen_mult,
            )
        cand["scores"] = cand_scores

    def normalized_rank(values: list[float]) -> list[float]:
        if not values:
            return []
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        denom = max(1, len(values) - 1)
        for r, i in enumerate(order):
            out[i] = r / denom
        return out

    for name in scorers:
        vals = [float(c["scores"][name]) for c in candidates]  # type: ignore[index]
        ranks = normalized_rank(vals)
        for cand, rk in zip(candidates, ranks, strict=True):
            cand.setdefault("rank", {})[name] = rk

    for cand in candidates:
        ranks = cand["rank"]
        assert isinstance(ranks, dict)
        rank_vals = [float(ranks[name]) for name in scorers]
        cand["rank_min"] = float(min(rank_vals))
        cand["rank_mean"] = float(sum(rank_vals) / len(rank_vals))
        cand["rank_combo"] = float(0.7 * cand["rank_mean"] + 0.3 * cand["rank_min"])
        cand["rank_base"] = float(ranks["base"])
        cand["rank_pk0"] = float(ranks["pk0"])

    def _write_qc(selected_cands: list[dict[str, object]]) -> None:
        try:
            final_pairs: list[frozenset[tuple[int, int]]] = []
            final_structs: list[str] = []
            n_selected_refined = 0
            for cand in selected_cands:
                final_structs.append(str(cand["struct"]))
                if bool(cand.get("is_refined")):
                    n_selected_refined += 1
                ps = cand.get("pairs")
                if isinstance(ps, frozenset):
                    final_pairs.append(ps)
                else:
                    final_pairs.append(normalized_pairs(str(cand["struct"])))

            counts = [len(ps) for ps in final_pairs]
            union: set[tuple[int, int]] = set()
            for ps in final_pairs:
                union |= set(ps)

            top1_pairs = final_pairs[0] if final_pairs else frozenset()
            jacc_to_top1 = (
                [jaccard_distance(ps, top1_pairs) for ps in final_pairs[1:]] if final_pairs else []
            )
            pairwise: list[float] = []
            for i in range(len(final_pairs)):
                for j in range(i + 1, len(final_pairs)):
                    pairwise.append(jaccard_distance(final_pairs[i], final_pairs[j]))

            seed_kind_counts: dict[str, int] = {}
            selected_seed_kind_counts: dict[str, int] = {}
            for cand in candidates:
                if not bool(cand.get("is_seed")):
                    continue
                kind = str(cand.get("seed_kind", "unknown"))
                seed_kind_counts[kind] = seed_kind_counts.get(kind, 0) + 1
            for cand in selected_cands:
                if not bool(cand.get("is_seed")):
                    continue
                kind = str(cand.get("seed_kind", "unknown"))
                selected_seed_kind_counts[kind] = selected_seed_kind_counts.get(kind, 0) + 1

            t_done = time.perf_counter()
            qc = {
                "inputs": {
                    "rfam_id": str(rfam_id),
                    "cov_present": bool(cov is not None and cov.exists()),
                    "consensus_pairs_n": int(consensus_pairs_n),
                    "fasta_had_nonstandard_bases": bool(fasta_had_nonstandard_bases),
                    "seed_boost": bool(seed_boost),
                },
                "sequence_length": int(L_fasta),
                "length_adaptive": bool(length_adaptive),
                "length_adaptive_scale": float(len_scale),
                "effective": {
                    "n_samples": int(n_samples_eff),
                    "burn_in": int(burn_in_eff),
                    "thin": int(thin_eff),
                    "max_scaffolds": int(max_scaffolds_eff),
                    "max_samples_per_scaffold": int(max_samples_per_scaffold_eff),
                    "refine_max_seconds": None
                    if refine_max_seconds_eff is None
                    else float(refine_max_seconds_eff),
                    "refine_max_structures": int(refine_max_structures_eff),
                    "refine_max_regions": int(refine_max_regions_eff),
                    "refine_max_seeds": int(refine_max_seeds_eff),
                    "refine_max_solutions": int(refine_max_solutions_eff),
                    "refine_kissing_candidates": int(refine_kissing_candidates_eff),
                    "include_unfixed_sampling": bool(include_unfixed_sampling),
                },
                "counts": {
                    "n_refined": int(len(refined_structs)),
                    "n_sampled_fixed": int(n_sampled_fixed),
                    "n_sampled_unfixed": int(n_sampled_unfixed),
                    "n_sampled_total": int(len(sampled_structs)),
                    "n_candidates": int(len(candidates)),
                    "n_selected": int(len(selected_cands)),
                    "n_selected_refined": int(n_selected_refined),
                },
                "seeds": {
                    "n_seed_candidates": int(sum(seed_kind_counts.values())),
                    "seed_kind_counts": seed_kind_counts,
                    "n_selected_seed": int(sum(selected_seed_kind_counts.values())),
                    "selected_seed_kind_counts": selected_seed_kind_counts,
                },
                "selected_pairs": {
                    "min": int(min(counts) if counts else 0),
                    "median": float(sorted(counts)[len(counts) // 2] if counts else 0),
                    "max": int(max(counts) if counts else 0),
                    "mean": float(sum(counts) / len(counts) if counts else 0.0),
                    "union_pairs": int(len(union)),
                },
                "diversity": {
                    "mean_jaccard_to_top1": float(
                        sum(jacc_to_top1) / len(jacc_to_top1) if jacc_to_top1 else 0.0
                    ),
                    "mean_pairwise_jaccard": float(sum(pairwise) / len(pairwise) if pairwise else 0.0),
                },
                "timing_seconds": {
                    "total": float(t_done - t0),
                    "refine": float(t_refine1 - t_refine0),
                    "sample_fixed": float(t_sample_fixed1 - t_sample_fixed0),
                    "sample_unfixed": None
                    if t_sample_unfixed0 is None or t_sample_unfixed1 is None
                    else float(t_sample_unfixed1 - t_sample_unfixed0),
                },
            }
            (out_dir / "qc.json").write_text(json.dumps(qc, indent=2) + "\n")
        except Exception:
            return

    # --- Multi-objective top-K selection (MMR / coverage@K proxy) ---
    # Even when we have fewer than K candidates, the ordering still matters for prefix (@50/@100)
    # oracle metrics.  So we always run the selection+ordering logic, but cap K to |candidates|.
    K_requested = int(top_k)
    K = min(K_requested, len(candidates))

    selected: list[dict[str, object]] = []
    selected_structs: set[str] = set()
    selected_pairs: list[frozenset[tuple[int, int]]] = []

    def add(cand: dict[str, object]) -> None:
        struct = str(cand["struct"])
        pairs = cand["pairs"]
        assert isinstance(pairs, frozenset)
        selected.append(cand)
        selected_structs.add(struct)
        selected_pairs.append(pairs)

    refined_all = [c for c in candidates if bool(c.get("is_refined"))]

    # Baseline behavior relied heavily on the refined DB ordering. In practice, the best-of-100
    # winner is often an early refined candidate (even when our proxy scorers disagree), so
    # preserve a large prefix of refined structures when covariation is absent.
    nseq = len(seqs)
    refined_prefix = 0
    reserved_allsub = max(0, min(int(force_allsub_output), int(K)))

    def _seed_order(c: dict[str, object]) -> int:
        try:
            return int(c.get("seed_order", 10**9))
        except Exception:
            return 10**9

    def _pick_evenly(pool: list[dict[str, object]], k: int) -> list[dict[str, object]]:
        if k <= 0 or not pool:
            return []
        if k >= len(pool):
            return pool
        if k == 1:
            return [pool[len(pool) // 2]]
        idxs = sorted({int(round(i * (len(pool) - 1) / (k - 1))) for i in range(k)})
        return [pool[i] for i in idxs]

    def _alternate_ends(items: list[dict[str, object]]) -> list[dict[str, object]]:
        out: list[dict[str, object]] = []
        i = 0
        j = len(items) - 1
        while i <= j:
            out.append(items[i])
            i += 1
            if i <= j:
                out.append(items[j])
                j -= 1
        return out

    allsub_sorted: list[dict[str, object]] = []
    if reserved_allsub > 0:
        allsub_pool = [c for c in candidates if str(c.get("seed_kind", "")) == "allsub"]
        allsub_sorted = sorted(allsub_pool, key=lambda c: (_seed_order(c), str(c["struct"])))

    if refined_all:
        if rfam_id == "no_rfam_hit":
            # For the fallback-only regime, refined candidates are local edits of an MFE scaffold
            # and can miss the true global fold; reserve substantial budget for diverse thermo seeds.
            refined_prefix = min(len(refined_all), min(K, max(60, int(round(0.70 * K)))))
        elif nseq <= 1:
            if len_scale > 0:
                frac = 0.80 - 0.40 * len_scale  # 0.80 @ <=150nt -> 0.40 @ 300nt
                refined_prefix = min(len(refined_all), min(K, max(20, int(round(frac * K)))))
            else:
                refined_prefix = min(len(refined_all), min(K, max(80, int(round(0.80 * K)))))
        else:
            refined_prefix = min(len(refined_all), min(K, max(20, int(round(0.30 * K)))))

    if reserved_allsub > 0:
        refined_prefix = min(refined_prefix, max(0, int(K) - int(reserved_allsub)))

    # If we generated many helix seeds (very sparse scaffold), cap refined-prefix to leave room.
    # This helps AU-rich / low-information targets where the true structure is a small helix that
    # the refined ordering may never reach.
    if refined_prefix > 0 and rfam_id != "no_rfam_hit" and consensus_pairs_n <= 6:
        helix_seed_n = sum(1 for c in candidates if str(c.get("seed_kind")) == "helix")
        if helix_seed_n > 0:
            reserve_other = 10  # keep some budget for non-seed MMR picks
            refined_prefix = min(
                refined_prefix,
                max(20, int(K) - int(helix_seed_n) - int(reserve_other)),
            )

    # Improve @1: pick the best scorer-ensemble candidate as the first structure.
    if candidates:
        best_overall = min(candidates, key=lambda c: (float(c["rank_combo"]), str(c["struct"])))
        add(best_overall)

    # Prefix-aware selection: ensure some AllSub candidates appear very early so prefixes (@50/@100)
    # contain thermo suboptimals even when refined-order is strongly front-loaded.
    if reserved_allsub > 0 and allsub_sorted:
        early_budget = min(K, 100)
        early_quota = min(reserved_allsub, int(round(0.25 * early_budget)))
        if K >= 50 and early_quota > 0:
            early_quota = max(10, early_quota)
        early_quota = min(int(early_quota), max(0, int(K) - len(selected)))
        if early_quota > 0:
            span = min(len(allsub_sorted), max(500, int(round(3.0 * float(K)))))
            pool = allsub_sorted[:span]
            head_n = min(len(allsub_sorted), min(10, int(early_quota)))
            probe = list(allsub_sorted[:head_n])
            if cov_empty:
                # Some weak-cov targets have their best AllSub structure slightly beyond the
                # very first few MFE-adjacent candidates (often around seed_order ~80–120).
                # Seed an additional early band so @200 doesn't miss these cases.
                band_start = min(len(pool), 80)
                band_end = min(len(pool), 140)
                if band_end > band_start:
                    probe += pool[band_start:band_end]
            remain = int(early_quota) - len(probe)
            if remain > 0:
                probe_n = min(len(pool), max(remain, remain * 3))
                probe += _alternate_ends(_pick_evenly(pool, probe_n))
            added = 0
            for cand in probe + pool + allsub_sorted:
                if added >= int(early_quota) or len(selected) >= K:
                    break
                if str(cand["struct"]) in selected_structs:
                    continue
                add(cand)
                added += 1

    # Adjust refined-prefix to ensure we can still satisfy the AllSub quota after prefix seeding.
    if refined_prefix > 0 and reserved_allsub > 0 and allsub_sorted:
        allsub_selected = sum(1 for c in selected if str(c.get("seed_kind", "")) == "allsub")
        need_allsub = max(0, int(reserved_allsub) - int(allsub_selected))
        refined_prefix = min(refined_prefix, max(0, int(K) - int(need_allsub) - len(selected)))

    if refined_prefix > 0:
        refined_pick = sorted(
            refined_all,
            key=lambda c: (
                float(c["rank_combo"]),
                int(c.get("refined_order", 10**9)),
                str(c["struct"]),
            ),
        )
        for cand in refined_pick[:refined_prefix]:
            if len(selected) >= K:
                break
            if str(cand["struct"]) in selected_structs:
                continue
            add(cand)

    def _topn(
        pool: list[dict[str, object]],
        key: str,
        n: int,
        reverse: bool = False,
    ) -> list[dict[str, object]]:
        return sorted(
            pool,
            key=lambda c: (
                float(c[key]) * (-1.0 if reverse else 1.0),
                str(c["struct"]),
            ),
        )[: max(0, int(n))]

    must_seed: list[dict[str, object]] = []
    if refined_all:
        must_seed += _topn(refined_all, "rank_base", 2)
        must_seed += _topn(refined_all, "rank_pk0", 2)
        must_seed += _topn(refined_all, "rank_min", 2)
        must_seed += _topn(refined_all, "rank_combo", 2)
        must_seed += _topn(refined_all, "support_sum", 2, reverse=True)

    for cand in must_seed:
        if len(selected) >= K:
            break
        if str(cand["struct"]) in selected_structs:
            continue
        add(cand)

    # Prefix-aware: pull a few high-scoring sampled candidates early so @200 can benefit from MCMC
    # rather than being dominated by refined/scaffold seeds.
    sampled_all = [
        c for c in candidates if (not bool(c.get("is_refined"))) and (not bool(c.get("is_seed")))
    ]
    if sampled_all:
        sampled_must: list[dict[str, object]] = []
        sampled_must += _topn(sampled_all, "rank_combo", 4)
        sampled_must += _topn(sampled_all, "rank_min", 4)
        sampled_must += sorted(
            sampled_all,
            key=lambda c: (float(c["support_sum"]), float(c["rank_combo"]), str(c["struct"])),
            reverse=True,
        )[:4]
        pk_sampled = [c for c in sampled_all if int(c.get("pk_pairs_count", 0)) > 0]
        if pk_sampled:
            sampled_must += _topn(pk_sampled, "rank_combo", 4)
            pk_sampled.sort(
                key=lambda c: (
                    int(c.get("pk_pairs_count", 0)),
                    int(c.get("total_crossings", 0)),
                    float(c.get("support_sum", 0.0)),
                ),
                reverse=True,
            )
            sampled_must += pk_sampled[:4]
        for cand in sampled_must:
            if len(selected) >= K:
                break
            if str(cand["struct"]) in selected_structs:
                continue
            add(cand)

    # If configured, ensure we include a large slice of AllSub-like candidates in the final top-K.
    # This is important for scaffold backends (LF/EF) where suboptimals can contain the best-of-100
    # structure, and where our proxy energy scorers might otherwise over-prioritize refined ordering.
    if reserved_allsub > 0 and allsub_sorted:
        allsub_selected = sum(1 for c in selected if str(c.get("seed_kind", "")) == "allsub")
        need = min(
            max(0, int(reserved_allsub) - int(allsub_selected)),
            max(0, int(K) - len(selected)),
        )
        if need > 0:
            # Spread picks across a wide span of the AllSub stream so we include both MFE-adjacent
            # and mid-energy alternatives (often important for tertiary-stabilized folds).
            span = min(len(allsub_sorted), max(500, int(round(6.0 * float(K)))))
            pool = allsub_sorted[:span]
            added = 0
            # First, add a small evenly-spaced probe set so mid-range suboptimals can land within
            # early prefixes (@200) instead of being pushed to the tail of the AllSub quota.
            probe_n = min(int(need), 24)
            if cov_empty:
                # Weak-evidence heuristic: reserve a small contiguous window around mid-to-high
                # energies (where the oracle-best structure can appear) rather than only sampling
                # sparse representatives. This is intentionally seed-order based (not proxy-score
                # based) so we keep candidates even when the proxy scorer misranks them.
                probe: list[dict[str, object]] = []
                seen_probe: set[str] = set()

                def add_probe(cand: dict[str, object]) -> None:
                    s = str(cand["struct"])
                    if s in seen_probe:
                        return
                    seen_probe.add(s)
                    probe.append(cand)

                window = min(len(pool), max(12, int(round(0.12 * float(need)))))
                # Include a small early-energy window (≈3% of the stream) because some targets'
                # best AllSub structure is not in the very first few seeds, but still near the
                # front of the ordering.
                centers = [0.03, 0.10, 0.50, 0.75, 0.90]
                for q in centers:
                    if len(probe) >= int(need):
                        break
                    if not pool:
                        break
                    c = int(round(q * float(len(pool) - 1)))
                    start = max(0, int(c) - (window // 2))
                    end = min(len(pool), start + window)
                    start = max(0, end - window)
                    for cand in pool[start:end]:
                        add_probe(cand)

                remaining_need = max(0, int(need) - len(probe))
                if remaining_need > 0:
                    fill_n = min(len(pool), max(remaining_need, remaining_need * 3))
                    fill = _pick_evenly(pool, int(fill_n))
                    for cand in _alternate_ends(fill):
                        if len(probe) >= int(need):
                            break
                        add_probe(cand)

                probe = probe[: int(need)]
            else:
                probe = _pick_evenly(pool, probe_n)

            for cand in (probe if cov_empty else _alternate_ends(probe)):
                if added >= need or len(selected) >= K:
                    break
                if str(cand["struct"]) in selected_structs:
                    continue
                add(cand)
                added += 1
            for cand in pool + allsub_sorted:
                if added >= need or len(selected) >= K:
                    break
                if str(cand["struct"]) in selected_structs:
                    continue
                add(cand)
                added += 1

    # Guardrail: keep some explicitly-injected seeds for weak-evidence cases.
    seed_all = [c for c in candidates if bool(c.get("is_seed"))]
    if seed_all and (seed_boost or rfam_id == "RF00005"):
        seed_must: list[dict[str, object]] = []

        if rfam_id == "no_rfam_hit":
            # Empirically, for no-Rfam targets the best AllSub structure is often not in the
            # first few (MFE-adjacent) candidates.  Pull a broader band so the top-100 contains
            # mid-energy alternatives that can match tertiary-constraint folds.
            allsub_pool = [c for c in seed_all if str(c.get("seed_kind")) == "allsub"]

            def _seed_order(c: dict[str, object]) -> int:
                try:
                    return int(c.get("seed_order", 10**9))
                except Exception:
                    return 10**9

            allsub_sorted = sorted(allsub_pool, key=lambda c: (_seed_order(c), str(c["struct"])))
            mid = [c for c in allsub_sorted if 80 <= _seed_order(c) <= 200]

            def _pick_evenly(pool: list[dict[str, object]], k: int) -> list[dict[str, object]]:
                if k <= 0 or not pool:
                    return []
                if k >= len(pool):
                    return pool
                if k == 1:
                    return [pool[len(pool) // 2]]
                idxs = sorted({int(round(i * (len(pool) - 1) / (k - 1))) for i in range(k)})
                return [pool[i] for i in idxs]

            seed_quota = min(len(allsub_sorted), max(20, int(round(0.30 * K))))
            pick_from = mid if len(mid) >= seed_quota else allsub_sorted
            seed_must += _pick_evenly(pick_from, seed_quota)
        else:
            # For extremely sparse scaffolds, explicitly include helix seeds (cheap, high recall).
            # These are generated from canonical complementarity scanning and can rescue cases
            # where both CaCoFold and thermo tools miss the correct small helix.
            if consensus_pairs_n <= 6:
                helix_pool = [c for c in seed_all if str(c.get("seed_kind")) == "helix"]
                if helix_pool:
                    seed_must += sorted(
                        helix_pool,
                        key=lambda c: (float(c["rank_combo"]), str(c["struct"])),
                    )
        if rfam_id != "no_rfam_hit" and seed_boost:
            # Prefer AllSub seeds when available, but fall back to any seed kind.
            seed_pool = [c for c in seed_all if str(c.get("seed_kind")) == "allsub"] or seed_all
            seed_must += _topn(seed_pool, "rank_base", 2)
            seed_must += _topn(seed_pool, "rank_pk0", 2)
            seed_must += _topn(seed_pool, "rank_min", 2)
            seed_must += _topn(seed_pool, "rank_combo", 2)
            seed_must += _topn(seed_pool, "support_sum", 2, reverse=True)
        if rfam_id == "RF00005":
            term_pool = [c for c in seed_all if str(c.get("seed_kind")) == "terminal"]
            if term_pool:
                seed_must += _topn(term_pool, "rank_combo", 4)

        for cand in seed_must:
            if len(selected) >= K:
                break
            if str(cand["struct"]) in selected_structs:
                continue
            add(cand)

    remaining = max(0, K - len(selected))
    q_energy = int(round(0.50 * remaining))
    q_pk = int(round(0.20 * remaining))
    q_evidence = int(round(0.10 * remaining))
    q_diverse = max(0, remaining - q_energy - q_pk - q_evidence)

    def mmr_pick(
        pool: list[dict[str, object]],
        score_key: str,
        score_weight: float,
        diversity_weight: float,
    ) -> dict[str, object] | None:
        # Back-compat helper (kept for small K). For large K, recomputing min-distance
        # to all selected pairs inside the inner loop is prohibitively expensive.
        best: dict[str, object] | None = None
        best_u: float | None = None
        best_tie: tuple[float, str] | None = None
        for cand in pool:
            struct = str(cand["struct"])
            if struct in selected_structs:
                continue
            pairs = cand["pairs"]
            assert isinstance(pairs, frozenset)
            base = 1.0 - float(cand[score_key])
            min_dist = (
                min(jaccard_distance(pairs, sp) for sp in selected_pairs) if selected_pairs else 1.0
            )
            u = score_weight * base + diversity_weight * min_dist
            tie = (float(cand[score_key]), struct)
            if best_u is None or u > best_u or (u == best_u and tie < (best_tie or tie)):
                best = cand
                best_u = u
                best_tie = tie
        return best

    def select_bucket(
        pool: list[dict[str, object]],
        k: int,
        score_key: str,
        score_weight: float,
        diversity_weight: float,
    ) -> None:
        # For large K, naive MMR is O(k * |pool| * |selected|) because it recomputes the
        # min-distance to the selected set for every candidate at every selection step.
        #
        # Use an incremental approximation: maintain each candidate's current min-distance to
        # the selected set and update it after each addition. This reduces the complexity to
        # O(k * |pool|) jaccard computations (plus an O(|pool| * |selected|) initialization).
        k = int(k)
        if k <= 0 or not pool:
            return

        # Heuristic threshold: once K grows beyond 150, the incremental path is much faster.
        # (For K<=100, the original behavior is already OK and keeps exact tie-breaking.)
        if K <= 150:
            for _ in range(k):
                if len(selected) >= K:
                    return
                nxt = mmr_pick(
                    pool,
                    score_key=score_key,
                    score_weight=score_weight,
                    diversity_weight=diversity_weight,
                )
                if nxt is None:
                    return
                add(nxt)
            return

        structs = [str(c["struct"]) for c in pool]
        pairs_list: list[frozenset[tuple[int, int]]] = []
        scores: list[float] = []
        for c in pool:
            ps = c["pairs"]
            assert isinstance(ps, frozenset)
            pairs_list.append(ps)
            scores.append(float(c[score_key]))

        # Initialize min distances to the current selected set.
        min_dist = [1.0] * len(pool)
        if selected_pairs:
            for i, ps in enumerate(pairs_list):
                md = 1.0
                for sp in selected_pairs:
                    d = jaccard_distance(ps, sp)
                    if d < md:
                        md = d
                        if md <= 0.0:
                            break
                min_dist[i] = md

        for _ in range(k):
            if len(selected) >= K:
                return

            best_i: int | None = None
            best_u: float | None = None
            best_tie: tuple[float, str] | None = None
            for i, cand in enumerate(pool):
                struct = structs[i]
                if struct in selected_structs:
                    continue
                base = 1.0 - scores[i]
                u = score_weight * base + diversity_weight * min_dist[i]
                tie = (scores[i], struct)
                if best_u is None or u > best_u or (u == best_u and tie < (best_tie or tie)):
                    best_i = i
                    best_u = u
                    best_tie = tie

            if best_i is None:
                return

            nxt = pool[best_i]
            add(nxt)
            new_pairs = pairs_list[best_i]

            # Update min-distances in one pass.
            for i in range(len(pool)):
                struct = structs[i]
                if struct in selected_structs:
                    continue
                d = jaccard_distance(pairs_list[i], new_pairs)
                if d < min_dist[i]:
                    min_dist[i] = d

    # Bucket 1: energy-ish (ensemble rank combo), mild diversity.
    # Anchor on the original/base scorer to avoid dropping high-quality candidates.
    pool_limit = (
        int(round(float(max(2000, K * 20)) * (1.0 + 2.0 * len_scale)))
        if len_scale > 0
        else max(2000, K * 20)
    )
    base_pool = sorted(candidates, key=lambda c: (float(c["rank_base"]), str(c["struct"])))
    base_pool = base_pool[:pool_limit]
    select_bucket(base_pool, q_energy, score_key="rank_base", score_weight=1.0, diversity_weight=0.25)

    # Bucket 2: PK-heavy (prefer high PK load, scored by pk0-rank), moderate diversity.
    pk_pool = [c for c in candidates if int(c["pk_pairs_count"]) > 0]
    pk_pool.sort(
        key=lambda c: (
            int(c["pk_pairs_count"]),
            int(c["total_crossings"]),
            1.0 - float(c["rank_pk0"]),
        ),
        reverse=True,
    )
    pk_pool = pk_pool[:pool_limit]
    select_bucket(pk_pool, q_pk, score_key="rank_pk0", score_weight=0.9, diversity_weight=0.6)

    # Bucket 3: evidence extremes (support_sum), mild diversity.
    ev_pool = sorted(candidates, key=lambda c: float(c["support_sum"]), reverse=True)
    ev_pool = ev_pool[:pool_limit]
    # Convert support_sum to a rank in [0,1] for MMR.
    ev_vals = [-float(c["support_sum"]) for c in ev_pool]
    ev_ranks = normalized_rank(ev_vals)
    for cand, rk in zip(ev_pool, ev_ranks, strict=True):
        cand["_ev_rank"] = rk
    select_bucket(ev_pool, q_evidence, score_key="_ev_rank", score_weight=0.9, diversity_weight=0.3)

    # Bucket 4: fill remaining with strong diversity.
    # Use rank_min so candidates that are good under *any* scorer can make it in.
    select_bucket(base_pool, q_diverse, score_key="rank_min", score_weight=0.7, diversity_weight=0.8)

    # Guardrail: ensure some refined coverage (helps when the scorer ensemble is noisy).
    refined_quota = min(len(refined_all), max(10, K // 5))
    refined_selected = sum(1 for c in selected if bool(c.get("is_refined")))
    if refined_selected < refined_quota:
        refined_pool = sorted(refined_all, key=lambda c: (float(c["rank_base"]), str(c["struct"])))
        select_bucket(
            refined_pool,
            refined_quota - refined_selected,
            score_key="rank_base",
            score_weight=1.0,
            diversity_weight=0.5,
        )

    # If we still have slack (e.g. all buckets exhausted early), append remaining candidates in
    # rank_combo order so we always emit exactly K outputs (or all candidates when |candidates|<K).
    if len(selected) < K:
        for cand in sorted(candidates, key=lambda c: (float(c["rank_combo"]), str(c["struct"]))):
            if len(selected) >= K:
                break
            if str(cand["struct"]) in selected_structs:
                continue
            add(cand)

    _write_qc(selected[:K])
    return [str(c["struct"]) for c in selected[:K]]


def main() -> None:
    parser = argparse.ArgumentParser(description="CaCoFold + MCMC predictor")
    parser.add_argument("--cm-db", required=True, help="Path to Rfam CM database")
    parser.add_argument("--rscape", default="rscape", help="Path to R-scape binary")
    parser.add_argument("--infernal-bin", default=None, help="Path to Infernal bin dir")
    parser.add_argument("--allsub-exe", required=True, help="Path to RNAstructure AllSub")
    parser.add_argument("--duplex-exe", default=None, help="Path to RNAstructure DuplexFold")
    parser.add_argument("--fold-exe", default=None, help="Path to RNAstructure Fold (MFE fallback)")
    parser.add_argument("--partition-exe", default=None, help="Path to RNAstructure Partition")
    parser.add_argument(
        "--probplot-exe", default=None, help="Path to RNAstructure ProbabilityPlot"
    )
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--burn-in", type=int, default=500)
    parser.add_argument("--thin", type=int, default=10)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--min-loop-sep",
        type=int,
        default=0,
        help="Minimum separation between paired indices (0 disables this constraint).",
    )
    parser.add_argument(
        "--pk-alpha",
        type=float,
        default=0.5,
        help="Pseudoknot penalty coefficient (lower makes PK more likely).",
    )
    parser.add_argument(
        "--pair-penalty",
        type=float,
        default=None,
        help="Optional per-base-pair cost added to the scoring energy.",
    )
    parser.add_argument(
        "--pair-penalty-scale",
        type=float,
        default=0.25,
        help="If --pair-penalty unset, set it to this times the median weight for the target.",
    )
    parser.add_argument(
        "--pair-penalty-mode",
        choices=["legacy", "length_aware"],
        default="legacy",
        help="How to set the default per-pair penalty when --pair-penalty is unset.",
    )
    parser.add_argument("--pair-penalty-c0", type=float, default=0.10)
    parser.add_argument("--pair-penalty-c1", type=float, default=0.50)
    parser.add_argument("--pair-penalty-min", type=float, default=None)
    parser.add_argument("--pair-penalty-max", type=float, default=None)
    parser.add_argument("--cov-mode", default="logE_power")
    parser.add_argument("--cov-alpha", type=float, default=3.0)
    parser.add_argument("--cov-min-power", type=float, default=0.1)
    parser.add_argument("--cov-forbid-negative", action="store_true")
    parser.add_argument(
        "--weight-calibration-method",
        choices=["none", "robust_z"],
        default="none",
    )
    parser.add_argument("--weight-calibration-zmax", type=float, default=3.0)
    parser.add_argument("--weight-alpha-core", type=float, default=1.0)
    parser.add_argument("--weight-alpha-alt", type=float, default=1.0)
    parser.add_argument("--weight-alpha-cov", type=float, default=1.0)
    parser.add_argument("--weight-alpha-thermo", type=float, default=1.0)
    parser.add_argument("--thermo-mode", choices=["allsub", "pf", "off"], default="allsub")
    parser.add_argument("--thermo-weight", type=float, default=1.0)
    parser.add_argument("--thermo-max-structures", type=int, default=50)
    parser.add_argument("--thermo-min-count", type=int, default=2)
    parser.add_argument("--thermo-min-prob", type=float, default=0.001)
    parser.add_argument("--thermo-log-eps", type=float, default=1e-6)
    parser.add_argument("--stem-start-penalty-scale", type=float, default=0.0)
    parser.add_argument("--stem-len1-penalty-scale", type=float, default=0.0)
    parser.add_argument("--stem-len2-penalty-scale", type=float, default=0.0)
    parser.add_argument("--stem-log-reward-scale", type=float, default=0.0)
    parser.add_argument("--stem-support-quantile", type=float, default=0.5)
    parser.add_argument("--refine-max-structures", type=int, default=200)
    parser.add_argument("--refine-min-unpaired", type=int, default=15)
    parser.add_argument("--refine-end-mask-step", type=int, default=5)
    parser.add_argument("--refine-max-end-mask-len", type=int, default=40)
    parser.add_argument("--refine-max-helices-sequential", type=int, default=20)
    parser.add_argument("--refine-max-helices-pairwise", type=int, default=10)
    parser.add_argument("--refine-max-regions", type=int, default=10)
    parser.add_argument("--refine-max-seeds", type=int, default=20)
    parser.add_argument("--refine-max-solutions", type=int, default=200)
    parser.add_argument("--refine-kissing-candidates", type=int, default=200)
    parser.add_argument(
        "--refine-max-seconds",
        type=float,
        default=None,
        help="Optional wall-clock timeout (seconds) for refining a single scaffold.",
    )
    parser.add_argument(
        "--max-scaffolds",
        type=int,
        default=20,
        help="Cap the number of refined scaffolds used for MCMC sampling.",
    )
    parser.add_argument(
        "--max-samples-per-scaffold",
        type=int,
        default=200,
        help="Cap the number of MCMC samples per scaffold (before thinning).",
    )
    parser.add_argument(
        "--length-adaptive",
        action="store_true",
        help="Enable length-aware compute allocation and selection tweaks for longer RNAs.",
    )
    parser.add_argument(
        "--include-unfixed-sampling",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Optionally run an extra sampling pass with fix_scaffold_pairs=False and merge candidates. "
            "Defaults to enabled when --length-adaptive is set."
        ),
    )
    parser.add_argument(
        "--reuse-cacofold-root",
        default=None,
        help=(
            "Optional directory containing per-target CaCoFold outputs to reuse (avoids rerunning cmscan/R-scape). "
            "Expected layout: {reuse_root}/{safe_id}/ with cmscan.tblout and *.cacofold.sto/*.cov."
        ),
    )
    parser.add_argument(
        "--inject-allsub-scaffolds",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Inject AllSub-like suboptimals into the scaffold pool before MCMC sampling (default: true).",
    )
    parser.add_argument(
        "--inject-allsub-scaffolds-max",
        type=int,
        default=25,
        help="Max number of injected AllSub-like scaffolds per target (default: 25).",
    )
    parser.add_argument(
        "--inject-allsub-timeout-s",
        type=float,
        default=10.0,
        help="Timeout (seconds) for generating injected AllSub-like scaffolds (default: 10s).",
    )
    parser.add_argument(
        "--force-allsub-output",
        type=int,
        default=0,
        help=(
            "Force inclusion of this many AllSub-like candidates in the final top-K (default: 0). "
            "Useful for LF/EF scaffold backends where best-of-100 often comes from suboptimals."
        ),
    )
    parser.add_argument("fasta", help="Input FASTA")
    parser.add_argument("outdir", help="Output dir")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[4]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ensure_rnastructure_datapath(exe_hint=Path(args.allsub_exe))

    raw_seq = read_fasta_sequence(Path(args.fasta))
    fasta_for_tools, sanitized_changed = write_sanitized_fasta(fasta_in=Path(args.fasta), out_dir=out_dir)
    # For seed-boosting heuristics, we only care whether we had to replace illegal characters
    # (e.g. '&') to make downstream tools happy. Ambiguous bases like 'X' are preserved.
    fasta_had_nonstandard_bases = bool(sanitized_changed)

    sto: Path | None = None
    cov: Path | None = None
    seq_id = parse_fasta_id(fasta_for_tools)
    rfam_id = "unknown"

    reused = False
    if args.reuse_cacofold_root:
        reuse_root = Path(str(args.reuse_cacofold_root)).expanduser().resolve()
        reuse_dir = reuse_root / out_dir.name
        if reuse_dir.is_dir():
            sto = find_cacofold_sto(reuse_dir)
            cov = find_cacofold_cov(reuse_dir)
            tbl = reuse_dir / "cmscan.tblout"
            model = parse_tblout_top_hit(tbl) if tbl.exists() else None
            if model:
                rfam_id = model
            if sto is not None:
                reused = True

    if not reused:
        try:
            sto, cov, seq_id, rfam_id = run_cacofold(
                fasta=fasta_for_tools,
                out_dir=out_dir,
                cm_db=Path(args.cm_db),
                rscape=args.rscape,
                infernal_bin=Path(args.infernal_bin) if args.infernal_bin else None,
            )
        except SystemExit as exc:
            if "No Rfam hit found" not in str(exc):
                raise
            if args.fold_exe is None:
                raise
            sto, seq_id = build_fallback_sto(
                fasta=fasta_for_tools,
                out_dir=out_dir,
                fold_exe=Path(args.fold_exe),
            )
            cov = None
            rfam_id = "no_rfam_hit"
    if sto is None:
        if args.fold_exe is None:
            raise SystemExit("No CaCoFold .sto produced and --fold-exe missing")
        sto, seq_id = build_fallback_sto(
            fasta=fasta_for_tools,
            out_dir=out_dir,
            fold_exe=Path(args.fold_exe),
        )
        cov = None
    else:
        # If CaCoFold alignment drops residues, fall back to MFE to keep lengths consistent.
        from src.lib import sample_cacofold_structures as scs

        fasta_seq = read_fasta_sequence(fasta_for_tools)
        seqs, _ = scs.parse_stockholm_single(str(sto))
        if seqs:
            sto_seq = "".join(ch for ch in next(iter(seqs.values())) if ch not in scs.GAP_CHARS)
            if len(sto_seq) != len(fasta_seq):
                sys.stderr.write(
                    "[WARN] CaCoFold sequence length mismatch; falling back to MFE scaffold.\n"
                )
                if args.fold_exe is None:
                    raise SystemExit("Length mismatch and --fold-exe missing")
                sto, seq_id = build_fallback_sto(
                    fasta=fasta_for_tools,
                    out_dir=out_dir,
                    fold_exe=Path(args.fold_exe),
                )
                cov = None

    structures = build_topk_predictions(
        sto=sto,
        cov=cov,
        seq_name=seq_id,
        rfam_id=rfam_id,
        fasta_had_nonstandard_bases=fasta_had_nonstandard_bases,
        allsub_exe=Path(args.allsub_exe),
        duplex_exe=Path(args.duplex_exe) if args.duplex_exe else None,
        fold_exe=Path(args.fold_exe) if args.fold_exe else None,
        partition_exe=Path(args.partition_exe) if args.partition_exe else None,
        probplot_exe=Path(args.probplot_exe) if args.probplot_exe else None,
        fasta=fasta_for_tools,
        out_dir=out_dir,
        top_k=args.top_k,
        n_samples=args.n_samples,
        burn_in=args.burn_in,
        thin=args.thin,
        beta=args.beta,
        seed=args.seed,
        min_loop_sep=args.min_loop_sep,
        pk_alpha=args.pk_alpha,
        pair_penalty=args.pair_penalty,
        pair_penalty_scale=args.pair_penalty_scale,
        cov_mode=args.cov_mode,
        cov_alpha=args.cov_alpha,
        cov_min_power=args.cov_min_power,
        cov_forbid_negative=args.cov_forbid_negative,
        weight_calibration_method=args.weight_calibration_method,
        weight_calibration_zmax=args.weight_calibration_zmax,
        weight_alpha_core=args.weight_alpha_core,
        weight_alpha_alt=args.weight_alpha_alt,
        weight_alpha_cov=args.weight_alpha_cov,
        weight_alpha_thermo=args.weight_alpha_thermo,
        thermo_mode=args.thermo_mode,
        thermo_weight=args.thermo_weight,
        thermo_max_structures=args.thermo_max_structures,
        thermo_min_count=args.thermo_min_count,
        thermo_min_prob=args.thermo_min_prob,
        thermo_log_eps=args.thermo_log_eps,
        stem_start_penalty_scale=args.stem_start_penalty_scale,
        stem_len1_penalty_scale=args.stem_len1_penalty_scale,
        stem_len2_penalty_scale=args.stem_len2_penalty_scale,
        stem_log_reward_scale=args.stem_log_reward_scale,
        stem_support_quantile=args.stem_support_quantile,
        pair_penalty_mode=args.pair_penalty_mode,
        pair_penalty_c0=args.pair_penalty_c0,
        pair_penalty_c1=args.pair_penalty_c1,
        pair_penalty_min=args.pair_penalty_min,
        pair_penalty_max=args.pair_penalty_max,
        refine_max_structures=args.refine_max_structures,
        refine_min_unpaired=args.refine_min_unpaired,
        refine_end_mask_step=args.refine_end_mask_step,
        refine_max_end_mask_len=args.refine_max_end_mask_len,
        refine_max_helices_sequential=args.refine_max_helices_sequential,
        refine_max_helices_pairwise=args.refine_max_helices_pairwise,
        refine_max_regions=args.refine_max_regions,
        refine_max_seeds=args.refine_max_seeds,
        refine_max_solutions=args.refine_max_solutions,
        refine_kissing_candidates=args.refine_kissing_candidates,
        refine_max_seconds=args.refine_max_seconds,
        max_scaffolds=args.max_scaffolds,
        max_samples_per_scaffold=args.max_samples_per_scaffold,
        length_adaptive=bool(args.length_adaptive),
        include_unfixed_sampling=(
            bool(args.length_adaptive)
            if args.include_unfixed_sampling is None
            else bool(args.include_unfixed_sampling)
        ),
        inject_allsub_scaffolds=bool(args.inject_allsub_scaffolds),
        inject_allsub_scaffolds_max=int(args.inject_allsub_scaffolds_max),
        inject_allsub_timeout_s=float(args.inject_allsub_timeout_s)
        if args.inject_allsub_timeout_s is not None
        else None,
        force_allsub_output=int(args.force_allsub_output),
    )

    # Write predictions.db
    # Use ungapped seq from CaCoFold consensus.
    from src.lib import refine_unpaired_regions as rup

    try:
        seq = rup.read_ungapped_seq_from_sto(sto, seq_id)
    except Exception:
        seq = read_fasta_sequence(fasta_for_tools)
    out_db = out_dir / "predictions.db"
    out_db.write_text("\n".join([seq] + structures) + "\n")


if __name__ == "__main__":
    main()
