#!/usr/bin/env python3
"""
High-level CaCoFold → sampling → refinement → Rosetta pipelines.

This module wires together:

- src.lib.sample_cacofold_structures
- src.lib.refine_unpaired_regions
- src.lib.filter_db_for_rosetta

into pre-defined pipelines:

    legacy:
        get_consensus_db → sample_pk → refine → rosetta_ready

    refine_first:
        get_consensus_db → refine → sample_pk → ensure_valid → rosetta_ready

All steps are available as Python functions so they can be reused from
other code, and there is also a small CLI (see `main()` below).
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Set

# Import the existing building blocks
from src.lib import sample_cacofold_structures as scs
from src.lib import refine_unpaired_regions as rup
from src.lib import filter_db_for_rosetta as fdr


# ---------------------------------------------------------------------------
# Config dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SampleConfig:
    n_samples: int = 1000
    burn_in: int = 1000
    thin: int = 10
    min_loop_sep: int = 1
    beta: float = 1.0
    seed: Optional[int] = None
    pk_alpha: float = 2.0
    pk_filter_frac: float = 0.2
    pk_filter_max_cross_per_pair: int = 2
    pk_filter_max_total_cross: int = 30
    cov_mode: str = "off"              # "off", "score", "power", ...
    cov_alpha: float = 1.0
    cov_min_power: float = 0.0
    cov_forbid_negative: bool = False
    cov_negative_E: float = 1.0


@dataclass
class RefineConfig:
    min_unpaired: int = 15
    min_terminal_unpaired: int = 15
    temperature: Optional[float] = None
    fold_extra_args: Sequence[str] = ()
    duplex_extra_args: Sequence[str] = ()


@dataclass
class EnsureValidConfig:
    max_total_crossings: int = 30
    max_crossings_per_pair: int = 2
    max_pk_fraction: float = 0.2  # fraction of positions that can be involved in PKs


@dataclass
class RosettaConfig:
    hamming_threshold: int = 0
    hamming_frac: float = 0.05
    max_structures: Optional[int] = None


@dataclass
class IOConfig:
    work_dir: Path
    consensus_db: Path
    refined_db: Path
    sampled_db: Path
    valid_db: Path
    rosetta_db: Path
    summary: Optional[Path]


@dataclass
class PipelineConfig:
    sto: Path
    cov: Optional[Path]
    seq_name: str
    fold_exe: Path
    duplex_exe: Optional[Path]
    sample: SampleConfig
    refine: RefineConfig
    ensure_valid: EnsureValidConfig
    rosetta: RosettaConfig
    io: IOConfig


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _resolve(path_str: Optional[str], base: Path) -> Optional[Path]:
    if path_str is None:
        return None
    p = Path(path_str)
    if not p.is_absolute():
        p = base / p
    return p


def load_pipeline_config(path: Path) -> PipelineConfig:
    """
    Load a config from JSON or YAML.

    YAML requires PyYAML; JSON uses builtin json.
    """
    text = path.read_text()
    base = path.parent

    if path.suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore
        except ImportError as e:  # pragma: no cover - only hit at runtime
            raise SystemExit(
                "PyYAML is required for YAML configs. "
                "Install with `pip install pyyaml` or use a JSON config instead."
            ) from e
        raw = yaml.safe_load(text)
    else:
        raw = json.loads(text)

    sto = _resolve(raw["sto"], base)
    cov = _resolve(raw.get("cov"), base)
    seq_name = raw["seq_name"]

    fold_exe = Path(raw["fold_exe"])
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

    # Default filenames encode pipeline order:
    #   00: consensus
    #   01: refined
    #   02: sampled
    #   03: valid
    #   04: Rosetta-ready
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

    # Basic validation
    if sto is None or not sto.is_file():
        raise FileNotFoundError(f"Stockholm file not found: {sto}")
    if cov is not None and not cov.is_file():
        raise FileNotFoundError(f"Covariation file not found: {cov}")
    if not fold_exe.is_file():
        raise FileNotFoundError(f"RNAstructure Fold executable not found: {fold_exe}")
    if duplex_exe is not None and not duplex_exe.is_file():
        raise FileNotFoundError(f"RNAstructure DuplexFold executable not found: {duplex_exe}")

    io.work_dir.mkdir(parents=True, exist_ok=True)

    return PipelineConfig(
        sto=sto,
        cov=cov,
        seq_name=seq_name,
        fold_exe=fold_exe,
        duplex_exe=duplex_exe,
        sample=sample,
        refine=refine,
        ensure_valid=ensure,
        rosetta=rosetta,
        io=io,
    )


# ---------------------------------------------------------------------------
# Step 1: get_consensus_db
# ---------------------------------------------------------------------------

def get_consensus_db(cfg: PipelineConfig) -> Path:
    """
    Project the CaCoFold SS_cons track for `cfg.seq_name` onto the ungapped
    sequence and write a one-structure .db file containing that consensus.

    Also optionally writes a human-readable summary (consensus + PK layers)
    if `cfg.io.summary` is set.
    """
    sto_path = cfg.sto
    seq_name = cfg.seq_name

    sys.stderr.write(f"[PIPELINE] Reading Stockholm file: {sto_path}\n")
    seqs, ss_tracks = scs.parse_stockholm_single(str(sto_path))

    if seq_name not in seqs:
        raise KeyError(f"Sequence {seq_name!r} not found in {sto_path}")

    # Ungapped sequence
    ungapped_seq = "".join(ch for ch in seqs[seq_name] if ch not in scs.GAP_CHARS)
    L = len(ungapped_seq)

    # Consensus projection (may fail if SS_cons missing)
    try:
        L_cons, ss_cons_str, _cons_pairs = scs.project_ss_cons_to_sequence(
            seq_name, seqs, ss_tracks
        )
        if L_cons != L:
            raise ValueError(
                f"Consensus length {L_cons} != ungapped length {L} for {seq_name}"
            )
    except Exception as e:
        sys.stderr.write(f"[WARN] Could not project SS_cons for {seq_name}: {e}\n")
        ss_cons_str = "." * L
    # Candidate pairs for summary file (weights not super important here)
    try:
        cov = None
        if cfg.cov is not None:
            try:
                # Build alignment→sequence map for this sequence
                aligned_seq = seqs[seq_name]
                aln2seq, _L = scs.aln_to_seq_map(aligned_seq)

                cov = scs.load_cov_stats(str(cfg.cov), aln2seq)
                sys.stderr.write(
                    f"[PIPELINE] get_consensus_db: loaded cov stats for "
                    f"{len(cov)} pairs from {cfg.cov}.\n"
                )
            except Exception as e:
                sys.stderr.write(
                    f"[WARN] get_consensus_db: failed to read cov file {cfg.cov}: {e}\n"
                    "Proceeding without covariation weights in summary.\n"
                )
                cov = None

        _, candidate_pairs = scs.extract_candidate_pairs(
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

    except Exception as e:
        sys.stderr.write(f"[WARN] Could not extract candidate pairs: {e}\n")
        candidate_pairs = []

    # --- NEW: enforce minimum spacing between helices on the consensus ---
    ss_cons_clean = enforce_min_helix_spacing(ss_cons_str, min_dots_between_helices=2)
    if ss_cons_clean != ss_cons_str:
        sys.stderr.write(
            "[HelixSpacing] Adjusted consensus to enforce min 2 dots between stems:\n"
            f"  before: {ss_cons_str}\n"
            f"  after : {ss_cons_clean}\n"
        )
        ss_cons_str = ss_cons_clean

    # --- NEW: prune non-complementary (non-canonical) base pairs on consensus ---
    try:
        ss_cons_pruned, removed = fdr.prune_noncomplementary_pairs(
            ss_cons_str,
            ungapped_seq,
        )
    except Exception as e:
        sys.stderr.write(
            f"[WARN] get_consensus_db: failed to prune non-complementary pairs: {e}\n"
        )
    else:
        if removed > 0:
            sys.stderr.write(
                f"[PIPELINE] get_consensus_db: pruned {removed} non-complementary pairs from consensus.\n"
            )
            ss_cons_str = ss_cons_pruned
        else:
            sys.stderr.write(
                "[PIPELINE] get_consensus_db: no non-complementary pairs detected in consensus.\n"
            )

    # Optional summary file
    if cfg.io.summary is not None:
        cfg.io.summary.parent.mkdir(parents=True, exist_ok=True)
        scs.write_summary_file(
            summary_path=str(cfg.io.summary),
            seq_name=seq_name,
            ungapped_seq=ungapped_seq,
            L=L,
            ss_cons_str=ss_cons_str,
            candidate_pairs=candidate_pairs,
        )

    # Write consensus .db:
    #   line 1: ungapped sequence
    #   line 2: consensus structure
    cfg.io.consensus_db.parent.mkdir(parents=True, exist_ok=True)
    with cfg.io.consensus_db.open("w") as fh:
        fh.write(ungapped_seq + "\n")
        fh.write(ss_cons_str + "\n")

    sys.stderr.write(
        f"[PIPELINE] Wrote consensus DB with 1 structure to {cfg.io.consensus_db}\n"
    )
    return cfg.io.consensus_db

# ---------------------------------------------------------------------------
# Step 2: refine (using RNAstructure)
# ---------------------------------------------------------------------------

def _check_balanced_parentheses_local(struct: str) -> bool:
    stack: List[int] = []
    for ch in struct:
        if ch == "(":
            stack.append(1)
        elif ch == ")":
            if not stack:
                return False
            stack.pop()
    return not stack


def refine_db(cfg: PipelineConfig, db_in: Path, db_out: Path) -> Path:
    """
    Refine structures in `db_in` by folding long INTERNAL unpaired regions
    (Fold) and optionally terminal ends (DuplexFold).

    If DuplexFold returns multiple alternative duplexes for a given structure,
    we keep *all* of them and write all refined structures to `db_out`.
    """
    full_seq = rup.read_ungapped_seq_from_sto(cfg.sto, cfg.seq_name)
    structs = rup.read_db_structures(db_in)

    sys.stderr.write(
        f"[PIPELINE] Refining {len(structs)} structures from {db_in} → {db_out}\n"
    )

    fold_exe = cfg.fold_exe
    duplex_exe = cfg.duplex_exe
    min_unpaired = cfg.refine.min_unpaired
    min_terminal = cfg.refine.min_terminal_unpaired
    temperature = cfg.refine.temperature
    fold_extra = list(cfg.refine.fold_extra_args)
    duplex_extra = list(cfg.refine.duplex_extra_args)

    # Caches keyed by (sequence, temperature, extra_args ...)
    fold_cache: Dict[Tuple[str, float | None, Tuple[str, ...]], str] = {}
    duplex_cache: Dict[
        Tuple[str, str, float | None, Tuple[str, ...]],
        List[List[Tuple[int, int]]],
    ] = {}

    refined_structs: List[str] = []

    for idx, s in enumerate(structs):
        # How long are the unpaired 5'/3' ends for this structure?
        len_5, len_3 = rup.find_terminal_unpaired_ends(s, min_terminal)

        # 1) INTERNAL runs via Fold
        refined = rup.refine_structure(
            struct=s,
            full_seq=full_seq,
            fold_exe=fold_exe,
            min_unpaired_len=min_unpaired,
            temperature=temperature,
            extra_args=fold_extra,
            cache=fold_cache,
        )

        # 2) Terminal 5'/3' refinement with DuplexFold
        #    Only if we actually have sufficiently long ends.
        if (
            duplex_exe is not None
            and min_terminal > 0
            and len_5 > 0
            and len_3 > 0
        ):
            refined_list = rup.refine_terminal_ends_with_duplex(
                struct=refined,
                full_seq=full_seq,
                len_5=len_5,
                len_3=len_3,
                duplex_exe=duplex_exe,
                temperature=temperature,
                extra_args=duplex_extra,
                cache=duplex_cache,
            )
        else:
            # No DuplexFold refinement → just keep the Fold-refined structure
            refined_list = [refined]

        # Sanity check + collect all alternatives
        for alt_idx, r in enumerate(refined_list):
            if not _check_balanced_parentheses_local(r):
                sys.stderr.write(
                    f"[WARN] Unbalanced parentheses in refined structure {idx}"
                    + (f".{alt_idx}" if len(refined_list) > 1 else "")
                    + f":\n{r}\n"
                )
            refined_structs.append(r)

    db_out.parent.mkdir(parents=True, exist_ok=True)
    with db_out.open("w") as fh:
        fh.write(full_seq + "\n")
        for rs in refined_structs:
            fh.write(rs + "\n")

    sys.stderr.write(
        f"[PIPELINE] Wrote {len(refined_structs)} refined structures to {db_out}\n"
    )
    return db_out


# ---------------------------------------------------------------------------
# Step 3: sample_pk (Metropolis sampling over candidate pairs)
# ---------------------------------------------------------------------------

def _parse_parentheses_pairs(struct: str) -> Set[Tuple[int, int]]:
    """
    Parse only standard parentheses '()' from a dot-bracket string
    into a set of 0-based index pairs (i, j) with i < j.

    This is used to lock in the refined consensus base pairs; PK layers
    and exotic brackets are handled separately by the sampler.
    """
    stack: List[int] = []
    pairs: Set[Tuple[int, int]] = set()
    for idx, ch in enumerate(struct):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                # Unbalanced right paren; we just skip it.
                continue
            i = stack.pop()
            if i < idx:
                pairs.add((i, idx))
    return pairs

def _pairs_from_struct(struct: str) -> Set[Tuple[int, int]]:
    """
    Parse a dot-bracket string into (i, j) pairs, ignoring pair type.
    Uses the same bracket logic as filter_db_for_rosetta._parse_pairs().
    """
    pairs_with_type = fdr._parse_pairs(struct)  # type: ignore[attr-defined]
    return {(i, j) for (i, j, _op) in pairs_with_type}

def _combine_scaffold_and_pk(scaffold: str, pk_pairs: Set[Tuple[int, int]], L: int) -> str:
    """
    Merge a refined scaffold structure with a set of PK pairs and return a
    fully consistent PK-annotated string.

    Any scaffold base pair that shares an endpoint with a PK pair is dropped
    so that no residue participates in more than one pair. This guarantees
    that we never leave behind a "half-helix" like

        ..((......)))...
        ..(([.....)))..]

    where the original partner ')' should have been converted to '.'.
    """
    # All pairs present in the scaffold (standard + PK-like brackets)
    scaffold_pairs = _pairs_from_struct(scaffold)

    # Positions occupied by any sampled PK pair
    pk_positions: Set[int] = set()
    for i, j in pk_pairs:
        pk_positions.add(i)
        pk_positions.add(j)

    # Drop scaffold pairs that conflict with PK endpoints
    pruned_scaffold_pairs = {
        (i, j) for (i, j) in scaffold_pairs
        if i not in pk_positions and j not in pk_positions
    }

    combined_pairs = pruned_scaffold_pairs.union(pk_pairs)
    return scs.pairs_to_pk_string(sorted(combined_pairs), L)

def sample_pk(
    cfg: PipelineConfig,
    refined_db: Optional[Path],
    out_db: Path,
) -> Path:
    """
    Sample pseudoknotted secondary structures using the CaCoFold-derived
    candidate pairs and Metropolis sampler from sample_cacofold_structures.py.

    If `refined_db` is provided, we treat *each* refined structure in it as a
    separate consensus scaffold:

        final_pairs = base_pairs_from_scaffold ∪ sampled_pk_pairs

    and generate a full set of PK samples for each scaffold.
    """
    sto_path = cfg.sto
    seq_name = cfg.seq_name

    sys.stderr.write(f"[PIPELINE] Sampling PK structures for {seq_name}\n")

    # Parse alignment + SS_cons tracks
    seqs, ss_tracks = scs.parse_stockholm_single(str(sto_path))

    # Ungapped sequence for this RNA (header for .db)
    ungapped_seq = "".join(
        ch for ch in seqs[seq_name] if ch not in scs.GAP_CHARS
    )
    # Optional covariation file
    cov = None
    if cfg.cov is not None:
        try:
            aligned_seq = seqs[seq_name]
            aln2seq, _L = scs.aln_to_seq_map(aligned_seq)

            cov = scs.load_cov_stats(str(cfg.cov), aln2seq)
            sys.stderr.write(
                f"[PIPELINE] sample_pk: loaded cov stats for {len(cov)} pairs "
                f"from {cfg.cov}.\n"
            )
        except Exception as e:
            sys.stderr.write(
                f"[WARN] sample_pk: failed to read cov file {cfg.cov}: {e}\n"
                "Proceeding without covariation weights.\n"
            )
            cov = None


    # Candidate pairs from CaCoFold summary (+ cov weights if available)
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

    sys.stderr.write(
        f"[PIPELINE] sample_pk: extracted {len(candidate_pairs)} candidate pair(s) from CaCoFold.\n"
    )

    # Load refined scaffolds, if any
    refined_structs: List[str] = []
    if refined_db is not None and refined_db.is_file():
        refined_structs = rup.read_db_structures(refined_db)
        if not refined_structs:
            sys.stderr.write(
                "[PIPELINE] sample_pk: refined_db is present but empty; "
                "no base pairs will be locked.\n"
            )
        else:
            sys.stderr.write(
                f"[PIPELINE] sample_pk: using {len(refined_structs)} refined "
                "consensus structure(s) as scaffolds.\n"
            )
            # Sanity: all scaffolds must match CaCoFold length
            for idx, refined in enumerate(refined_structs):
                if len(refined) != L:
                    raise ValueError(
                        f"Refined consensus #{idx} length ({len(refined)}) "
                        f"!= CaCoFold length ({L})."
                    )
    else:
        sys.stderr.write(
            "[PIPELINE] sample_pk: no refined_db provided; sampling PKs without "
            "a locked consensus.\n"
        )
    print(refined_structs)

    out_db.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with out_db.open("w") as fh:
        fh.write(ungapped_seq + "\n")

        if refined_structs:
            # For each refined scaffold, treat its base pairs as a starting
            # consensus and overlay PK pairs sampled from CaCoFold.
            for s_idx, refined in enumerate(refined_structs):
                refined = refined.strip()
                if not refined:
                    continue

                scaffold_pairs = _pairs_from_struct(refined)

                sys.stderr.write(
                    "[PIPELINE] sample_pk: refined scaffold provided; "
                    "locking its base pairs and sampling PKs only in unpaired regions.\n"
                )

                # NEW: lock scaffold base pairs.
                # Do not allow candidate PK pairs that reuse any residue
                # already paired in the refined scaffold. This preserves the
                # refined helices exactly; PKs can only decorate previously
                # unpaired positions.
                scaffold_positions: Set[int] = set()
                for i, j in scaffold_pairs:
                    scaffold_positions.add(i)
                    scaffold_positions.add(j)

                filtered_candidates = [
                    (i, j, w)
                    for (i, j, w) in candidate_pairs
                    if i not in scaffold_positions and j not in scaffold_positions
                ]

                sys.stderr.write(
                    f"[PIPELINE] sample_pk: scaffold {s_idx}: "
                    f"{len(filtered_candidates)} candidate pair(s) available for PK sampling "
                    f"(out of {len(candidate_pairs)} total; "
                    f"{len(scaffold_positions)} positions locked by scaffold).\n"
                )

                # Metropolis sampling of PK pairs for this scaffold
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
                )

                for pk_set in pk_samples:
                    # Merge refined scaffold + sampled PK pairs safely
                    pk = _combine_scaffold_and_pk(refined, pk_set, L)
                    #pk_clean = enforce_min_helix_spacing(
                    #    pk, min_dots_between_helices=2
                    #)
                    fh.write(pk + "\n")
                    written += 1
    #    else:
    #        # No refined consensus: unconstrained PK sampling on full candidate pool
    #        filtered_candidates = candidate_pairs
    #        sys.stderr.write(
    #            f"[PIPELINE] sample_pk: no refined consensus provided; using "
    #            f"{len(filtered_candidates)} candidate pair(s) for unconstrained PK sampling.\n"
    #        )

    #        pk_samples = scs.sample_matchings(
    #            L=L,
    #            candidate_pairs=filtered_candidates,
    #            n_samples=cfg.sample.n_samples,
    #            burn_in=cfg.sample.burn_in,
    #            thin=cfg.sample.thin,
    #            min_loop_sep=cfg.sample.min_loop_sep,
    #            beta=cfg.sample.beta,
    #            seed=cfg.sample.seed,
    #            pk_alpha=cfg.sample.pk_alpha,
    #            pk_filter_frac=cfg.sample.pk_filter_frac,
    #            pk_filter_max_cross_per_pair=cfg.sample.pk_filter_max_cross_per_pair,
    #            pk_filter_max_total_cross=cfg.sample.pk_filter_max_total_cross,
    #        )

    #        for pk_set in pk_samples:
    #            pk = scs.pairs_to_pk_string(sorted(pk_set), L)
    #            pk_clean = enforce_min_helix_spacing(pk, min_dots_between_helices=2)
    #            fh.write(pk_clean + "\n")
    #            written += 1

    sys.stderr.write(
        f"[PIPELINE] sample_pk: wrote {written} sampled structure(s) to {out_db}\n"
    )
    return out_db


# ---------------------------------------------------------------------------
# Step 4: ensure_valid (heuristic plausibility filter)
# ---------------------------------------------------------------------------
def enforce_min_helix_spacing(struct: str, min_dots_between_helices: int = 2) -> str:
    """
    Enforce two constraints on the standard parentheses secondary structure:

    1. Remove *isolated*, top-level one-bp helices "()" (no content inside).
       These are patterns like ....().... at depth 1.

    2. Enforce that between any closing ')' and the next opening '(' there are
       at least `min_dots_between_helices` '.' characters. If a helix starts
       too close to the previous one, we strip that *second* helix by turning
       all of its '(' and ')' into '.'.

    This only touches '()' pairs; PK brackets ([ ], { }, < >, letters) are
    left alone.
    """
    chars = list(struct)
    n = len(chars)

    # ------------------------------------------------------------------
    # Pass 1: remove isolated top-level "()" helices
    # ------------------------------------------------------------------
    stack: List[Tuple[int, int]] = []
    depth = 0
    to_strip_indices: set[int] = set()

    for idx, ch in enumerate(chars):
        if ch == "(":
            depth += 1
            # store the index and depth at which this '(' opened
            stack.append((idx, depth))
        elif ch == ")":
            if not stack:
                # unbalanced, just decrease depth defensively
                depth = max(depth - 1, 0)
                continue
            open_idx, open_depth = stack.pop()
            # close this level
            depth -= 1

            # Top-level single pair: "()" with no content inside
            if open_depth == 1 and idx == open_idx + 1:
                to_strip_indices.add(open_idx)
                to_strip_indices.add(idx)

    if to_strip_indices:
        sys.stderr.write(
            "[HelixSpacing] Removing isolated top-level '()' helices at positions: "
            + ", ".join(map(str, sorted(to_strip_indices))) + "\n"
        )
        for i in to_strip_indices:
            if 0 <= i < n and chars[i] in "()":
                chars[i] = "."

    # ------------------------------------------------------------------
    # Pass 2: enforce minimum dots between successive helices
    # ------------------------------------------------------------------
    i = 0
    while i < n:
        if chars[i] == ")":
            # Look ahead for the next '(' and count dots between them
            j = i + 1
            dot_count = 0
            while j < n and chars[j] != "(":
                if chars[j] == ".":
                    dot_count += 1
                j += 1

            # If we found a '(' and the dot run is too short, strip the helix
            if j < n and dot_count < min_dots_between_helices:
                # Find matching ')' for the '(' at position j
                stack_depth = 1
                k = j + 1
                while k < n and stack_depth > 0:
                    if chars[k] == "(":
                        stack_depth += 1
                    elif chars[k] == ")":
                        stack_depth -= 1
                    k += 1

                # Replace that helix's parentheses with dots
                sys.stderr.write(
                    f"[HelixSpacing] Stripped helix starting at {j} due to "
                    f"insufficient dots between stems (found {dot_count}, "
                    f"required {min_dots_between_helices}).\n"
                )
                for t in range(j, min(k, n)):
                    if chars[t] in "()":
                        chars[t] = "."

                # Continue scanning from end of stripped helix
                i = k
                continue

        i += 1

    return "".join(chars)



def _count_crossings(pairs: List[Tuple[int, int]]) -> Tuple[int, int]:
    """
    Return (total_crossings, max_cross_per_pair) for a list of (i, j) pairs.
    """
    total = 0
    max_per = 0
    for idx, (a, b) in enumerate(pairs):
        if a > b:
            a, b = b, a
        crossings_for_this = 0
        for (c, d) in pairs[idx + 1 :]:
            if c > d:
                c, d = d, c
            if (a < c < b < d) or (c < a < d < b):
                total += 1
                crossings_for_this += 1
        if crossings_for_this > max_per:
            max_per = crossings_for_this
    return total, max_per


def _pk_fraction(struct: str, pairs: List[Tuple[int, int]]) -> float:
    """
    Very rough measure: fraction of positions that participate in *any* pair,
    used as a proxy for how "densely" pseudoknotted things are.
    """
    pos = set()
    for i, j in pairs:
        pos.add(i)
        pos.add(j)
    return len(pos) / max(len(struct), 1)


def _is_structure_plausible(struct: str, cfg: EnsureValidConfig) -> bool:
    if not struct:
        return False

    # Only allow reasonable characters (., brackets, letters)
    allowed = set(".()[]{}<>") | set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    if any(ch not in allowed for ch in struct):
        return False

    try:
        pairs_with_type = fdr._parse_pairs(struct)  # type: ignore[attr-defined]
    except Exception:
        return False

    pairs = [(i, j) for (i, j, _op) in pairs_with_type]
    if not pairs:
        # allow completely unpaired structures
        return True

    total_cross, max_cross = _count_crossings(pairs)
    if total_cross > cfg.max_total_crossings:
        return False
    if max_cross > cfg.max_crossings_per_pair:
        return False

    if _pk_fraction(struct, pairs) > cfg.max_pk_fraction:
        return False

    return True


def ensure_valid(cfg: PipelineConfig, db_in: Path, db_out: Path) -> Path:
    """
    Filter the .db ensemble in `db_in` to keep only plausible structures based
    on coarse pseudoknot / crossing heuristics.
    """
    structs = fdr.read_db_structures(db_in)
    sys.stderr.write(
        f"[PIPELINE] ensure_valid: filtering {len(structs)} structures from {db_in}\n"
    )

    kept = [s for s in structs if _is_structure_plausible(s, cfg.ensure_valid)]

    if not kept:
        sys.stderr.write(
            "[PIPELINE] ensure_valid rejected all structures; "
            "falling back to original ensemble.\n"
        )
        kept = structs

    
    db_out.parent.mkdir(parents=True, exist_ok=True)

    # Sequence from Stockholm (for header)
    seq = fdr.read_ungapped_seq_from_sto(cfg.sto, cfg.seq_name)
    if len(seq) != len(kept[0]):
        raise ValueError(
            f"Sequence length from {cfg.sto} ({len(seq)}) does not match "
            f"structure length ({len(kept[0])}) in ensure_valid."
        )

    fdr.write_db_structures(db_out, kept, seq=seq)
    sys.stderr.write(
        f"[PIPELINE] ensure_valid kept {len(kept)} structures → {db_out}\n"
    )

    return db_out


# ---------------------------------------------------------------------------
# Step 5: make_rosetta_ready
# ---------------------------------------------------------------------------

def make_rosetta_ready(cfg: PipelineConfig, db_in: Path, db_out: Path) -> Path:
    """
    Deduplicate / cluster structures, convert to Rosetta notation, and
    prune non-complementary base pairs using the sequence from the
    Stockholm file, so Rosetta only sees canonical/wobble pairs.
    """
    structs = fdr.read_db_structures(db_in)
    sys.stderr.write(
        f"[PIPELINE] Rosetta-ready: starting from {len(structs)} structures\n"
    )

    # 1) Exact deduplication
    deduped = fdr.deduplicate(structs)
    sys.stderr.write(
        f"[PIPELINE] After exact deduplication: {len(deduped)} structures\n"
    )

    if not deduped:
        raise ValueError("No structures remain after deduplication")

    # 2) Optional Hamming-based merging
    thr = cfg.rosetta.hamming_threshold
    if thr < 0:
        thr = 0
    if thr > 0:
        L = len(deduped[0])
        thr = max(thr, int(round(cfg.rosetta.hamming_frac * L)))
        merged = fdr.merge_by_hamming(deduped, thr)  # type: ignore[attr-defined]
        sys.stderr.write(
            f"[PIPELINE] After Hamming merge (threshold={thr}): {len(merged)} structures\n"
        )
    else:
        merged = deduped

    # 3) Optional cap on number of structures
    if cfg.rosetta.max_structures is not None and cfg.rosetta.max_structures > 0:
        merged = merged[: cfg.rosetta.max_structures]
        sys.stderr.write(
            f"[PIPELINE] After max_structures cap: {len(merged)} structures\n"
        )

    # 4) Convert to Rosetta notation, skipping invalid structures
    converted: List[str] = []
    num_skipped = 0
    for idx, s in enumerate(merged):
        try:
            converted.append(fdr.convert_to_rosetta_notation(s))
        except Exception as e:
            sys.stderr.write(
                f"[WARN] Skipping structure {idx} during Rosetta conversion: {e}\n"
            )
            num_skipped += 1

    sys.stderr.write(
        f"[PIPELINE] Converted {len(converted)} structures to Rosetta; "
        f"skipped {num_skipped}.\n"
    )

    if not converted:
        raise ValueError("No structures remained after Rosetta conversion")

    # 5) Prune non-complementary pairs so Rosetta only sees canonical/wobble
    seq = fdr.read_ungapped_seq_from_sto(cfg.sto, cfg.seq_name)
    L_seq = len(seq)
    if len(converted[0]) != L_seq:
        raise ValueError(
            f"Sequence length from {cfg.sto} ({L_seq}) does not match "
            f"structure length ({len(converted[0])})."
        )

    pruned: List[str] = []
    total_removed = 0
    for s in converted:
        s_new, removed = fdr.prune_noncomplementary_pairs(s, seq)
        pruned.append(s_new)
        total_removed += removed

    if total_removed > 0:
        sys.stderr.write(
            f"[PIPELINE] Pruned {total_removed} non-complementary pairs across ensemble.\n"
        )
    else:
        sys.stderr.write("[PIPELINE] No non-complementary pairs detected.\n")

    db_out.parent.mkdir(parents=True, exist_ok=True)
    fdr.write_db_structures(db_out, pruned, seq=seq)
    sys.stderr.write(
        f"[PIPELINE] Wrote {len(pruned)} Rosetta-compatible structures to {db_out}\n"
    )

    return db_out


# ---------------------------------------------------------------------------
# High-level orchestrator
# ---------------------------------------------------------------------------

def run_pipeline(cfg: PipelineConfig, mode: str = "refine_first") -> Path:
    """
    Run one of the pre-wired pipelines.

    Modes
    -----
    legacy:
        get_consensus_db → sample_pk → refine → rosetta_ready

    refine_first:
        get_consensus_db → refine → sample_pk → ensure_valid → rosetta_ready
    """
    mode = mode.lower()
    if mode not in {"legacy", "refine_first"}:
        raise ValueError(f"Unknown mode: {mode!r}")

    # Step 0: consensus (always)
    consensus_db = get_consensus_db(cfg)

    if mode == "legacy":
        sampled_db = sample_pk(cfg, refined_db=None, out_db=cfg.io.sampled_db)
        refined_db = refine_db(cfg, db_in=sampled_db, db_out=cfg.io.refined_db)
        return make_rosetta_ready(cfg, db_in=refined_db, db_out=cfg.io.rosetta_db)

    # refine_first pipeline
    refined_consensus_db = refine_db(cfg, db_in=consensus_db, db_out=cfg.io.refined_db)
    #ensure_valid(cfg, refined_consensus_db, cfg.io.refined_db)
    sampled_db = sample_pk(cfg, refined_db=refined_consensus_db, out_db=cfg.io.sampled_db)
    #ensure_valid(cfg, sampled_db, cfg.io.sampled_db)
    valid_db = ensure_valid(cfg, db_in=sampled_db, db_out=cfg.io.valid_db)
    ensure_valid(cfg, valid_db, cfg.io.valid_db)
    return make_rosetta_ready(cfg, db_in=valid_db, db_out=cfg.io.rosetta_db)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Pre-wired CaCoFold → sampling → refinement → Rosetta pipelines.",
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="Path to a JSON or YAML config file describing inputs / parameters.",
    )
    parser.add_argument(
        "--mode",
        choices=["legacy", "refine_first"],
        default="refine_first",
        help="Which pipeline wiring to run.",
    )
    args = parser.parse_args(argv)

    cfg_path = Path(args.config).resolve()
    cfg = load_pipeline_config(cfg_path)

    final_db = run_pipeline(cfg, mode=args.mode)
    sys.stderr.write(f"[PIPELINE] Final Rosetta-ready DB: {final_db}\n")


if __name__ == "__main__":
    main()

