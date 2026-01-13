#!/usr/bin/env python3
"""
Sample alternative pseudoknotted secondary structures from a CaCoFold/R-scape
Stockholm (.sto) file, for a chosen sequence, and write them to a .db file.

Output .db format:
    line 1: ungapped RNA sequence
    line 2+: one PK-annotated dot-bracket structure per line
"""

import argparse
import json
import math
import multiprocessing as mp
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

from .energy import EnergyParams
from .moves import MoveType
from .sampler import SamplerConfig, sample_matchings_new
from .segments import build_candidate_segments

# Alignment gaps
GAP_CHARS = set("-._~")

# WUSS bracket pairs + Extended PK support (a-z)
# Matching the hierarchy used in refine_unpaired_regions
OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
# Add a-z -> A-Z
for i in range(26):
    c_lower = chr(ord("a") + i)
    c_upper = chr(ord("A") + i)
    OPEN_TO_CLOSE[c_lower] = c_upper

CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}

# Canonical base pairs (including GU wobble)
CANONICAL_BASE_PAIRS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),  # wobble
}


def is_canonical_pair(b1: str, b2: str) -> bool:
    """Return True if (b1, b2) is a canonical Watson–Crick or GU wobble pair."""
    b1u = b1.upper()
    b2u = b2.upper()
    return (b1u, b2u) in CANONICAL_BASE_PAIRS


# ------------------ Parsing Stockholm + WUSS ------------------ #


def parse_stockholm_single(path: str) -> tuple[dict[str, str], dict[str, str]]:
    """
    Parse a SINGLE Stockholm alignment.

    Returns
    -------
    sequences : dict
        name -> aligned sequence (with gaps)
    ss_tracks : dict
        tag -> full-length SS_cons-like track string (aligned columns)
    """
    sequences: dict[str, list[str]] = {}
    ss_tracks: dict[str, list[str]] = {}

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not line:
                continue
            if line.startswith("# STOCKHOLM"):
                continue
            if line.startswith("//"):
                break

            # Per-column markup (WUSS etc.)
            if line.startswith("#=GC "):
                parts = line.split(maxsplit=2)
                if len(parts) < 3:
                    continue
                tag, s = parts[1], parts[2]
                ss_tracks.setdefault(tag, []).append(s)
                continue

            # Ignore other comment lines
            if line.startswith("#"):
                continue

            # Sequence line: "name  aligned-sequence"
            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                continue
            name, s = parts[0], parts[1].strip()
            sequences.setdefault(name, []).append(s)

    seqs_joined = {name: "".join(chunks) for name, chunks in sequences.items()}
    ss_joined = {tag: "".join(chunks) for tag, chunks in ss_tracks.items()}

    lengths = {len(v) for v in seqs_joined.values()}
    lengths.update(len(v) for v in ss_joined.values())
    if len(lengths) != 1:
        raise ValueError(f"Alignment length mismatch in {path}: lengths={lengths}")

    return seqs_joined, ss_joined


def ss_tag_order(tag: str) -> tuple[int, int, str]:
    """
    For ordering SS_cons, SS_cons_1, SS_cons_2, ...
    """
    if tag == "SS_cons":
        return (0, 0, tag)
    if tag.startswith("SS_cons_"):
        try:
            n = int(tag.split("_")[-1])
        except ValueError:
            n = 9999
        return (1, n, tag)
    return (2, 9999, tag)


def pairs_from_track(track: str) -> list[tuple[int, int]]:
    """
    Given a WUSS-like track string (for all alignment columns),
    return a list of (i, j) pairs in alignment coordinates (0-based).

    Uses bracket types: <>, (), [], {}, and a-z.
    """
    stacks = {op: [] for op in OPEN_TO_CLOSE}
    result: list[tuple[int, int]] = []

    for i, ch in enumerate(track):
        if ch in OPEN_TO_CLOSE:
            stacks[ch].append(i)
        elif ch in CLOSE_TO_OPEN:
            op = CLOSE_TO_OPEN[ch]
            if stacks[op]:
                j = stacks[op].pop()
                result.append((j, i))

    result.sort()
    return result


def aln_to_seq_map(aligned_seq: str) -> tuple[dict[int, int], int]:
    """
    Build a map alignment_index -> sequence_index (ungapped)
    for one sequence. Returns (aln2seq, L), where:
      - aln2seq[i] = seq_index or -1 if gap
      - L = length of ungapped sequence
    """
    aln2seq: dict[int, int] = {}
    pos = 0
    for i, ch in enumerate(aligned_seq):
        if ch in GAP_CHARS:
            aln2seq[i] = -1
        else:
            aln2seq[i] = pos
            pos += 1
    return aln2seq, pos


def map_pairs_to_seq(
    pairs_aln: list[tuple[int, int]],
    aln2seq: dict[int, int],
) -> list[tuple[int, int]]:
    """
    Map pairs from alignment coordinates to sequence coordinates,
    dropping any that hit a gap in this sequence.
    """
    result: list[tuple[int, int]] = []
    for i, j in pairs_aln:
        si = aln2seq.get(i, -1)
        sj = aln2seq.get(j, -1)
        if si < 0 or sj < 0:
            continue
        if si > sj:
            si, sj = sj, si
        result.append((si, sj))
    result.sort()
    return result


# ------------------ Covariation stats + weights ------------------ #


@dataclass
class CovRecord:
    in_cacofold: bool
    score: float
    evalue: float
    pvalue: float
    substitutions: int
    power: float


def load_cov_stats(
    cov_path: str,
    aln2seq: dict[int, int],
) -> dict[tuple[int, int], CovRecord]:
    """
    Load a CaCoFold/R-scape .cov file and map pairs into sequence coordinates.

    Returns
    -------
    cov : dict
        (i, j) -> CovRecord
        in ungapped sequence coordinates for this sequence.
    """
    cov: dict[tuple[int, int], CovRecord] = {}

    with open(cov_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("-"):
                continue

            parts = line.split()
            # There are 2 layouts:
            # 9 fields:  in_CaCoFold in_given left right score E p subs power
            # 8 fields:  in_CaCoFold          left right score E p subs power
            if len(parts) == 9:
                in_cacofold_flag, _in_given = parts[0], parts[1]
                left_pos, right_pos = int(parts[2]), int(parts[3])
                score = float(parts[4])
                evalue = float(parts[5])
                pvalue = float(parts[6])
                subs = int(parts[7])
                power = float(parts[8])
            elif len(parts) == 8:
                in_cacofold_flag = parts[0]
                left_pos, right_pos = int(parts[1]), int(parts[2])
                score = float(parts[3])
                evalue = float(parts[4])
                pvalue = float(parts[5])
                subs = int(parts[6])
                power = float(parts[7])
            else:
                # unexpected format
                continue

            i_aln = left_pos - 1
            j_aln = right_pos - 1
            si = aln2seq.get(i_aln, -1)
            sj = aln2seq.get(j_aln, -1)
            if si < 0 or sj < 0:
                continue
            if si > sj:
                si, sj = sj, si

            key = (si, sj)
            rec = CovRecord(
                in_cacofold=(in_cacofold_flag == "*"),
                score=score,
                evalue=evalue,
                pvalue=pvalue,
                substitutions=subs,
                power=power,
            )

            # If duplicates show up, keep the strongest (smallest E-value)
            if key in cov:
                if rec.evalue < cov[key].evalue:
                    cov[key] = rec
            else:
                cov[key] = rec

    return cov


def cov_weight(rec: CovRecord, mode: str = "logE_power") -> float:
    """
    Convert a CovRecord into a non-negative scalar weight.
    """
    if mode == "logE":
        # Only positive evidence: E < 1 → positive; otherwise 0
        if rec.evalue <= 0:
            return 0.0
        return max(0.0, -math.log10(rec.evalue))

    if mode == "logE_power":
        if rec.evalue <= 0:
            return 0.0
        return max(0.0, -math.log10(rec.evalue)) * rec.power

    if mode == "score":
        return max(0.0, rec.score)

    if mode == "indicator":
        # Just "has covariation" vs not
        return 1.0

    return 0.0


# ------------------ Candidate pairs + weights ------------------ #


def extract_candidate_pair_components(
    seq_name: str,
    seqs: dict[str, str],
    ss_tracks: dict[str, str],
    w_core: float = 2.0,
    w_alt: float = 1.0,
    cov: dict[tuple[int, int], CovRecord] | None = None,
    cov_mode: str = "off",
    cov_alpha: float = 1.0,
    cov_min_power: float = 0.0,
    cov_forbid_negative: bool = False,
    cov_negative_E: float = 1.0,
) -> tuple[int, str, dict[str, dict[tuple[int, int], float]]]:
    if seq_name not in seqs:
        raise KeyError(f"Sequence '{seq_name}' not in alignment.")

    aligned_seq = seqs[seq_name]
    aln2seq, L = aln_to_seq_map(aligned_seq)
    ungapped_seq = "".join(ch for ch in aligned_seq if ch not in GAP_CHARS)
    tags = [t for t in ss_tracks if t.startswith("SS_cons")]

    if not tags:
        raise ValueError("No SS_cons* tracks found in the Stockholm file.")

    tags = sorted(tags, key=ss_tag_order)

    core: dict[tuple[int, int], float] = {}
    alt: dict[tuple[int, int], float] = {}
    cov_comp: dict[tuple[int, int], float] = {}
    dropped_noncanonical = 0

    for tag in tags:
        track = ss_tracks[tag]
        pairs_aln = pairs_from_track(track)
        pairs_seq = map_pairs_to_seq(pairs_aln, aln2seq)
        for i, j in pairs_seq:
            if i > j:
                i, j = j, i

            # Enforce canonical WC/GU pairing at the *candidate* level
            b_i = ungapped_seq[i]
            b_j = ungapped_seq[j]
            if not is_canonical_pair(b_i, b_j):
                dropped_noncanonical += 1
                continue

            key = (i, j)

            # Base layer-based weight
            layer_w = w_core if tag == "SS_cons" else w_alt
            if tag == "SS_cons":
                core[key] = core.get(key, 0.0) + layer_w
            else:
                alt[key] = alt.get(key, 0.0) + layer_w

            # Optional covariation component
            rec: CovRecord | None = None
            if cov is not None and cov_mode != "off":
                rec = cov.get(key)

                if rec is not None:
                    # Optionally filter on power / negative evidence
                    if rec.power < cov_min_power:
                        if cov_forbid_negative:
                            continue

                    # Negative evidence: high power but weak covariation
                    if (
                        cov_forbid_negative
                        and rec.power >= cov_min_power
                        and rec.evalue >= cov_negative_E
                    ):
                        continue

                    cov_w = cov_weight(rec, mode=cov_mode)
                    cov_comp[key] = cov_comp.get(key, 0.0) + cov_w

    components = {"core": core, "alt": alt, "cov": cov_comp}
    if dropped_noncanonical > 0:
        sys.stderr.write(
            f"[CaCoFoldSample] Dropped {dropped_noncanonical} "
            f"non-canonical candidate pair(s) for {seq_name} based on sequence.\n"
        )
    return L, ungapped_seq, components


def extract_candidate_pairs(
    seq_name: str,
    seqs: dict[str, str],
    ss_tracks: dict[str, str],
    w_core: float = 2.0,
    w_alt: float = 1.0,
    cov: dict[tuple[int, int], CovRecord] | None = None,
    cov_mode: str = "off",
    cov_alpha: float = 1.0,
    cov_min_power: float = 0.0,
    cov_forbid_negative: bool = False,
    cov_negative_E: float = 1.0,
) -> tuple[int, list[tuple[int, int, float]]]:
    """
    From SS_cons, SS_cons_1, SS_cons_2, ... build:
      - L: length of ungapped sequence
      - candidate_pairs: list of (i, j, weight)

    Base layer weight heuristic:
      For each layer where a pair (i,j) appears:
        - if layer == SS_cons: weight += w_core
        - else: weight += w_alt

    If cov and cov_mode != "off", augment layer weights with covariation-based
    weights:

      layer_w = w_core or w_alt
      cov_w   = cov_weight(rec, mode=cov_mode)

      combined:
        w = layer_w + cov_alpha * cov_w

    and optionally:
      - filter on power (< cov_min_power)
      - drop pairs with strong negative evidence
        (power >= cov_min_power and evalue >= cov_negative_E) if
        cov_forbid_negative is True.
    """
    L, _seq, components = extract_candidate_pair_components(
        seq_name=seq_name,
        seqs=seqs,
        ss_tracks=ss_tracks,
        w_core=w_core,
        w_alt=w_alt,
        cov=cov,
        cov_mode=cov_mode,
        cov_alpha=cov_alpha,
        cov_min_power=cov_min_power,
        cov_forbid_negative=cov_forbid_negative,
        cov_negative_E=cov_negative_E,
    )

    core = components.get("core", {})
    alt = components.get("alt", {})
    cov_comp = components.get("cov", {})

    pair_scores: dict[tuple[int, int], float] = {}
    for key in set(core) | set(alt) | set(cov_comp):
        pair_scores[key] = (
            core.get(key, 0.0) + alt.get(key, 0.0) + cov_alpha * cov_comp.get(key, 0.0)
        )

    candidate_pairs = [(i, j, w) for (i, j), w in sorted(pair_scores.items())]
    return L, candidate_pairs


# ------------------ Consensus projection + PK summary ----------- #


def project_ss_cons_to_sequence(
    seq_name: str,
    seqs: dict[str, str],
    ss_tracks: dict[str, str],
) -> tuple[int, str, list[tuple[int, int]]]:
    """
    Project the primary SS_cons track onto the ungapped sequence.

    Returns
    -------
    L : int
        Ungapped sequence length
    consensus_str : str
        Dot-bracket string using only '(' and ')'
    pairs_seq : list of (i, j)
        Base pairs in sequence coordinates
    """
    if seq_name not in seqs:
        raise KeyError(f"Sequence '{seq_name}' not in alignment.")

    aligned_seq = seqs[seq_name]
    aln2seq, L = aln_to_seq_map(aligned_seq)

    if "SS_cons" not in ss_tracks:
        # Fallback: no primary SS_cons found
        return L, "." * L, []

    track = ss_tracks["SS_cons"]
    pairs_aln = pairs_from_track(track)
    pairs_seq = map_pairs_to_seq(pairs_aln, aln2seq)

    chars = ["."] * L
    for i, j in pairs_seq:
        if i > j:
            i, j = j, i
        chars[i] = "("
        chars[j] = ")"

    return L, "".join(chars), pairs_seq


# ------------------ Pseudoknot layering + PK string -------------- #


def pairs_cross(p: tuple[int, int], q: tuple[int, int]) -> bool:
    """
    Return True if arcs (i, j) and (k, l) cross (pseudoknot)
    in the standard 1D arc diagram sense.
    """
    i, j = p
    k, l = q
    if i > j:
        i, j = j, i
    if k > l:
        k, l = l, k
    return (i < k < j < l) or (k < i < l < j)


def pairs_to_layers(pairs: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """
    Given a list of pairs (i, j) (sequence coords, i<j),
    partition them into layers so that within each layer
    there are no crossing arcs.
    """
    layers: list[list[tuple[int, int]]] = []

    for p in sorted(pairs):
        placed = False
        for layer in layers:
            if any(pairs_cross(p, q) for q in layer):
                continue
            layer.append(p)
            placed = True
            break
        if not placed:
            layers.append([p])

    return layers


def compute_pk_stats_for_set(
    pairs: set[tuple[int, int]],
) -> tuple[set[tuple[int, int]], dict[tuple[int, int], int], int, int]:
    """
    Given a set of pairs, compute:

      - pk_pairs: pairs involved in >= 1 crossing
      - cross_counts: pair -> # of crossings it participates in
      - total_crossings: total number of pair-pair crossings
      - max_cross: max #crossings over all pairs
    """
    pairs_list = sorted(pairs)
    cross_counts: dict[tuple[int, int], int] = dict.fromkeys(pairs_list, 0)
    n = len(pairs_list)

    for i in range(n):
        p = pairs_list[i]
        for j in range(i + 1, n):
            q = pairs_list[j]
            if pairs_cross(p, q):
                cross_counts[p] += 1
                cross_counts[q] += 1

    pk_pairs = {p for p, c in cross_counts.items() if c > 0}
    total_crossings = sum(cross_counts.values()) // 2
    max_cross = max(cross_counts.values()) if cross_counts else 0

    return pk_pairs, cross_counts, total_crossings, max_cross


def prune_pseudoknots(
    pair_set: set[tuple[int, int]],
    pair_weights: dict[tuple[int, int], float],
    L: int,
    pk_filter_frac: float,
    pk_filter_max_cross_per_pair: int,
    pk_filter_max_total_cross: int,
    fixed_pairs: set[tuple[int, int]] = None,
) -> set[tuple[int, int]]:
    """
    Given a matching (pair_set), greedily remove the "worst" pseudoknotted
    pairs until it satisfies the PK thresholds:

      len(pk_pairs) <= pk_filter_frac * L          (if pk_filter_frac > 0)
      max_crossings_per_pair <= pk_filter_max_cross_per_pair (if >= 0)
      total_crossings <= pk_filter_max_total_cross (if >= 0)

    "Worst" = highest crossing count, tie-broken by lowest weight.
    Returns a NEW set; does not modify the input set.

    If fixed_pairs is provided, these pairs are NEVER removed, even if they
    cause violations.
    """
    if (
        (pk_filter_frac is None or pk_filter_frac <= 0)
        and (pk_filter_max_cross_per_pair is None or pk_filter_max_cross_per_pair < 0)
        and (pk_filter_max_total_cross is None or pk_filter_max_total_cross < 0)
    ):
        return set(pair_set)

    current = set(pair_set)
    if fixed_pairs is None:
        fixed_pairs = set()

    def violates(pk_pairs, max_cross, total_crossings) -> bool:
        if pk_filter_frac is not None and pk_filter_frac > 0:
            if len(pk_pairs) > pk_filter_frac * float(L):
                return True
        if pk_filter_max_cross_per_pair is not None and pk_filter_max_cross_per_pair >= 0:
            if max_cross > pk_filter_max_cross_per_pair:
                return True
        if pk_filter_max_total_cross is not None and pk_filter_max_total_cross >= 0:
            if total_crossings > pk_filter_max_total_cross:
                return True
        return False

    while current:
        pk_pairs, cross_counts, total_crossings, max_cross = compute_pk_stats_for_set(current)

        # If we now satisfy thresholds, stop pruning
        if not violates(pk_pairs, max_cross, total_crossings):
            break

        # If thresholds still violated but there are no pk_pairs, bail out
        if not pk_pairs:
            break

        # Filter out fixed pairs from consideration for removal
        removable_pk_pairs = [p for p in pk_pairs if p not in fixed_pairs]

        if not removable_pk_pairs:
            # We can't remove anything else to satisfy constraints because all
            # problematic pairs are fixed. We must accept the violation.
            break

        # Remove the "worst" pair: most crossings, then lowest weight
        worst_p = None
        worst_score = None  # higher is "worse"
        for p in removable_pk_pairs:
            c = cross_counts.get(p, 0)
            w = pair_weights.get(p, 0.0)
            score = (c, -w)  # more crossings, lower weight → more likely to remove
            if worst_score is None or score > worst_score:
                worst_score = score
                worst_p = p

        if worst_p is None:
            break
        current.remove(worst_p)

    return current


def pairs_to_pk_string(pairs: list[tuple[int, int]], L: int) -> str:
    """
    Convert a set/list of pairs (i, j) into a pseudoknot-annotated string.

    Uses Rosetta-safe hierarchy:
    - Layer 0 -> '(' and ')'
    - Layer 1 -> '[' and ']'
    - Layer 2 -> '{' and '}'
    - Layers 3+ -> 'a/A', 'b/B', ...

    Assumes 'pairs' is a matching: no residue appears in more than one pair.
    """
    # Sanity check: Ensure pairs form a valid matching (no shared indices)
    indices = [idx for pair in pairs for idx in pair]
    if len(indices) != len(set(indices)):
        raise ValueError("Invalid matching: overlapping pairs detected in structure.")

    chars = ["."] * L
    layers = pairs_to_layers(pairs)

    rosetta_brackets = [("(", ")"), ("[", "]"), ("{", "}")]
    # Add a-z for deeper layers
    for i in range(26):
        rosetta_brackets.append((chr(ord("a") + i), chr(ord("A") + i)))

    for layer_idx, layer in enumerate(layers):
        if layer_idx < len(rosetta_brackets):
            open_char, close_char = rosetta_brackets[layer_idx]
        else:
            # Fallback
            open_char, close_char = "{", "}"

        for i, j in layer:
            chars[i] = open_char
            chars[j] = close_char

    return "".join(chars)


def write_summary_file(
    summary_path: str,
    seq_name: str,
    ungapped_seq: str,
    L: int,
    ss_cons_str: str,
    candidate_pairs: list[tuple[int, int, float]],
) -> None:
    """
    Write a human-readable summary of:

      - Ungapped sequence
      - Projected consensus secondary structure (SS_cons)
      - PK layers derived from all SS_cons* tracks
    """
    # Collect unique (i, j) pairs from candidate_pairs
    unique_pairs = sorted({(i, j) if i < j else (j, i) for (i, j, _w) in candidate_pairs})

    # Combined PK string
    # Candidates may conflict (overlap), which is fine for a pool but not for a single structure.
    # We catch the validation error to prevent pipeline crashes.
    try:
        pk_all = pairs_to_pk_string(unique_pairs, L)
    except ValueError:
        pk_all = (
            "(Cannot represent candidates as a single string due to overlapping/conflicting pairs)"
        )

    # Per-layer strings
    layers = pairs_to_layers(unique_pairs)

    def layer_string(layer_idx: int, layer_pairs: list[tuple[int, int]]) -> str:
        chars = ["."] * L

        # Check for overlaps within this layer and PRUNE them for display
        # instead of erroring out.
        clean_pairs = []
        occupied = set()

        # Sort pairs (e.g. by first index) to be deterministic
        for p in sorted(layer_pairs):
            i, j = p
            if i > j:
                i, j = j, i

            if i not in occupied and j not in occupied:
                occupied.add(i)
                occupied.add(j)
                clean_pairs.append((i, j))

        op, cl = "(", ")"
        for i, j in clean_pairs:
            chars[i] = op
            chars[j] = cl

        return "".join(chars)

    with open(summary_path, "w") as fh:
        fh.write(f"# CaCoFold summary for {seq_name}\n")
        fh.write(f"# Length (ungapped): {L}\n\n")

        fh.write(">sequence\n")
        fh.write(ungapped_seq + "\n\n")

        fh.write(">consensus_projected\n")
        fh.write(ss_cons_str + "\n\n")

        fh.write(">pk_all_layers\n")
        fh.write(pk_all + "\n\n")

        for idx, layer in enumerate(layers):
            fh.write(f">pk_layer_{idx}\n")
            fh.write(layer_string(idx, layer) + "\n")


# ------------------ Metropolis sampler on matchings -------------- #


def sample_matchings(
    L: int,
    candidate_pairs: list[tuple[int, int, float]],
    n_samples: int,
    burn_in: int = 1000,
    thin: int = 10,
    min_loop_sep: int = 1,  # minimum unpaired residues inside a hairpin
    beta: float = 1.0,
    seed: int | None = None,
    pk_alpha: float = 2.0,
    pk_filter_frac: float = 0.2,
    pk_filter_max_cross_per_pair: int = 2,
    pk_filter_max_total_cross: int = 30,
    pk_depth_limit: int | None = None,
    scaffold_pairs: set[tuple[int, int]] | None = None,
) -> list[set[tuple[int, int]]]:
    """
    Metropolis sampler over matchings (sets of non-overlapping base pairs).

    scaffold_pairs: A set of pairs that MUST be present in every sample.
                    These pairs are fixed and cannot be removed.
                    Sampling only explores adding/removing *other* candidates.

    pk_depth_limit: If None, no limit.
                    If set, we interpret it as allowing this many *additional*
                    crossing layers on top of the scaffold's inherent complexity.
                    Effective Limit = Scaffold Layers + 1.
    """
    if seed is not None:
        random.seed(seed)

    # Initialize state with scaffold
    pair_set: set[tuple[int, int]] = set()
    partners = [-1] * L

    fixed_pairs = set()
    if scaffold_pairs:
        for i, j in scaffold_pairs:
            if i > j:
                i, j = j, i
            fixed_pairs.add((i, j))
            pair_set.add((i, j))
            partners[i] = j
            partners[j] = i

    M = len(candidate_pairs)

    # If there are no candidate pairs to sample from, simply return the scaffold.
    if M == 0:
        current_solution = set(fixed_pairs) if fixed_pairs else set()
        return [current_solution for _ in range(n_samples)]

    # Determine dynamic depth limit based on scaffold
    # If scaffold has N layers, we allow N + 1 layers total.
    effective_depth_limit = 9999
    if pk_depth_limit is not None or True:  # Force logic: allow 1 extra layer
        # Calculate scaffold layers
        if not fixed_pairs:
            scaffold_layers = 0
        else:
            scaffold_layers = len(pairs_to_layers(list(fixed_pairs)))

        # We allow 1 extra layer on top of the scaffold
        # So total layers allowed = scaffold_layers + 1
        effective_depth_limit = scaffold_layers + 1

    score = 0.0  # Score of *variable* pairs (delta)

    # Precompute weight lookup for pruning
    pair_weights: dict[tuple[int, int], float] = {}
    for i, j, w in candidate_pairs:
        if i > j:
            i, j = j, i
        pair_weights[(i, j)] = w

    samples: list[set[tuple[int, int]]] = []
    total_steps = burn_in + n_samples * thin

    for step in range(total_steps):
        idx = random.randrange(M)
        i, j, w = candidate_pairs[idx]

        if i > j:
            i, j = j, i

        # Hairpin length constraint
        if (j - i - 1) < min_loop_sep:
            continue

        pair = (i, j)

        # If pair is fixed (scaffold), we cannot toggle it.
        if pair in fixed_pairs:
            continue

        if pair in pair_set:
            # Propose removal
            delta = -w
            if delta >= 0 or random.random() < math.exp(beta * delta):
                pair_set.remove(pair)
                partners[i] = -1
                partners[j] = -1
                score += delta

        else:
            # Propose addition only if both residues currently unpaired
            if partners[i] == -1 and partners[j] == -1:
                # Disallow immediate ')(' adjacency (local constraint)
                invalid = False
                if i - 1 >= 0 and partners[i - 1] != -1 and partners[i - 1] < (i - 1):
                    invalid = True
                if (
                    not invalid
                    and j + 1 < L
                    and partners[j + 1] != -1
                    and partners[j + 1] > (j + 1)
                ):
                    invalid = True

                if invalid:
                    continue

                # Check PK depth limit (relative to scaffold)
                temp_set = list(pair_set)
                temp_set.append(pair)
                layers = pairs_to_layers(temp_set)
                if len(layers) > effective_depth_limit:
                    continue

                delta = w
                if delta >= 0:
                    accept_prob = 1.0
                else:
                    accept_prob = math.exp(beta * delta)

                # PK-specific penalty
                if pk_alpha > 0.0 and L > 0:
                    is_crossing = any(pairs_cross(pair, q) for q in pair_set)
                    if is_crossing:
                        pk_pairs, _, _, _ = compute_pk_stats_for_set(pair_set)
                        pk_load = len(pk_pairs) / float(L)
                        accept_prob *= math.exp(-pk_alpha * pk_load)

                if random.random() < accept_prob:
                    pair_set.add(pair)
                    partners[i] = j
                    partners[j] = i
                    score += delta

        # Record samples
        if step >= burn_in and (step - burn_in) % thin == 0:
            # Prune logic: ensure we don't prune fixed pairs
            # Pass explicit fixed_pairs set to prune_pseudoknots
            pruned = prune_pseudoknots(
                set(pair_set),
                pair_weights=pair_weights,
                L=L,
                pk_filter_frac=pk_filter_frac,
                pk_filter_max_cross_per_pair=pk_filter_max_cross_per_pair,
                pk_filter_max_total_cross=pk_filter_max_total_cross,
                fixed_pairs=fixed_pairs,
            )
            samples.append(pruned)

    return samples


def sample_matchings_modular(
    L: int,
    candidate_pairs: list[tuple[int, int, float]],
    n_samples: int,
    burn_in: int = 1000,
    thin: int = 10,
    min_hairpin: int = 3,
    beta: float = 1.0,
    seed: int | None = None,
    pk_alpha: float = 2.0,
    pk_gamma: float = 0.0,
    lonely_penalty: float = 0.0,
    scaffold_pairs: set[tuple[int, int]] | None = None,
    use_segment_moves: bool = True,
    segment_birth_prob: float = 0.2,
    segment_death_prob: float = 0.2,
    segment_swap_prob: float = 0.1,
    diagnostics_file: str | None = None,
) -> list[set[tuple[int, int]]]:
    """
    Sample matchings using the new modular MCMC sampler.

    This uses the refactored sampler with:
    - Proper delta energy computation with incremental cache updates
    - Segment moves (birth/death/swap) for better mixing
    - Guaranteed detailed balance through Hastings ratios
    - Optional diagnostics output

    Args:
        L: Sequence length
        candidate_pairs: List of (i, j, weight) candidates
        n_samples: Number of samples to collect
        burn_in: Number of burn-in steps
        thin: Thinning interval
        min_hairpin: Minimum hairpin loop size
        beta: Inverse temperature
        seed: Random seed
        pk_alpha: Pseudoknot penalty per crossing pair
        pk_gamma: Pseudoknot penalty for max crossings
        lonely_penalty: Penalty for isolated base pairs
        scaffold_pairs: Fixed pairs that must be in every sample
        use_segment_moves: Whether to enable segment moves
        segment_birth_prob: Probability of segment birth move
        segment_death_prob: Probability of segment death move
        segment_swap_prob: Probability of segment swap move
        diagnostics_file: Optional path to write diagnostics JSON

    Returns:
        List of sampled matchings
    """
    # Build energy parameters
    energy_params = EnergyParams(
        beta=beta,
        pk_alpha=pk_alpha,
        pk_gamma=pk_gamma,
        lonely_penalty=lonely_penalty,
    )

    # Build move probabilities
    toggle_prob = 1.0 - segment_birth_prob - segment_death_prob - segment_swap_prob
    move_probs: dict[MoveType, float] | None = None

    if use_segment_moves:
        move_probs = {
            MoveType.TOGGLE: max(0.0, toggle_prob),
            MoveType.SEGMENT_BIRTH: segment_birth_prob,
            MoveType.SEGMENT_DEATH: segment_death_prob,
            MoveType.SEGMENT_SWAP: segment_swap_prob,
        }

    # Build sampler config
    config = SamplerConfig(
        n_samples=n_samples,
        burn_in=burn_in,
        thin=thin,
        seed=seed,
        energy_params=energy_params,
        move_probs=move_probs,
        min_hairpin=min_hairpin,
    )

    # Build segment index if using segment moves
    segment_index = None
    if use_segment_moves:
        all_pairs = [(i, j) for i, j, _ in candidate_pairs]
        segment_index = build_candidate_segments(all_pairs, min_len=2, max_len=10)

    # Run sampler
    samples, diagnostics = sample_matchings_new(
        length=L,
        candidate_pairs=candidate_pairs,
        config=config,
        scaffold_pairs=scaffold_pairs,
        segment_index=segment_index,
    )

    # Write diagnostics if requested
    if diagnostics_file:
        summary = diagnostics.summary()
        # Add move-by-move acceptance rates
        summary["move_details"] = {
            mt.name: {
                "proposals": diagnostics.total_proposals.get(mt, 0),
                "accepted": diagnostics.accepted.get(mt, 0),
                "rate": diagnostics.acceptance_rate(mt),
            }
            for mt in MoveType
            if diagnostics.total_proposals.get(mt, 0) > 0
        }
        with open(diagnostics_file, "w") as f:
            json.dump(summary, f, indent=2)
        sys.stderr.write(f"[INFO] Wrote diagnostics to {diagnostics_file}\n")

    return samples


def _generate_beta_ladder(
    center_beta: float,
    n_chains: int,
    min_factor: float = 0.5,
    max_factor: float = 2.0,
) -> list[float]:
    if n_chains <= 1:
        return [center_beta]
    if center_beta <= 0:
        raise ValueError("center_beta must be > 0")
    if min_factor <= 0 or max_factor <= 0:
        raise ValueError("beta factors must be > 0")
    if min_factor > max_factor:
        min_factor, max_factor = max_factor, min_factor

    beta_min = center_beta * min_factor
    beta_max = center_beta * max_factor

    # Geometric spacing is generally better than linear for MH acceptance.
    ratio = (beta_max / beta_min) ** (1.0 / (n_chains - 1))
    return [beta_min * (ratio**k) for k in range(n_chains)]


def _scaled_move_probs_for_beta(
    base_probs: dict[MoveType, float],
    base_beta: float,
    beta: float,
    min_scale: float = 0.5,
    max_scale: float = 2.0,
) -> dict[MoveType, float]:
    """
    Heuristic: hotter chains (smaller beta) get more global moves.
    We scale non-toggle moves by s = clamp(base_beta / beta).
    """
    if beta <= 0 or base_beta <= 0:
        return dict(base_probs)

    scale = base_beta / beta
    if scale < min_scale:
        scale = min_scale
    if scale > max_scale:
        scale = max_scale

    toggle = base_probs.get(MoveType.TOGGLE, 0.0)
    out: dict[MoveType, float] = {MoveType.TOGGLE: 0.0}
    other_sum = 0.0
    for mt, p in base_probs.items():
        if mt == MoveType.TOGGLE:
            continue
        sp = max(0.0, p * scale)
        out[mt] = sp
        other_sum += sp

    # Keep toggle as the remainder, preserving normalization.
    # If other moves blow up, renormalize all.
    total = other_sum + max(0.0, toggle)
    if total <= 0:
        return {MoveType.TOGGLE: 1.0}

    # First, set toggle to remainder if possible.
    rem = 1.0 - other_sum
    if rem >= 0:
        out[MoveType.TOGGLE] = rem
        # Normalize small numeric drift.
        s = sum(out.values())
        if s > 0:
            for k in out:
                out[k] /= s
        return out

    # Otherwise renormalize everything.
    for k in out:
        out[k] /= other_sum
    out[MoveType.TOGGLE] = 0.0
    return out


def _partners_from_pairs(pair_set: set[tuple[int, int]], L: int) -> list[int]:
    partners = [-1] * L
    for i, j in pair_set:
        if i > j:
            i, j = j, i
        partners[i] = j
        partners[j] = i
    return partners


def _hamming_partner_distance(a: list[int], b: list[int]) -> int:
    return sum(x != y for x, y in zip(a, b))


def _approx_median_pairwise_hamming(partners_list: list[list[int]], max_pairs: int = 2000) -> float:
    if len(partners_list) < 2:
        return 0.0
    rng = random.Random(0)
    n = len(partners_list)
    dists: list[int] = []
    draws = min(max_pairs, n * (n - 1) // 2)
    for _ in range(draws):
        i = rng.randrange(n)
        j = rng.randrange(n - 1)
        if j >= i:
            j += 1
        dists.append(_hamming_partner_distance(partners_list[i], partners_list[j]))
    dists.sort()
    mid = len(dists) // 2
    return float(dists[mid])


def _binary_entropy(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))


def sample_matchings_multibeta_modular(
    L: int,
    candidate_pairs: list[tuple[int, int, float]],
    n_samples: int,
    burn_in: int = 1000,
    thin: int = 10,
    min_hairpin: int = 3,
    beta: float = 1.0,
    betas: list[float] | None = None,
    n_beta_chains: int = 4,
    beta_min_factor: float = 0.5,
    beta_max_factor: float = 2.0,
    seed: int | None = None,
    pk_alpha: float = 2.0,
    pk_gamma: float = 0.0,
    lonely_penalty: float = 0.0,
    scaffold_pairs: set[tuple[int, int]] | None = None,
    initial_pairs: set[tuple[int, int]] | None = None,
    use_segment_moves: bool = True,
    segment_birth_prob: float = 0.2,
    segment_death_prob: float = 0.2,
    segment_swap_prob: float = 0.1,
    diversity_dist_frac: float = 0.05,
    pk_filter_frac: float = 0.2,
    pk_filter_max_cross_per_pair: int = 2,
    pk_filter_max_total_cross: int = 30,
    delayed_accept: bool = True,
    parallel_chains: bool = True,
    max_workers: int = 8,
    use_beam_seeds: bool = True,
    beam_width: int = 64,
    beam_expansions_per_state: int = 24,
    beam_max_segments: int = 2000,
    beam_max_depth: int = 6,
    diagnostics_file: str | None = None,
    label: str | None = None,
) -> tuple[list[set[tuple[int, int]]], dict]:
    """
    Multi-beta strategy (no exchange): run short independent chains at different betas,
    then deduplicate + diversity-filter the union under a fixed output budget.

    `scaffold_pairs` are fixed (always-present) constraints.
    `initial_pairs` are only used to seed the chain state and are NOT fixed.
    """
    if betas is None or len(betas) == 0:
        betas = _generate_beta_ladder(
            center_beta=beta,
            n_chains=n_beta_chains,
            min_factor=beta_min_factor,
            max_factor=beta_max_factor,
        )
    betas = [float(b) for b in betas if b is not None]
    if not betas:
        betas = [beta]

    # Distribute samples (and burn-in) across chains so total compute stays similar.
    if n_samples <= 0:
        return [], {"label": label, "L": L, "n_requested": n_samples, "n_returned": 0}

    n_chains = min(len(betas), n_samples)
    betas = betas[:n_chains]
    base_burn_in = max(0, int(burn_in))
    burn_in_each = max(50, base_burn_in // n_chains) if n_chains > 1 else base_burn_in

    base = n_samples // n_chains
    rem = n_samples % n_chains
    n_samples_each = [base + (1 if k < rem else 0) for k in range(n_chains)]

    # Build energy parameters base (beta per-chain below).
    base_energy_params = EnergyParams(
        beta=beta,
        pk_alpha=pk_alpha,
        pk_gamma=pk_gamma,
        lonely_penalty=lonely_penalty,
        hairpin_min=min_hairpin,
    )

    # Base move probabilities (temperature-scaled per chain).
    toggle_prob = 1.0 - segment_birth_prob - segment_death_prob - segment_swap_prob
    base_move_probs: dict[MoveType, float] = {
        MoveType.TOGGLE: max(0.0, toggle_prob),
        MoveType.SEGMENT_BIRTH: max(0.0, segment_birth_prob),
        MoveType.SEGMENT_DEATH: max(0.0, segment_death_prob),
        MoveType.SEGMENT_SWAP: max(0.0, segment_swap_prob),
    }

    # Precompute segment index once (shared across chains) only if needed.
    segment_index = None
    if use_segment_moves or use_beam_seeds:
        all_pairs = [(i, j) for i, j, _ in candidate_pairs]
        segment_index = build_candidate_segments(all_pairs, min_len=2, max_len=10)

    # Weight lookup used for energy and PK pruning postprocess.
    pair_weights: dict[tuple[int, int], float] = {}
    for i, j, w in candidate_pairs:
        if i > j:
            i, j = j, i
        pair_weights[(i, j)] = w

    # Optional beam/constructive seeding to give each chain a distinct macro-state.
    beam_seeds: list[set[tuple[int, int]]] = []
    if use_beam_seeds and segment_index is not None and segment_index.all_segments_list:
        beam_seeds = beam_seed_segments(
            L=L,
            segment_index=segment_index,
            pair_weights=pair_weights,
            scaffold_pairs=set(scaffold_pairs or set()),
            min_hairpin=min_hairpin,
            n_seeds=n_chains,
            diversity_dist_frac=diversity_dist_frac,
            beam_width=beam_width,
            expansions_per_state=beam_expansions_per_state,
            max_segments=beam_max_segments,
            max_depth=beam_max_depth,
        )

    t0 = time.perf_counter()
    all_samples_with_energy: list[
        tuple[set[tuple[int, int]], float, float]
    ] = []  # (pairs, E, beta)
    per_chain_reports: list[dict] = []

    tasks: list[dict] = []
    for chain_idx, (b, ns) in enumerate(zip(betas, n_samples_each, strict=False)):
        if ns <= 0:
            continue
        chain_seed = int(seed) + 100_000 * chain_idx if seed is not None else None
        init_pairs = None
        if initial_pairs is not None:
            init_pairs = set(initial_pairs)
        elif beam_seeds:
            init_pairs = beam_seeds[chain_idx % len(beam_seeds)]
        tasks.append(
            {
                "chain_idx": chain_idx,
                "beta": float(b),
                "n_samples": int(ns),
                "burn_in": int(burn_in_each),
                "thin": int(thin),
                "seed": chain_seed,
                "move_probs": _scaled_move_probs_for_beta(base_move_probs, base_beta=beta, beta=b),
                "initial_pairs": init_pairs,
            }
        )

    # Run chains (parallel or sequential).
    chain_results: list[dict] = []
    if parallel_chains and len(tasks) > 1 and max_workers > 1:
        try:
            ctx = mp.get_context("fork")
        except ValueError:
            ctx = None

        if ctx is None:
            chain_results = [
                _run_chain_task(
                    task,
                    L=L,
                    candidate_pairs=candidate_pairs,
                    scaffold_pairs=scaffold_pairs,
                    segment_index=segment_index if use_segment_moves else None,
                    pair_weights=pair_weights,
                    base_energy_params=base_energy_params,
                    delayed_accept=delayed_accept,
                )
                for task in tasks
            ]
        else:
            # Set globals for forked workers to avoid pickling large data per task.
            _CHAIN_GLOBALS["L"] = L
            _CHAIN_GLOBALS["candidate_pairs"] = candidate_pairs
            _CHAIN_GLOBALS["scaffold_pairs"] = scaffold_pairs
            _CHAIN_GLOBALS["segment_index"] = segment_index if use_segment_moves else None
            _CHAIN_GLOBALS["pair_weights"] = pair_weights
            _CHAIN_GLOBALS["base_energy_params"] = base_energy_params
            _CHAIN_GLOBALS["delayed_accept"] = delayed_accept

            try:
                with ProcessPoolExecutor(
                    max_workers=min(int(max_workers), len(tasks)),
                    mp_context=ctx,
                ) as ex:
                    futures = [ex.submit(_run_chain_task_global, task) for task in tasks]
                    for fut in as_completed(futures):
                        chain_results.append(fut.result())
            except (PermissionError, OSError) as exc:
                # Some environments disallow process-based parallelism (e.g. blocked semaphores).
                # Fall back to sequential chain execution.
                sys.stderr.write(
                    f"[WARN] Chain parallelism unavailable ({exc}); running chains sequentially.\n"
                )
                chain_results = [
                    _run_chain_task(
                        task,
                        L=L,
                        candidate_pairs=candidate_pairs,
                        scaffold_pairs=scaffold_pairs,
                        segment_index=segment_index if use_segment_moves else None,
                        pair_weights=pair_weights,
                        base_energy_params=base_energy_params,
                        delayed_accept=delayed_accept,
                    )
                    for task in tasks
                ]
    else:
        chain_results = [
            _run_chain_task(
                task,
                L=L,
                candidate_pairs=candidate_pairs,
                scaffold_pairs=scaffold_pairs,
                segment_index=segment_index if use_segment_moves else None,
                pair_weights=pair_weights,
                base_energy_params=base_energy_params,
                delayed_accept=delayed_accept,
            )
            for task in tasks
        ]

    chain_results.sort(key=lambda r: r["chain_idx"])
    for res in chain_results:
        b = res["beta"]
        for s, e in zip(res["samples"], res["energies"], strict=False):
            all_samples_with_energy.append((s, float(e), float(b)))
        per_chain_reports.append(res["report"])

    # Deduplicate by structure (pair set), keeping best (lowest) energy.
    best_by_struct: dict[frozenset[tuple[int, int]], tuple[float, float, set[tuple[int, int]]]] = {}
    for s, e, b in all_samples_with_energy:
        fs = frozenset(s)
        prev = best_by_struct.get(fs)
        if prev is None or e < prev[0]:
            best_by_struct[fs] = (e, b, s)

    uniques = [(e, b, s) for (e, b, s) in best_by_struct.values()]
    uniques.sort(key=lambda x: x[0])  # by energy

    # Diversity filtering on partner vectors.
    threshold = max(1, int(round(diversity_dist_frac * L)))
    kept: list[set[tuple[int, int]]] = []
    kept_partners: list[list[int]] = []

    for e, b, s in uniques:
        p = _partners_from_pairs(s, L)
        if all(_hamming_partner_distance(p, q) > threshold for q in kept_partners):
            kept.append(s)
            kept_partners.append(p)
            if len(kept) >= n_samples:
                break

    # If diversity filter is too strict, fill by best energy.
    if len(kept) < n_samples:
        kept_fs = {frozenset(s) for s in kept}
        for e, b, s in uniques:
            if len(kept) >= n_samples:
                break
            if frozenset(s) in kept_fs:
                continue
            kept.append(s)
            kept_partners.append(_partners_from_pairs(s, L))

    # Postprocess: enforce PK thresholds only on final output (avoid O(n^2) in the loop).
    out: list[set[tuple[int, int]]] = []
    pk_stats_out: list[dict] = []
    for s in kept:
        pruned = prune_pseudoknots(
            set(s),
            pair_weights=pair_weights,
            L=L,
            pk_filter_frac=pk_filter_frac,
            pk_filter_max_cross_per_pair=pk_filter_max_cross_per_pair,
            pk_filter_max_total_cross=pk_filter_max_total_cross,
            fixed_pairs=set(scaffold_pairs or set()),
        )
        out.append(pruned)
        pk_pairs, _cc, total_cross, max_cross = compute_pk_stats_for_set(pruned)
        pk_stats_out.append(
            {
                "pk_pairs": len(pk_pairs),
                "total_crossings": total_cross,
                "max_crossings": max_cross,
            }
        )

    elapsed = time.perf_counter() - t0

    # Diversity summary
    out_partners = [_partners_from_pairs(s, L) for s in out]
    median_hamming = _approx_median_pairwise_hamming(out_partners)

    # Pair entropy on observed pairs (in final output)
    edge_counts: dict[tuple[int, int], int] = {}
    for s in out:
        for i, j in s:
            if i > j:
                i, j = j, i
            edge_counts[(i, j)] = edge_counts.get((i, j), 0) + 1

    k = max(1, len(out))
    mean_edge_entropy = 0.0
    if edge_counts:
        mean_edge_entropy = sum(_binary_entropy(c / k) for c in edge_counts.values()) / len(
            edge_counts
        )

    report = {
        "label": label,
        "L": L,
        "n_requested": n_samples,
        "n_returned": len(out),
        "n_total_raw": len(all_samples_with_energy),
        "n_unique_raw": len(uniques),
        "betas": betas,
        "beam_seeds_used": bool(beam_seeds),
        "burn_in_each": burn_in_each,
        "thin": thin,
        "wall_seconds": elapsed,
        "steps_per_chain": burn_in_each + thin * max(n_samples_each) if n_samples_each else 0,
        "diversity": {
            "hamming_threshold": threshold,
            "median_pairwise_hamming_approx": median_hamming,
            "median_pairwise_hamming_frac_approx": median_hamming / float(L) if L else 0.0,
            "mean_observed_edge_entropy_nats": mean_edge_entropy,
        },
        "pk_out": pk_stats_out,
        "chains": per_chain_reports,
    }

    if diagnostics_file:
        with open(diagnostics_file, "w") as f:
            json.dump(report, f, indent=2)
        sys.stderr.write(f"[INFO] Wrote sampler diagnostics to {diagnostics_file}\n")

    return out, report


_CHAIN_GLOBALS: dict[str, object] = {}


def _run_chain_task(
    task: dict,
    *,
    L: int,
    candidate_pairs: list[tuple[int, int, float]],
    scaffold_pairs: set[tuple[int, int]] | None,
    segment_index,
    pair_weights: dict[tuple[int, int], float],
    base_energy_params: EnergyParams,
    delayed_accept: bool,
) -> dict:
    chain_idx = int(task["chain_idx"])
    b = float(task["beta"])
    ns = int(task["n_samples"])
    burn_in = int(task["burn_in"])
    thin = int(task["thin"])
    chain_seed = task.get("seed")
    move_probs = task.get("move_probs")
    initial_pairs = task.get("initial_pairs")

    energy_params = EnergyParams(**base_energy_params.__dict__)
    energy_params.beta = b

    cfg = SamplerConfig(
        n_samples=ns,
        burn_in=burn_in,
        thin=thin,
        seed=chain_seed,
        energy_params=energy_params,
        move_probs=move_probs,
        min_hairpin=base_energy_params.hairpin_min,
        validate_structures=False,
        delayed_accept=delayed_accept,
    )

    samples, diagnostics = sample_matchings_new(
        length=L,
        candidate_pairs=candidate_pairs,
        config=cfg,
        scaffold_pairs=scaffold_pairs,
        initial_pairs=initial_pairs,
        segment_index=segment_index,
        weights=pair_weights,
    )

    return {
        "chain_idx": chain_idx,
        "beta": b,
        "samples": samples,
        "energies": list(diagnostics.energy_trace),
        "report": {
            "chain_idx": chain_idx,
            "beta": b,
            "n_samples": len(samples),
            "burn_in": burn_in,
            "thin": thin,
            "seed": chain_seed,
            "initial_pairs_size": len(initial_pairs) if initial_pairs else 0,
            "summary": diagnostics.summary(),
            "move_probs": {mt.name: float(move_probs.get(mt, 0.0)) for mt in MoveType}
            if isinstance(move_probs, dict)
            else None,
        },
    }


def _run_chain_task_global(task: dict) -> dict:
    return _run_chain_task(
        task,
        L=_CHAIN_GLOBALS["L"],
        candidate_pairs=_CHAIN_GLOBALS["candidate_pairs"],
        scaffold_pairs=_CHAIN_GLOBALS["scaffold_pairs"],
        segment_index=_CHAIN_GLOBALS["segment_index"],
        pair_weights=_CHAIN_GLOBALS["pair_weights"],
        base_energy_params=_CHAIN_GLOBALS["base_energy_params"],
        delayed_accept=_CHAIN_GLOBALS["delayed_accept"],
    )


def beam_seed_segments(
    *,
    L: int,
    segment_index,
    pair_weights: dict[tuple[int, int], float],
    scaffold_pairs: set[tuple[int, int]],
    min_hairpin: int,
    n_seeds: int,
    diversity_dist_frac: float,
    beam_width: int,
    expansions_per_state: int,
    max_segments: int,
    max_depth: int,
) -> list[set[tuple[int, int]]]:
    """
    Constructive beam search over helical segments to generate diverse initial states.

    Returns a list of initial (non-scaffold) pair sets.
    """
    segments = segment_index.all_segments_list
    if not segments or n_seeds <= 0:
        return []

    seg_pool = []
    for seg in segments:
        score = 0.0
        ok = True
        for i, j in seg.pairs():
            if (j - i - 1) < min_hairpin:
                ok = False
                break
            score += pair_weights.get((i, j) if i < j else (j, i), 0.0)
        if ok and score > 0:
            seg_pool.append((score, seg))

    seg_pool.sort(key=lambda x: x[0], reverse=True)
    seg_pool = seg_pool[: max(1, int(max_segments))]
    if not seg_pool:
        return []

    seg_pairs = [set(seg.pairs()) for _score, seg in seg_pool]
    seg_positions = [seg.positions() for _score, seg in seg_pool]
    seg_scores = [float(score) for score, _seg in seg_pool]

    scaffold_occ = {p for pair in scaffold_pairs for p in pair}
    base_pairs = set(scaffold_pairs)

    # Beam state: (score, pairs_set, occupied_positions_set)
    beam: list[tuple[float, set[tuple[int, int]], set[int]]] = [(0.0, base_pairs, scaffold_occ)]

    for _depth in range(max_depth):
        candidates: list[tuple[float, set[tuple[int, int]], set[int]]] = []
        for score, pairs_set, occ in beam:
            added = 0
            for idx in range(len(seg_pool)):
                if added >= expansions_per_state:
                    break
                if seg_positions[idx] & occ:
                    continue
                new_pairs = pairs_set | seg_pairs[idx]
                new_occ = occ | seg_positions[idx]
                candidates.append((score + seg_scores[idx], new_pairs, new_occ))
                added += 1
        if not candidates:
            break
        # Keep top beam_width candidates (plus current beam to allow early stop)
        merged = beam + candidates
        merged.sort(key=lambda x: x[0], reverse=True)
        beam = merged[: max(1, int(beam_width))]

    # Extract best candidates and enforce diversity.
    beam.sort(key=lambda x: x[0], reverse=True)
    threshold = max(1, int(round(diversity_dist_frac * L)))
    seeds: list[set[tuple[int, int]]] = []
    seed_partners: list[list[int]] = []

    for _score, pairs_set, _occ in beam:
        init = set(pairs_set) - scaffold_pairs
        partners = _partners_from_pairs(pairs_set, L)
        if all(_hamming_partner_distance(partners, p) > threshold for p in seed_partners):
            seeds.append(init)
            seed_partners.append(partners)
            if len(seeds) >= n_seeds:
                break

    return seeds


# ------------------ CLI ------------------------------------------ #


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Sample alternative pseudoknotted structures from a "
            "CaCoFold .sto for a given sequence and write to a .db file."
        )
    )
    parser.add_argument("sto", help="CaCoFold/R-scape Stockholm (.sto) file")
    parser.add_argument(
        "--seq",
        required=True,
        help="Sequence name to project/summarize (must match a sequence ID in the .sto).",
    )
    parser.add_argument(
        "--list-seqs",
        action="store_true",
        help="List available sequence names and exit.",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=10,
        help="Number of structures to sample (default: 10).",
    )
    parser.add_argument(
        "--burn-in",
        type=int,
        default=2000,
        help="Number of burn-in steps (default: 2000).",
    )
    parser.add_argument(
        "--thin",
        type=int,
        default=20,
        help="Thinning interval (keep 1 in every N steps after burn-in; default: 20).",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="Inverse temperature for Metropolis acceptance (default: 1.0).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed (optional).",
    )
    parser.add_argument(
        "--out-db",
        required=True,
        help="Output .db file to write sequence + sampled structures.",
    )
    parser.add_argument(
        "--min-loop-sep",
        type=int,
        default=1,
        help=(
            "Minimum number of unpaired residues between the two ends "
            "of a base pair (default: 1, which forbids '()')."
        ),
    )
    parser.add_argument(
        "--summary",
        default=None,
        help=(
            "Optional summary file. If provided, writes the projected "
            "consensus structure (SS_cons) and PK layers derived from "
            "all SS_cons* tracks."
        ),
    )
    # ---- Covariation options ----
    parser.add_argument(
        "--cov-file",
        default=None,
        help="Optional R-scape/CaCoFold .cov file to use for covariation weights.",
    )
    parser.add_argument(
        "--cov-mode",
        choices=["off", "logE", "logE_power", "score", "indicator"],
        default="off",
        help="How to convert covariation stats into pair weights (default: off).",
    )
    parser.add_argument(
        "--cov-alpha",
        type=float,
        default=1.0,
        help="Scale factor for covariation weight when combining with layer weights.",
    )
    parser.add_argument(
        "--cov-min-power",
        type=float,
        default=0.0,
        help="Drop or downweight pairs with covariation power below this threshold.",
    )
    parser.add_argument(
        "--cov-forbid-negative",
        action="store_true",
        help=(
            "If set, drop pairs with high power but weak covariation "
            "(large E-value) and optionally pairs with very low power."
        ),
    )
    # ---- Pseudoknot penalty / filtering options ----
    parser.add_argument(
        "--pk-alpha",
        type=float,
        default=2.0,
        help=(
            "Strength of pseudoknot penalty. If > 0, then whenever a proposed "
            "pair crosses any existing pair, its Metropolis acceptance "
            "probability is multiplied by exp(-pk_alpha * pk_load), where "
            "pk_load = (#pk_pairs)/L using the current matching."
        ),
    )
    parser.add_argument(
        "--pk-filter-frac",
        type=float,
        default=0.2,
        help=(
            "Maximum fraction of pseudoknotted pairs allowed before a sampled "
            "structure is rejected. The condition is len(pk_pairs) <= "
            "pk_filter_frac * L (default: 0.2 ~ L/5). Set <= 0 to disable "
            "this filter."
        ),
    )
    parser.add_argument(
        "--pk-filter-max-cross_per_pair",
        type=int,
        default=2,
        help=(
            "Maximum number of crossings allowed for any single pair in a "
            "sampled structure (default: 2). Set < 0 to disable this filter."
        ),
    )
    parser.add_argument(
        "--pk-filter-max-total-cross",
        type=int,
        default=30,
        help=(
            "Maximum total number of pair-pair crossings allowed in a "
            "sampled structure (default: 30). Set < 0 to disable this filter."
        ),
    )
    parser.add_argument(
        "--pk-depth-limit",
        type=int,
        default=None,
        help="Limit the pseudoknot depth during sampling (e.g. 1 = 1 crossing layer).",
    )
    # ---- New modular sampler options ----
    sampler_group = parser.add_mutually_exclusive_group()
    sampler_group.add_argument(
        "--legacy-sampler",
        action="store_true",
        help="Use the legacy sampler (deprecated; not detailed-balance-correct).",
    )
    sampler_group.add_argument(
        "--use-new-sampler",
        action="store_true",
        help=(
            "Use the new modular MCMC sampler (default). This flag is deprecated "
            "and kept for compatibility."
        ),
    )
    parser.add_argument(
        "--pk-gamma",
        type=float,
        default=0.0,
        help="Pseudoknot penalty based on max crossings per pair (new sampler only).",
    )
    parser.add_argument(
        "--lonely-penalty",
        type=float,
        default=0.0,
        help="Penalty for isolated (lonely) base pairs (new sampler only).",
    )
    parser.add_argument(
        "--no-segment-moves",
        action="store_true",
        help="Disable segment birth/death/swap moves (new sampler only).",
    )
    parser.add_argument(
        "--segment-birth-prob",
        type=float,
        default=0.2,
        help="Probability of segment birth move (new sampler only, default: 0.2).",
    )
    parser.add_argument(
        "--segment-death-prob",
        type=float,
        default=0.2,
        help="Probability of segment death move (new sampler only, default: 0.2).",
    )
    parser.add_argument(
        "--segment-swap-prob",
        type=float,
        default=0.1,
        help="Probability of segment swap move (new sampler only, default: 0.1).",
    )
    parser.add_argument(
        "--diagnostics",
        default=None,
        help="Path to write sampler diagnostics JSON (new sampler only).",
    )

    args = parser.parse_args()

    try:
        seqs, ss_tracks = parse_stockholm_single(args.sto)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to parse Stockholm file: {e}\n")
        sys.exit(1)

    if args.list_seqs:
        for name in sorted(seqs.keys()):
            print(name)
        return

    if args.seq not in seqs:
        sys.stderr.write(f"[ERROR] Sequence '{args.seq}' not found in {args.sto}\n")
        sys.exit(1)

    # Alignment → sequence map for this sequence, used both for cov and candidates
    aligned_seq = seqs[args.seq]
    aln2seq, _L_aln = aln_to_seq_map(aligned_seq)

    cov_dict: dict[tuple[int, int], CovRecord] | None = None
    if args.cov_file is not None and args.cov_mode != "off":
        try:
            cov_dict = load_cov_stats(args.cov_file, aln2seq)
            sys.stderr.write(
                f"[INFO] Loaded covariation stats for {len(cov_dict)} pairs from {args.cov_file}\n"
            )
        except Exception as e:
            sys.stderr.write(f"[WARN] Failed to load cov stats from {args.cov_file}: {e}\n")
            cov_dict = None

    try:
        L, candidate_pairs = extract_candidate_pairs(
            seq_name=args.seq,
            seqs=seqs,
            ss_tracks=ss_tracks,
            w_core=2.0,
            w_alt=1.0,
            cov=cov_dict,
            cov_mode=args.cov_mode,
            cov_alpha=args.cov_alpha,
            cov_min_power=args.cov_min_power,
            cov_forbid_negative=args.cov_forbid_negative,
        )
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to extract candidate pairs: {e}\n")
        sys.exit(1)

    if not candidate_pairs:
        sys.stderr.write("[ERROR] No candidate base pairs found for this sequence.\n")
        sys.exit(1)

    sys.stderr.write(
        f"[INFO] Sequence: {args.seq}, length (ungapped): {L}\n"
        f"[INFO] Candidate pairs: {len(candidate_pairs)}\n"
    )

    # Ungapped sequence (no FASTA headers in the .db file)
    ungapped_seq = "".join(ch for ch in seqs[args.seq] if ch not in GAP_CHARS)

    # Project the primary SS_cons track for the consensus secondary structure
    try:
        L_cons, ss_cons_str, _cons_pairs = project_ss_cons_to_sequence(args.seq, seqs, ss_tracks)
        if L_cons != L:
            sys.stderr.write(
                f"[WARN] Ungapped length from SS_cons ({L_cons}) != "
                f"length from candidate pairs ({L}); using {L}.\n"
            )
    except Exception as e:
        sys.stderr.write(f"[WARN] Could not project SS_cons: {e}\n")
        ss_cons_str = "." * L

    # Optional summary file: consensus + PK layers
    if args.summary is not None:
        write_summary_file(
            summary_path=args.summary,
            seq_name=args.seq,
            ungapped_seq=ungapped_seq,
            L=L,
            ss_cons_str=ss_cons_str,
            candidate_pairs=candidate_pairs,
        )
        sys.stderr.write(f"[INFO] Wrote CaCoFold summary to {args.summary}\n")

    use_new_sampler = True
    if args.legacy_sampler:
        use_new_sampler = False
    elif args.use_new_sampler:
        use_new_sampler = True

    # Sample structures
    if use_new_sampler:
        sys.stderr.write("[INFO] Using new modular MCMC sampler\n")
        samples = sample_matchings_modular(
            L=L,
            candidate_pairs=candidate_pairs,
            n_samples=args.n_samples,
            burn_in=args.burn_in,
            thin=args.thin,
            min_hairpin=args.min_loop_sep,
            beta=args.beta,
            seed=args.seed,
            pk_alpha=args.pk_alpha,
            pk_gamma=args.pk_gamma,
            lonely_penalty=args.lonely_penalty,
            use_segment_moves=not args.no_segment_moves,
            segment_birth_prob=args.segment_birth_prob,
            segment_death_prob=args.segment_death_prob,
            segment_swap_prob=args.segment_swap_prob,
            diagnostics_file=args.diagnostics,
        )
    else:
        sys.stderr.write("[WARN] Using legacy sampler (detailed balance not guaranteed)\n")
        samples = sample_matchings(
            L=L,
            candidate_pairs=candidate_pairs,
            n_samples=args.n_samples,
            burn_in=args.burn_in,
            thin=args.thin,
            min_loop_sep=args.min_loop_sep,
            beta=args.beta,
            seed=args.seed,
            pk_alpha=args.pk_alpha,
            pk_filter_frac=args.pk_filter_frac,
            pk_filter_max_cross_per_pair=args.pk_filter_max_cross_per_pair,
            pk_filter_max_total_cross=args.pk_filter_max_total_cross,
            pk_depth_limit=args.pk_depth_limit,
            # scaffold_pairs argument is not exposed in CLI directly but used by pipeline code calling `sample_matchings`
        )

    # Write to .db file:
    #   line 1: sequence
    #   line 2+: one PK dot-bracket string per sample
    with open(args.out_db, "w") as out_f:
        out_f.write(ungapped_seq + "\n")
        for pair_set in samples:
            pk = pairs_to_pk_string(sorted(pair_set), L)
            out_f.write(pk + "\n")

    sys.stderr.write(
        f"[INFO] Wrote {len(samples)} structures to {args.out_db} (requested {args.n_samples}).\n"
    )


if __name__ == "__main__":
    main()
