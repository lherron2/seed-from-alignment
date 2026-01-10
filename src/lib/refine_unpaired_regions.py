#!/usr/bin/env python3
"""
Refine long unpaired regions in CaCoFold-sampled structures using an enumerative combinatorial approach.

Refinement Tracks:
    0. Consensus Masking (New Layer):
       - Before refinement begins, we generate variants of the input consensus structure.
       - We identify all helices in the input.
       - We generate scaffolds by masking (removing) these helices:
            a) Sequential Masking: Remove 1 helix at a time (top N by length).
            b) Pairwise Masking: Remove pairs of helices (top M by length).
            c) End Masking Grid: Systematically unpair 5' and 3' ends in increments of 5nt
               to force exploration of alternative terminal structures.
       - This forces the pipeline to re-explore regions that were "locked" by the consensus,
         allowing for alternative folds in those specific areas.

    1. Terminal Initialization (Seed Mode):
       - For every Masked Scaffold, we explicitly calculate 5' and 3' end interactions
         using dense, incremental windowing to identify potential global closing stems.
       - Seeds include the Masked Scaffold + these new 5'-3' pairs.

    2. Independent Combinatorial Assembly (Enhanced for Kissing Loops):
       - Input: Unpaired regions in the (Masked + Seeded) scaffold.
       - Candidates: Local AllSub vs. Interaction DuplexFold.
       - Logic: Recursively explore compatible combinations.

    3. Sequential Hierarchical Assembly:
       - Hierarchical approach: Force local folds first, then interact loops.

    4. Post-Process Kissing Loop Discovery:
       - Intercept top candidates.
       - Check unpaired apical loops for complementarity.
       - Add inter-loop (kissing) interactions if topologically valid.

    - Validation & Filtering:
       - STRICT CANONICAL PAIR CHECK.
       - ROSETTA FORMAT.
       - Topological complexity limits.
"""

from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

# Re-use logic from sample_cacofold_structures where appropriate
from src.lib.sample_cacofold_structures import pairs_from_track

UNPAIRED_CHAR = "."
GAP_CHARS = set("-._~")
CANONICAL_PAIRS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}

# --- Default Constants (used if not overridden by config) ---
DEFAULT_END_MASK_STEP = 5  # Increment size for unpairing ends
DEFAULT_MAX_END_MASK_LEN = 40  # Max length to unpair from either end during grid search

DEFAULT_MAX_REGIONS_TO_REFINE = (
    30  # Safety cap: Only refine the N longest runs to prevent explosion
)
DEFAULT_MAX_SOLUTIONS = 2000  # Applied PER TRACK (so up to 2x this total before filtering)
DEFAULT_KISSING_LOOP_CANDIDATES = 1000  # Max structures to check for kissing loops

DEFAULT_MAX_HELICES_SEQUENTIAL = 20  # Consider top N helices for single removal
DEFAULT_MAX_HELICES_PAIRWISE = 10  # Consider top N helices for pairwise removal (N*(N-1)/2 combos)
DEFAULT_MAX_SEEDS = 50  # Maximum number of 5'-3' initialization seeds to process per mask

# Scoring Constants
DEFAULT_SCAFFOLD_PAIR_ENERGY = -1.5  # Estimated kcal/mol per scaffold base pair
DEFAULT_WEIGHT_L0 = 1.0  # Nested: Full stability
DEFAULT_WEIGHT_L1 = 0.6  # PK Layer 1: Reduced stability
DEFAULT_WEIGHT_L2_PLUS = -1.0  # PK Layer 2+: Penalty (flips sign of stability)

# --- Graph Coloring for Layers (Local Definition for Safety) ---


def pairs_cross(p: tuple[int, int], q: tuple[int, int]) -> bool:
    """Return True if arcs (i, j) and (k, l) cross."""
    i, j = sorted(p)
    k, l = sorted(q)
    return (i < k < j < l) or (k < i < l < j)


def pairs_to_layers(pairs: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """Partition pairs into layers such that no pairs in a layer cross."""
    layers: list[list[tuple[int, int]]] = []
    for p in sorted(pairs):
        placed = False
        for layer in layers:
            if not any(pairs_cross(p, q) for q in layer):
                layer.append(p)
                placed = True
                break
        if not placed:
            layers.append([p])
    return layers


def calculate_weighted_score(
    all_pairs: set[tuple[int, int]],
    refined_energy: float,
    scaffold_pairs: set[tuple[int, int]],
    scaffold_pair_energy: float = DEFAULT_SCAFFOLD_PAIR_ENERGY,
    weight_l0: float = DEFAULT_WEIGHT_L0,
    weight_l1: float = DEFAULT_WEIGHT_L1,
    weight_l2_plus: float = DEFAULT_WEIGHT_L2_PLUS,
) -> float:
    """
    Calculate a topology-aware weighted score.

    Formula:
      Score = Sum( Weight(Layer(p)) * Energy(p) )

    Where:
      - Energy(p) for refined pairs is (refined_energy / N_refined).
      - Energy(p) for scaffold pairs is scaffold_pair_energy.
      - Weight depends on topological layer (0, 1, 2+).
    """
    if not all_pairs:
        return 0.0

    # 1. Determine global layers for ALL pairs together
    sorted_pairs = sorted(all_pairs)
    layers = pairs_to_layers(sorted_pairs)

    # Map pair -> layer index
    pair_to_layer = {}
    for idx, layer in enumerate(layers):
        for p in layer:
            pair_to_layer[p] = idx

    # 2. Identify refined pairs
    refined_pairs = all_pairs - scaffold_pairs
    n_refined = len(refined_pairs)

    avg_refined_energy = 0.0
    if n_refined > 0:
        avg_refined_energy = refined_energy / n_refined

    score = 0.0

    # 3. Sum weighted energy for Refined Pairs
    for p in refined_pairs:
        layer_idx = pair_to_layer.get(p, 0)
        w = weight_l0
        if layer_idx == 1:
            w = weight_l1
        elif layer_idx >= 2:
            w = weight_l2_plus

        score += w * avg_refined_energy

    # 4. Sum weighted energy for Scaffold Pairs
    for p in scaffold_pairs:
        # Only count if present in current structure (masking might have removed some)
        if p in all_pairs:
            layer_idx = pair_to_layer.get(p, 0)
            w = weight_l0
            if layer_idx == 1:
                w = weight_l1
            elif layer_idx >= 2:
                w = weight_l2_plus

            score += w * scaffold_pair_energy

    return score


def pairs_to_rosetta_string(pairs: list[tuple[int, int]], L: int) -> str:
    """
    Convert pairs to a dot-bracket string using Rosetta-safe hierarchy.
    """
    indices = [idx for pair in pairs for idx in pair]
    if len(indices) != len(set(indices)):
        raise ValueError("Invalid matching: overlapping pairs detected in refined structure.")

    chars = ["."] * L
    layers = pairs_to_layers(pairs)

    brackets = [("(", ")"), ("[", "]"), ("{", "}")]
    for i in range(26):
        brackets.append((chr(ord("a") + i), chr(ord("A") + i)))

    for layer_idx, layer in enumerate(layers):
        if layer_idx < len(brackets):
            op, cl = brackets[layer_idx]
        else:
            op, cl = "{", "}"  # Fallback

        for i, j in layer:
            chars[i] = op
            chars[j] = cl
    return "".join(chars)


# --- Validators ---


def is_canonical(b1: str, b2: str) -> bool:
    """Check if two bases form a canonical Watson-Crick or Wobble pair."""
    return (b1.upper(), b2.upper()) in CANONICAL_PAIRS


def validate_and_filter_pairs(pairs: set[tuple[int, int]], seq: str) -> set[tuple[int, int]]:
    """Return only pairs that are valid indices and canonical matches."""
    valid = set()
    L = len(seq)
    for i, j in pairs:
        if 0 <= i < L and 0 <= j < L:
            if is_canonical(seq[i], seq[j]):
                valid.add((i, j))
    return valid


def remove_isolated_pairs(pairs: set[tuple[int, int]]) -> set[tuple[int, int]]:
    """Iteratively remove pairs that have no stacking neighbor (isolated)."""
    current_pairs = pairs.copy()
    while True:
        to_remove = set()
        for i, j in current_pairs:
            p_down = (i - 1, j + 1)
            p_up = (i + 1, j - 1)
            has_down = (p_down in current_pairs) or ((p_down[1], p_down[0]) in current_pairs)
            has_up = (p_up in current_pairs) or ((p_up[1], p_up[0]) in current_pairs)
            if not (has_down or has_up):
                to_remove.add((i, j))
        if not to_remove:
            break
        current_pairs -= to_remove
    return current_pairs


# --- I/O ---


def read_ungapped_seq_from_sto(path: Path, seq_name: str | None = None) -> str:
    sequences: dict[str, list[str]] = {}
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                continue
            name, s = parts[0], parts[1].strip()
            sequences.setdefault(name, []).append(s)

    if not sequences:
        raise ValueError(f"No sequences found in Stockholm file {path}")

    if seq_name is None:
        name = next(iter(sequences))
    else:
        if seq_name not in sequences:
            raise ValueError(f"Sequence '{seq_name}' not found.")
        name = seq_name

    aligned = "".join(sequences[name])
    ungapped = "".join(ch for ch in aligned if ch not in GAP_CHARS)
    return ungapped.upper().replace("T", "U")


def read_fasta_sequence(path: Path, seq_name: str | None = None) -> str:
    seqs: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    seqs[current_name] = "".join(chunks)
                header = line[1:].strip()
                current_name = header.split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_name is not None:
        seqs[current_name] = "".join(chunks)

    if not seqs:
        raise ValueError(f"No sequences found in {path}")

    if seq_name is None:
        name = next(iter(seqs))
    else:
        if seq_name not in seqs:
            raise ValueError(f"Sequence '{seq_name}' not found.")
        name = seq_name

    return seqs[name].upper().replace("T", "U")


def read_db_structures(path: Path) -> list[str]:
    structs: list[str] = []
    sequence: str | None = None
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if sequence is None:
                if all(ch not in ".()[]{}<>" for ch in s):
                    sequence = s
                    continue
            structs.append(s)
    if not structs:
        raise ValueError(f"No structures found in {path}")
    return structs


def find_unpaired_runs(struct: str, min_len: int) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    n = len(struct)
    i = 0
    while i < n:
        if struct[i] == UNPAIRED_CHAR:
            start = i
            while i < n and struct[i] == UNPAIRED_CHAR:
                i += 1
            if (i - start) >= min_len:
                runs.append((start, i))
        else:
            i += 1
    return runs


def parse_energy_from_header(header: str) -> float:
    """Extract free energy from CT file header line (standard RNAstructure format)."""
    # Example: "  30  -12.3  ENERGY = -12.3  seq_name"
    # Sometimes: "  30  ENERGY = -12.3  seq_name"
    # Or just "  30  seq_name" (no energy)
    upper = header.upper()
    if "ENERGY =" in upper:
        try:
            # "... ENERGY = -12.3 ..."
            val_str = upper.split("ENERGY =")[1].strip().split()[0]
            return float(val_str)
        except (IndexError, ValueError):
            pass

    # Fallback: try second token if it looks like a float
    parts = header.split()
    if len(parts) >= 2:
        try:
            return float(parts[1])
        except ValueError:
            pass

    return 0.0


def parse_ct_file(ct_path: Path) -> list[tuple[str, float]]:
    """Parse CT file into (structure, energy) tuples."""
    structures: list[tuple[str, float]] = []
    try:
        with ct_path.open() as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
    except FileNotFoundError:
        return []

    if not lines:
        return []

    idx = 0
    while idx < len(lines):
        header = lines[idx]
        energy = parse_energy_from_header(header)
        parts = header.split()
        if not parts:
            idx += 1
            continue
        try:
            seq_len = int(parts[0])
        except ValueError:
            idx += 1
            continue
        if idx + seq_len >= len(lines):
            break
        db_chars = ["."] * seq_len
        for offset in range(1, seq_len + 1):
            line = lines[idx + offset]
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                i = int(fields[0]) - 1
                j = int(fields[4]) - 1
            except ValueError:
                continue
            if j > i:
                if 0 <= i < seq_len and 0 <= j < seq_len:
                    db_chars[i] = "("
                    db_chars[j] = ")"
        structures.append(("".join(db_chars), energy))
        idx += seq_len + 1
    return structures


# --- RNAstructure Wrappers ---


def call_rnastructure_allsub(
    allsub_exe: Path,
    subseq: str,
    temperature: float | None = None,
    absolute_energy: float | None = None,
    percent_energy: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> list[tuple[str, float]]:
    extra_args = list(extra_args) if extra_args is not None else []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        seq_path = tmpdir_path / "subseq.fa"
        ct_path = tmpdir_path / "output.ct"

        with seq_path.open("w") as fh:
            fh.write(">subseq\n")
            fh.write(subseq + "\n")

        cmd: list[str] = [str(allsub_exe), str(seq_path), str(ct_path)]
        if temperature is not None:
            cmd.extend(["-t", str(temperature)])
        if absolute_energy is not None:
            cmd.extend(["-a", str(absolute_energy)])
        if percent_energy is not None:
            cmd.extend(["-p", str(percent_energy)])
        cmd.extend(extra_args)

        try:
            subprocess.run(cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError:
            return [("." * len(subseq), 0.0)]

        if not ct_path.is_file():
            return [("." * len(subseq), 0.0)]
        structures = parse_ct_file(ct_path)
        if not structures:
            return [("." * len(subseq), 0.0)]
        return structures


def call_rnastructure_duplexfold(
    duplex_exe: Path,
    seq5: str,
    seq3: str,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> list[tuple[list[tuple[int, int]], float]]:
    extra_args = list(extra_args) if extra_args is not None else []
    if not any(arg.startswith("-m") for arg in extra_args):
        extra_args.extend(["-m", "100"])

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        seq5_path = tmpdir_path / "seq5.fa"
        seq3_path = tmpdir_path / "seq3.fa"
        ct_path = tmpdir_path / "duplex.ct"

        with seq5_path.open("w") as fh:
            fh.write(">seq5\n")
            fh.write(seq5 + "\n")
        with seq3_path.open("w") as fh:
            fh.write(">seq3\n")
            fh.write(seq3 + "\n")

        cmd: list[str] = [str(duplex_exe), str(seq5_path), str(seq3_path), str(ct_path)]
        if temperature is not None:
            cmd += ["-t", str(temperature)]
        cmd += extra_args

        try:
            subprocess.run(cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError:
            return []

        structures = parse_ct_file(ct_path)
        if not structures:
            return []

        seq_len_ct = len(structures[0][0])
        results = []
        n1 = len(seq5)
        n2 = len(seq3)
        linker_len = max(0, seq_len_ct - (n1 + n2))

        for db, energy in structures:
            pairs = []
            stack = []
            for i, ch in enumerate(db):
                if ch == "(":
                    stack.append(i)
                elif ch == ")":
                    if stack:
                        j = stack.pop()
                        u, v = sorted((j, i))
                        if u < n1 and v >= (n1 + linker_len):
                            i5_local = u + 1
                            j3_local = v - (n1 + linker_len) + 1
                            if i5_local <= n1 and j3_local <= n2:
                                pairs.append((i5_local, j3_local))
            results.append((pairs, energy))
        return results


def _pairs_from_struct(struct: str) -> set[tuple[int, int]]:
    return set(pairs_from_track(struct))


def _struct_string_to_pairs(s: str, offset: int = 0) -> set[tuple[int, int]]:
    pairs = set()
    stack = []
    for i, ch in enumerate(s):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                pairs.add((offset + j, offset + i))
    return pairs


def get_loop_interaction_options(
    sub_runs,
    full_seq,
    duplex_exe,
    temperature,
    extra_args,
    duplex_cache,
    L,
    dense_window: bool = False,
) -> dict[tuple[int, int], list[tuple[set[tuple[int, int]], float]]]:
    options = {}

    def get_len_steps(length: int) -> list[int]:
        if dense_window:
            return list(range(4, length + 1))

        if length <= 8:
            return list(range(4, length + 1))
        stride = 4 if length < 30 else 6
        steps = list(range(4, length, stride))
        steps.append(length)
        return sorted(set(steps))

    for r_idx_a in range(len(sub_runs)):
        for r_idx_b in range(r_idx_a + 1, len(sub_runs)):
            sa, ea = sub_runs[r_idx_a]
            sb, eb = sub_runs[r_idx_b]

            dist = sb - ea
            if dist < 3:
                continue

            len_a_total = ea - sa
            len_b_total = eb - sb
            if len_a_total < 4 or len_b_total < 4:
                continue

            steps_a = get_len_steps(len_a_total)
            steps_b = get_len_steps(len_b_total)

            all_valid_int_sets = []

            for len_a in steps_a:
                seq_a = full_seq[sa : sa + len_a]
                for len_b in steps_b:
                    seq_b = full_seq[eb - len_b : eb]

                    dkey = (seq_a, seq_b, temperature, tuple(extra_args or []))
                    if dkey in duplex_cache:
                        d_res = duplex_cache[dkey]
                    else:
                        d_res = call_rnastructure_duplexfold(
                            duplex_exe, seq_a, seq_b, temperature, extra_args
                        )
                        duplex_cache[dkey] = d_res

                    if d_res:
                        for pair_list, energy in d_res:
                            if pair_list:
                                g_pairs = set()
                                for i_loc, j_loc in pair_list:
                                    gi = sa + i_loc - 1
                                    gj = (eb - len_b) + j_loc - 1
                                    if gi > gj:
                                        gi, gj = gj, gi
                                    g_pairs.add((gi, gj))
                                canon_pairs = validate_and_filter_pairs(g_pairs, full_seq)
                                if canon_pairs:
                                    all_valid_int_sets.append((canon_pairs, energy))

            if all_valid_int_sets:
                # Deduplicate based on pair set, keeping best energy
                seen_sets = {}
                for p_set, en in all_valid_int_sets:
                    fs = frozenset(p_set)
                    if fs not in seen_sets or en < seen_sets[fs][1]:
                        seen_sets[fs] = (set(fs), en)

                options[(r_idx_a, r_idx_b)] = list(seen_sets.values())
    return options


def merge_disjoint_options(
    option_sets: list[tuple[set[tuple[int, int]], float]], max_combo=50
) -> list[tuple[set[tuple[int, int]], float]]:
    combined = list(option_sets)
    n = len(option_sets)
    limit = min(n, 20)
    new_combos = []

    # Store combined as (frozenset, energy) to check dupes
    seen_combos = {frozenset(p): e for p, e in combined}

    for i in range(limit):
        for j in range(i + 1, limit):
            s1, e1 = option_sets[i]
            s2, e2 = option_sets[j]
            idx1 = {x for p in s1 for x in p}
            idx2 = {x for p in s2 for x in p}
            if idx1.isdisjoint(idx2):
                union_set = s1.union(s2)
                union_energy = e1 + e2
                fs = frozenset(union_set)

                if fs not in seen_combos or union_energy < seen_combos[fs]:
                    seen_combos[fs] = union_energy
                    new_combos.append((union_set, union_energy))
                    if len(new_combos) >= max_combo:
                        break
        if len(new_combos) >= max_combo:
            break

    # Convert back to list format
    # Note: we don't strictly enforce uniqueness here, but it helps
    combined.extend(new_combos)
    return combined


def get_5p3p_seeds(
    scaffold_pairs: set[tuple[int, int]],
    full_seq: str,
    duplex_exe: Path,
    temperature: float | None,
    extra_args: Sequence[str] | None,
    duplex_cache: dict,
    L: int,
    max_seeds: int = DEFAULT_MAX_SEEDS,
) -> list[tuple[set[tuple[int, int]], float]]:
    """Explicitly find 5'-3' interaction candidates to seed the refinement."""
    paired_mask = [False] * L
    for i, j in scaffold_pairs:
        if 0 <= i < L:
            paired_mask[i] = True
        if 0 <= j < L:
            paired_mask[j] = True

    start_5p = 0
    end_5p = 0
    while end_5p < L and not paired_mask[end_5p]:
        end_5p += 1

    start_3p = L
    end_3p = L
    while start_3p > 0 and not paired_mask[start_3p - 1]:
        start_3p -= 1

    # Overlap/Intersection Logic
    if end_5p >= start_3p:
        # Overlapping regions. Bisect into two domains with a 3nt gap for the linker constraint.
        total_range_len = end_3p - start_5p
        if total_range_len < 10:
            return []

        mid = start_5p + (total_range_len // 2)
        new_end_5p = mid - 1
        new_start_3p = mid + 2

        if (new_end_5p - start_5p) < 4 or (end_3p - new_start_3p) < 4:
            return []

        term_runs = [(start_5p, new_end_5p), (new_start_3p, end_3p)]

    elif (end_5p - start_5p) < 3 or (end_3p - start_3p) < 3:
        return []
    else:
        term_runs = [(start_5p, end_5p), (start_3p, end_3p)]

    opts_map = get_loop_interaction_options(
        term_runs, full_seq, duplex_exe, temperature, extra_args, duplex_cache, L, dense_window=True
    )

    seeds = []
    if (0, 1) in opts_map:
        candidates = opts_map[(0, 1)]
        scaffold_indices = {x for p in scaffold_pairs for x in p}
        for cand, energy in candidates:
            cand_indices = {x for p in cand for x in p}
            if cand_indices.isdisjoint(scaffold_indices):
                seeds.append((cand, energy))

    seeds.sort(key=lambda x: len(x[0]), reverse=True)
    return seeds[:max_seeds]


# --- Consensus Masking Logic ---


def identify_helices(pairs: set[tuple[int, int]]) -> list[set[tuple[int, int]]]:
    """Group pairs into contiguous helices (stems)."""
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    helices = []

    if not sorted_pairs:
        return []

    current_helix = {sorted_pairs[0]}

    for k in range(1, len(sorted_pairs)):
        prev = sorted_pairs[k - 1]
        curr = sorted_pairs[k]

        # Check for stacking: (i, j) stacks on (i-1, j+1)
        if curr[0] == prev[0] + 1 and curr[1] == prev[1] - 1:
            current_helix.add(curr)
        else:
            helices.append(current_helix)
            current_helix = {curr}

    if current_helix:
        helices.append(current_helix)

    return helices


def generate_masked_scaffolds(
    initial_pairs: set[tuple[int, int]],
    seq_len: int,
    end_mask_step: int = DEFAULT_END_MASK_STEP,
    max_end_mask_len: int = DEFAULT_MAX_END_MASK_LEN,
    max_helices_sequential: int = DEFAULT_MAX_HELICES_SEQUENTIAL,
    max_helices_pairwise: int = DEFAULT_MAX_HELICES_PAIRWISE,
) -> list[set[tuple[int, int]]]:
    """
    Generate variants of the scaffold by systematically removing helices.
    Includes:
      1. Original scaffold.
      2. Scaffolds with 1 helix removed (Sequential).
      3. Scaffolds with 2 helices removed (Pairwise).
      4. End Masking Grid: Grid search removing 5' and 3' ends in steps.
    """
    variants = [initial_pairs]

    helices = identify_helices(initial_pairs)
    # Sort helices by length (descending) to target main structural elements
    helices.sort(key=len, reverse=True)

    # 1. Sequential Masking (Top N)
    limit_seq = min(len(helices), max_helices_sequential)
    for i in range(limit_seq):
        h = helices[i]
        # Create new scaffold without this helix
        new_scaffold = initial_pairs - h
        if new_scaffold not in variants:
            variants.append(new_scaffold)

    # 2. Pairwise Masking (Top M)
    limit_pair = min(len(helices), max_helices_pairwise)
    for i in range(limit_pair):
        for j in range(i + 1, limit_pair):
            h1 = helices[i]
            h2 = helices[j]
            new_scaffold = initial_pairs - h1 - h2
            if new_scaffold not in variants:
                variants.append(new_scaffold)

    # 3. End Masking Grid Search
    # Incrementally unpair 5' and 3' ends to force exploration of terminal structures
    # Loop 0, 5, 10, ... MAX_END_MASK_LEN
    end_steps_5 = list(range(0, max_end_mask_len + 1, end_mask_step))
    end_steps_3 = list(range(0, max_end_mask_len + 1, end_mask_step))

    for len_5 in end_steps_5:
        for len_3 in end_steps_3:
            if len_5 == 0 and len_3 == 0:
                continue

            # Safety check: don't unpair the whole molecule
            if (len_5 + len_3) >= seq_len:
                continue

            end_masked_scaffold = set()
            for i, j in initial_pairs:
                # Keep pair only if both bases are outside the mask windows
                # i and j are 0-indexed.
                # 5' mask removes 0 to len_5-1.
                # 3' mask removes (L - len_3) to L-1.
                cutoff_3 = seq_len - len_3

                in_5_region = (i < len_5) or (j < len_5)
                in_3_region = (i >= cutoff_3) or (j >= cutoff_3)

                if not in_5_region and not in_3_region:
                    end_masked_scaffold.add((i, j))

            if end_masked_scaffold not in variants:
                variants.append(end_masked_scaffold)

    return variants


# --- Kissing Loop Discovery ---


def find_hairpin_loops(pairs: set[tuple[int, int]], L: int) -> list[tuple[int, int]]:
    """
    Identify unpaired 'hairpin' loops.
    A hairpin loop is a region (i+1...j-1) that is completely unpaired,
    where (i, j) is a closing pair.
    """
    paired_pos = {p for pair in pairs for p in pair}
    hairpins = []

    # Iterate all pairs to see if they close an empty region
    for i, j in pairs:
        # Normalize i < j
        u, v = sorted((i, j))
        if (v - u) < 2:
            continue

        is_empty = True
        for k in range(u + 1, v):
            if k in paired_pos:
                is_empty = False
                break

        if is_empty:
            hairpins.append((u + 1, v))  # The loop region itself (start, end-exclusive)

    return hairpins


def check_complementarity_heuristic(seq1: str, seq2: str, min_len: int = 3) -> bool:
    """
    Fast check if two sequences have any contiguous complementary substring
    of at least min_len (handling standard WC+GU).
    """
    # Simple sliding window check
    # Check if a k-mer in seq1 (reversed) is comp to seq2
    # Because kissing loops are antiparallel, we check if seq1 is reverse-complementary to seq2
    len1, len2 = len(seq1), len(seq2)
    if len1 < min_len or len2 < min_len:
        return False

    for i in range(len1 - min_len + 1):
        kmer = seq1[i : i + min_len]
        # Generate target complement
        target = []
        for b in kmer:
            b = b.upper()
            if b == "A":
                target.append("U")
            elif b == "U":
                target.append("A")  # Matches A or G
            elif b == "G":
                target.append("C")  # Matches C or U
            elif b == "C":
                target.append("G")
            else:
                target.append(".")  # Wildcard?

        # Since it's antiparallel, an occurrence of 'target' in seq2 (read 5'-3')
        # implies a match.
        # But wobbles (G-U) make exact string matching hard.
        # Let's brute force checking validity instead of string search.

        for j in range(len2 - min_len + 1):
            sub2 = seq2[j : j + min_len]
            # Check pairing kmer (i) vs sub2 (j) antiparallel
            # kmer[0] pairs with sub2[end], etc.
            match = True
            for k in range(min_len):
                if not is_canonical(kmer[k], sub2[min_len - 1 - k]):
                    match = False
                    break
            if match:
                return True
    return False


def detect_and_add_kissing_loops(
    candidates: list[tuple[str, float]],
    full_seq: str,
    duplex_exe: Path,
    temperature: float | None,
    extra_args: Sequence[str] | None,
    duplex_cache: dict,
    kissing_loop_candidates: int = DEFAULT_KISSING_LOOP_CANDIDATES,
    scaffold_pair_energy: float = DEFAULT_SCAFFOLD_PAIR_ENERGY,
    weight_l0: float = DEFAULT_WEIGHT_L0,
    weight_l1: float = DEFAULT_WEIGHT_L1,
    weight_l2_plus: float = DEFAULT_WEIGHT_L2_PLUS,
) -> list[tuple[str, float]]:
    """
    Intercept the top N candidates, scan for interacting hairpin loops,
    validate with DuplexFold, and add new structures with updated scores.
    """
    # 1. Sort and slice
    sorted_cands = sorted(candidates, key=lambda x: x[1])
    top_candidates = sorted_cands[:kissing_loop_candidates]

    new_structures = []
    L = len(full_seq)

    sys.stderr.write(
        f"[REFINE] Checking top {len(top_candidates)} structures for kissing loops...\n"
    )

    for struct_str, base_score in top_candidates:
        pairs = _pairs_from_struct(struct_str)
        # Find hairpin loops (pure unpaired regions closed by a pair)
        loops = find_hairpin_loops(pairs, L)

        if len(loops) < 2:
            continue

        # Check pairs of loops
        for l_idx_a in range(len(loops)):
            for l_idx_b in range(l_idx_a + 1, len(loops)):
                start_a, end_a = loops[l_idx_a]
                start_b, end_b = loops[l_idx_b]

                # Check distance (sterics)
                if abs(start_a - start_b) < 4:
                    continue

                seq_a = full_seq[start_a:end_a]
                seq_b = full_seq[start_b:end_b]

                # Heuristic Filter
                if not check_complementarity_heuristic(seq_a, seq_b, min_len=3):
                    continue

                # Run DuplexFold
                dkey = (seq_a, seq_b, temperature, tuple(extra_args or []))
                if dkey in duplex_cache:
                    d_res = duplex_cache[dkey]
                else:
                    d_res = call_rnastructure_duplexfold(
                        duplex_exe, seq_a, seq_b, temperature, extra_args
                    )
                    duplex_cache[dkey] = d_res

                if not d_res:
                    continue

                for pair_list, d_energy in d_res:
                    if not pair_list:
                        continue

                    # Convert local loop coords to global
                    new_pairs = set()
                    for i_loc, j_loc in pair_list:
                        # DuplexFold returns 1-based index relative to input subseqs
                        gi = start_a + i_loc - 1
                        gj = start_b + j_loc - 1
                        if gi > gj:
                            gi, gj = gj, gi
                        new_pairs.add((gi, gj))

                    # Validate
                    canon_pairs = validate_and_filter_pairs(new_pairs, full_seq)
                    if not canon_pairs:
                        continue

                    # Ensure no overlap with existing structure (should be safe by def of hairpin loop, but double check)
                    existing_indices = {x for p in pairs for x in p}
                    cand_indices = {x for p in canon_pairs for x in p}
                    if not existing_indices.isdisjoint(cand_indices):
                        continue

                    combined_set = pairs.union(canon_pairs)

                    # Check complexity
                    layers = pairs_to_layers(list(combined_set))
                    if len(layers) > 3:
                        continue

                    new_str = pairs_to_rosetta_string(list(combined_set), L)

                    # --- SCORING ---
                    # We accept the 'base_score' as the current refined score.
                    # We add the new interaction score, weighted by its layer depth.
                    # This is an approximation since we don't have the original raw energy components separately.

                    # Map new pairs to layers
                    pair_to_layer = {}
                    for lay_idx, lay in enumerate(layers):
                        for p in lay:
                            pair_to_layer[p] = lay_idx

                    # Calculate additive score term
                    avg_duplex_energy = d_energy / len(canon_pairs)
                    added_score = 0.0
                    for p in canon_pairs:
                        layer = pair_to_layer.get(p, 0)
                        w = weight_l0
                        if layer == 1:
                            w = weight_l1
                        elif layer >= 2:
                            w = weight_l2_plus
                        added_score += w * avg_duplex_energy

                    new_total_score = base_score + added_score
                    new_structures.append((new_str, new_total_score))

    sys.stderr.write(f"[REFINE] Found {len(new_structures)} kissing-loop variants.\n")

    # Combine and sort
    all_results = candidates + new_structures
    all_results.sort(key=lambda x: x[1])

    # Deduplicate by structure string
    seen = set()
    unique = []
    for s, score in all_results:
        if s not in seen:
            seen.add(s)
            unique.append((s, score))

    return unique


# --- Main Refinement Logic ---


def _refine_structure_impl(
    struct: str,
    full_seq: str,
    allsub_exe: Path,
    duplex_exe: Path,
    min_unpaired_len: int,
    temperature: float | None = None,
    absolute_energy: float | None = None,
    percent_energy: float | None = None,
    extra_args: Sequence[str] | None = None,
    allsub_cache: dict = None,
    duplex_cache: dict = None,
    # Scriptable parameters
    end_mask_step: int = DEFAULT_END_MASK_STEP,
    max_end_mask_len: int = DEFAULT_MAX_END_MASK_LEN,
    max_helices_sequential: int = DEFAULT_MAX_HELICES_SEQUENTIAL,
    max_helices_pairwise: int = DEFAULT_MAX_HELICES_PAIRWISE,
    max_seeds: int = DEFAULT_MAX_SEEDS,
    max_regions_to_refine: int = DEFAULT_MAX_REGIONS_TO_REFINE,
    max_solutions: int = DEFAULT_MAX_SOLUTIONS,
    kissing_loop_candidates: int = DEFAULT_KISSING_LOOP_CANDIDATES,
    scaffold_pair_energy: float = DEFAULT_SCAFFOLD_PAIR_ENERGY,
    weight_l0: float = DEFAULT_WEIGHT_L0,
    weight_l1: float = DEFAULT_WEIGHT_L1,
    weight_l2_plus: float = DEFAULT_WEIGHT_L2_PLUS,
) -> list[tuple[str, float]]:
    L = len(full_seq)
    if len(struct) != L:
        raise ValueError("Structure/Sequence length mismatch")

    if allsub_cache is None:
        allsub_cache = {}
    if duplex_cache is None:
        duplex_cache = {}

    # 1. Parse and SANITIZE the initial input scaffold.
    raw_scaffold = _pairs_from_struct(struct)
    initial_scaffold = validate_and_filter_pairs(raw_scaffold, full_seq)

    # 2. GENERATE MASKED VARIANTS
    #    Create alternative consensus baselines by removing helices.
    masked_scaffolds = generate_masked_scaffolds(
        initial_scaffold,
        L,
        end_mask_step=end_mask_step,
        max_end_mask_len=max_end_mask_len,
        max_helices_sequential=max_helices_sequential,
        max_helices_pairwise=max_helices_pairwise,
    )
    sys.stderr.write(f"[REFINE] Generated {len(masked_scaffolds)} masked variants for structure.\n")

    final_unique_structs: list[tuple[str, float]] = []
    seen_global: dict[str, float] = {}

    # Iterate over every Masked Variant (Consensus perturbation)
    for mask_idx, base_scaffold in enumerate(masked_scaffolds):
        sys.stderr.write(
            f"[REFINE] Processing Mask {mask_idx + 1}/{len(masked_scaffolds)} (Pairs: {len(base_scaffold)})...\n"
        )

        # 3. GENERATE SEED SCAFFOLDS (Initialization) for THIS mask
        #    Explicitly calculate 5'-3' interactions.

        # Tuple: (scaffold_set, added_energy)
        # Base scaffold has 0.0 refined energy
        scaffolds_to_process = [(base_scaffold, 0.0)]
        seed_interactions = get_5p3p_seeds(
            base_scaffold,
            full_seq,
            duplex_exe,
            temperature,
            extra_args,
            duplex_cache,
            L,
            max_seeds=max_seeds,
        )
        sys.stderr.write(
            f"[REFINE]   Generated {len(seed_interactions)} end-to-end seeds for Mask {mask_idx + 1}.\n"
        )

        for seed, energy in seed_interactions:
            new_scaffold = base_scaffold.union(seed)
            scaffolds_to_process.append((new_scaffold, energy))

        # 4. Process each Seeded Scaffold independently
        for seed_idx, (current_scaffold, base_energy) in enumerate(scaffolds_to_process):
            # Recalculate runs based on the current scaffold
            current_struct_str = pairs_to_rosetta_string(list(current_scaffold), L)
            runs = find_unpaired_runs(current_struct_str, min_unpaired_len)

            # --- HEURISTIC: Limit to top N largest runs ---
            if len(runs) > max_regions_to_refine:
                sorted_runs = sorted(runs, key=lambda r: (r[1] - r[0]), reverse=True)
                runs = sorted_runs[:max_regions_to_refine]
                runs.sort(key=lambda r: r[0])

            sys.stderr.write(
                f"[REFINE]     Mask {mask_idx + 1}, Seed {seed_idx + 1}: Runs to refine = {len(runs)}\n"
            )

            if not runs:
                s_str = pairs_to_rosetta_string(list(current_scaffold), L)

                # Apply new scoring logic
                s_adj_score = calculate_weighted_score(
                    all_pairs=current_scaffold,
                    refined_energy=base_energy,
                    scaffold_pairs=base_scaffold,
                    scaffold_pair_energy=scaffold_pair_energy,
                    weight_l0=weight_l0,
                    weight_l1=weight_l1,
                    weight_l2_plus=weight_l2_plus,
                )

                if s_str not in seen_global or s_adj_score < seen_global[s_str]:
                    seen_global[s_str] = s_adj_score
                    final_unique_structs.append((s_str, s_adj_score))
                continue

            # --- Pre-computation for LOCAL Options ---
            local_options: dict[int, list[tuple[set[tuple[int, int]], float]]] = {}
            for idx, (start, end) in enumerate(runs):
                subseq = full_seq[start:end]
                key = (
                    subseq,
                    temperature,
                    absolute_energy,
                    percent_energy,
                    tuple(extra_args or []),
                )

                if key in allsub_cache:
                    sub_structs = allsub_cache[key]
                else:
                    sub_structs = call_rnastructure_allsub(
                        allsub_exe, subseq, temperature, absolute_energy, percent_energy, extra_args
                    )
                    allsub_cache[key] = sub_structs

                option_sets = []
                for s_str, energy in sub_structs:
                    p_sub = _struct_string_to_pairs(s_str, offset=start)
                    p_clean = validate_and_filter_pairs(p_sub, full_seq)
                    option_sets.append((p_clean, energy))

                if not option_sets:
                    option_sets = [(set(), 0.0)]
                option_sets = merge_disjoint_options(option_sets)
                local_options[idx] = option_sets

            # --- Pre-computation for INTERACTION Options ---
            interaction_options: dict[
                tuple[int, int], list[tuple[set[tuple[int, int]], float]]
            ] = {}
            if len(runs) >= 2:
                interaction_options = get_loop_interaction_options(
                    runs, full_seq, duplex_exe, temperature, extra_args, duplex_cache, L
                )

            # (pair_set, current_score)
            track1_results: list[tuple[set[tuple[int, int]], float]] = []
            track2_results: list[tuple[set[tuple[int, int]], float]] = []

            # TRACK 1: Enumerative Combinatorial
            sys.stderr.write(
                f"[REFINE]     Mask {mask_idx + 1}, Seed {seed_idx + 1}: Starting Track 1 (Combinatorial)...\n"
            )
            num_regions = len(runs)

            def _backtrack_track1(
                idx: int,
                covered: set[int],
                current_pairs: set[tuple[int, int]],
                current_score: float,
            ):
                if len(track1_results) >= max_solutions:
                    return
                if idx >= num_regions:
                    track1_results.append((current_pairs, current_score))
                    return
                if idx in covered:
                    _backtrack_track1(idx + 1, covered, current_pairs, current_score)
                    return

                # Interaction
                for k in range(idx + 1, num_regions):
                    if k in covered:
                        continue
                    opts = interaction_options.get((idx, k))
                    if opts:
                        new_covered_int = covered.copy()
                        new_covered_int.add(idx)
                        new_covered_int.add(k)
                        for int_set, int_energy in opts:
                            if len(track1_results) >= max_solutions:
                                return
                            used_indices = {x for p in int_set for x in p}

                            # Compatible local
                            compatible_idx = [
                                l
                                for l in local_options[idx]
                                if {x for p in l[0] for x in p}.isdisjoint(used_indices)
                            ] or [(set(), 0.0)]
                            compatible_k = [
                                l
                                for l in local_options[k]
                                if {x for p in l[0] for x in p}.isdisjoint(used_indices)
                            ] or [(set(), 0.0)]

                            for c_i, e_i in compatible_idx:
                                for c_k, e_k in compatible_k:
                                    if len(track1_results) >= max_solutions:
                                        return
                                    combined_set = int_set.union(c_i).union(c_k)
                                    _backtrack_track1(
                                        idx + 1,
                                        new_covered_int,
                                        current_pairs.union(combined_set),
                                        current_score + int_energy + e_i + e_k,
                                    )

                # Local
                new_covered_local = covered.copy()
                new_covered_local.add(idx)
                for loc_set, loc_energy in local_options[idx]:
                    if len(track1_results) >= max_solutions:
                        return
                    _backtrack_track1(
                        idx + 1,
                        new_covered_local,
                        current_pairs.union(loc_set),
                        current_score + loc_energy,
                    )

            _backtrack_track1(0, set(), current_scaffold, base_energy)
            sys.stderr.write(
                f"[REFINE]     Mask {mask_idx + 1}, Seed {seed_idx + 1}: Finished Track 1. Found {len(track1_results)} solutions.\n"
            )

            # TRACK 2: Enumerative Sequential
            sys.stderr.write(
                f"[REFINE]     Mask {mask_idx + 1}, Seed {seed_idx + 1}: Starting Track 2 (Sequential)...\n"
            )
            all_local_lists = [local_options[i] for i in range(num_regions)]

            # Using Cartesian product
            for local_combo in itertools.product(*all_local_lists):
                if len(track2_results) >= max_solutions:
                    break

                # Unpack (set, energy) tuples
                base_set = current_scaffold.copy()
                base_score = base_energy
                for p_set, en in local_combo:
                    base_set.update(p_set)
                    base_score += en

                paired_mask = [False] * L
                for i, j in base_set:
                    if 0 <= i < L:
                        paired_mask[i] = True
                    if 0 <= j < L:
                        paired_mask[j] = True

                sub_runs = []
                i = 0
                while i < L:
                    if not paired_mask[i]:
                        start_i = i
                        while i < L and not paired_mask[i]:
                            i += 1
                        if (i - start_i) >= 3:
                            sub_runs.append((start_i, i))
                    else:
                        i += 1

                if len(sub_runs) > max_regions_to_refine:
                    sorted_sub = sorted(sub_runs, key=lambda r: (r[1] - r[0]), reverse=True)
                    sub_runs = sorted_sub[:max_regions_to_refine]
                    sub_runs.sort(key=lambda r: r[0])

                if len(sub_runs) < 2:
                    track2_results.append((base_set, base_score))
                    continue

                loop_int_options = get_loop_interaction_options(
                    sub_runs, full_seq, duplex_exe, temperature, extra_args, duplex_cache, L
                )
                loop_local_options = {}
                for l_idx, (start, end) in enumerate(sub_runs):
                    subseq = full_seq[start:end]
                    key = (
                        subseq,
                        temperature,
                        absolute_energy,
                        percent_energy,
                        tuple(extra_args or []),
                    )
                    if key in allsub_cache:
                        sub_structs = allsub_cache[key]
                    else:
                        sub_structs = call_rnastructure_allsub(
                            allsub_exe,
                            subseq,
                            temperature,
                            absolute_energy,
                            percent_energy,
                            extra_args,
                        )
                        allsub_cache[key] = sub_structs
                    l_opt_sets = []
                    for s_str, energy in sub_structs:
                        p_sub = _struct_string_to_pairs(s_str, offset=start)
                        p_clean = validate_and_filter_pairs(p_sub, full_seq)
                        l_opt_sets.append((p_clean, energy))
                    if not l_opt_sets:
                        l_opt_sets = [(set(), 0.0)]
                    loop_local_options[l_idx] = l_opt_sets

                num_loops = len(sub_runs)

                def _backtrack_track2_loops(
                    l_idx: int, l_covered: set[int], l_pairs: set[tuple[int, int]], l_score: float
                ):
                    if len(track2_results) >= max_solutions:
                        return
                    if l_idx >= num_loops:
                        track2_results.append((l_pairs, l_score))
                        return
                    if l_idx in l_covered:
                        _backtrack_track2_loops(l_idx + 1, l_covered, l_pairs, l_score)
                        return

                    for k in range(l_idx + 1, num_loops):
                        if k in l_covered:
                            continue
                        opts = loop_int_options.get((l_idx, k))
                        if opts:
                            new_cov_pair = l_covered.copy()
                            new_cov_pair.add(l_idx)
                            new_cov_pair.add(k)
                            for int_set, int_en in opts:
                                if len(track2_results) >= max_solutions:
                                    return
                                _backtrack_track2_loops(
                                    l_idx + 1,
                                    new_cov_pair,
                                    l_pairs.union(int_set),
                                    l_score + int_en,
                                )

                    new_cov_local = l_covered.copy()
                    new_cov_local.add(l_idx)
                    for l_set, l_en in loop_local_options[l_idx]:
                        if len(track2_results) >= max_solutions:
                            return
                        _backtrack_track2_loops(
                            l_idx + 1, new_cov_local, l_pairs.union(l_set), l_score + l_en
                        )

                _backtrack_track2_loops(0, set(), base_set, base_score)

            sys.stderr.write(
                f"[REFINE]     Mask {mask_idx + 1}, Seed {seed_idx + 1}: Finished Track 2. Found {len(track2_results)} solutions.\n"
            )

            pair_sets = track1_results + track2_results
            for p_set, p_score in pair_sets:
                cleaned_set = remove_isolated_pairs(p_set)

                flattened = []
                valid_struct = True
                for i, j in cleaned_set:
                    if i < 0 or i >= L or j < 0 or j >= L:
                        valid_struct = False
                        break
                    flattened.append(i)
                    flattened.append(j)
                if not valid_struct:
                    continue
                if len(flattened) != len(set(flattened)):
                    continue

                layers = pairs_to_layers(list(cleaned_set))
                if len(layers) > 3:
                    continue

                pk_str = pairs_to_rosetta_string(sorted(cleaned_set), L)

                # Apply new scoring logic
                adj_score = calculate_weighted_score(
                    all_pairs=cleaned_set,
                    refined_energy=p_score,
                    scaffold_pairs=base_scaffold,
                    scaffold_pair_energy=scaffold_pair_energy,
                    weight_l0=weight_l0,
                    weight_l1=weight_l1,
                    weight_l2_plus=weight_l2_plus,
                )

                if pk_str not in seen_global or adj_score < seen_global[pk_str]:
                    seen_global[pk_str] = adj_score
                    final_unique_structs.append((pk_str, adj_score))

    # --- STEP 4: Post-Process Kissing Loop Discovery ---
    # Intercept top candidates and look for interacting hairpin loops
    # limiting to top N to save runtime.

    enhanced_structs = detect_and_add_kissing_loops(
        candidates=final_unique_structs,
        full_seq=full_seq,
        duplex_exe=duplex_exe,
        temperature=temperature,
        extra_args=extra_args,
        duplex_cache=duplex_cache,
        kissing_loop_candidates=kissing_loop_candidates,
        scaffold_pair_energy=scaffold_pair_energy,
        weight_l0=weight_l0,
        weight_l1=weight_l1,
        weight_l2_plus=weight_l2_plus,
    )

    return enhanced_structs


def refine_structure(
    struct: str,
    full_seq: str,
    allsub_exe: Path,
    duplex_exe: Path,
    min_unpaired_len: int,
    temperature: float | None = None,
    absolute_energy: float | None = None,
    percent_energy: float | None = None,
    extra_args: Sequence[str] | None = None,
    allsub_cache: dict = None,
    duplex_cache: dict = None,
    return_all: bool = False,
    # Scriptable parameters
    end_mask_step: int = DEFAULT_END_MASK_STEP,
    max_end_mask_len: int = DEFAULT_MAX_END_MASK_LEN,
    max_helices_sequential: int = DEFAULT_MAX_HELICES_SEQUENTIAL,
    max_helices_pairwise: int = DEFAULT_MAX_HELICES_PAIRWISE,
    max_seeds: int = DEFAULT_MAX_SEEDS,
    max_regions_to_refine: int = DEFAULT_MAX_REGIONS_TO_REFINE,
    max_solutions: int = DEFAULT_MAX_SOLUTIONS,
    kissing_loop_candidates: int = DEFAULT_KISSING_LOOP_CANDIDATES,
    scaffold_pair_energy: float = DEFAULT_SCAFFOLD_PAIR_ENERGY,
    weight_l0: float = DEFAULT_WEIGHT_L0,
    weight_l1: float = DEFAULT_WEIGHT_L1,
    weight_l2_plus: float = DEFAULT_WEIGHT_L2_PLUS,
) -> list[tuple[str, float]]:
    return _refine_structure_impl(
        struct=struct,
        full_seq=full_seq,
        allsub_exe=allsub_exe,
        duplex_exe=duplex_exe,
        min_unpaired_len=min_unpaired_len,
        temperature=temperature,
        absolute_energy=absolute_energy,
        percent_energy=percent_energy,
        extra_args=extra_args,
        allsub_cache=allsub_cache,
        duplex_cache=duplex_cache,
        end_mask_step=end_mask_step,
        max_end_mask_len=max_end_mask_len,
        max_helices_sequential=max_helices_sequential,
        max_helices_pairwise=max_helices_pairwise,
        max_seeds=max_seeds,
        max_regions_to_refine=max_regions_to_refine,
        max_solutions=max_solutions,
        kissing_loop_candidates=kissing_loop_candidates,
        scaffold_pair_energy=scaffold_pair_energy,
        weight_l0=weight_l0,
        weight_l1=weight_l1,
        weight_l2_plus=weight_l2_plus,
    )


def filter_and_rank_structures(
    structs: list[tuple[str, float]], max_output: int = 1000, diversity_dist_frac: float = 0.05
) -> list[str]:
    # 1. Deduplicate keeping best score
    best_scores = {}
    for s, score in structs:
        if s not in best_scores or score < best_scores[s]:
            best_scores[s] = score

    # Convert to list and SORT by score (ascending: lower energy is better)
    unique = sorted(best_scores.items(), key=lambda x: x[1])

    if not unique:
        return []

    L = len(unique[0][0])
    threshold = max(1, int(round(diversity_dist_frac * L)))

    kept = []

    def hamming_dist(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    for s, score in unique:
        if len(kept) >= max_output:
            break
        is_far_enough = True
        for k in kept:
            if hamming_dist(s, k) <= threshold:
                is_far_enough = False
                break
        if is_far_enough:
            kept.append(s)

    return kept


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Refine unpaired regions using AllSub and DuplexFold."
    )
    parser.add_argument("--sto", help="Stockholm file for sequence.")
    parser.add_argument("--seq-file", help="FASTA/SEQ file for sequence.")
    parser.add_argument("--seq-name", default=None)
    parser.add_argument("--db-in", required=True)
    parser.add_argument(
        "--db-out", required=True, help="Path for filtered structures (default output)."
    )
    parser.add_argument(
        "--db-out-all",
        default=None,
        help="Optional path to save ALL refined structures (before filtering).",
    )

    parser.add_argument("--allsub-exe", default="$PROJECT/repos/RNAstructure/exe/AllSub")
    parser.add_argument("--duplex-exe", default="$PROJECT/repos/RNAstructure/exe/DuplexFold")

    parser.add_argument("--min-unpaired", type=int, default=4)

    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--allsub-abs", type=float, default=None, help="-a flag for AllSub")
    parser.add_argument("--allsub-pct", type=float, default=None, help="-p flag for AllSub")

    parser.add_argument(
        "--max-structures",
        type=int,
        default=1000,
        help="Maximum number of refined structures to output (default: 1000).",
    )
    parser.add_argument(
        "--diversity-dist",
        type=float,
        default=0.05,
        help="Minimum Hamming distance fraction (0.0-1.0) for diversity clustering.",
    )

    parser.add_argument("--extra-arg", action="append", default=[])

    # Scriptable parameters exposed to CLI
    parser.add_argument("--end-mask-step", type=int, default=DEFAULT_END_MASK_STEP)
    parser.add_argument("--max-end-mask-len", type=int, default=DEFAULT_MAX_END_MASK_LEN)
    parser.add_argument(
        "--max-helices-sequential", type=int, default=DEFAULT_MAX_HELICES_SEQUENTIAL
    )
    parser.add_argument("--max-helices-pairwise", type=int, default=DEFAULT_MAX_HELICES_PAIRWISE)
    parser.add_argument("--max-seeds", type=int, default=DEFAULT_MAX_SEEDS)
    parser.add_argument("--max-regions-to-refine", type=int, default=DEFAULT_MAX_REGIONS_TO_REFINE)
    parser.add_argument("--max-solutions", type=int, default=DEFAULT_MAX_SOLUTIONS)
    parser.add_argument(
        "--kissing-loop-candidates", type=int, default=DEFAULT_KISSING_LOOP_CANDIDATES
    )
    parser.add_argument("--scaffold-pair-energy", type=float, default=DEFAULT_SCAFFOLD_PAIR_ENERGY)
    parser.add_argument("--weight-l0", type=float, default=DEFAULT_WEIGHT_L0)
    parser.add_argument("--weight-l1", type=float, default=DEFAULT_WEIGHT_L1)
    parser.add_argument("--weight-l2-plus", type=float, default=DEFAULT_WEIGHT_L2_PLUS)

    args = parser.parse_args(argv)

    allsub_exe = Path(os.path.expandvars(args.allsub_exe))
    duplex_exe = Path(os.path.expandvars(args.duplex_exe))

    if not allsub_exe.is_file():
        sys.exit(f"[ERROR] AllSub not found: {allsub_exe}")
    if not duplex_exe.is_file():
        sys.exit(f"[ERROR] DuplexFold not found: {duplex_exe}")

    if args.sto:
        full_seq = read_ungapped_seq_from_sto(Path(args.sto), args.seq_name)
    elif args.seq_file:
        full_seq = read_fasta_sequence(Path(args.seq_file), args.seq_name)
    else:
        sys.exit("[ERROR] Must provide --sto or --seq-file")

    structs = read_db_structures(Path(args.db_in))

    allsub_cache = {}
    duplex_cache = {}

    collect_all_mode = args.db_out_all is not None
    all_structs_for_output: list[tuple[str, float]] = []

    for idx, s in enumerate(structs):
        sys.stderr.write(
            f"[INFO] Structure {idx + 1}/{len(structs)}: Generating consensus variants...\n"
        )
        variants = refine_structure(
            struct=s,
            full_seq=full_seq,
            allsub_exe=allsub_exe,
            duplex_exe=duplex_exe,
            min_unpaired_len=args.min_unpaired,
            temperature=args.temperature,
            absolute_energy=args.allsub_abs,
            percent_energy=args.allsub_pct,
            extra_args=args.extra_arg,
            allsub_cache=allsub_cache,
            duplex_cache=duplex_cache,
            return_all=collect_all_mode,
            # Pass through new args
            end_mask_step=args.end_mask_step,
            max_end_mask_len=args.max_end_mask_len,
            max_helices_sequential=args.max_helices_sequential,
            max_helices_pairwise=args.max_helices_pairwise,
            max_seeds=args.max_seeds,
            max_regions_to_refine=args.max_regions_to_refine,
            max_solutions=args.max_solutions,
            kissing_loop_candidates=args.kissing_loop_candidates,
            scaffold_pair_energy=args.scaffold_pair_energy,
            weight_l0=args.weight_l0,
            weight_l1=args.weight_l1,
            weight_l2_plus=args.weight_l2_plus,
        )
        all_structs_for_output.extend(variants)
        sys.stderr.write(f"[INFO] Structure {idx + 1} yielded {len(variants)} models.\n")

    if args.db_out_all:
        sys.stderr.write(
            f"[INFO] Writing ALL {len(all_structs_for_output)} structures to {args.db_out_all}\n"
        )

        # Sort by score for the "all" output as well
        sorted_all = sorted(all_structs_for_output, key=lambda x: x[1])

        with Path(args.db_out_all).open("w") as f:
            f.write(full_seq + "\n")
            seen_full = set()
            for rs, score in sorted_all:
                if rs not in seen_full:
                    # Optional: Could write score as comment if DB format supported it,
                    # but spec says read-only DBs. Just writing structure.
                    f.write(rs + "\n")
                    seen_full.add(rs)

    final_structs = filter_and_rank_structures(
        all_structs_for_output,
        max_output=args.max_structures,
        diversity_dist_frac=args.diversity_dist,
    )

    with Path(args.db_out).open("w") as f:
        f.write(full_seq + "\n")
        for rs in final_structs:
            f.write(rs + "\n")

    sys.stderr.write(f"[INFO] Wrote {len(final_structs)} filtered structures to {args.db_out}\n")


if __name__ == "__main__":
    main()
