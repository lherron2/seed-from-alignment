"""
Refine long unpaired regions in CaCoFold-sampled structures using an enumerative combinatorial approach.

Refinement Tracks:
    1. Independent Combinatorial Assembly (Enhanced for Kissing Loops):
       - Input: Original unpaired regions (runs).
       - Candidates: 
            a) Local folding of run i (AllSub).
            b) Interaction between run i and run j (DuplexFold).
       - Kissing Loop Logic: When selecting an interaction (b), the algorithm now actively
         looks for compatible local structures (a) that can exist simultaneously (disjoint indices),
         allowing for stem-loops that also participate in tertiary contacts.
       - Algorithm: Recursive backtracking to generate ALL valid assignments.
       - PRIORITY: Interactions are explored BEFORE local folds to bias results toward
         complex topologies (e.g. simultaneous PKs).
       - MERGING: Disjoint local structures for the same region are combined.

    2. Sequential Hierarchical Assembly:
       - Step A: Force all original regions to fold locally. Enumerate all combinations of AllSub structures.
       - Step B: For each combination (Base Structure):
            - Identify REMAINING unpaired segments (loops/bulges), min_len=4.
            - Filter pairs of loops by linker length (>= 3nt).
            - Run DuplexFold on valid pairs.
            - Run AllSub on loops to find missed local helices.
            - Algorithm: Recursive backtracking to enumerate ALL valid sets of non-overlapping
              loop-loop interactions AND loop local folds.
       - Result: Base Structure + Selected Loop Interactions/Folds.

    - Validation & Filtering:
       - STRICT CANONICAL PAIR CHECK: All pairs must be WC or GU.
       - ROSETTA FORMAT: Use () [] {} aA bB... to prevent a...a pairing errors.
       - Remove isolated base pairs.
       - Discard structures with > 2 layers of crossings (total layers > 3).
"""

from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import sys
import tempfile
import string
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional, Set
from itertools import combinations

# Re-use logic from sample_cacofold_structures where appropriate
from src.lib.sample_cacofold_structures import pairs_from_track

UNPAIRED_CHAR = "."
GAP_CHARS = set("-._~")
CANONICAL_PAIRS = {
    ("A", "U"), ("U", "A"),
    ("G", "C"), ("C", "G"),
    ("G", "U"), ("U", "G"),
}
UNPAIR_ENDS_LEN = 25       # Length of 5'/3' ends to force unpaired in that specific variant
MAX_REGIONS_TO_REFINE = 50  # Safety cap: Only refine the N longest runs to prevent explosion
MAX_SOLUTIONS = 20000      # Applied PER TRACK (so up to 2x this total before filtering)
MAX_SCAFFOLD_VARIANTS = 50 # Limit how many consensus masks we try

# --- Graph Coloring for Layers (Local Definition for Safety) ---

def pairs_cross(p: Tuple[int, int], q: Tuple[int, int]) -> bool:
    """Return True if arcs (i, j) and (k, l) cross."""
    i, j = sorted(p)
    k, l = sorted(q)
    return (i < k < j < l) or (k < i < l < j)

def pairs_to_layers(pairs: List[Tuple[int, int]]) -> List[List[Tuple[int, int]]]:
    """Partition pairs into layers such that no pairs in a layer cross."""
    layers: List[List[Tuple[int, int]]] = []
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

def pairs_to_rosetta_string(pairs: List[Tuple[int, int]], L: int) -> str:
    """
    Convert pairs to a dot-bracket string using Rosetta-safe hierarchy:
    Layer 0: ()
    Layer 1: []
    Layer 2: {}
    Layer 3+: aA, bB ...
    """
    # Sanity check: Ensure pairs form a valid matching (no shared indices)
    indices = [idx for pair in pairs for idx in pair]
    if len(indices) != len(set(indices)):
        raise ValueError("Invalid matching: overlapping pairs detected in refined structure.")

    chars = ["."] * L
    layers = pairs_to_layers(pairs)
    
    brackets = [("(", ")"), ("[", "]"), ("{", "}")]
    # Add a-z/A-Z for deeper layers
    for i in range(26):
        brackets.append((chr(ord('a')+i), chr(ord('A')+i)))
        
    for layer_idx, layer in enumerate(layers):
        if layer_idx < len(brackets):
            op, cl = brackets[layer_idx]
        else:
            op, cl = "{", "}" # Fallback
            
        for i, j in layer:
            chars[i] = op
            chars[j] = cl
    return "".join(chars)

# --- Validators ---

def is_canonical(b1: str, b2: str) -> bool:
    """Check if two bases form a canonical Watson-Crick or Wobble pair."""
    return (b1.upper(), b2.upper()) in CANONICAL_PAIRS

def validate_and_filter_pairs(pairs: Set[Tuple[int, int]], seq: str) -> Set[Tuple[int, int]]:
    """Return only pairs that are valid indices and canonical matches."""
    valid = set()
    L = len(seq)
    for i, j in pairs:
        if 0 <= i < L and 0 <= j < L:
            if is_canonical(seq[i], seq[j]):
                valid.add((i, j))
    return valid

def remove_isolated_pairs(pairs: Set[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    """Iteratively remove pairs that have no stacking neighbor (isolated)."""
    current_pairs = pairs.copy()
    while True:
        to_remove = set()
        for i, j in current_pairs:
            # Check for stack: (i-1, j+1) or (i+1, j-1)
            p_down = (i - 1, j + 1)
            p_up = (i + 1, j - 1)
            
            # Helper to check existence regardless of order
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
    """Read the ungapped RNA sequence for ``seq_name`` from a Stockholm file."""
    sequences: Dict[str, List[str]] = {}
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
            available = ", ".join(sorted(sequences.keys()))
            raise ValueError(f"Sequence '{seq_name}' not found. Available: {available}")
        name = seq_name

    aligned = "".join(sequences[name])
    ungapped = "".join(ch for ch in aligned if ch not in GAP_CHARS)
    return ungapped.upper().replace("T", "U")


def read_fasta_sequence(path: Path, seq_name: str | None = None) -> str:
    """Read the first (or named) sequence from FASTA/SEQ and normalise to RNA."""
    seqs: Dict[str, str] = {}
    current_name: str | None = None
    chunks: List[str] = []

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


def read_db_structures(path: Path) -> List[str]:
    """Read structures from a .db file."""
    structs: List[str] = []
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


def find_unpaired_runs(struct: str, min_len: int) -> List[Tuple[int, int]]:
    """Find contiguous runs of '.' in the structure with length >= min_len."""
    runs: List[Tuple[int, int]] = []
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


def parse_ct_file(ct_path: Path) -> List[str]:
    """Parse an RNAstructure CT file containing one or more structures."""
    structures: List[str] = []
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
            sys.stderr.write(f"[WARN] CT file truncated at structure header: {header}\n")
            break
            
        db_chars = ["."] * seq_len
        
        for offset in range(1, seq_len + 1):
            line = lines[idx + offset]
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                i = int(fields[0])
                j = int(fields[4])
            except ValueError:
                continue
            
            i -= 1
            j -= 1
            
            if j > i:
                if 0 <= i < seq_len and 0 <= j < seq_len:
                    db_chars[i] = "("
                    db_chars[j] = ")"

        structures.append("".join(db_chars))
        idx += (seq_len + 1)

    return structures


# --- RNAstructure Wrappers ---

def call_rnastructure_allsub(
    allsub_exe: Path,
    subseq: str,
    temperature: float | None = None,
    absolute_energy: float | None = None,
    percent_energy: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> List[str]:
    """Call RNAstructure AllSub to find suboptimal structures."""
    extra_args = list(extra_args) if extra_args is not None else []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        seq_path = tmpdir_path / "subseq.fa"
        ct_path = tmpdir_path / "output.ct"

        with seq_path.open("w") as fh:
            fh.write(">subseq\n")
            fh.write(subseq + "\n")

        cmd: List[str] = [str(allsub_exe), str(seq_path), str(ct_path)]
        if temperature is not None:
            cmd.extend(["-t", str(temperature)]) 
        if absolute_energy is not None:
            cmd.extend(["-a", str(absolute_energy)])
        if percent_energy is not None:
            cmd.extend(["-p", str(percent_energy)])
        cmd.extend(extra_args)

        try:
            subprocess.run(cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                f"[ERROR] RNAstructure AllSub failed (exit {e.returncode})\n"
                f"Command: {' '.join(cmd)}\n{e.stderr}\n"
            )
            return ["." * len(subseq)]

        if not ct_path.is_file():
             return ["." * len(subseq)]

        structures = parse_ct_file(ct_path)
        if not structures:
            return ["." * len(subseq)]
        return structures


def call_rnastructure_duplexfold(
    duplex_exe: Path,
    seq5: str,
    seq3: str,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> List[List[Tuple[int, int]]]:
    """
    Call RNAstructure DuplexFold on two sequences.
    Returns list of pair sets. Pairs are (i_local_seq5, j_local_seq3) 1-based.
    Includes bounds checking.
    """
    extra_args = list(extra_args) if extra_args is not None else []

    # Ensure exploration of suboptimals to catch transient interactions (like 5'-3' kisses)
    # Only add -m 100 if user hasn't specified -m in extra_args (allows user control)
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

        cmd: List[str] = [str(duplex_exe), str(seq5_path), str(seq3_path), str(ct_path)]
        if temperature is not None:
            cmd += ["-t", str(temperature)]
        cmd += extra_args

        try:
            subprocess.run(cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            return []

        structures = parse_ct_file(ct_path)
        if not structures:
            return []

        # DuplexFold usually inserts a 3-nt linker (e.g. 'III') between sequences in the CT file.
        # We detect this dynamically.
        seq_len_ct = len(structures[0])
        all_pair_sets = []
        n1 = len(seq5)
        n2 = len(seq3)
        linker_len = max(0, seq_len_ct - (n1 + n2))
        
        for db in structures:
            pairs = []
            stack = []
            for i, ch in enumerate(db):
                if ch == '(':
                    stack.append(i)
                elif ch == ')':
                    if stack:
                        j = stack.pop()
                        u, v = sorted((j, i)) 
                        # u must be in seq5 (0..n1-1)
                        # v must be in seq3 (n1+linker..end)
                        if u < n1 and v >= (n1 + linker_len):
                            i5_local = u + 1
                            j3_local = v - (n1 + linker_len) + 1
                            if i5_local <= n1 and j3_local <= n2:
                                pairs.append((i5_local, j3_local))
            all_pair_sets.append(pairs)
        return all_pair_sets


def _pairs_from_struct(struct: str) -> Set[Tuple[int, int]]:
    return set(pairs_from_track(struct))


def _struct_string_to_pairs(s: str, offset: int = 0) -> Set[Tuple[int, int]]:
    pairs = set()
    stack = []
    for i, ch in enumerate(s):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
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
    L
) -> Dict[Tuple[int, int], List[Set[Tuple[int, int]]]]:
    """Helper to compute pairwise DuplexFold options for a set of loops."""
    options = {}
    for r_idx_a in range(len(sub_runs)):
        for r_idx_b in range(r_idx_a + 1, len(sub_runs)):
            sa, ea = sub_runs[r_idx_a]
            sb, eb = sub_runs[r_idx_b]
            
            # LINKER FILTER:
            # Distance between end of A and start of B
            dist = sb - ea
            if dist < 3:
                continue
            
            seq_a = full_seq[sa:ea]
            seq_b = full_seq[sb:eb]
            
            dkey = (seq_a, seq_b, temperature, tuple(extra_args or []))
            if dkey in duplex_cache:
                d_res = duplex_cache[dkey]
            else:
                d_res = call_rnastructure_duplexfold(
                    duplex_exe, seq_a, seq_b, temperature, extra_args
                )
                duplex_cache[dkey] = d_res
            
            if d_res:
                valid_int_sets = []
                for pair_list in d_res:
                    if pair_list:
                        g_pairs = set()
                        for i_loc, j_loc in pair_list:
                            gi = sa + i_loc - 1
                            gj = sb + j_loc - 1
                            if gi > gj: gi, gj = gj, gi
                            # No strict bounds check needed here as derived from valid slices,
                            # but filtering happens later.
                            g_pairs.add((gi, gj))
                        
                        # Validate CANONICAL immediately
                        canon_pairs = validate_and_filter_pairs(g_pairs, full_seq)
                        if canon_pairs:
                            valid_int_sets.append(canon_pairs)

                if valid_int_sets:
                    options[(r_idx_a, r_idx_b)] = valid_int_sets
    return options

def merge_disjoint_options(option_sets: List[Set[Tuple[int, int]]], max_combo=50) -> List[Set[Tuple[int, int]]]:
    """
    If AllSub returns multiple suboptimal structures for a single region,
    some might be disjoint (non-overlapping). We merge them here so we don't
    miss the combination of 'Helix A' and 'Helix B' if AllSub only returned them separately.
    """
    combined = list(option_sets)
    n = len(option_sets)
    # Simple N^2 pass to merge disjoint pairs
    # Limit to top options to prevent explosion
    limit = min(n, 20)
    
    new_combos = []
    for i in range(limit):
        for j in range(i + 1, limit):
            s1 = option_sets[i]
            s2 = option_sets[j]
            
            # Check overlap
            # Indices involved
            idx1 = {x for p in s1 for x in p}
            idx2 = {x for p in s2 for x in p}
            
            if idx1.isdisjoint(idx2):
                union_set = s1.union(s2)
                # Check if this union is already in the list
                if union_set not in combined and union_set not in new_combos:
                    new_combos.append(union_set)
                    if len(new_combos) >= max_combo: break
        if len(new_combos) >= max_combo: break
    
    combined.extend(new_combos)
    return combined

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
    allsub_cache: Dict = None,
    duplex_cache: Dict = None,
) -> List[str]:
    """
    Refine structure using two distinct tracks with heuristic filters.
    Core implementation.
    """
    L = len(full_seq)
    if len(struct) != L:
        raise ValueError("Structure/Sequence length mismatch")

    if allsub_cache is None: allsub_cache = {}
    if duplex_cache is None: duplex_cache = {}

    runs = find_unpaired_runs(struct, min_unpaired_len)
    
    # --- HEURISTIC: Limit to top N largest runs ---
    # Combinatorial complexity is exponential in len(runs).
    # If we have too many, sort by length and keep only the largest ones.
    if len(runs) > MAX_REGIONS_TO_REFINE:
        # Sort by length (descending)
        sorted_runs = sorted(runs, key=lambda r: (r[1] - r[0]), reverse=True)
        # Keep top N
        runs = sorted_runs[:MAX_REGIONS_TO_REFINE]
        # Re-sort by position to ensure logic flows left-to-right
        runs.sort(key=lambda r: r[0])
    
    # 1. Parse and SANITIZE the scaffold. 
    raw_scaffold = _pairs_from_struct(struct)
    scaffold_pairs = validate_and_filter_pairs(raw_scaffold, full_seq)

    if not runs:
        # Just return sanitized scaffold
        return [pairs_to_rosetta_string(list(scaffold_pairs), L)]

    # --- Pre-computation for LOCAL Options (Original Runs) ---
    local_options: Dict[int, List[Set[Tuple[int, int]]]] = {}
    
    for idx, (start, end) in enumerate(runs):
        subseq = full_seq[start:end]
        key = (subseq, temperature, absolute_energy, percent_energy, tuple(extra_args or []))
        
        if key in allsub_cache:
            sub_structs = allsub_cache[key]
        else:
            sub_structs = call_rnastructure_allsub(
                allsub_exe, subseq, temperature, absolute_energy, percent_energy, extra_args
            )
            allsub_cache[key] = sub_structs
        
        option_sets = []
        for s_str in sub_structs:
            # Parse pairs from AllSub
            p_sub = _struct_string_to_pairs(s_str, offset=start)
            # Filter non-canonical
            p_clean = validate_and_filter_pairs(p_sub, full_seq)
            option_sets.append(p_clean)
        
        # Ensure at least one option (empty set)
        if not option_sets: option_sets = [set()]
        
        # MERGE DISJOINT OPTIONS: If "Helix A" and "Helix B" appear separately, combine them.
        option_sets = merge_disjoint_options(option_sets)
        
        local_options[idx] = option_sets

    # --- Pre-computation for INTERACTION Options (Original Runs) ---
    interaction_options: Dict[Tuple[int, int], List[Set[Tuple[int, int]]]] = {}

    if len(runs) >= 2:
        interaction_options = get_loop_interaction_options(
            runs, full_seq, duplex_exe, temperature, extra_args, duplex_cache, L
        )

    # Separate Lists for Tracks
    track1_results: List[Set[Tuple[int, int]]] = []
    track2_results: List[Set[Tuple[int, int]]] = []

    # =========================================================================
    # TRACK 1: Enumerative Combinatorial (Independent)
    # =========================================================================
    num_regions = len(runs)
    
    def _backtrack_track1(idx: int, covered: Set[int], current_pairs: Set[Tuple[int, int]]):
        if len(track1_results) >= MAX_SOLUTIONS:
            return

        if idx >= num_regions:
            track1_results.append(current_pairs)
            return

        if idx in covered:
            _backtrack_track1(idx + 1, covered, current_pairs)
            return

        # STRATEGY CHANGE: Prioritize Interactions over Local Folds.
        # Choice A: Interaction with future region k (Hybrid: Kiss + Compatible Local)
        for k in range(idx + 1, num_regions):
            if k in covered: continue
            
            opts = interaction_options.get((idx, k))
            if opts:
                new_covered_int = covered.copy()
                new_covered_int.add(idx)
                new_covered_int.add(k)
                for int_set in opts:
                    if len(track1_results) >= MAX_SOLUTIONS: return
                    
                    # KISSING LOOP ENHANCEMENT (Hybrid Logic)
                    used_indices = {x for p in int_set for x in p}
                    
                    # Find compatible local options for idx
                    compatible_idx = []
                    for l_set in local_options[idx]:
                        l_indices = {x for p in l_set for x in p}
                        if l_indices.isdisjoint(used_indices):
                            compatible_idx.append(l_set)
                    if not compatible_idx: compatible_idx = [set()]
                    
                    # Find compatible local options for k
                    compatible_k = []
                    for l_set in local_options[k]:
                        l_indices = {x for p in l_set for x in p}
                        if l_indices.isdisjoint(used_indices):
                            compatible_k.append(l_set)
                    if not compatible_k: compatible_k = [set()]

                    # Iterate combinations
                    for c_i in compatible_idx:
                        for c_k in compatible_k:
                            if len(track1_results) >= MAX_SOLUTIONS: return
                            
                            combined_set = int_set.union(c_i).union(c_k)
                            _backtrack_track1(idx + 1, new_covered_int, current_pairs.union(combined_set))

        # Choice B: Local Fold (or empty if local_options has empty set)
        new_covered_local = covered.copy()
        new_covered_local.add(idx)
        for loc_set in local_options[idx]:
            if len(track1_results) >= MAX_SOLUTIONS: return
            _backtrack_track1(idx + 1, new_covered_local, current_pairs.union(loc_set))

    _backtrack_track1(0, set(), scaffold_pairs)

    # =========================================================================
    # TRACK 2: Enumerative Sequential (Hierarchical)
    # =========================================================================
    
    # 1. Generate all combinations of strictly local folds
    all_local_lists = [local_options[i] for i in range(num_regions)]
    
    for local_combo in itertools.product(*all_local_lists):
        if len(track2_results) >= MAX_SOLUTIONS: break

        base_set = scaffold_pairs.copy()
        for p_set in local_combo:
            base_set.update(p_set)
        
        # 2. Identify remaining unpaired regions (LOOPS)
        paired_mask = [False] * L
        for i, j in base_set:
            if 0 <= i < L: paired_mask[i] = True
            if 0 <= j < L: paired_mask[j] = True
            
        sub_runs = []
        i = 0
        while i < L:
            if not paired_mask[i]:
                start_i = i
                while i < L and not paired_mask[i]:
                    i += 1
                # Use min_len=4 per request
                if (i - start_i) >= 3:
                    sub_runs.append((start_i, i))
            else:
                i += 1
        
        # Limit sub-runs here as well for safety
        if len(sub_runs) > MAX_REGIONS_TO_REFINE:
             sorted_sub = sorted(sub_runs, key=lambda r: (r[1] - r[0]), reverse=True)
             sub_runs = sorted_sub[:MAX_REGIONS_TO_REFINE]
             sub_runs.sort(key=lambda r: r[0])

        if len(sub_runs) < 2:
            track2_results.append(base_set)
            continue

        # PRECOMPUTE OPTIONS FOR LOOPS: Interaction AND Local
        # A. Interaction (DuplexFold)
        loop_int_options = get_loop_interaction_options(
            sub_runs, full_seq, duplex_exe, temperature, extra_args, duplex_cache, L
        )
        
        # B. Local Fold (AllSub) - Recalculate for these specific loop sequences
        loop_local_options = {}
        for l_idx, (start, end) in enumerate(sub_runs):
            subseq = full_seq[start:end]
            key = (subseq, temperature, absolute_energy, percent_energy, tuple(extra_args or []))
            
            if key in allsub_cache:
                sub_structs = allsub_cache[key]
            else:
                sub_structs = call_rnastructure_allsub(
                    allsub_exe, subseq, temperature, absolute_energy, percent_energy, extra_args
                )
                allsub_cache[key] = sub_structs
            
            l_opt_sets = []
            for s_str in sub_structs:
                p_sub = _struct_string_to_pairs(s_str, offset=start)
                p_clean = validate_and_filter_pairs(p_sub, full_seq)
                l_opt_sets.append(p_clean)
            if not l_opt_sets: l_opt_sets = [set()]
            loop_local_options[l_idx] = l_opt_sets

        num_loops = len(sub_runs)
        
        def _backtrack_track2_loops(l_idx: int, l_covered: Set[int], l_pairs: Set[Tuple[int, int]]):
            if len(track2_results) >= MAX_SOLUTIONS: return

            if l_idx >= num_loops:
                track2_results.append(l_pairs)
                return
            
            if l_idx in l_covered:
                _backtrack_track2_loops(l_idx + 1, l_covered, l_pairs)
                return
            
            # Choice A: Loop interacts with future loop k
            for k in range(l_idx + 1, num_loops):
                if k in l_covered: continue
                
                opts = loop_int_options.get((l_idx, k))
                if opts:
                    new_cov_pair = l_covered.copy()
                    new_cov_pair.add(l_idx)
                    new_cov_pair.add(k)
                    for int_set in opts:
                        if len(track2_results) >= MAX_SOLUTIONS: return
                        _backtrack_track2_loops(l_idx + 1, new_cov_pair, l_pairs.union(int_set))

            # Choice B: Loop folds locally (New Feature!)
            new_cov_local = l_covered.copy()
            new_cov_local.add(l_idx)
            for l_set in loop_local_options[l_idx]:
                if len(track2_results) >= MAX_SOLUTIONS: return
                _backtrack_track2_loops(l_idx + 1, new_cov_local, l_pairs.union(l_set))

            # Choice C: Loop remains unpaired (Implicitly covered if local_options includes empty set, but ensures safety)
            # Actually, we should only do this if we haven't exhausted options, or if local options are non-empty but we want to skip.
            # Just to be safe and thorough:
            # _backtrack_track2_loops(l_idx + 1, new_cov_local, l_pairs) 
            # (Skipped because local_options usually includes empty set or we rely on loop_local_options)
        
        _backtrack_track2_loops(0, set(), base_set)

    # =========================================================================
    # Final Output Generation & Filtering
    # =========================================================================
    
    # Merge both tracks
    final_pair_sets = track1_results + track2_results

    unique_structs = []
    seen = set()

    for p_set in final_pair_sets:
        # 1. Filter: Remove Isolated Pairs
        cleaned_set = remove_isolated_pairs(p_set)
        
        # Check for conflicts
        flattened = []
        valid_struct = True
        for (i, j) in cleaned_set:
            if i < 0 or i >= L or j < 0 or j >= L:
                valid_struct = False; break
            flattened.append(i)
            flattened.append(j)
        if not valid_struct: continue
        if len(flattened) != len(set(flattened)): continue 

        # 2. Filter: Topological Complexity
        # Discard if total layers > 3 (1 nested + 2 crossing)
        layers = pairs_to_layers(list(cleaned_set))
        if len(layers) > 3:
            continue
        
        # Generate Rosetta-safe string (using [], {} before aA)
        pk_str = pairs_to_rosetta_string(sorted(list(cleaned_set)), L)
        
        if pk_str not in seen:
            seen.add(pk_str)
            unique_structs.append(pk_str)

    return unique_structs


def get_helical_stacks(pairs: Set[Tuple[int, int]]) -> List[Set[Tuple[int, int]]]:
    """
    Group pairs into continuous helical stacks.
    A pair (i, j) is connected to (i+1, j-1) in the same stack.
    """
    sorted_pairs = sorted(list(pairs))
    stacks = []
    visited = set()
    
    for p in sorted_pairs:
        if p in visited:
            continue
        
        current_stack = {p}
        visited.add(p)
        
        # Expand inwards: (i+1, j-1)
        curr = p
        while True:
            next_p = (curr[0] + 1, curr[1] - 1)
            if next_p in pairs:
                current_stack.add(next_p)
                visited.add(next_p)
                curr = next_p
            else:
                break
        
        stacks.append(current_stack)
        
    return stacks


def generate_scaffold_variants(struct: str, L: int) -> List[Tuple[str, str]]:
    """
    Generate distinct variants of the consensus structure by systematically
    masking structural domains (helical stacks).

    Returns list of (description, structure_string).
    """
    variants: List[Tuple[str, str]] = []
    
    # 1. Base Structure
    variants.append(("Base", struct))

    # 2. 5'/3' End Unpairing (Existing Strategy)
    # Parse generic to handle PKs, filter, then reconstruct simple string for masking
    pairs_all = pairs_from_track(struct)
    
    # Generate end-unpaired variant
    valid_pairs_ends = []
    for i, j in pairs_all:
        u, v = sorted((i, j))
        if u >= UNPAIR_ENDS_LEN and v < (L - UNPAIR_ENDS_LEN):
            valid_pairs_ends.append((u, v))
            
    struct_ends_unpaired = pairs_to_rosetta_string(valid_pairs_ends, L)
    if struct_ends_unpaired != struct:
        variants.append(("Unpair 5'/3' Ends", struct_ends_unpaired))
    
    # 3. Domain Masking (Ambitious Strategy)
    # Identify Stacks
    all_pairs_set = set(pairs_all)
    stacks = get_helical_stacks(all_pairs_set)
    
    # Filter for significant stacks (length >= 2) to avoid noise
    significant_stacks = [s for s in stacks if len(s) >= 2]
    
    # Strategy A: Mask each stack individually
    for idx, stack in enumerate(significant_stacks):
        remaining_pairs = all_pairs_set - stack
        var_str = pairs_to_rosetta_string(list(remaining_pairs), L)
        variants.append((f"Masked Helix #{idx+1}", var_str))
        
    # Strategy B: Mask pairs of stacks (if we have multiple)
    # Limit combinatorial explosion: Only do this if we have > 2 and < 8 stacks
    if 2 < len(significant_stacks) <= 20:
        combo_count = 0
        for stack_pair in combinations(significant_stacks, 2):
            if combo_count > 20: break # Safety limit
            combined_mask = stack_pair[0].union(stack_pair[1])
            remaining_pairs = all_pairs_set - combined_mask
            var_str = pairs_to_rosetta_string(list(remaining_pairs), L)
            variants.append((f"Masked 2 Helices", var_str))
            combo_count += 1

    # Deduplicate by string content (keeping first occurrence description)
    unique_variants = []
    seen_strs = set()
    for desc, s_str in variants:
        if s_str not in seen_strs:
            seen_strs.add(s_str)
            unique_variants.append((desc, s_str))
            
    # Cap total variants to prevent pipeline timeouts
    if len(unique_variants) > MAX_SCAFFOLD_VARIANTS:
        unique_variants = unique_variants[:MAX_SCAFFOLD_VARIANTS]
        
    return unique_variants

# =============================================================================
# Ranking and Selection Logic
# =============================================================================

def calculate_structure_score(struct: str) -> float:
    """
    Score a secondary structure based on fraction paired and knot complexity.
    
    Formula:
        Score = Sum(Value(pair)) / L
        
    Where Value(pair) depends on its topological layer:
        Layer 0 (Nested):  1.0  (Maximizes standard stability)
        Layer 1 (Simple):  0.8  (Good, but effectively penalized vs nested)
        Layer 2+ (Deep):  -1.0  (Heavily penalized artifact/hairball)
    """
    L = len(struct)
    if L == 0: return 0.0
    
    pairs = pairs_from_track(struct)
    layers = pairs_to_layers(pairs)
    
    score_sum = 0.0
    
    for layer_idx, layer in enumerate(layers):
        # Count residues involved (2 per pair)
        n_residues = len(layer) * 2
        
        if layer_idx == 0:
            score_sum += n_residues * 1.0
        elif layer_idx == 1:
            score_sum += n_residues * 0.6
        else:
            score_sum += n_residues * -1.0
            
    return score_sum / L


def hamming_dist(s1: str, s2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def filter_and_rank_structures(
    structs: List[str], 
    max_output: int = 100,
    diversity_dist_frac: float = 0.05
) -> List[str]:
    """
    Select the best structures using a Greedy Diversity Clustering approach.
    """
    if not structs: return []
    L = len(structs[0])
    min_dist = int(L * diversity_dist_frac)
    
    # 1. Score and Sort
    scored = []
    for s in structs:
        score = calculate_structure_score(s)
        scored.append((score, s))
        
    scored.sort(key=lambda x: x[0], reverse=True)
    
    selected = []
    
    # 2. Greedy Selection
    for score, candidate in scored:
        if len(selected) >= max_output:
            break
        
        is_distinct = True
        for sel in selected:
            if hamming_dist(candidate, sel) < min_dist:
                is_distinct = False
                break
        
        if is_distinct:
            selected.append(candidate)
            
    return selected

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
    allsub_cache: Dict = None,
    duplex_cache: Dict = None,
    return_all: bool = False, # If False, returns filtered set.
) -> List[str]:
    """
    Wrapper that generates multiple consensus variants (by masking domains)
    and refines each one.
    
    By default (return_all=False), this returns a FILTERED list (top 100),
    to prevent pipeline explosion.
    """
    L = len(full_seq)
    
    # 1. Generate Consensus Variants
    scaffold_candidates = generate_scaffold_variants(struct, L)
    
    all_refined = []
    
    for idx, (desc, variant_struct) in enumerate(scaffold_candidates):
        sys.stderr.write(
            f"   > Variant {idx+1}/{len(scaffold_candidates)} [{desc}]: "
            f"Refining...\r"
        )
        results = _refine_structure_impl(
            variant_struct, full_seq, allsub_exe, duplex_exe, min_unpaired_len,
            temperature, absolute_energy, percent_energy, extra_args,
            allsub_cache, duplex_cache
        )
        all_refined.extend(results)
    
    sys.stderr.write("\n")
    
    # 2. Deduplicate
    unique_structs = []
    seen = set()
    for s in all_refined:
        if s not in seen:
            seen.add(s)
            unique_structs.append(s)
            
    if return_all:
        return unique_structs
        
    # 3. Filter by Default
    filtered = filter_and_rank_structures(unique_structs, max_output=1000, diversity_dist_frac=0.0)
    return filtered


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Refine unpaired regions using AllSub and DuplexFold.")
    parser.add_argument("--sto", help="Stockholm file for sequence.")
    parser.add_argument("--seq-file", help="FASTA/SEQ file for sequence.")
    parser.add_argument("--seq-name", default=None)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--db-out", required=True, help="Path for filtered structures (default output).")
    parser.add_argument("--db-out-all", default=None, help="Optional path to save ALL refined structures (before filtering).")
    
    parser.add_argument("--allsub-exe", default="$PROJECT/repos/RNAstructure/exe/AllSub")
    parser.add_argument("--duplex-exe", default="$PROJECT/repos/RNAstructure/exe/DuplexFold")
    
    parser.add_argument("--min-unpaired", type=int, default=4)
    
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--allsub-abs", type=float, default=None, help="-a flag for AllSub")
    parser.add_argument("--allsub-pct", type=float, default=None, help="-p flag for AllSub")
    
    parser.add_argument("--max-structures", type=int, default=1000, 
                        help="Maximum number of refined structures to output (default: 100).")
    parser.add_argument("--diversity-dist", type=float, default=0.05,
                        help="Minimum Hamming distance fraction (0.0-1.0) for diversity clustering.")
    
    parser.add_argument("--extra-arg", action="append", default=[])

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
    
    # We collect everything first if user wants "all", or just filtered?
    # To support --db-out-all, we must ask refine_structure to return all, then filter here.
    
    collect_all_mode = (args.db_out_all is not None)
    
    all_structs_for_output = [] # This will hold the big list if requested, or the small list if not.
    
    # If explicit full output requested, we ask refine_structure for everything.
    # If not, we let it filter internally to save memory/time passing lists around.
    
    for idx, s in enumerate(structs):
        sys.stderr.write(f"[INFO] Structure {idx+1}/{len(structs)}: Generating consensus variants...\n")
        
        # Call refine_structure
        # If we need the full list (collect_all_mode), pass return_all=True
        # Otherwise pass False (default) to get pre-filtered list.
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
            return_all=collect_all_mode
        )
        all_structs_for_output.extend(variants)
        sys.stderr.write(f"[INFO] Structure {idx+1} yielded {len(variants)} models.\n")
    
    # 1. Handle "Full" Output
    if args.db_out_all:
        sys.stderr.write(f"[INFO] Writing ALL {len(all_structs_for_output)} structures to {args.db_out_all}\n")
        with Path(args.db_out_all).open("w") as f:
            f.write(full_seq + "\n")
            seen_full = set()
            for rs in all_structs_for_output:
                if rs not in seen_full:
                    f.write(rs + "\n")
                    seen_full.add(rs)
    
    # 2. Filter for "Main" Output
    # If we collected all, we must filter now.
    # If we collected filtered (default), we just write it (maybe re-filter globally if multiple input structs).
    
    if collect_all_mode:
        # We have the big list, filter it down.
        final_structs = filter_and_rank_structures(
            all_structs_for_output, 
            max_output=args.max_structures, 
            diversity_dist_frac=args.diversity_dist
        )
    else:
        # We have the per-structure filtered lists combined. 
        # If input db had multiple structures, we might still have > max_structures total.
        # So we filter again globally to be safe.
        final_structs = filter_and_rank_structures(
            all_structs_for_output, 
            max_output=args.max_structures, 
            diversity_dist_frac=args.diversity_dist
        )
    
    with Path(args.db_out).open("w") as f:
        f.write(full_seq + "\n")
        for rs in final_structs:
            f.write(rs + "\n")
    
    sys.stderr.write(f"[INFO] Wrote {len(final_structs)} filtered structures to {args.db_out}\n")

if __name__ == "__main__":
    main()
