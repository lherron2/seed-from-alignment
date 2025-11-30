#!/usr/bin/env python3
"""
Filter/cluster an ensemble of secondary structures in .db format
and convert them to Rosetta-compatible dot-parens notation.

Features
--------
1. Read a .db file with one dot-bracket structure per line.
2. Remove exact duplicates (always).
3. Optional greedy merging based on Hamming distance:
   - For each structure, if it is within a threshold of any already-kept
     representative, it is discarded (merged).
4. Convert the remaining structures to Rosetta notation:

   Rosetta pair types (in order of assignment):
       1. ()  (main nested structure)
       2. []
       3. {}
       4. aA
       5. bB
       6. cC
       7. dD
       ...

   Mapping rule:
       - Keep () as ().
       - For all other CaCoFold pair types ([], {}, <>, aA, bB, cC, ...),
         assign them to the next available Rosetta pair type in that order,
         based on *order of first appearance* in the structure.

So, for example, if a structure uses () + aA + bB + cC (and no [] or {}),
you get:
    () -> ()
    aA -> []
    bB -> {}
    cC -> aA
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Iterable
import sys


# -------------------------------------------------------------
# Sequence + base-pair helpers
# -------------------------------------------------------------

def read_fasta_sequence(path: Path) -> str:
    """
    Read a single sequence from a FASTA file and return it as an
    uppercase string with no whitespace.
    """
    seq_lines: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    if not seq_lines:
        raise ValueError(f"No sequence found in FASTA file {path}")
    return "".join(seq_lines).upper()


def is_complementary(base1: str, base2: str, allow_wobble: bool = True) -> bool:
    """
    Return True if the two nucleotides form a canonical Watson-Crick
    pair (AU/UA/GC/CG), optionally allowing GU wobble (GU/UG).
    Any non-ACGU nucleotide is treated as non-complementary.
    """
    b1 = base1.upper()
    b2 = base2.upper()

    if b1 not in "ACGU" or b2 not in "ACGU":
        return False

    pair = b1 + b2
    canonical = {"AU", "UA", "GC", "CG"}

    if allow_wobble:
        canonical |= {"GU", "UG"}

    return pair in canonical


def prune_noncomplementary_pairs(struct: str, seq: str) -> tuple[str, int]:
    """
    Given a dot-bracket structure (Rosetta notation) and its sequence,
    remove any base pairs that are not between complementary nucleotides
    by setting those positions to '.'.

    Returns:
        (new_struct, num_pairs_removed)
    """
    if len(struct) != len(seq):
        raise ValueError(
            f"Structure length ({len(struct)}) and sequence length ({len(seq)}) do not match"
        )

    pairs = _parse_pairs(struct)
    chars = list(struct)
    removed = 0

    for i, j, _op in pairs:
        if not is_complementary(seq[i], seq[j]):
            if chars[i] != "." or chars[j] != ".":
                chars[i] = "."
                chars[j] = "."
                removed += 1

    return "".join(chars), removed


# -------------------------------------------------------------
# I/O helpers
# -------------------------------------------------------------

def read_db_structures(path: Path) -> List[str]:
    """Read one structure per non-empty line from a .db file."""
    structs: List[str] = []
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            structs.append(s)
    if not structs:
        raise ValueError(f"No structures found in {path}")
    # sanity: check consistent lengths
    lengths = {len(s) for s in structs}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent structure lengths in {path}: {lengths}")
    return structs


def write_db_structures(path: Path, structs: Iterable[str]) -> None:
    with path.open("w") as fh:
        for s in structs:
            fh.write(s + "\n")


# -------------------------------------------------------------
# Filtering / clustering
# -------------------------------------------------------------

def hamming_leq(a: str, b: str, threshold: int) -> bool:
    """
    Return True if Hamming(a, b) <= threshold, False otherwise.

    Optimized for small thresholds:
        - Early-exit as soon as mismatches exceed threshold.
    """
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")

    mismatches = 0
    for c1, c2 in zip(a, b):
        if c1 != c2:
            mismatches += 1
            if mismatches > threshold:
                return False
    return True


def count_paired(s: str) -> int:
    """
    Count the number of paired positions (i.e., non-dot) in the structure.
    This is a cheap prefilter: if |paired(a) - paired(b)| > threshold,
    they cannot be within Hamming distance <= threshold.
    """
    return sum(ch != "." for ch in s)


def merge_by_hamming(
    structs: list[str],
    threshold: int,
) -> list[str]:
    """
    Greedy diversity filter based on Hamming distance with optimizations:

    - Always keeps the first structure.
    - For each subsequent structure s:
        * First do a cheap filter on paired count vs each kept rep:
              |paired(s) - paired(rep)| <= threshold
          otherwise skip rep entirely.
        * If it passes that, use hamming_leq(s, rep, threshold) which
          early-exits once mismatches > threshold.
        * If s is within threshold of ANY rep, s is discarded.
        * Otherwise s starts a new cluster (rep).

    threshold <= 0 => no merging (returns input as-is).
    """
    if threshold <= 0:
        return structs

    reps: list[str] = []
    rep_paired: list[int] = []

    for s in structs:
        paired_s = count_paired(s)
        keep = True

        for r, paired_r in zip(reps, rep_paired):
            # Cheap prefilter on number of paired positions
            if abs(paired_s - paired_r) > threshold:
                continue

            # Full Hamming check with early exit
            if hamming_leq(s, r, threshold):
                keep = False
                break

        if keep:
            reps.append(s)
            rep_paired.append(paired_s)

    return reps


def deduplicate(structs: List[str]) -> List[str]:
    """Preserve order while removing exact duplicates."""
    seen = set()
    out: List[str] = []
    for s in structs:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# -------------------------------------------------------------
# Parsing CaCoFold-style pairs and mapping to Rosetta
# -------------------------------------------------------------

def _parse_pairs(struct: str) -> List[Tuple[int, int, str]]:
    """
    Parse a generic CaCoFold-style dot-bracket string into pairs.

    Supports:
        - () [] {} <>   (usual bracket types)
        - aA, bB, cC, ... (letters: open = lowercase, close = uppercase)

    Returns:
        List of (i, j, open_char) pairs, with i < j.
    """
    # open -> close mapping
    OPEN_TO_CLOSE: Dict[str, str] = {
        "(": ")",
        "[": "]",
        "{": "}",
        "<": ">",
        # letters a-z use A-Z as closers
        **{chr(ord("a") + i): chr(ord("A") + i) for i in range(26)},
    }

    CLOSE_TO_OPEN: Dict[str, str] = {v: k for k, v in OPEN_TO_CLOSE.items()}

    # stacks for open characters
    stacks: Dict[str, List[int]] = {op: [] for op in OPEN_TO_CLOSE}
    pairs: List[Tuple[int, int, str]] = []

    for idx, ch in enumerate(struct):
        if ch in OPEN_TO_CLOSE:
            stacks[ch].append(idx)
        elif ch in CLOSE_TO_OPEN:
            op = CLOSE_TO_OPEN[ch]
            st = stacks.get(op)
            if st is None or not st:
                raise ValueError(
                    f"Unmatched closing '{ch}' at position {idx} in structure:\n{struct}"
                )
            i = st.pop()
            if i > idx:
                i, idx = idx, i
            pairs.append((i, idx, op))
        elif ch == ".":
            continue
        else:
            raise ValueError(
                f"Unexpected character '{ch}' in structure. "
                f"Add it to OPEN_TO_CLOSE/CLOSE_TO_OPEN if it is a valid pair symbol."
            )

    # Check unmatched opens
    for op, st in stacks.items():
        if st:
            raise ValueError(
                f"Unmatched opening '{op}' at positions {st} in structure:\n{struct}"
            )

    return pairs


def convert_to_rosetta_notation(struct: str) -> str:
    """
    Convert a generic CaCoFold dot-bracket string to Rosetta-compatible notation.

    Rosetta pair types (in order of assignment):
        1. ()       (main)
        2. []       (first pseudoknot)
        3. {}       (second pseudoknot)
        4. aA       (third)
        5. bB       (fourth)
        6. cC       (fifth)
        ...
    Mapping rule:
        - Keep () -> ().
        - For all other CaCoFold pair types ([], {}, <>, aA, bB, cC, ...),
          assign them to the next available Rosetta pair type in that order,
          based on *order of first appearance* in the structure.
    """
    pairs = _parse_pairs(struct)
    n = len(struct)

    # Rosetta pair types in order
    rosetta_types: List[Tuple[str, str]] = [
        ("(", ")"),
        ("[", "]"),
        ("{", "}"),
    ] + [(chr(ord("a") + i), chr(ord("A") + i)) for i in range(26)]

    # Determine order of first appearance of each pair type
    first_pos: Dict[str, int] = {}
    for i, j, op in pairs:
        pos = i  # where this pair opens
        if op not in first_pos or pos < first_pos[op]:
            first_pos[op] = pos

    # Sort pair types by first appearance
    types_in_order = sorted(first_pos.keys(), key=lambda op: first_pos[op])

    # Build mapping from CaCoFold pair type -> Rosetta pair type
    mapping: Dict[str, Tuple[str, str]] = {}

    # First: handle main '()' explicitly
    if "(" in types_in_order:
        mapping["("] = rosetta_types[0]  # () -> ()
        # Remaining types, excluding '('
        other_types = [t for t in types_in_order if t != "("]
        start_idx = 1  # start assigning from '[]'
    else:
        # No parentheses at all; just assign everything from the top
        other_types = types_in_order
        start_idx = 0

    # Assign other types sequentially to the remaining Rosetta pair types
    for k, op in enumerate(other_types):
        idx = start_idx + k
        if idx >= len(rosetta_types):
            raise RuntimeError(
                f"Ran out of Rosetta pair types (more than {len(rosetta_types)} "
                f"unique pair types encountered)."
            )
        mapping[op] = rosetta_types[idx]

    # Build the output structure
    new_struct = ["." for _ in range(n)]
    for i, j, op in pairs:
        open_out, close_out = mapping[op]
        new_struct[i] = open_out
        new_struct[j] = close_out

    return "".join(new_struct)


# -------------------------------------------------------------
# Main CLI
# -------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Filter/cluster .db secondary structures and convert them "
            "to Rosetta-compatible dot-parens notation."
        )
    )
    parser.add_argument(
        "--db-in",
        required=True,
        help=".db file with one dot-bracket structure per line.",
    )
    parser.add_argument(
        "--db-out",
        required=True,
        help="Output .db file with filtered + Rosetta-compatible structures.",
    )
    parser.add_argument(
        "--hamming-threshold",
        type=int,
        default=0,
        help=(
            "Hamming distance threshold for merging similar structures. "
            "0 => no merging (only exact dedup). Typical values: 2 or 3."
        ),
    )
    parser.add_argument(
        "--hamming-frac",
        type=float,
        default=0.05,
        help=(
            "If --hamming-threshold == 0, choose an adaptive threshold "
            "as round(hamming_frac * L), where L is the structure length. "
            "Set hamming_frac <= 0 to disable Hamming merging."
        ),
    )
    parser.add_argument(
        "--max-structures",
        type=int,
        default=None,
        help="Optional cap on number of structures to keep after filtering.",
    )
    parser.add_argument(
        "--seq-file",
        help=(
            "FASTA file with the RNA sequence corresponding to the "
            "structures (used to prune non-complementary base pairs)."
        ),
    )

    args = parser.parse_args()
    db_in = Path(args.db_in)
    db_out = Path(args.db_out)

    structs = read_db_structures(db_in)
    print(f"[INFO] Read {len(structs)} structures from {db_in}")

    # 1) Remove exact duplicates
    structs = deduplicate(structs)
    print(f"[INFO] After deduplication: {len(structs)} unique structures")

    # 2) Determine Hamming threshold (static or adaptive)
    L = len(structs[0])

    # Start from user-specified threshold
    threshold = args.hamming_threshold

    if threshold == 0:
        # Use adaptive threshold based on structure length
        if args.hamming_frac > 0.0:
            threshold = max(1, int(round(args.hamming_frac * L)))
            print(
                f"[INFO] Using adaptive Hamming threshold: "
                f"threshold = round({args.hamming_frac} * {L}) = {threshold}"
            )
        else:
            threshold = 0
            print("[INFO] Hamming merging disabled (hamming_frac <= 0 and threshold == 0)")
    else:
        print(f"[INFO] Using fixed Hamming threshold: {threshold}")

    # 3) Greedy Hamming-based merging (optional)
    if threshold > 0:
        structs = merge_by_hamming(structs, threshold)
        print(
            f"[INFO] After Hamming merge (threshold={threshold}): "
            f"{len(structs)} structures"
        )

    # 4) Optional cap on number of structures
    if args.max_structures is not None and len(structs) > args.max_structures:
        structs = structs[: args.max_structures]
        print(
            f"[INFO] Truncated to max_structures={args.max_structures}: "
            f"{len(structs)} structures"
        )

    converted = []
    num_skipped = 0
    for idx, s in enumerate(structs):
        try:
            converted.append(convert_to_rosetta_notation(s))
        except ValueError as e:
            # Log and skip any invalid structures (e.g., unbalanced parentheses)
            sys.stderr.write(
                f"[WARN] Skipping structure {idx} due to parsing error: {e}\n"
            )
            num_skipped += 1
            continue

    sys.stderr.write(
        f"[INFO] Successfully converted {len(converted)} structures; "
        f"skipped {num_skipped} invalid structures.\n"
    )

    # 6) Optionally prune non-complementary base pairs using the sequence
    if args.seq_file is not None:
        seq_path = Path(args.seq_file)
        seq = read_fasta_sequence(seq_path)
        if len(seq) != L:
            raise ValueError(
                f"Sequence length ({len(seq)}) from {seq_path} does not match "
                f"structure length ({L})"
            )

        pruned: List[str] = []
        total_removed = 0
        for s in converted:
            new_s, removed = prune_noncomplementary_pairs(s, seq)
            pruned.append(new_s)
            total_removed += removed

        converted = pruned
        if total_removed > 0:
            print(
                f"[INFO] Removed {total_removed} non-complementary base pairs "
                f"across all structures based on sequence in {seq_path}"
            )
        else:
            print("[INFO] No non-complementary base pairs detected.")

    write_db_structures(db_out, converted)
    print(f"[INFO] Wrote {len(converted)} Rosetta-compatible structures to {db_out}")


if __name__ == "__main__":
    main()

