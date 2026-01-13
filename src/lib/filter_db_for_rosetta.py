"""
Filter/cluster an ensemble of secondary structures in .db format
and convert them to Rosetta-compatible dot-parens notation.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

GAP_CHARS = set("-._~")

# Copying logic from sample_cacofold_structures to ensure clean layering
OPEN_TO_CLOSE_LAYER = {"(": ")", "[": "]", "{": "}"}
CLOSE_TO_OPEN_LAYER = {v: k for k, v in OPEN_TO_CLOSE_LAYER.items()}


def pairs_cross(p: tuple[int, int], q: tuple[int, int]) -> bool:
    """
    Return True if arcs (i, j) and (k, l) cross (pseudoknot).
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
    Given a list of pairs (i, j), partition them into layers so that
    within each layer there are no crossing arcs.
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


def read_ungapped_seq_from_sto(path: Path, seq_name: str | None = None) -> str:
    """
    Read the ungapped RNA sequence for `seq_name` from a Stockholm (.sto) file.
    """
    sequences: dict[str, list[str]] = {}

    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("//"):
                break

            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                continue
            name, s = parts[0], parts[1].strip()
            sequences.setdefault(name, []).append(s)

    if not sequences:
        raise ValueError(f"No sequences found in Stockholm file {path}")

    if seq_name is None:
        name = next(iter(sequences.keys()))
    else:
        if seq_name not in sequences:
            available = ", ".join(sorted(sequences.keys()))
            raise ValueError(f"Sequence '{seq_name}' not found in {path}. Available: {available}")
        name = seq_name

    aligned = "".join(sequences[name])
    ungapped = "".join(ch for ch in aligned if ch not in GAP_CHARS)
    return ungapped.upper().replace("T", "U")


# -------------------------------------------------------------
# Sequence + base-pair helpers
# -------------------------------------------------------------


def read_fasta_sequence(path: Path) -> str:
    """
    Read a single sequence from a FASTA file and return it as an
    uppercase RNA string (T converted to U) with no whitespace.
    """
    seq_lines: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)

    if not seq_lines:
        raise ValueError(f"No sequence found in FASTA file {path}")

    seq = "".join(seq_lines).upper()
    return seq.replace("T", "U")


def is_complementary(base1: str, base2: str, allow_wobble: bool = True) -> bool:
    """
    Heuristic complementarity check mimicking Rosetta's expectations.
    """

    def _norm(b: str) -> str:
        b = b.upper()
        if b == "T":
            b = "U"
        return b

    b1 = _norm(base1)
    b2 = _norm(base2)

    # Strict check when both bases are unambiguous
    if b1 in "ACGU" and b2 in "ACGU":
        pair = b1 + b2
        canonical = {"AU", "UA", "GC", "CG"}
        if allow_wobble:
            canonical |= {"GU", "UG"}
        return pair in canonical

    # At least one ambiguous / unknown base: be permissive here
    return True


def prune_noncomplementary_pairs(
    struct: str, seq: str = None, allow_wobble: bool = True
) -> tuple[str, int]:
    """
    Given a dot-bracket structure (Rosetta notation) and its sequence,
    remove any base pairs that are not between complementary nucleotides.
    """
    if struct is None:
        raise ValueError("Structure passed to prune_noncomplementary_pairs is None.")
    if seq is None:
        raise ValueError("Sequence passed to prune_noncomplementary_pairs is None.")

    if len(struct) != len(seq):
        raise ValueError(
            f"Structure length ({len(struct)}) and sequence length ({len(seq)}) do not match."
        )

    pairs = _parse_pairs(struct)
    chars = list(struct)
    removed = 0

    for i, j, op in pairs:
        b1 = seq[i]
        b2 = seq[j]
        if not is_complementary(b1, b2, allow_wobble=allow_wobble):
            chars[i] = "."
            chars[j] = "."
            removed += 1

    return "".join(chars), removed


# -------------------------------------------------------------
# I/O helpers
# -------------------------------------------------------------


def read_db_structures(path: Path) -> list[str]:
    """
    Read structures from a .db file.
    """
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


def write_db_structures(
    path: Path,
    structs: Iterable[str],
    seq: str | None = None,
) -> None:
    """
    Write a .db file.
    """
    with path.open("w") as fh:
        if seq is not None:
            fh.write(seq + "\n")
        for s in structs:
            fh.write(s + "\n")


# -------------------------------------------------------------
# Filtering / clustering
# -------------------------------------------------------------


def hamming_leq(a: str, b: str, threshold: int) -> bool:
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
    return sum(ch != "." for ch in s)


def merge_by_hamming(
    structs: list[str],
    threshold: int,
) -> list[str]:
    if threshold <= 0:
        return structs

    reps: list[str] = []
    rep_paired: list[int] = []

    for s in structs:
        paired_s = count_paired(s)
        keep = True

        for r, paired_r in zip(reps, rep_paired):
            if abs(paired_s - paired_r) > threshold:
                continue
            if hamming_leq(s, r, threshold):
                keep = False
                break

        if keep:
            reps.append(s)
            rep_paired.append(paired_s)

    return reps


def deduplicate(structs: list[str]) -> list[str]:
    seen = set()
    out: list[str] = []
    for s in structs:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# -------------------------------------------------------------
# Parsing CaCoFold-style pairs and mapping to Rosetta
# -------------------------------------------------------------


def _parse_pairs(struct: str) -> list[tuple[int, int, str]]:
    """
    Parse a generic CaCoFold-style dot-bracket string into pairs.
    """
    OPEN_TO_CLOSE: dict[str, str] = {
        "(": ")",
        "[": "]",
        "{": "}",
        "<": ">",
        **{chr(ord("a") + i): chr(ord("A") + i) for i in range(26)},
    }

    CLOSE_TO_OPEN: dict[str, str] = {v: k for k, v in OPEN_TO_CLOSE.items()}

    stacks: dict[str, list[int]] = {op: [] for op in OPEN_TO_CLOSE}
    pairs: list[tuple[int, int, str]] = []

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
            raise ValueError(f"Unexpected character '{ch}' in structure.")

    for op, st in stacks.items():
        if st:
            raise ValueError(f"Unmatched opening '{op}' at positions {st} in structure:\n{struct}")

    return pairs


def convert_to_rosetta_notation(struct: str) -> str:
    """
    Convert a generic CaCoFold dot-bracket string to Rosetta-compatible notation.
    """
    parsed = _parse_pairs(struct)
    pairs = [(i, j) for i, j, _ in parsed]

    # Sanity check: Ensure input structure didn't contain overlapping pairs
    indices = [idx for pair in pairs for idx in pair]
    if len(indices) != len(set(indices)):
        raise ValueError("Invalid matching: overlapping pairs detected in input structure.")

    n = len(struct)
    layers = pairs_to_layers(pairs)

    new_chars = ["."] * n

    rosetta_brackets = [("(", ")"), ("[", "]"), ("{", "}")]
    rosetta_brackets += [(chr(ord("a") + i), chr(ord("A") + i)) for i in range(26)]

    for layer_idx, layer_pairs in enumerate(layers):
        if layer_idx < len(rosetta_brackets):
            op, cl = rosetta_brackets[layer_idx]
        else:
            op, cl = "{", "}"

        for i, j in layer_pairs:
            new_chars[i] = op
            new_chars[j] = cl

    return "".join(new_chars)


# -------------------------------------------------------------
# Main CLI
# -------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Filter/cluster .db secondary structures and convert them "
            "to Rosetta-compatible dot-parens notation."
        )
    )
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--db-out", required=True)
    parser.add_argument("--hamming-threshold", type=int, default=0)
    parser.add_argument("--hamming-frac", type=float, default=0.05)
    parser.add_argument("--max-structures", type=int, default=None)
    parser.add_argument("--sto", help="CaCoFold/R-scape Stockholm (.sto) file")
    parser.add_argument("--seq-name", help="Name of sequence in .sto file.")

    args = parser.parse_args(argv)
    db_in = Path(args.db_in)
    db_out = Path(args.db_out)

    structs = read_db_structures(db_in)
    print(f"[INFO] Read {len(structs)} structures from {db_in}")

    structs = deduplicate(structs)
    L = len(structs[0])
    threshold = args.hamming_threshold

    if threshold == 0:
        if args.hamming_frac > 0.0:
            threshold = max(1, int(round(args.hamming_frac * L)))

    if threshold > 0:
        structs = merge_by_hamming(structs, threshold)

    if args.max_structures is not None and len(structs) > args.max_structures:
        structs = structs[: args.max_structures]

    sto_path = Path(args.sto)
    if not sto_path.is_file():
        raise FileNotFoundError(f"Stockholm file not found at {sto_path}")

    seq = read_ungapped_seq_from_sto(sto_path, args.seq_name)

    converted = []
    for s in structs:
        try:
            rn = convert_to_rosetta_notation(s)
            ps, _ = prune_noncomplementary_pairs(rn, seq)
            converted.append(ps)
        except ValueError as e:
            sys.stderr.write(f"[WARN] Dropped invalid structure: {e}\n")
            continue

    write_db_structures(db_out, converted, seq=seq)
    print(f"[INFO] Wrote {len(converted)} Rosetta-compatible structures to {db_out}")


if __name__ == "__main__":
    main()
