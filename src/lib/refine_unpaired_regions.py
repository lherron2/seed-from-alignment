"""
Refine long unpaired regions in CaCoFold-sampled structures using RNAstructure Fold,
and optionally refine unpaired 5'/3' terminal ends jointly with RNAstructure DuplexFold.

Workflow:
    - Input:
        * A FASTA/SEQ file with the full ungapped RNA sequence.
        * A .db file containing one dot-bracket structure per line
          (as produced by sample_cacofold_structures.py).
    - For each structure:
        * Find contiguous runs of '.' longer than a given threshold (default: 15).
          For INTERNAL runs (not touching either terminus):
              - Take the corresponding subsequence of the RNA.
              - Call RNAstructure Fold on that subsequence to get a local dot-bracket structure.
              - Splice that local structure back into the global structure, replacing the '.' run.
        * Additionally, if BOTH the 5' end and 3' end have long unpaired runs (by default,
          length >= min_unpaired, or min_terminal_unpaired if provided):
              - Extract the 5' and 3' terminal subsequences.
              - Call RNAstructure DuplexFold on those two subsequences.
              - Map the resulting duplex base pairs back into the global structure, so that
                5' positions receive '(' and 3' positions receive ')'.
    - Write refined structures to a new .db file.

Notes:
    - This assumes:
        * '.' = unpaired
        * Any other character indicates an occupied/base-paired position
          (e.g., '()', 'aA', 'bB', etc. from your pseudoknot layers).
    - RNAstructure Fold is called via the command line, e.g.:
        $PROJECT/repos/RNAstructure/exe/Fold seq.fa - -k -mfe -q
    - RNAstructure DuplexFold is called via:
        $PROJECT/repos/RNAstructure/exe/DuplexFold seq5.fa seq3.fa duplex.ct
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

UNPAIRED_CHAR = "."
# For parsing RNAstructure Fold bracket output:
ALLOWED_DB_CHARS = set(".()")
# Alignment gaps (same as in sample_cacofold_structures.py)
GAP_CHARS = set("-._~")


def read_ungapped_seq_from_sto(path: Path, seq_name: str | None = None) -> str:
    """Read the ungapped RNA sequence for ``seq_name`` from a Stockholm file.

    The returned sequence is upper-case and any DNA ``T`` is converted to
    RNA ``U`` so that all downstream steps see a consistent alphabet.
    If ``seq_name`` is None, the first sequence in the file is used.
    """
    sequences: Dict[str, List[str]] = {}
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not line:
                continue
            if line.startswith("#"):
                # comments, #=GC, etc.
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
        name = next(iter(sequences))
    else:
        if seq_name not in sequences:
            available = ", ".join(sorted(sequences.keys()))
            raise ValueError(
                f"Sequence '{seq_name}' not found in {path}. "
                f"Available: {available}"
            )
        name = seq_name

    aligned = "".join(sequences[name])
    ungapped = "".join(ch for ch in aligned if ch not in GAP_CHARS)
    return ungapped.upper().replace("T", "U")

def read_fasta_sequence(path: Path, seq_name: str | None = None) -> str:
    """Read the first (or named) sequence from FASTA/SEQ and normalise to RNA.

    The returned sequence is upper-case and any DNA ``T`` is converted to
    RNA ``U``.
    """
    seqs: Dict[str, str] = {}
    current_name: str | None = None
    chunks: List[str] = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
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
        # just take the first
        name = next(iter(seqs))
    else:
        if seq_name not in seqs:
            raise ValueError(
                f"Sequence '{seq_name}' not found in {path}. "
                f"Available: {', '.join(seqs.keys())}"
            )
        name = seq_name

    return seqs[name].upper().replace("T", "U")


def read_db_structures(path: Path) -> List[str]:
    """
    Read structures from a .db file.

    Supports two layouts:

      1. Legacy: one structure per non-empty line.
      2. New:    first non-empty line is the ungapped sequence
                 (no dot/bracket chars), remaining non-empty
                 lines are structures.

    Returns
    -------
    structs : list of dot-bracket strings
        The structures only (sequence header, if present,
        is ignored here).
    """
    structs: List[str] = []
    sequence: str | None = None

    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue

            if sequence is None:
                # Heuristic: treat as sequence if it has no dot/bracket chars
                if all(ch not in ".()[]{}<>" for ch in s):
                    sequence = s
                    continue

            structs.append(s)

    if not structs:
        raise ValueError(f"No structures found in {path}")

    lengths = {len(s) for s in structs}
    if len(lengths) != 1:
        raise ValueError(
            f"Structures in {path} have inconsistent lengths: {lengths}"
        )

    if sequence is not None and len(sequence) not in lengths:
        raise ValueError(
            f"Sequence length in {path} ({len(sequence)}) does not match "
            f"structure length(s) {lengths}"
        )

    return structs



def find_unpaired_runs(struct: str, min_len: int) -> List[Tuple[int, int]]:
    """
    Find contiguous runs of '.' in the structure with length >= min_len.

    Returns a list of (start, end) indices in Python slice convention:
        struct[start:end] is the unpaired region, length = end - start.
    """
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

def find_terminal_unpaired_ends(struct: str, min_len: int) -> Tuple[int, int]:
    """
    Find unpaired runs at the 5' and 3' ends.

    Returns (len_5prime, len_3prime). If both ends are non-zero but
    shorter than `min_len`, or if the entire sequence is a single unpaired
    run, returns (0, 0).

    Note: we treat an end as "terminal" if it has at least one '.', and we
    only *gate* DuplexFold on whether at least one end is >= min_len. This
    makes it possible to refine asymmetric cases like len_5=9, len_3=17
    with min_len=15.
    """
    n = len(struct)

    # 5' end
    len_5 = 0
    for i in range(n):
        if struct[i] == UNPAIRED_CHAR:
            len_5 += 1
        else:
            break

    # 3' end
    len_3 = 0
    for j in range(n - 1, -1, -1):
        if struct[j] == UNPAIRED_CHAR:
            len_3 += 1
        else:
            break

    # If everything is unpaired, don't treat this as "two ends"
    if len_5 == n or len_3 == n:
        return 0, 0

    # Require that both ends are present (non-zero), but only one of them
    # needs to exceed the threshold. This is friendlier for asymmetric
    # ends like 9/17 nt with min_len=15.
    if len_5 > 0 and len_3 > 0 and (len_5 >= min_len or len_3 >= min_len):
        return len_5, len_3
    else:
        return 0, 0


def call_rnastructure_fold(
    fold_exe: Path,
    subseq: str,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> str:
    """
    Call RNAstructure Fold on a subsequence and return the dot-bracket structure.

    Uses:
        Fold subseq.fa - -k -mfe -q [extra args]

    Assumes that the bracket output has a line with only '().' of len
    len(subseq). We take the *last* such line as the structure.
    """
    extra_args = list(extra_args) if extra_args is not None else []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        seq_path = tmpdir_path / "subseq.fa"

        with seq_path.open("w") as fh:
            fh.write(">subseq\n")
            fh.write(subseq + "\n")

        cmd: List[str] = [
            str(fold_exe),
            str(seq_path),
            "-",     # write DBN to stdout (since we use -k)
            "-k",    # dot-bracket output
            "-mfe",  # minimum free energy only (fast, single structure)
            "-q",    # quiet
        ]
        if temperature is not None:
            cmd.extend(["-T", str(temperature)])
        cmd.extend(extra_args)

        sys.stderr.write(
            f"[Fold] subseq_len={len(subseq)}; running: {' '.join(cmd)}\n"
        )

        try:
            result = subprocess.run(
                cmd,
                check=True,
                text=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                f"[ERROR] RNAstructure Fold failed (exit {e.returncode})\n"
                f"Command: {' '.join(cmd)}\n"
                f"STDOUT:\n{e.stdout}\n"
                f"STDERR:\n{e.stderr}\n"
            )
            raise

        lines = [ln.strip() for ln in result.stdout.splitlines() if ln.strip()]
        # Heuristic: structure is the last line whose length matches subseq
        # and whose chars are subset of .()
        for line in reversed(lines):
            if len(line) == len(subseq) and set(line) <= ALLOWED_DB_CHARS:
                return line

        raise RuntimeError(
            f"Could not parse Fold output for subseq of length {len(subseq)}. "
            f"Output lines were:\n" + "\n".join(lines)
        )


def call_rnastructure_duplexfold(
    duplex_exe: Path,
    seq5: str,
    seq3: str,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> List[List[Tuple[int, int]]]:
    """
    Call RNAstructure DuplexFold on two terminal subsequences.

    Returns
    -------
    all_pairs : list of list of (i5, j3)
        Outer list indexes the CT structure (1, 2, ...).
        Inner list contains inter-strand base pairs in *local* 1-based
        coordinates:

            i5 ∈ [1, len(seq5)]
            j3 ∈ [1, len(seq3)]

    Notes
    -----
    DuplexFold writes multiple CT structures by repeating:

        header: "<N>  ENERGY = ..."
        N lines: CT entries (indices 1..N)

    The base identities and positions are the same for each block;
    only the pairing column changes. We parse each block separately.
    """
    extra_args = list(extra_args) if extra_args is not None else []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        seq5_path = tmpdir_path / "seq5.fa"
        seq3_path = tmpdir_path / "seq3.fa"
        ct_path = tmpdir_path / "duplex.ct"

        # Write sequences
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

        sys.stderr.write(f"[DuplexFold] Running: {' '.join(cmd)}\n")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                text=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                f"[ERROR] RNAstructure DuplexFold failed (exit {e.returncode})\n"
                f"Command: {' '.join(cmd)}\n"
                f"STDOUT:\n{e.stdout}\n"
                f"STDERR:\n{e.stderr}\n"
            )
            raise

        try:
            with ct_path.open() as fh:
                raw_lines = [ln.rstrip("\n") for ln in fh]
        except FileNotFoundError:
            raise RuntimeError(
                f"DuplexFold did not produce CT file at {ct_path}. "
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}\n"
            )

        # Strip leading/trailing blank lines but keep internal structure
        lines = [ln for ln in raw_lines if ln.strip()]
        if not lines:
            sys.stderr.write("[DuplexFold] CT file empty; no pairs returned.\n")
            return []

        n1 = len(seq5)
        n2 = len(seq3)
        valid_bases = set("ACGUT")  # treat others (e.g. I) as non-nucleotide

        all_pair_sets: List[List[Tuple[int, int]]] = []
        i = 0

        while i < len(lines):
            header = lines[i].strip()
            if not header:
                i += 1
                continue

            header_fields = header.split()
            try:
                total_nt = int(header_fields[0])
            except Exception:
                # Not a CT header; stop parsing further
                sys.stderr.write(
                    f"[DuplexFold] Stopping CT parse at non-header line: {header!r}\n"
                )
                break

            # Ensure we have enough lines for this structure
            if i + 1 + total_nt > len(lines):
                sys.stderr.write(
                    "[DuplexFold] Incomplete CT block at end of file; "
                    "ignoring partial structure.\n"
                )
                break

            body = lines[i + 1 : i + 1 + total_nt]

            # --- Build CT index -> (strand, local_index) for this block --- #
            ct_map: Dict[int, Tuple[int, int] | None] = {}
            real_count = 0

            for ln in body:
                fields = ln.split()
                if len(fields) < 2:
                    continue
                try:
                    idx = int(fields[0])
                except ValueError:
                    continue
                base = fields[1].upper()

                if base in valid_bases:
                    real_count += 1
                    if real_count <= n1:
                        ct_map[idx] = (1, real_count)  # 5' strand
                    elif real_count <= n1 + n2:
                        ct_map[idx] = (2, real_count - n1)  # 3' strand
                    else:
                        ct_map[idx] = None
                else:
                    # 'I' or any other placeholder; not part of either RNA
                    ct_map[idx] = None

            if real_count != n1 + n2:
                sys.stderr.write(
                    f"[DuplexFold] WARNING: CT block at line {i+1} has {real_count} "
                    f"real bases, but len(seq5)+len(seq3)={n1+n2}. "
                    "Proceeding with mapped subset.\n"
                )

            # --- Extract inter-strand pairs for this structure --- #
            pair_set: set[Tuple[int, int]] = set()

            for ln in body:
                fields = ln.split()
                if len(fields) < 5:
                    continue
                try:
                    idx_i = int(fields[0])
                    idx_j = int(fields[4])
                except ValueError:
                    continue

                if idx_j <= 0:
                    continue

                info_i = ct_map.get(idx_i)
                info_j = ct_map.get(idx_j)
                if not info_i or not info_j:
                    continue

                strand_i, pos_i = info_i
                strand_j, pos_j = info_j
                if strand_i == strand_j:
                    continue

                if strand_i == 1 and strand_j == 2:
                    i5_local, j3_local = pos_i, pos_j
                elif strand_i == 2 and strand_j == 1:
                    i5_local, j3_local = pos_j, pos_i
                else:
                    continue

                if 1 <= i5_local <= n1 and 1 <= j3_local <= n2:
                    pair_set.add((i5_local, j3_local))

            all_pair_sets.append(sorted(pair_set))
            i += 1 + total_nt  # jump to next CT block

        sys.stderr.write(
            f"[DuplexFold] Parsed {len(all_pair_sets)} structure(s) from CT, "
            "each with local 5'/3' indices.\n"
        )
        return all_pair_sets


def refine_structure(
    struct: str,
    full_seq: str,
    fold_exe: Path,
    min_unpaired_len: int,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
    cache: Dict[Tuple[str, float | None, Tuple[str, ...]], str] | None = None,
) -> str:
    """
    Refine a single dot-bracket structure string by folding long INTERNAL unpaired regions.

    Terminal runs at the 5' or 3' ends are skipped here and are handled separately
    (optionally) via DuplexFold.
    """
    if len(struct) != len(full_seq):
        raise ValueError(
            f"Structure length ({len(struct)}) does not match sequence length ({len(full_seq)})"
        )

    if cache is None:
        cache = {}

    runs = find_unpaired_runs(struct, min_unpaired_len)
    if not runs:
        return struct

    new_struct = list(struct)
    n = len(struct)

    for start, end in runs:
        # Skip terminal runs here; DuplexFold handles them separately.
        if start == 0 and end < n:
            continue
        if start > 0 and end == n:
            continue

        subseq = full_seq[start:end]
        extra_args_tuple = tuple(extra_args) if extra_args is not None else tuple()
        key = (subseq, temperature, extra_args_tuple)

        sys.stderr.write(
            f"[Refine] Internal run [{start}:{end}] (len={end-start}); "
            f"subseq_len={len(subseq)}; cache_hit={key in cache}\n"
        )

        if key in cache:
            local_db = cache[key]
        else:
            local_db = call_rnastructure_fold(
                fold_exe=fold_exe,
                subseq=subseq,
                temperature=temperature,
                extra_args=extra_args,
            )
            cache[key] = local_db

        if len(local_db) != (end - start):
            raise RuntimeError(
                f"RNAstructure Fold output length {len(local_db)} != region length {end - start}"
            )

        # Splice local structure back in:
        # - If local_db[k] == '.', leave global position unchanged ('.').
        # - If local_db[k] in '()', write that character to the global structure.
        for offset, ch in enumerate(local_db):
            global_pos = start + offset
            if ch == UNPAIRED_CHAR:
                continue
            elif ch in "()":
                new_struct[global_pos] = ch
            else:
                raise RuntimeError(
                    f"Unexpected character {ch!r} in local DBN from RNAstructure Fold."
                )

    return "".join(new_struct)


def refine_terminal_ends_with_duplex(
    struct: str,
    full_seq: str,
    len_5: int,
    len_3: int,
    duplex_exe: Path,
    cache: Dict[Tuple[str, str, float | None, Tuple[str, ...]], List[List[Tuple[int, int]]]],
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
) -> List[str]:
    """
    Refine terminal unpaired regions using DuplexFold.

    Returns multiple structures, one per DuplexFold CT structure.
    If DuplexFold returns K alternative duplexes, the output list
    will have length K (or 1 if DuplexFold produced nothing).
    """
    n = len(struct)
    if len(full_seq) != n:
        raise ValueError("full_seq and struct must have same length")

    seq5 = full_seq[:len_5]
    seq3 = full_seq[n - len_3 :]

    extra_args_tuple = tuple(extra_args) if extra_args is not None else tuple()
    key = (seq5, seq3, temperature, extra_args_tuple)

    if key in cache:
        all_pair_sets = cache[key]
        sys.stderr.write(
            f"[DuplexRefine] Using cached DuplexFold result for len_5={len_5}, "
            f"len_3={len_3}: {len(all_pair_sets)} structure(s).\n"
        )
    else:
        all_pair_sets = call_rnastructure_duplexfold(
            duplex_exe=duplex_exe,
            seq5=seq5,
            seq3=seq3,
            temperature=temperature,
            extra_args=extra_args,
        )
        cache[key] = all_pair_sets

    if not all_pair_sets:
        # DuplexFold gave nothing; just return the original structure as-is
        sys.stderr.write(
            "[DuplexRefine] DuplexFold returned no structures; "
            "leaving terminal ends unmodified.\n"
        )
        return [struct]

    n1 = len(seq5)
    n2 = len(seq3)
    if n1 != len_5 or n2 != len_3:
        raise RuntimeError(
            f"Inconsistent terminal lengths: len_5={len_5}, len_3={len_3}, "
            f"but seq5={n1}, seq3={n2}"
        )

    def _is_canonical_pair(b5: str, b3: str) -> bool:
        b5 = b5.upper()
        b3 = b3.upper()
        return (b5, b3) in {
            ("A", "U"), ("U", "A"),
            ("G", "C"), ("C", "G"),
            ("G", "U"), ("U", "G"),  # wobble
        }

    results: List[str] = []

    for pair_set in all_pair_sets:
        new_struct = list(struct)

        for i5_local, j3_local in pair_set:
            # Local 1-based → global 0-based
            idx5_global = i5_local - 1
            idx3_global = n - n2 + (j3_local - 1)

            if not (0 <= idx5_global < n and 0 <= idx3_global < n):
                sys.stderr.write(
                    f"[WARN] Skipping out-of-bounds DuplexFold mapping: "
                    f"i5_local={i5_local}, j3_local={j3_local}, "
                    f"idx5_global={idx5_global}, idx3_global={idx3_global}, n={n}\n"
                )
                continue

            if (
                new_struct[idx5_global] == UNPAIRED_CHAR
                and new_struct[idx3_global] == UNPAIRED_CHAR
            ):
                b5 = full_seq[idx5_global]
                b3 = full_seq[idx3_global]
                if not _is_canonical_pair(b5, b3):
                    sys.stderr.write(
                        f"[WARN] Skipping non-canonical DuplexFold pair "
                        f"{b5}{idx5_global+1}–{b3}{idx3_global+1}\n"
                    )
                    continue

                new_struct[idx5_global] = "("
                new_struct[idx3_global] = ")"

        results.append("".join(new_struct))

    return results

def _check_balanced_parentheses(struct: str) -> bool:
    stack = []
    for i, ch in enumerate(struct):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                return False
            stack.pop()
    return not stack

def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Refine long unpaired regions in CaCoFold-sampled .db structures "
            "by locally folding internal regions with RNAstructure Fold and, "
            "optionally, refining unpaired 5'/3' terminal ends jointly with RNAstructure DuplexFold."
        )
    )
    parser.add_argument(
        "--sto",
        required=True,
        help="FASTA/SEQ file with the full ungapped sequence for this RNA.",
    )
    parser.add_argument(
        "--seq-name",
        default=None,
        help=(
            "Optional sequence name (FASTA header, up to first whitespace). "
            "If omitted, the first sequence in the file is used."
        ),
    )
    parser.add_argument(
        "--db-in",
        required=True,
        help=".db file from sample_cacofold_structures.py (one structure per line).",
    )
    parser.add_argument(
        "--db-out",
        required=True,
        help="Output .db file with refined structures.",
    )
    parser.add_argument(
        "--fold-exe",
        default="$PROJECT/repos/RNAstructure/exe/Fold",
        help="Path to RNAstructure Fold executable (default: $PROJECT/repos/RNAstructure/exe/Fold).",
    )
    parser.add_argument(
        "--duplex-exe",
        default="$PROJECT/repos/RNAstructure/exe/DuplexFold",
        help="Path to RNAstructure DuplexFold executable (default: $PROJECT/repos/RNAstructure/exe/DuplexFold).",
    )
    parser.add_argument(
        "--min-unpaired",
        type=int,
        default=4,
        help="Minimum length of an unpaired run ('.') to refine internally with Fold (default: 15).",
    )
    parser.add_argument(
        "--min-terminal-unpaired",
        type=int,
        default=2,
        help=(
            "Minimum length of unpaired 5' AND 3' ends ('.') to refine jointly with DuplexFold. "
            "If omitted, defaults to --min-unpaired."
        ),
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=None,
        help=(
            "Temperature in Kelvin for RNAstructure (-T). "
            "If omitted, uses RNAstructure defaults for both Fold and DuplexFold."
        ),
    )
    parser.add_argument(
        "--fold-extra-arg",
        action="append",
        default=[],
        help=(
            "Extra argument(s) to pass directly to Fold (e.g. --fold-extra-arg=-d "
            "to use DNA parameters). Can be specified multiple times."
        ),
    )
    parser.add_argument(
        "--duplex-extra-arg",
        action="append",
        default=[],
        help=(
            "Extra argument(s) to pass directly to DuplexFold. "
            "Can be specified multiple times."
        ),
    )

    args = parser.parse_args(argv)

    fold_exe_str = os.path.expandvars(args.fold_exe)
    fold_exe = Path(fold_exe_str)
    if not fold_exe.is_file():
        sys.stderr.write(f"[ERROR] Fold executable not found at: {fold_exe}\n")
        sys.exit(1)

    duplex_exe_str = os.path.expandvars(args.duplex_exe)
    duplex_exe = Path(duplex_exe_str)
    if not duplex_exe.is_file():
        sys.stderr.write(f"[ERROR] DuplexFold executable not found at: {duplex_exe}\n")
        sys.exit(1)

    db_in = Path(args.db_in)
    db_out = Path(args.db_out)

    # Decide where to get the sequence from:
    # 1) If --sto is provided, read ungapped sequence from the Stockholm file.
    # 2) Else, fall back to --seq-file FASTA/SEQ.
    if args.sto is not None:
        sto_path = Path(args.sto)
        if not sto_path.is_file():
            sys.stderr.write(f"[ERROR] Stockholm file not found at: {sto_path}\n")
            sys.exit(1)
        full_seq = read_ungapped_seq_from_sto(sto_path, args.seq_name)
    else:
        if args.seq_file is None:
            sys.stderr.write(
                "[ERROR] You must provide either --sto or --seq-file.\n"
            )
            sys.exit(1)
        seq_file = Path(args.seq_file)
        if not seq_file.is_file():
            sys.stderr.write(f"[ERROR] Sequence file not found at: {seq_file}\n")
            sys.exit(1)
        full_seq = read_fasta_sequence(seq_file, args.seq_name)

    structs = read_db_structures(db_in)

    if len(full_seq) != len(structs[0]):
        raise ValueError(
            f"Sequence length ({len(full_seq)}) != structure length ({len(structs[0])}). "
            "Make sure you are using the ungapped sequence corresponding to these structures."
        )

    if len(full_seq) != len(structs[0]):
        raise ValueError(
            f"Sequence length ({len(full_seq)}) != structure length ({len(structs[0])}). "
            "Make sure you are using the ungapped sequence corresponding to these structures."
        )

    min_terminal = args.min_terminal_unpaired or args.min_unpaired

    sys.stderr.write(
        f"[INFO] Sequence length          : {len(full_seq)}\n"
        f"[INFO] # of input structures     : {len(structs)}\n"
        f"[INFO] min internal unpaired len : {args.min_unpaired}\n"
        f"[INFO] min terminal unpaired len : {min_terminal}\n"
    )

    fold_cache: Dict[Tuple[str, float | None, Tuple[str, ...]], str] = {}
    duplex_cache: Dict[Tuple[str, str, float | None, Tuple[str, ...]], List[Tuple[int, int]]] = {}
    refined_structs: List[str] = []

    for idx, s in enumerate(structs):
        # 1) Internal refinement with Fold (excluding terminal runs)
        internal_runs = find_unpaired_runs(s, args.min_unpaired)
        internal_has_any = False
        n = len(s)
        for start, end in internal_runs:
            if start == 0 and end < n:
                continue
            if start > 0 and end == n:
                continue
            internal_has_any = True
            break

        if internal_has_any:
            sys.stderr.write(
                f"[INFO] Structure {idx}: refining internal long unpaired regions with Fold.\n"
            )
            refined = refine_structure(
                struct=s,
                full_seq=full_seq,
                fold_exe=fold_exe,
                min_unpaired_len=args.min_unpaired,
                temperature=310,
                extra_args=args.fold_extra_arg,
                cache=fold_cache,
            )
        else:
            refined = s

        # 2) Terminal 5'/3' refinement with DuplexFold
        len_5, len_3 = find_terminal_unpaired_ends(refined, min_terminal)
        if len_5 and len_3:
            sys.stderr.write(
                f"[INFO] Structure {idx}: refining 5'({len_5}) and 3'({len_3}) unpaired ends with DuplexFold.\n"
            )
            variants = refine_terminal_ends_with_duplex(
                struct=refined,
                full_seq=full_seq,
		len_5=len_5,
		len_3=len_3,
                duplex_exe=duplex_exe,
                temperature=310,
                extra_args=args.duplex_extra_arg,
                cache=duplex_cache,
            )

        refined_structs.append(variants)

    for idx, rs in enumerate(refined_structs):
        if not _check_balanced_parentheses(rs):
            sys.stderr.write(
                f"[WARN] Unbalanced parentheses in refined structure {idx}:\n{rs}\n"
            )


    with db_out.open("w") as out_f:
        out_f.write(full_seq + "\n")
        for rs in refined_structs:
            out_f.write(rs + "\n")

    sys.stderr.write(
        f"[INFO] Wrote {len(refined_structs)} refined structures to {db_out}\n"
    )


if __name__ == "__main__":
    main()

