#!/usr/bin/env python3
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


def read_fasta_sequence(path: Path, seq_name: str | None = None) -> str:
    """Read the first sequence (or the one matching seq_name) from a FASTA/SEQ file."""
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
        return seqs[name].upper()

    if seq_name not in seqs:
        raise ValueError(
            f"Sequence '{seq_name}' not found in {path}. "
            f"Available: {', '.join(seqs.keys())}"
        )
    return seqs[seq_name].upper()


def read_db_structures(path: Path) -> List[str]:
    """Read a .db file that has ONE structure per non-empty line."""
    structs: List[str] = []
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            structs.append(s)
    if not structs:
        raise ValueError(f"No structures found in {path}")
    # sanity: all same length
    lengths = {len(s) for s in structs}
    if len(lengths) != 1:
        raise ValueError(f"Structures in {path} have inconsistent lengths: {lengths}")
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

    Returns (len_5prime, len_3prime). If either is < min_len, or if the entire
    sequence is a single unpaired run, returns (0, 0).
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

    if len_5 >= min_len and len_3 >= min_len:
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

    Assumes that the bracket output has a line with only '().' of length len(subseq),
    and returns the last such line.
    """
    extra_args = list(extra_args) if extra_args is not None else []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        seq_path = tmpdir_path / "subseq.fa"

        # Write subsequence as a one-record FASTA
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
) -> List[Tuple[int, int]]:
    """
    Call RNAstructure DuplexFold on two terminal subsequences.

    Returns a list of base pairs as (i, j) in 1-based duplex coordinates, where
    indices run from 1..len(seq5)+len(seq3), with seq5 first, then seq3.

    We only keep INTER-STRAND pairs (one index <= len(seq5), the other > len(seq5)).
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

        cmd: List[str] = [
            str(duplex_exe),
            str(seq5_path),
            str(seq3_path),
            str(ct_path),
        ]
        if temperature is not None:
            cmd.extend(["-T", str(temperature)])
        cmd.extend(extra_args)

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

        # Parse the CT file for base pairs.
        try:
            with ct_path.open() as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
        except FileNotFoundError:
            raise RuntimeError(
                f"DuplexFold did not produce CT file at {ct_path}. "
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}\n"
            )

        if not lines:
            return []

        # First line is header: "<N> <energy> ..."
        try:
            header_fields = lines[0].split()
            total_nt = int(header_fields[0])
        except Exception:
            raise RuntimeError(
                f"Could not parse CT header line: {lines[0]!r} "
                f"in file {ct_path}"
            )

        n1 = len(seq5)
        n2 = len(seq3)
        if total_nt != n1 + n2:
            sys.stderr.write(
                f"[WARN] CT length {total_nt} != len(seq5)+len(seq3)={n1+n2}. "
                "Proceeding anyway.\n"
            )

        pairs: List[Tuple[int, int]] = []
        max_valid = n1 + n2  # we only trust indices that map into our two strands

        for line in lines[1:]:
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                i = int(fields[0])
                j = int(fields[4])
            except ValueError:
                continue

            # Ignore unpaired or obviously bad entries
            if i <= 0 or j <= 0:
                continue
            if i > total_nt or j > total_nt:
                continue

            if i > j:
                i, j = j, i

            # Only keep inter-strand pairs whose indices also fit into our seq lengths
            if i <= n1 < j <= max_valid:
                pairs.append((i, j))

        return pairs

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

    all_runs = find_unpaired_runs(struct, min_unpaired_len)
    if not all_runs:
        return struct

    n = len(struct)
    internal_runs: List[Tuple[int, int]] = []
    for start, end in all_runs:
        # Entire sequence unpaired: treat as an internal run.
        if start == 0 and end == n:
            internal_runs.append((start, end))
        # 5' terminal run: start == 0, end < n  -> skip
        elif start == 0 and end < n:
            continue
        # 3' terminal run: start > 0, end == n  -> skip
        elif start > 0 and end == n:
            continue
        else:
            internal_runs.append((start, end))

    if not internal_runs:
        # Nothing internal to refine
        return struct

    new_struct = list(struct)
    extra_args_tuple = tuple(extra_args) if extra_args is not None else tuple()

    for start, end in internal_runs:
        subseq = full_seq[start:end]
        key = (subseq, temperature, extra_args_tuple)

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
    duplex_exe: Path,
    min_terminal_unpaired_len: int,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
    cache: Dict[Tuple[str, str, float | None, Tuple[str, ...]], List[Tuple[int, int]]] | None = None,
) -> str:
    """
    If BOTH the 5' and 3' ends have long unpaired runs, refine them jointly using DuplexFold.

    The 5' terminal region (positions [0:len_5)) and 3' terminal region
    (positions [n-len_3:n)) are extracted as separate strands and passed to DuplexFold.
    Inter-strand base pairs are then mapped back into the global dot-bracket structure,
    assigning '(' to 5' positions and ')' to 3' positions.
    """
    if len(struct) != len(full_seq):
        raise ValueError(
            f"Structure length ({len(struct)}) does not match sequence length ({len(full_seq)})"
        )

    if cache is None:
        cache = {}

    n = len(struct)
    len_5, len_3 = find_terminal_unpaired_ends(struct, min_terminal_unpaired_len)
    if len_5 == 0 or len_3 == 0:
        # Nothing to do
        return struct

    seq5 = full_seq[:len_5]
    seq3 = full_seq[n - len_3 :]

    extra_args_tuple = tuple(extra_args) if extra_args is not None else tuple()
    key = (seq5, seq3, temperature, extra_args_tuple)

    if key in cache:
        pairs = cache[key]
    else:
        pairs = call_rnastructure_duplexfold(
            duplex_exe=duplex_exe,
            seq5=seq5,
            seq3=seq3,
            temperature=temperature,
            extra_args=extra_args,
        )
        cache[key] = pairs

    if not pairs:
        return struct

    new_struct = list(struct)
    n1 = len(seq5)
    n2 = len(seq3)
    if n1 != len_5 or n2 != len_3:
        raise RuntimeError(
            f"Inconsistent terminal lengths: len_5={len_5}, len_3={len_3}, "
            f"but seq5={n1}, seq3={n2}"
        )

    # Map duplex indices back to global indices
    for i, j in pairs:
        # Inter-strand pairs were filtered in call_rnastructure_duplexfold;
        # still handle both orientations just in case.
        if i <= n1 < j:
            idx5_local = i - 1
            idx3_local = j - n1 - 1
        elif j <= n1 < i:
            idx5_local = j - 1
            idx3_local = i - n1 - 1
        else:
            # Shouldn't happen, but skip intra-strand pairs if present
            continue

        # Extra safety: make sure the local indices are in bounds for the
        # terminal segments we actually extracted.
        if not (0 <= idx5_local < n1 and 0 <= idx3_local < n2):
            # You could optionally log here if you want to see how often this happens:
            sys.stderr.write(
                f"[WARN] Skipping inconsistent DuplexFold pair (i={i}, j={j}) "
                f"-> idx5_local={idx5_local}, idx3_local={idx3_local}, "
                f"n1={n1}, n2={n2}\n"
            )
            continue

        idx5_global = idx5_local
        idx3_global = n - n2 + idx3_local

        # Only write the pair if BOTH ends are unpaired; otherwise skip the pair entirely.
        if (
            new_struct[idx5_global] == UNPAIRED_CHAR
            and new_struct[idx3_global] == UNPAIRED_CHAR
        ):
            new_struct[idx5_global] = "("
            new_struct[idx3_global] = ")"
        # else: conflict with existing structure -> ignore this DuplexFold pair

    return "".join(new_struct)

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

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Refine long unpaired regions in CaCoFold-sampled .db structures "
            "by locally folding internal regions with RNAstructure Fold and, "
            "optionally, refining unpaired 5'/3' terminal ends jointly with RNAstructure DuplexFold."
        )
    )
    parser.add_argument(
        "--seq-file",
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

    args = parser.parse_args()

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

    seq_file = Path(args.seq_file)
    db_in = Path(args.db_in)
    db_out = Path(args.db_out)

    full_seq = read_fasta_sequence(seq_file, args.seq_name)
    structs = read_db_structures(db_in)

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
                temperature=args.temperature,
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
            refined = refine_terminal_ends_with_duplex(
                struct=refined,
                full_seq=full_seq,
                duplex_exe=duplex_exe,
                min_terminal_unpaired_len=min_terminal,
                temperature=args.temperature,
                extra_args=args.duplex_extra_arg,
                cache=duplex_cache,
            )

        refined_structs.append(refined)

    for idx, rs in enumerate(refined_structs):
        if not _check_balanced_parentheses(rs):
            sys.stderr.write(
                f"[WARN] Unbalanced parentheses in refined structure {idx}:\n{rs}\n"
            )


    with db_out.open("w") as out_f:
        for rs in refined_structs:
            out_f.write(rs + "\n")

    sys.stderr.write(f"[INFO] Wrote {len(refined_structs)} refined structures to {db_out}\n")


if __name__ == "__main__":
    main()

