"""RNAstructure AllSub wrapper (MFE + suboptimals) emitting `predictions.db`."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

from .fasta_utils import read_single_fasta, sanitize_rna_sequence_preserve_length


def ensure_datapath(*, exe_hint: Path | None = None) -> None:
    # Treat an explicitly exported (possibly empty) DATAPATH as "set".
    # This supports docker-wrapped RNAstructure executables, which provide their own in-container DATAPATH.
    if "DATAPATH" in os.environ:
        return
    # If we're using the docker-wrapped executables, do not auto-point DATAPATH to a host directory.
    # The wrapper will inject a valid in-container default when DATAPATH is empty/unset.
    if exe_hint is not None and "rnastructure_docker" in exe_hint.parts:
        os.environ["DATAPATH"] = ""
        return
    candidates: list[Path] = []
    if exe_hint is not None:
        candidates.append(exe_hint.parent.parent / "data_tables")
        candidates.append(exe_hint.parent / "data_tables")

    repo_root = Path(__file__).resolve().parents[4]
    candidates.append(repo_root / "RNAstructure" / "data_tables")

    for path in candidates:
        if path.is_dir():
            os.environ["DATAPATH"] = str(path)
            return


def ct_to_dotbracket(ct_text: str) -> list[str]:
    """Parse an RNAstructure CT file containing one or more structures."""
    lines = [ln.strip() for ln in ct_text.splitlines() if ln.strip()]
    out: list[str] = []
    idx = 0
    while idx < len(lines):
        header = lines[idx]
        parts = header.split()
        if not parts:
            idx += 1
            continue
        try:
            n = int(parts[0])
        except ValueError:
            break
        idx += 1
        pairs = [0] * n
        for _ in range(n):
            if idx >= len(lines):
                break
            cols = lines[idx].split()
            idx += 1
            if len(cols) < 5:
                continue
            try:
                i = int(cols[0])
                j = int(cols[4])
            except ValueError:
                continue
            if 1 <= i <= n:
                pairs[i - 1] = j
        db = ["."] * n
        for i, j in enumerate(pairs, start=1):
            if j > i and 1 <= j <= n:
                db[i - 1] = "("
                db[j - 1] = ")"
        out.append("".join(db))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RNAstructure AllSub and emit predictions.db")
    parser.add_argument("fasta", help="Input FASTA path (single sequence)")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("--k", type=int, default=100, help="Number of structures to emit (default: 100)")
    parser.add_argument(
        "--allsub-exe",
        type=Path,
        default=None,
        help="Path to RNAstructure AllSub (default: repo RNAstructure/exe/AllSub)",
    )
    parser.add_argument(
        "--fold-exe",
        type=Path,
        default=None,
        help="Path to RNAstructure Fold (fallback; default: repo RNAstructure/exe/Fold)",
    )
    parser.add_argument(
        "--abs",
        type=float,
        default=20.0,
        help="Absolute energy window for AllSub (-a, kcal/mol).",
    )
    parser.add_argument(
        "--max-allsub-len",
        type=int,
        default=None,
        help=(
            "If set, skip AllSub for sequences longer than this length and emit Fold(MFE) only. "
            "Useful for keeping runtime bounded on large benchmarks."
        ),
    )
    parser.add_argument(
        "--keep-ct",
        action="store_true",
        help="Keep the (potentially large) AllSub CT output in the per-target directory.",
    )
    parser.add_argument(
        "--max-seconds",
        type=float,
        default=60.0,
        help="Per-target wall clock budget for AllSub (default: 60s).",
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[4]
    allsub = args.allsub_exe or (repo_root / "RNAstructure" / "exe" / "AllSub")
    fold_exe = args.fold_exe or (repo_root / "RNAstructure" / "exe" / "Fold")

    rec = read_single_fasta(fasta_path.read_text())
    seq = sanitize_rna_sequence_preserve_length(rec.sequence)
    sanitized_fasta = out_dir / "input_sanitized.fa"
    sanitized_fasta.write_text(f">{rec.seq_id}\n{seq}\n")

    ensure_datapath(exe_hint=allsub)

    ct_path = out_dir / "allsub.ct"
    skipped_allsub = bool(args.max_allsub_len is not None and len(seq) > int(args.max_allsub_len))
    if skipped_allsub:
        out = "[SKIP] AllSub skipped due to --max-allsub-len\n"
        elapsed = 0.0
    else:
        t0 = time.time()
        cmd = [str(allsub), str(sanitized_fasta), str(ct_path), "-a", str(float(args.abs))]
        try:
            out = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                timeout=float(args.max_seconds) if args.max_seconds is not None else None,
            ).stdout
        except subprocess.TimeoutExpired:
            out = f"[TIMEOUT] AllSub exceeded {args.max_seconds} seconds\n"
        elapsed = time.time() - t0

    structs: list[str] = []
    if (not skipped_allsub) and ct_path.exists():
        try:
            structs = ct_to_dotbracket(ct_path.read_text())
        except Exception:
            structs = []
    if ct_path.exists() and not bool(args.keep_ct):
        try:
            ct_path.unlink()
        except Exception:
            pass

    used_fallback_fold = False
    if not structs:
        # AllSub can fail on some targets (or time out). Fall back to MFE (Fold) so the
        # baseline at least produces a valid top-1 structure.
        mfe_ct = out_dir / "mfe.ct"
        try:
            fold_out = subprocess.run(
                [str(fold_exe), str(sanitized_fasta), str(mfe_ct)],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                timeout=float(args.max_seconds) if args.max_seconds is not None else None,
            ).stdout
            _ = fold_out
            if mfe_ct.exists():
                structs = ct_to_dotbracket(mfe_ct.read_text())[:1]
                used_fallback_fold = bool(structs)
        except subprocess.TimeoutExpired:
            structs = []
        finally:
            if mfe_ct.exists():
                try:
                    mfe_ct.unlink()
                except Exception:
                    pass

    # Dedupe while preserving order and enforce length.
    seen: set[str] = set()
    filtered: list[str] = []
    for s in structs:
        if len(s) != len(seq):
            continue
        if s in seen:
            continue
        seen.add(s)
        filtered.append(s)
        if len(filtered) >= int(args.k):
            break

    if not filtered:
        filtered = ["." * len(seq)]
    while len(filtered) < int(args.k):
        filtered.append("." * len(seq))
    filtered = filtered[: int(args.k)]

    (out_dir / "predictions.db").write_text("\n".join([seq] + filtered) + "\n")
    (out_dir / "meta.json").write_text(
        json.dumps(
            {
                "method": "RNAstructure AllSub",
                "allsub_exe": str(allsub),
                "fold_exe": str(fold_exe),
                "abs": float(args.abs),
                "k": int(args.k),
                "max_allsub_len": int(args.max_allsub_len) if args.max_allsub_len is not None else None,
                "skipped_allsub": skipped_allsub,
                "elapsed_seconds": elapsed,
                "n_unique_structs": int(len(set(filtered))),
                "used_fallback_fold": used_fallback_fold,
                "stdout_head": out.splitlines()[:50],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
