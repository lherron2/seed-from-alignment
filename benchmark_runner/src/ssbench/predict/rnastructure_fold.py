"""RNAstructure Fold (MFE) wrapper emitting `predictions.db`.

This is a lightweight baseline: top-1 is the MFE structure, and the remaining
entries are padded to reach K structures (default: 100).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

from .fasta_utils import read_single_fasta, sanitize_rna_sequence_preserve_length
from .rnastructure_allsub import ct_to_dotbracket, ensure_datapath


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RNAstructure Fold and emit predictions.db")
    parser.add_argument("fasta", help="Input FASTA path (single sequence)")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("--k", type=int, default=100, help="Number of structures to emit (default: 100)")
    parser.add_argument(
        "--fold-exe",
        type=Path,
        default=None,
        help="Path to RNAstructure Fold (default: repo RNAstructure/exe/Fold)",
    )
    parser.add_argument(
        "--max-seconds",
        type=float,
        default=30.0,
        help="Per-target wall clock budget for Fold (default: 30s).",
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[4]
    fold_exe = args.fold_exe or (repo_root / "RNAstructure" / "exe" / "Fold")

    rec = read_single_fasta(fasta_path.read_text())
    seq = sanitize_rna_sequence_preserve_length(rec.sequence)
    sanitized_fasta = out_dir / "input_sanitized.fa"
    sanitized_fasta.write_text(f">{rec.seq_id}\n{seq}\n")

    ensure_datapath(exe_hint=fold_exe)

    ct_path = out_dir / "mfe.ct"
    t0 = time.time()
    try:
        out = subprocess.run(
            [str(fold_exe), str(sanitized_fasta), str(ct_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=float(args.max_seconds) if args.max_seconds is not None else None,
        ).stdout
    except subprocess.TimeoutExpired:
        out = f"[TIMEOUT] Fold exceeded {args.max_seconds} seconds\n"
    elapsed = time.time() - t0

    structs: list[str] = []
    if ct_path.exists():
        try:
            structs = ct_to_dotbracket(ct_path.read_text())[:1]
        except Exception:
            structs = []
        try:
            ct_path.unlink()
        except Exception:
            pass

    if not structs:
        structs = ["." * len(seq)]

    # Pad to K.
    while len(structs) < int(args.k):
        structs.append(structs[0])
    structs = structs[: int(args.k)]

    (out_dir / "predictions.db").write_text("\n".join([seq] + structs) + "\n")
    (out_dir / "meta.json").write_text(
        json.dumps(
            {
                "method": "RNAstructure Fold (MFE)",
                "fold_exe": str(fold_exe),
                "k": int(args.k),
                "elapsed_seconds": elapsed,
                "datapath": os.environ.get("DATAPATH", ""),
                "stdout_head": out.splitlines()[:50],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()

