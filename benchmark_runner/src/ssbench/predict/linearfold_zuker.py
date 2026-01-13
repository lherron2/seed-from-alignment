"""LinearFold wrapper that emits a top-1 structure plus Zuker suboptimals.

This is designed to plug into `ssbench.cli predict run` as a command-based predictor.
It reads a single-sequence FASTA and writes `predictions.db` in `out_dir`.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

from .fasta_utils import read_single_fasta, sanitize_rna_sequence_preserve_length


def _is_dotbracket(token: str, length: int) -> bool:
    if len(token) != length:
        return False
    allowed = set(".()[]{}<>")
    return all(ch in allowed for ch in token)


def _run_linearfold_stream(
    *,
    linearfold: Path,
    docker_image: str | None,
    seq_id: str,
    seq: str,
    k: int,
    beamsize: int,
    delta: float,
    max_seconds: float | None,
) -> tuple[list[str], str]:
    if docker_image:
        uid = os.getuid()
        gid = os.getgid()
        cmd = [
            "docker",
            "run",
            "--rm",
            "-i",
            "--user",
            f"{uid}:{gid}",
            docker_image,
            "--fasta",
            "-V",
            "-b",
            str(beamsize),
            "--zuker",
            "--delta",
            str(delta),
        ]
    else:
        cmd = [
            str(linearfold),
            "--fasta",
            "-V",
            "-b",
            str(beamsize),
            "--zuker",
            "--delta",
            str(delta),
        ]
    fasta_text = f">{seq_id}\n{seq}\n"

    start = time.time()
    env = os.environ.copy()
    env.setdefault("PYTHONWARNINGS", "ignore")
    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )
    assert proc.stdin is not None
    assert proc.stdout is not None
    proc.stdin.write(fasta_text)
    proc.stdin.close()

    structs: list[str] = []
    seen: set[str] = set()
    captured_lines: list[str] = []

    try:
        for line in proc.stdout:
            captured_lines.append(line)
            if max_seconds is not None and (time.time() - start) > max_seconds:
                break
            stripped = line.strip()
            if not stripped:
                continue
            token = stripped.split()[0]
            if not _is_dotbracket(token, len(seq)):
                continue
            if token in seen:
                continue
            seen.add(token)
            structs.append(token)
            if len(structs) >= k:
                break
    finally:
        try:
            if proc.poll() is None:
                proc.terminate()
        except Exception:
            pass
        try:
            proc.wait(timeout=2.0)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass
            try:
                proc.wait(timeout=2.0)
            except Exception:
                pass

    return structs, "".join(captured_lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run LinearFold (+ Zuker suboptimals) and emit predictions.db")
    parser.add_argument("fasta", help="Input FASTA path (single sequence)")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("--k", type=int, default=100, help="Number of structures to emit (default: 100)")
    parser.add_argument("--beamsize", type=int, default=100, help="LinearFold beam size (default: 100)")
    parser.add_argument("--delta", type=float, default=80.0, help="Zuker delta (kcal/mol) for suboptimals")
    parser.add_argument(
        "--linearfold",
        type=Path,
        default=None,
        help="Path to LinearFold wrapper script (default: repo external/LinearFold/linearfold)",
    )
    parser.add_argument(
        "--docker-image",
        default=None,
        help="Run LinearFold via Docker (image name/tag). If set, --linearfold is ignored.",
    )
    parser.add_argument(
        "--max-seconds",
        type=float,
        default=30.0,
        help="Per-target wall clock budget; stops collecting beyond this (default: 30s)",
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[4]
    linearfold = args.linearfold or (repo_root / "external" / "LinearFold" / "linearfold")
    docker_image = str(args.docker_image) if args.docker_image else os.environ.get("LINEARFOLD_DOCKER_IMAGE")

    rec = read_single_fasta(fasta_path.read_text())
    seq = sanitize_rna_sequence_preserve_length(rec.sequence)

    t0 = time.time()
    structs, raw_out = _run_linearfold_stream(
        linearfold=linearfold,
        docker_image=docker_image,
        seq_id=rec.seq_id,
        seq=seq,
        k=int(args.k),
        beamsize=int(args.beamsize),
        delta=float(args.delta),
        max_seconds=float(args.max_seconds) if args.max_seconds is not None else None,
    )
    elapsed = time.time() - t0

    if not structs:
        structs = ["." * len(seq)]
    while len(structs) < int(args.k):
        structs.append("." * len(seq))
    structs = structs[: int(args.k)]

    (out_dir / "predictions.db").write_text("\n".join([seq] + structs) + "\n")
    (out_dir / "meta.json").write_text(
        json.dumps(
            {
                "method": "LinearFold-V (Zuker suboptimals)",
                "linearfold": str(linearfold),
                "docker_image": docker_image or "",
                "beamsize": int(args.beamsize),
                "delta": float(args.delta),
                "k": int(args.k),
                "elapsed_seconds": elapsed,
                "n_unique_structs": int(len(set(structs))),
                "stdout_head": raw_out.splitlines()[:50],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
