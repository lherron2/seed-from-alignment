"""EternaFold wrapper that emits a top-1 MEA structure plus suboptimal samples.

Top-1: `contrafold predict` (MEA).
@100: top-1 + 99 samples from `contrafold sample`.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from .fasta_utils import read_single_fasta, sanitize_rna_sequence_preserve_length


def _is_dotbracket_line(line: str, length: int) -> bool:
    s = line.strip()
    if len(s) != length:
        return False
    allowed = set(".()")
    return all(ch in allowed for ch in s)


def _extract_first_struct(output: str, length: int) -> str | None:
    for line in output.splitlines():
        if _is_dotbracket_line(line, length):
            return line.strip()
    return None


def _extract_structs(output: str, length: int) -> list[str]:
    out: list[str] = []
    for line in output.splitlines():
        if _is_dotbracket_line(line, length):
            out.append(line.strip())
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run EternaFold and emit predictions.db")
    parser.add_argument("fasta", help="Input FASTA path (single sequence)")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("--k", type=int, default=100, help="Number of structures to emit (default: 100)")
    parser.add_argument(
        "--contrafold",
        type=Path,
        default=None,
        help="Path to EternaFold contrafold binary (default: repo external/EternaFold/src/contrafold)",
    )
    parser.add_argument(
        "--params",
        type=Path,
        default=None,
        help="Path to EternaFold params file (default: repo external/EternaFold/parameters/EternaFoldParams.v1)",
    )
    parser.add_argument(
        "--sample",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include stochastic samples for suboptimal structures (default: true).",
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[4]
    contrafold = args.contrafold or (repo_root / "external" / "EternaFold" / "src" / "contrafold")
    params = args.params or (repo_root / "external" / "EternaFold" / "parameters" / "EternaFoldParams.v1")

    rec = read_single_fasta(fasta_path.read_text())
    seq = sanitize_rna_sequence_preserve_length(rec.sequence)
    sanitized_fasta = out_dir / "input_sanitized.fa"
    sanitized_fasta.write_text(f">{rec.seq_id}\n{seq}\n")

    # Top-1 MEA prediction
    t0 = time.time()
    pred_out = subprocess.run(
        [str(contrafold), "predict", str(sanitized_fasta), "--params", str(params)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    ).stdout
    top1 = _extract_first_struct(pred_out, len(seq))

    # Suboptimal structures: stochastic sampling
    samples: list[str] = []
    if bool(args.sample):
        n_samples = max(0, int(args.k) - 1)
        sample_out = subprocess.run(
            [
                str(contrafold),
                "sample",
                str(sanitized_fasta),
                "--params",
                str(params),
                "--nsamples",
                str(n_samples),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        ).stdout
        samples = _extract_structs(sample_out, len(seq))[:n_samples]

    elapsed = time.time() - t0

    structs: list[str] = []
    if top1 is not None:
        structs.append(top1)
    structs.extend(samples)
    if not structs:
        structs = ["." * len(seq)]
    while len(structs) < int(args.k):
        structs.append("." * len(seq))
    structs = structs[: int(args.k)]

    (out_dir / "predictions.db").write_text("\n".join([seq] + structs) + "\n")
    (out_dir / "meta.json").write_text(
        json.dumps(
            {
                "method": "EternaFold (MEA + samples)",
                "contrafold": str(contrafold),
                "params": str(params),
                "k": int(args.k),
                "elapsed_seconds": elapsed,
                "n_unique_structs": int(len(set(structs))),
                "predict_stdout_head": pred_out.splitlines()[:50],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()

