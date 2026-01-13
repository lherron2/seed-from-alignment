"""Command-based predictor runner."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .parse_dotbracket import parse_dotbracket


def run_predictor(
    predictor_cmd: str,
    fasta_path: str | Path,
    out_dir: str | Path,
    k: int,
) -> list[str]:
    """Run predictor command template and collect dot-bracket outputs."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = predictor_cmd.format(fasta=str(fasta_path), outdir=str(out_dir))
    subprocess.run(cmd, shell=True, check=True)

    dotbrackets: list[str] = []
    # Prefer the predictor's final output if present.
    # Many predictors (including ours) also emit intermediate `00_*.db`, `01_*.db`, etc.
    preferred = out_dir / "predictions.db"
    if preferred.exists():
        _seq, structs = parse_dotbracket(preferred.read_text())
        return structs[:k]

    # Strategy A: any .db files (excluding common intermediate prefixes)
    for path in sorted(out_dir.glob("*.db")):
        if path.name.startswith(("00_", "01_", "02_", "03_", "04_")):
            continue
        _seq, structs = parse_dotbracket(path.read_text())
        dotbrackets.extend(structs)

    # Strategy B: any .dot or .txt files
    for path in sorted(out_dir.glob("*.dot")):
        _seq, structs = parse_dotbracket(path.read_text())
        dotbrackets.extend(structs)

    return dotbrackets[:k]
