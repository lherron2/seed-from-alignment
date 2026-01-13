"""
RNAstructure partition function utilities.
"""

from __future__ import annotations

import math
import subprocess
import tempfile
from collections.abc import Sequence
from pathlib import Path


def _write_fasta(seq: str, path: Path) -> None:
    path.write_text(">seq\n" + seq + "\n")


def run_partition(
    partition_exe: Path,
    fasta_path: Path,
    pfs_path: Path,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
    env: dict[str, str] | None = None,
) -> None:
    cmd = [str(partition_exe), str(fasta_path), str(pfs_path)]
    if temperature is not None:
        cmd += ["-t", str(temperature)]
    if extra_args:
        cmd += list(extra_args)
    subprocess.run(cmd, check=True, text=True, capture_output=True, env=env)


def run_probability_plot(
    probplot_exe: Path,
    pfs_path: Path,
    out_txt: Path,
    extra_args: Sequence[str] | None = None,
    env: dict[str, str] | None = None,
) -> None:
    cmd = [str(probplot_exe), str(pfs_path), str(out_txt)]
    if extra_args:
        cmd += list(extra_args)
    subprocess.run(cmd, check=True, text=True, capture_output=True, env=env)


def parse_probability_plot(text: str) -> dict[tuple[int, int], float]:
    pairs: dict[tuple[int, int], float] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", ";")):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            i = int(parts[0]) - 1
            j = int(parts[1]) - 1
            val = float(parts[2])
        except ValueError:
            continue
        if i < 0 or j < 0 or i == j:
            continue
        if 0.0 <= val <= 1.0:
            p = val
        elif val > 1.0:
            p = math.pow(10.0, -val)
        else:
            continue
        p = max(0.0, min(1.0, p))
        if i > j:
            i, j = j, i
        pairs[(i, j)] = p
    return pairs


def call_rnastructure_pf_probs(
    partition_exe: Path,
    probplot_exe: Path,
    seq: str,
    temperature: float | None = None,
    extra_args: Sequence[str] | None = None,
    env: dict[str, str] | None = None,
) -> dict[tuple[int, int], float]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        fasta_path = tmpdir_path / "seq.fa"
        pfs_path = tmpdir_path / "out.pfs"
        out_txt = tmpdir_path / "prob.txt"
        _write_fasta(seq, fasta_path)
        run_partition(
            partition_exe=partition_exe,
            fasta_path=fasta_path,
            pfs_path=pfs_path,
            temperature=temperature,
            extra_args=extra_args,
            env=env,
        )
        run_probability_plot(
            probplot_exe=probplot_exe,
            pfs_path=pfs_path,
            out_txt=out_txt,
            extra_args=None,
            env=env,
        )
        if not out_txt.is_file():
            return {}
        return parse_probability_plot(out_txt.read_text())
