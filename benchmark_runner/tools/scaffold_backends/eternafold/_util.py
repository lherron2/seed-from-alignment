from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _docker import docker_run_in_dir  # noqa: E402


def extract_db_lines(text: str, length: int) -> list[str]:
    out: list[str] = []
    for ln in text.splitlines():
        s = ln.strip()
        if len(s) != length:
            continue
        if any(ch not in ".()" for ch in s):
            continue
        out.append(s)
    return out


def _ensure_local_datapath(env: dict[str, str]) -> dict[str, str]:
    if env.get("DATAPATH"):
        return env
    repo_root = Path(__file__).resolve().parents[4]
    tables = repo_root / "RNAstructure" / "data_tables"
    if tables.is_dir():
        env = dict(env)
        env["DATAPATH"] = str(tables)
    return env


def score_with_efn2(*, rnastructure_image: str | None, work_dir: Path, ct_name: str) -> list[float]:
    if rnastructure_image:
        rc, out = docker_run_in_dir(
            image=rnastructure_image,
            work_dir=work_dir,
            args=["efn2", ct_name, "-", "-s", "-p", "-q"],
        )
        if rc != 0:
            return []
    else:
        repo_root = Path(__file__).resolve().parents[4]
        efn2 = Path(os.environ.get("RNASTRUCTURE_EFN2_BIN") or (repo_root / "RNAstructure" / "exe" / "efn2"))
        proc = subprocess.run(
            [str(efn2), str(Path(ct_name)), "-", "-s", "-p", "-q"],
            cwd=str(work_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            env=_ensure_local_datapath(os.environ.copy()),
        )
        if proc.returncode != 0:
            return []
        out = proc.stdout
    energies: list[float] = []
    for ln in out.splitlines():
        ln = ln.strip()
        if not ln.startswith("Structure:"):
            continue
        if "Energy =" not in ln:
            continue
        try:
            rhs = ln.split("Energy =")[1].strip()
            val = rhs.split()[0]
            energies.append(float(val))
        except Exception:
            continue
    return energies
