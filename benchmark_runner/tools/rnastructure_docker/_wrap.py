from __future__ import annotations

import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Mount:
    host_dir: Path
    container_dir: str


def _is_probable_path(token: str) -> bool:
    if token in {"-", ""}:
        return False
    if token.startswith("-"):
        return False
    # Avoid treating numeric CLI arguments (e.g. "-a 80.0") as paths.
    if re.fullmatch(r"[+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", token):
        return False
    # Keep this permissive: these wrappers are only used for RNAstructure tools
    # where file args come first.
    return any(sep in token for sep in ("/", "\\")) or "." in Path(token).name


def _resolve_mounts(argv: list[str]) -> tuple[list[Mount], dict[str, str]]:
    mounts: list[Mount] = []
    mapping: dict[str, str] = {}

    # Collect unique parent dirs for any file-like args that exist or look like paths.
    parents: list[Path] = []
    for tok in argv:
        if not _is_probable_path(tok):
            continue
        p = Path(tok)
        parent = p.parent if p.parent != Path("") else Path(".")
        parent = parent.resolve()
        if parent not in parents:
            parents.append(parent)

    for idx, host_dir in enumerate(parents):
        container_dir = f"/mnt/{idx}"
        mounts.append(Mount(host_dir=host_dir, container_dir=container_dir))

    # Map each path token to the first mount that contains it (or its parent for non-existent output paths).
    for tok in argv:
        if not _is_probable_path(tok):
            continue
        p = Path(tok)
        p_abs = p.resolve() if p.exists() else (p.parent.resolve() / p.name)
        chosen: Mount | None = None
        for m in mounts:
            try:
                p_abs.relative_to(m.host_dir)
            except ValueError:
                continue
            chosen = m
            rel = p_abs.relative_to(m.host_dir)
            mapping[tok] = str(Path(m.container_dir) / rel)
            break
        if chosen is None:
            # Fall back: mount the file's parent itself.
            host_dir = p_abs.parent
            container_dir = f"/mnt/{len(mounts)}"
            m = Mount(host_dir=host_dir, container_dir=container_dir)
            mounts.append(m)
            mapping[tok] = str(Path(container_dir) / p_abs.name)

    return mounts, mapping


def main(tool: str) -> None:
    image = os.environ.get("RNASTRUCTURE_DOCKER_IMAGE", "rnanneal/rnastructure:v6.4.0")
    uid = os.getuid()
    gid = os.getgid()

    tool_args = sys.argv[1:]
    mounts, mapping = _resolve_mounts(tool_args)

    translated: list[str] = []
    for tok in tool_args:
        translated.append(mapping.get(tok, tok))

    cmd = ["docker", "run", "--rm", "--user", f"{uid}:{gid}"]

    # Optional resource limits to avoid overwhelming constrained environments (e.g. WSL).
    cpus = os.environ.get("RNASTRUCTURE_DOCKER_CPUS")
    if cpus:
        cmd += ["--cpus", str(cpus)]
    memory = os.environ.get("RNASTRUCTURE_DOCKER_MEMORY")
    if memory:
        # Keep swap equal to memory by default to provide a hard ceiling.
        cmd += ["--memory", str(memory), "--memory-swap", str(memory)]
    pids_limit = os.environ.get("RNASTRUCTURE_DOCKER_PIDS_LIMIT")
    if pids_limit:
        cmd += ["--pids-limit", str(pids_limit)]
    for m in mounts:
        cmd += ["-v", f"{m.host_dir}:{m.container_dir}"]

    # Some RNAstructure binaries depend on DATAPATH; allow user override but provide a sane default
    # that works with the rnanneal/rnastructure image. Treat empty as unset.
    datapath = os.environ.get("DATAPATH")
    if not datapath:
        datapath = "/opt/RNAstructure/data_tables"
    cmd += ["-e", f"DATAPATH={datapath}"]

    cmd += [image, tool, *translated]
    res = subprocess.run(cmd, text=True)
    raise SystemExit(res.returncode)
