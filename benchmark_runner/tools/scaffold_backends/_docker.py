from __future__ import annotations

import os
import subprocess
from pathlib import Path


def docker_run_stream(*, image: str, args: list[str], stdin_text: str) -> tuple[int, str]:
    uid = os.getuid()
    gid = os.getgid()
    proc = subprocess.Popen(
        ["docker", "run", "--rm", "-i", "--user", f"{uid}:{gid}", image, *args],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert proc.stdin is not None
    assert proc.stdout is not None
    proc.stdin.write(stdin_text)
    proc.stdin.close()
    out = proc.stdout.read()
    rc = proc.wait()
    return rc, out


def docker_run_in_dir(*, image: str, work_dir: Path, args: list[str]) -> tuple[int, str]:
    uid = os.getuid()
    gid = os.getgid()
    work_dir = work_dir.resolve()
    proc = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "--user",
            f"{uid}:{gid}",
            "-v",
            f"{work_dir}:/work",
            "-w",
            "/work",
            image,
            *args,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return int(proc.returncode), proc.stdout

