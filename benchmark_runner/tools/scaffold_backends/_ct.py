from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CtStruct:
    struct: str
    energy: float
    name: str = "structure"


def pairs_from_db(db: str) -> list[tuple[int, int]]:
    stack: list[int] = []
    out: list[tuple[int, int]] = []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                continue
            j = stack.pop()
            out.append((j, i) if j < i else (i, j))
    return out


def write_ct(*, seq: str, structs: list[CtStruct], out_path: Path) -> None:
    n = len(seq)
    lines: list[str] = []
    for idx, item in enumerate(structs, start=1):
        energy = float(item.energy)
        lines.append(f"{n} ENERGY = {energy:.6f} {item.name} {idx}")
        pair_map = {i: 0 for i in range(n)}
        for i, j in pairs_from_db(item.struct):
            pair_map[i] = j + 1
            pair_map[j] = i + 1
        for i in range(n):
            base = seq[i]
            prev_i = i if i > 0 else 0
            next_i = i + 2 if i + 1 < n else 0
            pair_i = pair_map[i]
            lines.append(f"{i+1} {base} {prev_i} {next_i} {pair_i} {i+1}")
    out_path.write_text("\n".join(lines) + "\n")

