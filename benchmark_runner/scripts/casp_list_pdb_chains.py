#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ChainStats:
    pdb_file: str
    chain_id: str
    n_residues_atom: int
    residue_names: str


def _iter_atom_records(lines: list[str]) -> list[tuple[str, str, str, str]]:
    out: list[tuple[str, str, str, str]] = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 27:
            continue
        resname = line[17:20].strip()
        chain = line[21].strip() or "_"
        resseq = line[22:26].strip()
        icode = line[26].strip()
        resid = f"{resseq}{icode}"
        out.append((chain, resname, resid, line[:6].strip()))
    return out


def _chain_stats(pdb_path: Path) -> list[ChainStats]:
    lines = pdb_path.read_text(errors="ignore").splitlines()
    atoms = _iter_atom_records(lines)

    # Count residues based on ATOM records only (polymer-like heuristic).
    chain_to_res: dict[str, set[str]] = {}
    chain_to_names: dict[str, set[str]] = {}
    for chain, resname, resid, rectype in atoms:
        if rectype != "ATOM":
            continue
        chain_to_res.setdefault(chain, set()).add(resid)
        chain_to_names.setdefault(chain, set()).add(resname)

    stats: list[ChainStats] = []
    for chain_id in sorted(chain_to_res):
        names = sorted(chain_to_names.get(chain_id, set()))
        stats.append(
            ChainStats(
                pdb_file=pdb_path.name,
                chain_id=chain_id,
                n_residues_atom=len(chain_to_res[chain_id]),
                residue_names=";".join(names),
            )
        )
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description="List chains in PDB files (ATOM-based residue counts)")
    parser.add_argument("--pdb-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ChainStats] = []
    for pdb_path in sorted(pdb_dir.glob("*.pdb")):
        rows.extend(_chain_stats(pdb_path))

    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb_file", "chain_id", "n_residues_atom", "residue_names"])
        for r in rows:
            w.writerow([r.pdb_file, r.chain_id, r.n_residues_atom, r.residue_names])

    print(f"Wrote {len(rows)} chain rows to {out_path}")


if __name__ == "__main__":
    main()

