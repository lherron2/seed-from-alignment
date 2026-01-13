#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _extract_chain_pdb(src_pdb: Path, chain_id: str, out_pdb: Path) -> None:
    # Minimal PDB chain extraction:
    # - keep ATOM/HETATM/TER for matching chain_id (column 22)
    # - keep header-ish records as-is (MODEL/ENDMDL preserved by pass-through)
    keep_prefixes = (
        "HEADER",
        "TITLE ",
        "REMARK",
        "CRYST1",
        "ORIGX",
        "SCALE",
        "MTRIX",
        "MODEL ",
        "ENDMDL",
        "END",
    )
    kept: list[str] = []
    for line in src_pdb.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM", "TER")):
            if len(line) > 21 and (line[21].strip() or "_") == chain_id:
                kept.append(line)
            continue
        if line.startswith(keep_prefixes):
            kept.append(line)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(kept) + ("\n" if kept else ""))


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract per-target chain PDBs for ssbench truth build")
    parser.add_argument("--targets", required=True, help="CSV with columns: pdb_id,chain_id,pdb_path")
    parser.add_argument("--out-dir", required=True, help="Write {pdb_id}_{chain_id}.pdb here")
    args = parser.parse_args()

    targets_csv = Path(args.targets)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    skipped = 0
    with targets_csv.open(newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            pdb_id = str(row.get("pdb_id", "")).strip()
            chain_id = str(row.get("chain_id", "")).strip()
            pdb_path = Path(str(row.get("pdb_path", "")).strip())
            if not pdb_id or not chain_id or not pdb_path:
                skipped += 1
                continue
            if not pdb_path.exists():
                print(f"[WARN] missing pdb_path for {pdb_id}:{chain_id}: {pdb_path}")
                skipped += 1
                continue
            out_pdb = out_dir / f"{pdb_id.lower()}_{chain_id}.pdb"
            _extract_chain_pdb(pdb_path, chain_id, out_pdb)
            ok += 1

    print(f"Prepared chain PDBs: ok={ok} skipped={skipped} out_dir={out_dir}")


if __name__ == "__main__":
    main()

