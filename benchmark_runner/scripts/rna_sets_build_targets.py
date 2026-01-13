#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TargetRow:
    target_id: str
    pdb_id: str
    chain_id: str
    pdb_path: str
    nts_observed: int
    split: str
    bucket: str


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


def _chain_counts(pdb_path: Path) -> dict[str, int]:
    lines = pdb_path.read_text(errors="ignore").splitlines()
    atoms = _iter_atom_records(lines)
    chain_to_res: dict[str, set[str]] = {}
    for chain, _resname, resid, rectype in atoms:
        if rectype != "ATOM":
            continue
        chain_to_res.setdefault(chain, set()).add(resid)
    return {chain: len(resids) for chain, resids in chain_to_res.items()}


def _pick_reference_pdb(ref_dir: Path) -> Path | None:
    for name in ("native.pdb", "ref.pdb"):
        cand = ref_dir / name
        if cand.exists():
            return cand
    pdbs = sorted(ref_dir.glob("*.pdb"))
    return pdbs[0] if pdbs else None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build CASP/ribo targets from rna_sets (single-chain, <=max-nt)."
    )
    parser.add_argument("--rna-sets-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-nt", type=int, default=300)
    args = parser.parse_args()

    rna_sets = Path(args.rna_sets_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[TargetRow] = []
    skipped = 0

    for subset_dir in sorted(p for p in rna_sets.iterdir() if p.is_dir()):
        subset = subset_dir.name
        for target_dir in sorted(p for p in subset_dir.iterdir() if p.is_dir()):
            ref_dir = target_dir / "structures" / "references"
            if not ref_dir.exists():
                skipped += 1
                continue
            pdb_path = _pick_reference_pdb(ref_dir)
            if pdb_path is None:
                skipped += 1
                continue

            chain_counts = _chain_counts(pdb_path)
            if not chain_counts:
                skipped += 1
                continue
            if len(chain_counts) != 1:
                skipped += 1
                continue

            chain_id, nts = next(iter(chain_counts.items()))
            if nts > int(args.max_nt):
                skipped += 1
                continue

            target_id = f"{subset}_{target_dir.name}"
            pdb_id = target_dir.name
            rows.append(
                TargetRow(
                    target_id=target_id,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    pdb_path=str(pdb_path),
                    nts_observed=nts,
                    split=subset,
                    bucket=f"len_le_{args.max_nt}",
                )
            )

    if not rows:
        raise SystemExit("No single-chain targets found. Check rna_sets contents.")

    with out_path.open("w", newline="") as f:
        fieldnames = [
            "target_id",
            "pdb_id",
            "chain_id",
            "pdb_path",
            "nts_observed",
            "split",
            "bucket",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row.__dict__)

    print(f"Wrote {len(rows)} targets to {out_path} (skipped={skipped})")


if __name__ == "__main__":
    main()

