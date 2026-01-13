#!/usr/bin/env python3
"""Build a manifest from a FASTA file by matching sequences to PDB chains."""

from __future__ import annotations

import argparse
import csv
import difflib
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Entry:
    header: str
    seq: str


def read_fasta(path: Path) -> list[Entry]:
    entries: list[Entry] = []
    header = None
    seq_chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append(Entry(header=header, seq="".join(seq_chunks)))
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if header is not None:
        entries.append(Entry(header=header, seq="".join(seq_chunks)))
    return entries


def residue_to_base(residue, gemmi) -> str:
    try:
        tab = gemmi.find_tabulated_residue(residue.name)
        if not tab.is_nucleic_acid():
            return ""
        base = tab.one_letter_code.upper()
    except Exception:
        base = "N"
    if base == "T":
        base = "U"
    if base not in {"A", "C", "G", "U"}:
        base = "N"
    return base


def chain_sequence(chain, gemmi) -> str:
    seq = []
    for res in chain:
        base = residue_to_base(res, gemmi)
        if not base:
            continue
        seq.append(base)
    return "".join(seq)


def best_chain_match(structure, target_seq: str, min_ratio: float, gemmi):
    best = None
    best_ratio = 0.0
    for model in structure:
        for chain in model:
            seq = chain_sequence(chain, gemmi)
            if not seq:
                continue
            ratio = difflib.SequenceMatcher(a=target_seq, b=seq).ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best = (chain.name, seq, ratio)
    if best is None or best_ratio < min_ratio:
        return None
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description="Build manifest from FASTA and PDBs.")
    parser.add_argument("--fasta", required=True, help="Input FASTA")
    parser.add_argument("--out-manifest", required=True, help="Output manifest CSV")
    parser.add_argument("--pdb-dir", required=True, help="Directory for PDB/mmCIF files")
    parser.add_argument("--download", action="store_true", help="Download PDBs if missing")
    parser.add_argument("--min-identity", type=float, default=0.9, help="Min sequence identity ratio")
    parser.add_argument("--allow-unmatched", action="store_true", help="Skip entries with no chain match")
    parser.add_argument("--unmatched-out", default=None, help="Write unmatched headers to this file")
    args = parser.parse_args()

    try:
        import gemmi
    except ImportError as exc:
        raise SystemExit("gemmi is required for chain matching") from exc

    pdb_dir = Path(args.pdb_dir)
    pdb_dir.mkdir(parents=True, exist_ok=True)
    entries = read_fasta(Path(args.fasta))

    if args.download:
        from ssbench.dataset.rcsb_fetch import fetch_pdb

    def find_local_structure(pdb_id: str) -> Path | None:
        pdb_id_norm = pdb_id.strip()
        if not pdb_id_norm:
            return None
        candidates = []
        for stem in {pdb_id_norm.upper(), pdb_id_norm.lower()}:
            candidates.extend([pdb_dir / f"{stem}.pdb", pdb_dir / f"{stem}.cif"])
        for path in candidates:
            if path.exists():
                return path
        for path in pdb_dir.iterdir():
            if not path.is_file():
                continue
            if path.suffix.lower() not in {".pdb", ".cif"}:
                continue
            if path.stem.lower() == pdb_id_norm.lower():
                return path
        return None

    rows = []
    unmatched = []
    for entry in entries:
        parts = entry.header.split("|")
        if len(parts) < 2:
            raise SystemExit(f"Header missing PDB id: {entry.header}")
        target_id = parts[0]
        pdb_id = parts[1].upper()
        pdb_path = find_local_structure(pdb_id)
        if pdb_path is None and args.download:
            pdb_path = fetch_pdb(pdb_id, pdb_dir / f"{pdb_id}.pdb")
        if pdb_path is None or not pdb_path.exists():
            raise SystemExit(f"Missing structure file for {pdb_id}. Use --download.")

        structure = gemmi.read_structure(str(pdb_path))
        match = best_chain_match(structure, entry.seq, args.min_identity, gemmi)
        if match is None:
            msg = f"No chain match for {target_id} ({pdb_id})."
            if not args.allow_unmatched:
                raise SystemExit(msg)
            unmatched.append(entry.header)
            continue
        chain_id, _seq, ratio = match

        rows.append(
            {
                "target_id": target_id,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "ife_id": "",
                "ec_id": "",
                "rfam": "",
                "nts_observed": len(entry.seq),
                "pdb_resolution": "",
                "split": "train",
                "source_header": entry.header,
                "match_ratio": f"{ratio:.3f}",
            }
        )

    out_path = Path(args.out_manifest)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    if args.unmatched_out and unmatched:
        Path(args.unmatched_out).write_text("\n".join(unmatched) + "\n")


if __name__ == "__main__":
    main()
