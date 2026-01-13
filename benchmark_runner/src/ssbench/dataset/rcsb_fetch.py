"""Fetch PDB structures from RCSB."""

from __future__ import annotations

from pathlib import Path

import requests

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb}.cif"


def fetch_pdb(pdb_id: str, out_path: str | Path) -> Path:
    """Download a PDB or mmCIF file and return the path."""
    pdb_id_norm = pdb_id.upper()
    out_path = Path(out_path)
    url = RCSB_PDB_URL.format(pdb=pdb_id_norm)
    resp = requests.get(url, timeout=60)
    if resp.status_code == 404:
        url = RCSB_CIF_URL.format(pdb=pdb_id_norm)
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        out_path = out_path.with_suffix(".cif")
    else:
        resp.raise_for_status()
    out_path.write_text(resp.text)
    return out_path
