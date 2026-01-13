"""Extract a single chain from a PDB file."""

from __future__ import annotations

from pathlib import Path


def extract_chain(
    pdb_path: str | Path,
    chain_id: str,
    out_path: str | Path,
    *,
    model_index: int = 0,
) -> Path:
    """Extract a single RNA/DNA chain using gemmi and write to out_path.

    Notes:
        - We intentionally write a minimal PDB for downstream tools (Barnaba, etc).
        - We keep only nucleic-acid residues (including many modified residues), which avoids
          polluting the output with waters/ions/ligands that can confuse annotators.
        - For multi-model structures, we extract a single model (default: first model).
    """
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("gemmi is required for chain extraction") from exc

    structure = gemmi.read_structure(str(pdb_path))
    if len(structure) == 0:
        raise ValueError(f"No models found in structure: {pdb_path}")
    if model_index < 0 or model_index >= len(structure):
        raise ValueError(f"model_index {model_index} out of range for {pdb_path} (n_models={len(structure)})")

    model = structure[model_index]
    new_model = gemmi.Model(str(getattr(model, "num", model_index + 1)))

    src_chain = None
    for chain in model:
        if chain.name == chain_id:
            src_chain = chain
            break
    if src_chain is None:
        raise ValueError(f"Chain {chain_id!r} not found in {pdb_path}")

    # BGSU/FR3D chain/asym IDs can be multi-character (valid in mmCIF, invalid in PDB).
    # gemmi's PDB writer will emit malformed fixed-width PDB records if chain.name > 1,
    # which then breaks downstream parsers (Barnaba/mdtraj). Since we always extract a
    # single chain per file, we can safely normalize the output chain ID.
    chain_out = gemmi.Chain("A")

    kept = 0
    for res in src_chain:
        info = gemmi.find_tabulated_residue(res.name)
        if info.kind in (gemmi.ResidueKind.RNA, gemmi.ResidueKind.DNA):
            chain_out.add_residue(res.clone())
            kept += 1

    if kept == 0:
        raise ValueError(f"Chain {chain_id!r} in {pdb_path} contains no RNA/DNA residues")

    new_model.add_chain(chain_out)
    out = gemmi.Structure()
    out.add_model(new_model)
    out_path = Path(out_path)
    out.write_minimal_pdb(str(out_path))
    return out_path
