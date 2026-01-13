"""Sequence normalization utilities for ssbench.

The benchmark pipeline interacts with external tools that often expect a valid
RNA alphabet. Truth sequences extracted from 3D structures can include modified
residue codes that Barnaba renders as non-ACGU characters (e.g. "X"), or can
contain separators (e.g. "&") when discontinuous fragments are present.

These helpers provide conservative normalization and a best-effort
canonicalization for Barnaba-derived sequences.
"""

from __future__ import annotations

from typing import Any


def normalize_rna_sequence(seq: str) -> str:
    """Uppercase, remove whitespace, and convert T->U."""
    seq = "".join(str(seq or "").split())
    return seq.upper().replace("T", "U")


def sanitize_rna_sequence(seq: str, *, unknown: str = "N") -> str:
    """Replace any non-ACGU character with `unknown` (default: N)."""
    if len(unknown) != 1:
        raise ValueError("unknown must be a single character")
    seq = normalize_rna_sequence(seq)
    out: list[str] = []
    for ch in seq:
        out.append(ch if ch in "ACGU" else unknown)
    return "".join(out)


def canonicalize_barnaba_sequence(seq: str, res: Any | None) -> str:
    """Best-effort canonicalization of a Barnaba dot_bracket() sequence.

    Barnaba's `dot_bracket(pairings, res)` returns a sequence string aligned to its
    internal residue list `res`, but may contain non-ACGU letters for modified
    residues. When `res` is available, we map residue names to canonical bases
    via gemmi's residue table.

    If gemmi is unavailable or mapping fails for a position, we fall back to the
    normalized sequence character (and finally to "N").
    """
    seq_norm = normalize_rna_sequence(seq)

    if res is None or not isinstance(res, (list, tuple)):
        return sanitize_rna_sequence(seq_norm, unknown="N")

    try:
        import gemmi  # type: ignore
    except ImportError:
        return sanitize_rna_sequence(seq_norm, unknown="N")

    res_idx = 0
    out: list[str] = []
    for ch in seq_norm:
        if ch == "&":
            # Barnaba may insert '&' between discontinuous fragments.
            out.append("N")
            continue

        res_name = None
        if res_idx < len(res):
            rid = str(res[res_idx])
            res_idx += 1
            res_name = rid.split("_", 1)[0].strip()

        base: str | None = None
        if res_name:
            info = gemmi.find_tabulated_residue(res_name)
            code = (info.one_letter_code or "").strip()
            if code:
                candidate = code.upper().replace("T", "U")
                if candidate in "ACGU":
                    base = candidate

        if base is None:
            base = ch if ch in "ACGU" else "N"
        out.append(base)

    return "".join(out)

