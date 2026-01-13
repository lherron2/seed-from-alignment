"""Extract canonical base pairs from 3D structures using Barnaba."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..sequence import canonicalize_barnaba_sequence
from .helices import identify_helices
from .pk_detect import split_nested_pk

CANONICAL_LABELS = {"WCc", "WCt", "GUc", "GUt"}


@dataclass
class TruthRecord:
    """Truth structure extracted from 3D."""

    sequence: str
    canonical_pairs: list[tuple[int, int]]
    nested_pairs: list[tuple[int, int]]
    pk_pairs: list[tuple[int, int]]
    helices: list[list[tuple[int, int]]]


def _parse_pairings(pairings: Any) -> list[tuple[int, int, str]]:
    parsed: list[tuple[int, int, str]] = []
    for entry in pairings:
        if isinstance(entry, dict):
            i = int(entry.get("nt1", entry.get("i", -1)))
            j = int(entry.get("nt2", entry.get("j", -1)))
            label = str(entry.get("interaction", entry.get("type", "")))
            if i >= 0 and j >= 0:
                if i > j:
                    i, j = j, i
                parsed.append((i, j, label))
            continue

        if isinstance(entry, (list, tuple)) and len(entry) == 2 and isinstance(entry[0], list):
            pairs_list = entry[0]
            labels_list = entry[1] if isinstance(entry[1], list) else []
            for idx, pair in enumerate(pairs_list):
                if len(pair) < 2:
                    continue
                i = int(pair[0])
                j = int(pair[1])
                label = str(labels_list[idx]) if idx < len(labels_list) else ""
                if i > j:
                    i, j = j, i
                parsed.append((i, j, label))
            continue

        if isinstance(entry, (list, tuple)) and len(entry) >= 2:
            i = int(entry[0])
            j = int(entry[1])
            label = str(entry[2]) if len(entry) > 2 else ""
            if i > j:
                i, j = j, i
            parsed.append((i, j, label))
    return parsed


def build_truth(pdb_path: str | Path) -> TruthRecord:
    """Run Barnaba to extract canonical base pairs and derive truth record."""
    try:
        import barnaba as bb
    except ImportError as exc:
        raise RuntimeError("barnaba is required to build truth") from exc

    anno = bb.annotate(str(pdb_path))
    pairings: Any
    res = None
    if isinstance(anno, dict):
        pairings = anno.get("pairings", [])
        res = anno.get("res", None)
    elif isinstance(anno, (list, tuple)) and len(anno) >= 2:
        # Barnaba returns (stackings, pairings, res)
        pairings = anno[1]
        if len(anno) >= 3:
            res = anno[2]
    else:
        pairings = []

    parsed = _parse_pairings(pairings)

    canonical_pairs = [(i, j) for i, j, label in parsed if label in CANONICAL_LABELS]
    canonical_pairs = sorted(set(canonical_pairs))

    # Use Barnaba dot_bracket for sequence indexing
    dotbr, seq = bb.dot_bracket(pairings, res)
    if isinstance(seq, list):
        seq = "".join(seq)
    seq = canonicalize_barnaba_sequence(str(seq), res)

    nested_pairs, pk_pairs = split_nested_pk(canonical_pairs)
    helices = identify_helices(canonical_pairs)

    return TruthRecord(
        sequence=seq,
        canonical_pairs=canonical_pairs,
        nested_pairs=nested_pairs,
        pk_pairs=pk_pairs,
        helices=helices,
    )
