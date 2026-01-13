"""FASTA parsing and conservative RNA sanitization utilities for predictors."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FastaRecord:
    seq_id: str
    sequence: str


def read_single_fasta(text: str) -> FastaRecord:
    """Parse a single-sequence FASTA (or raw sequence) into (id, sequence)."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return FastaRecord(seq_id="sequence", sequence="")
    if lines[0].startswith(">"):
        seq_id = lines[0][1:].split()[0].strip() or "sequence"
        seq = "".join(lines[1:]).strip()
        return FastaRecord(seq_id=seq_id, sequence=seq)
    return FastaRecord(seq_id="sequence", sequence="".join(lines).strip())


def normalize_rna_sequence(seq: str) -> str:
    """Uppercase, remove whitespace, convert T->U."""
    seq = "".join(str(seq or "").split())
    return seq.upper().replace("T", "U")


def sanitize_rna_sequence_preserve_length(seq: str) -> str:
    """Return a tool-safe sequence while preserving length.

    Rules:
    - T->U
    - '&' -> 'N' (Barnaba fragment separator; rejected by many tools)
    - Any non-alphanumeric -> 'N'
    - Keep A/C/G/U and common ambiguity codes (X/N); other letters -> 'N'
    """
    seq = normalize_rna_sequence(seq)
    out: list[str] = []
    for ch in seq:
        if ch == "&":
            out.append("N")
            continue
        if not ch.isalnum():
            out.append("N")
            continue
        if ch in {"A", "C", "G", "U", "N", "X"}:
            out.append(ch)
            continue
        out.append("N")
    return "".join(out)

