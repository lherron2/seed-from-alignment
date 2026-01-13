"""Dot-bracket parsing utilities."""

from __future__ import annotations


def parse_dotbracket(text: str) -> tuple[str, list[str]]:
    """Parse dot-bracket content supporting 2-line or FASTA-like formats.

    Returns:
        Tuple of (sequence, list of dot-bracket strings).
    """
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return "", []

    if lines[0].startswith(">"):
        if len(lines) < 3:
            return "", []
        seq = lines[1]
        dotbrs = lines[2:]
        return seq, dotbrs

    if len(lines) >= 2:
        seq = lines[0]
        dotbrs = lines[1:]
        return seq, dotbrs

    return "", []


def pairs_from_dotbracket(struct: str) -> list[tuple[int, int]]:
    """Parse pairs from dot-bracket with multiple bracket alphabets."""
    open_to_close = {"(": ")", "[": "]", "{": "}", "<": ">"}
    for i in range(26):
        open_to_close[chr(ord("a") + i)] = chr(ord("A") + i)
    close_to_open = {v: k for k, v in open_to_close.items()}

    stacks = {op: [] for op in open_to_close}
    pairs: list[tuple[int, int]] = []
    for idx, ch in enumerate(struct):
        if ch in open_to_close:
            stacks[ch].append(idx)
        elif ch in close_to_open:
            op = close_to_open[ch]
            if stacks[op]:
                i = stacks[op].pop()
                pairs.append((i, idx))
    return pairs
