"""
Structure validation module for RNA secondary structures.

This module is the SINGLE SOURCE OF TRUTH for RNA structure validation.
Use it everywhere: unit tests, integration tests, pipeline post-checks,
CLI --validate mode.

See TASK_TRACKER.md Structure Validation Module Spec for detailed requirements.
"""

from __future__ import annotations

__all__ = [
    "validate_structure",
    "parse_structure_to_pairs",
    "pairs_to_structure",
]

# Bracket types for WUSS notation
OPEN_TO_CLOSE = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
}
# Add a-z -> A-Z
for i in range(26):
    OPEN_TO_CLOSE[chr(ord("a") + i)] = chr(ord("A") + i)

CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}

# Canonical base pairs
CANONICAL_PAIRS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),  # wobble
}


def parse_structure_to_pairs(struct: str) -> list[tuple[int, int]]:
    """Parse a dot-bracket structure into base pairs.

    Args:
        struct: Dot-bracket string (may include pseudoknot brackets)

    Returns:
        List of (i, j) pairs with i < j
    """
    stacks = {op: [] for op in OPEN_TO_CLOSE}
    pairs = []

    for i, ch in enumerate(struct):
        if ch in OPEN_TO_CLOSE:
            stacks[ch].append(i)
        elif ch in CLOSE_TO_OPEN:
            open_ch = CLOSE_TO_OPEN[ch]
            if stacks[open_ch]:
                j = stacks[open_ch].pop()
                pairs.append((j, i) if j < i else (i, j))

    return sorted(pairs)


def pairs_to_structure(pairs: list[tuple[int, int]], length: int) -> str:
    """Convert pairs to a dot-bracket structure.

    Uses standard bracket hierarchy: () [] {} <> aA bB ...

    Args:
        pairs: List of (i, j) pairs
        length: Sequence length

    Returns:
        Dot-bracket string
    """
    from .pk import pairs_to_layers

    chars = ["."] * length
    layers = pairs_to_layers(pairs)

    bracket_types = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
    for i in range(26):
        bracket_types.append((chr(ord("a") + i), chr(ord("A") + i)))

    for layer_idx, layer in enumerate(layers):
        if layer_idx < len(bracket_types):
            open_ch, close_ch = bracket_types[layer_idx]
        else:
            open_ch, close_ch = "{", "}"

        for i, j in layer:
            chars[i] = open_ch
            chars[j] = close_ch

    return "".join(chars)


def validate_structure(
    seq: str,
    struct: str,
    *,
    allow_pk: bool = True,
    min_hairpin: int = 3,
    allow_lonely_pairs: bool = False,
    max_pk_depth: int | None = None,
    canonical_only: bool = False,
    lowercase_forced_unpaired: bool = False,
) -> tuple[bool, list[str]]:
    """Validate an RNA structure against all invariants.

    This is the single source of truth for structure validation.

    Args:
        seq: RNA sequence
        struct: Dot-bracket structure
        allow_pk: Whether pseudoknots are allowed
        min_hairpin: Minimum unpaired bases in hairpin loops
        allow_lonely_pairs: Whether isolated pairs are allowed
        max_pk_depth: Maximum PK layer depth (None = unlimited)
        canonical_only: Only allow canonical base pairs (AU, GC, GU)
        lowercase_forced_unpaired: Lowercase bases must be unpaired

    Returns:
        Tuple of (is_valid, list_of_error_messages)
    """
    errors = []

    # 1. Length match
    if len(seq) != len(struct):
        errors.append(f"Length mismatch: seq={len(seq)}, struct={len(struct)}")
        return False, errors

    # 2. Parse structure and check balance
    pairs = []
    stacks = {op: [] for op in OPEN_TO_CLOSE}
    paired_positions = set()

    for i, ch in enumerate(struct):
        if ch in OPEN_TO_CLOSE:
            stacks[ch].append(i)
        elif ch in CLOSE_TO_OPEN:
            open_ch = CLOSE_TO_OPEN[ch]
            if not stacks[open_ch]:
                errors.append(f"Unmatched close bracket '{ch}' at position {i}")
            else:
                j = stacks[open_ch].pop()
                pair = (j, i) if j < i else (i, j)
                pairs.append(pair)

                # Check for double-pairing
                if j in paired_positions:
                    errors.append(f"Position {j} paired twice")
                if i in paired_positions:
                    errors.append(f"Position {i} paired twice")
                paired_positions.add(i)
                paired_positions.add(j)
        elif ch != ".":
            errors.append(f"Unknown character '{ch}' at position {i}")

    # Check for unmatched open brackets
    for open_ch, stack in stacks.items():
        for pos in stack:
            errors.append(f"Unmatched open bracket '{open_ch}' at position {pos}")

    if errors:
        return False, errors

    # 3. Check pairs are in range and ordered
    for i, j in pairs:
        if i < 0 or j >= len(seq):
            errors.append(f"Pair ({i}, {j}) out of range")
        if i >= j:
            errors.append(f"Pair ({i}, {j}) not properly ordered (i >= j)")

    # 4. Canonical pairing check
    if canonical_only:
        for i, j in pairs:
            b1 = seq[i].upper()
            b2 = seq[j].upper()
            if (b1, b2) not in CANONICAL_PAIRS:
                errors.append(f"Non-canonical pair ({b1}, {b2}) at ({i}, {j})")

    # 5. Hairpin length check
    if min_hairpin > 0:
        for i, j in pairs:
            loop_len = j - i - 1
            # Check if any positions in between are unpaired
            # This is a simple hairpin check - more complex loops handled elsewhere
            all_unpaired = all(k not in paired_positions for k in range(i + 1, j))
            if all_unpaired and loop_len < min_hairpin:
                errors.append(f"Hairpin loop at ({i}, {j}) too small: {loop_len} < {min_hairpin}")

    # 6. Lonely pair check
    if not allow_lonely_pairs:
        pair_set = set(pairs)
        for i, j in pairs:
            # A pair is lonely if neither adjacent stacking pair exists
            has_stack_inside = (i + 1, j - 1) in pair_set
            has_stack_outside = (i - 1, j + 1) in pair_set
            if not has_stack_inside and not has_stack_outside:
                errors.append(f"Lonely pair at ({i}, {j})")

    # 7. Pseudoknot checks
    if not allow_pk or max_pk_depth is not None:
        from .pk import crosses, pairs_to_layers

        # Check for crossings
        has_pk = False
        for idx1, p1 in enumerate(pairs):
            for p2 in pairs[idx1 + 1 :]:
                if crosses(p1, p2):
                    has_pk = True
                    if not allow_pk:
                        errors.append(f"Pseudoknot between ({p1}) and ({p2})")
                    break
            if has_pk and not allow_pk:
                break

        # Check depth limit
        if max_pk_depth is not None and has_pk:
            layers = pairs_to_layers(pairs)
            if len(layers) > max_pk_depth + 1:  # +1 because layer 0 is non-PK
                errors.append(f"PK depth {len(layers) - 1} exceeds limit {max_pk_depth}")

    # 8. Lowercase forced unpaired check
    if lowercase_forced_unpaired:
        for i, ch in enumerate(seq):
            if ch.islower() and i in paired_positions:
                errors.append(f"Lowercase base '{ch}' at position {i} is paired")

    return len(errors) == 0, errors
