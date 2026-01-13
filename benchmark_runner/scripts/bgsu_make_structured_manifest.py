#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import random
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class Bucket:
    lo: int
    hi: int

    def label(self) -> str:
        return f"{self.lo}-{self.hi}"

    def contains(self, length: int) -> bool:
        return self.lo <= length <= self.hi


DEFAULT_BUCKETS = [Bucket(30, 80), Bucket(81, 150), Bucket(151, 300)]


def parse_buckets(text: str | None) -> list[Bucket]:
    if not text:
        return list(DEFAULT_BUCKETS)
    buckets: list[Bucket] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" not in part:
            raise ValueError(f"Invalid bucket {part!r}; expected 'lo-hi'")
        lo_s, hi_s = part.split("-", 1)
        lo = int(lo_s)
        hi = int(hi_s)
        if lo > hi:
            lo, hi = hi, lo
        buckets.append(Bucket(lo, hi))
    if not buckets:
        return list(DEFAULT_BUCKETS)
    return buckets


def assign_bucket(length: int, buckets: list[Bucket]) -> str:
    for b in buckets:
        if b.contains(length):
            return b.label()
    return "other"


def load_truth_metrics(truth_path: Path) -> tuple[str, int, int]:
    data = json.loads(truth_path.read_text())
    seq = str(data.get("sequence", ""))
    pairs = data.get("canonical_pairs", [])
    return seq, len(seq), len(pairs)


def _cmscan_exe(infernal_bin: str | None) -> str:
    if infernal_bin:
        return str(Path(infernal_bin) / "cmscan")
    return "cmscan"


def _cmscan_top_hit(
    *,
    cmscan: str,
    rfam_cm: Path,
    fasta_path: Path,
    out_dir: Path,
) -> tuple[str, str, str] | None:
    tblout = out_dir / "cmscan.tblout"
    cmd = [cmscan, "--tblout", str(tblout), str(rfam_cm), str(fasta_path)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return None
    for line in tblout.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        target_name = parts[0].strip()
        accession = parts[1].strip()
        # tblout columns end with "description of target" (may contain spaces).
        desc = " ".join(parts[18:]).strip() if len(parts) > 18 else ""
        return target_name, accession, desc
    return None


def is_ribosomal_rna_hit(target_name: str, accession: str, description: str) -> bool:
    # Heuristic: rRNA models reliably contain "rRNA" in the model name (e.g., 5S_rRNA, SSU_rRNA_*).
    # Avoid over-filtering families whose description mentions "ribosomal" but are not rRNAs.
    _ = accession  # reserved for future allow/deny lists
    # Avoid false positives like "xrRNA" (exoribonuclease-resistant RNAs), which contain "rRNA"
    # as a substring but are not ribosomal RNAs.
    import re

    if re.search(r"(^|_)rRNA(_|$)", target_name):
        return True
    if "ribosomal RNA" in description:
        return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter a BGSU/FR3D representative-set manifest using truth-derived structure metrics.",
    )
    parser.add_argument("--manifest", required=True, help="Input manifest CSV (must include target_id).")
    parser.add_argument("--truth-dir", required=True, help="Directory of truth JSONs named {target_id}.json.")
    parser.add_argument("--out-manifest", required=True, help="Output manifest CSV.")
    parser.add_argument("--min-len", type=int, default=30)
    parser.add_argument("--max-len", type=int, default=300)
    parser.add_argument("--min-canonical-pairs", type=int, default=5)
    parser.add_argument(
        "--min-pair-density",
        type=float,
        default=0.0,
        help="Optional min (n_pairs / length) threshold.",
    )
    parser.add_argument(
        "--max-targets",
        type=int,
        default=None,
        help="If set, downsample to this many targets (roughly balanced across buckets).",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--buckets",
        default=None,
        help="Comma-separated length buckets like '30-80,81-150,151-300'.",
    )
    parser.add_argument(
        "--dedupe-ec",
        action="store_true",
        help="If ec_id exists, keep at most one target per ec_id (preferring higher n_pairs then length).",
    )
    parser.add_argument(
        "--exclude-ribosomal",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude targets whose best cmscan hit appears to be a ribosomal RNA (rRNA).",
    )
    parser.add_argument(
        "--rfam-cm",
        default=None,
        help="Path to Rfam.cm (required when --exclude-ribosomal).",
    )
    parser.add_argument(
        "--infernal-bin",
        default=None,
        help="Optional path to Infernal bin dir containing cmscan.",
    )
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    truth_dir = Path(args.truth_dir)
    out_path = Path(args.out_manifest)

    df = pd.read_csv(manifest_path)
    if "target_id" not in df.columns:
        raise SystemExit("manifest is missing required column: target_id")

    rfam_cm = None
    cmscan = None
    if bool(args.exclude_ribosomal):
        rfam_cm_val = args.rfam_cm or str(Path(__file__).resolve().parents[1] / "data" / "rfam" / "Rfam.cm")
        rfam_cm = Path(rfam_cm_val)
        if not rfam_cm.exists():
            raise SystemExit(f"--exclude-ribosomal enabled but Rfam.cm not found: {rfam_cm}")
        cmscan = _cmscan_exe(args.infernal_bin)
        if args.infernal_bin and not Path(cmscan).exists():
            raise SystemExit(f"--infernal-bin provided but cmscan not found: {cmscan}")

    rows = []
    missing = 0
    excluded_rrna = 0
    with tempfile.TemporaryDirectory(prefix="cmscan_rrna_") as tmp:
        tmpdir = Path(tmp)
        fasta_path = tmpdir / "q.fa"
        for _, row in df.iterrows():
            target_id = str(row["target_id"])
            truth_path = truth_dir / f"{target_id}.json"
            if not truth_path.exists():
                missing += 1
                continue
            seq, length, n_pairs = load_truth_metrics(truth_path)
            length = int(length)
            n_pairs = int(n_pairs)

            # Apply structure/length filters early to avoid running cmscan on targets
            # that will be filtered out anyway (cmscan can be expensive).
            if length < int(args.min_len) or length > int(args.max_len) or n_pairs < int(args.min_canonical_pairs):
                continue
            if float(args.min_pair_density) > 0.0:
                density = float(n_pairs) / float(max(1, length))
                if density < float(args.min_pair_density):
                    continue

            if bool(args.exclude_ribosomal):
                seq_u = seq.strip().upper().replace("T", "U")
                # cmscan expects a valid nucleotide alphabet; tolerate modified/unknown residues.
                seq_u = "".join(ch if ch in "ACGU" else "N" for ch in seq_u)
                if seq_u:
                    fasta_path.write_text(f">{target_id}\n{seq_u}\n")
                    hit = _cmscan_top_hit(
                        cmscan=str(cmscan),
                        rfam_cm=Path(rfam_cm),
                        fasta_path=fasta_path,
                        out_dir=tmpdir,
                    )
                    if hit is not None:
                        name, acc, desc = hit
                        if is_ribosomal_rna_hit(name, acc, desc):
                            excluded_rrna += 1
                            continue
            rows.append(
                {
                    **row.to_dict(),
                    "truth_len": length,
                    "truth_canonical_pairs": n_pairs,
                }
            )

    if not rows:
        raise SystemExit("No rows had matching truth JSONs. Check --truth-dir and target_id naming.")
    if missing:
        print(f"[WARN] Missing truth for {missing} targets (skipped).")
    if excluded_rrna:
        print(f"[INFO] Excluded {excluded_rrna} ribosomal RNA targets.")

    df2 = pd.DataFrame(rows)
    df2 = df2.copy()
    if float(args.min_pair_density) > 0:
        df2["truth_pair_density"] = df2["truth_canonical_pairs"] / df2["truth_len"].clip(lower=1)

    if args.dedupe_ec and "ec_id" in df2.columns:
        df2 = df2.sort_values(
            by=["truth_canonical_pairs", "truth_len", "target_id"],
            ascending=[False, False, True],
        )
        df2 = df2.drop_duplicates(subset=["ec_id"], keep="first")

    if df2.empty:
        raise SystemExit("No targets remain after filtering; try lowering thresholds.")

    buckets = parse_buckets(args.buckets)
    df2["bucket_truth_len"] = [assign_bucket(int(x), buckets) for x in df2["truth_len"].tolist()]

    if args.max_targets is not None and len(df2) > int(args.max_targets):
        rng = random.Random(int(args.seed))
        max_targets = int(args.max_targets)
        # Allocate roughly evenly across declared buckets; fill leftovers from remaining pool.
        bucket_labels = [b.label() for b in buckets]
        per_bucket = max_targets // max(1, len(bucket_labels))
        selected_rows = []
        remaining = []
        for label in bucket_labels + ["other"]:
            sub = df2[df2["bucket_truth_len"] == label]
            if sub.empty:
                continue
            idx = list(sub.index)
            rng.shuffle(idx)
            take = per_bucket if label != "other" else max_targets
            chosen = idx[:take]
            selected_rows.extend(chosen)
        selected_set = set(selected_rows)
        remaining_idx = [i for i in df2.index if i not in selected_set]
        rng.shuffle(remaining_idx)
        selected_rows = selected_rows[:max_targets]
        while len(selected_rows) < max_targets and remaining_idx:
            selected_rows.append(remaining_idx.pop())
        df2 = df2.loc[selected_rows].copy()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df2.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(df2)} rows)")


if __name__ == "__main__":
    main()
