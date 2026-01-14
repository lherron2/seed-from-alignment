# ribo-benchmark oracle top-K evaluation (F1 + MCC)

This documents how to reproduce the **per-RNA oracle best-of-K** results on
`rna_sets/ribo-benchmark` for `K ∈ {1,50,100,200,500}` and compare:

- RNAnneal-ss top-K pipeline output (`predictions.db`)
- EternaFold samples (`eternafold/ef.ss`)
- LinearFold samples (`linearfold-v/stochastic.ss`)

Oracle best-of-K means: for each `K`, score every structure in the first `K` outputs and take the
best (highest) F1/MCC.

## Install (fresh clone)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

If RNAstructure complains about missing thermodynamic tables:

```bash
export DATAPATH="$(pwd)/RNAstructure/data_tables"
```

## Run RNAnneal-ss on ribo-benchmark (top-k 500)

From repo root:

```bash
out_root="out/ribo_benchmark_rnanneal_ss_top500_opt6"
mkdir -p "$out_root"
for d in rna_sets/ribo-benchmark/*; do
  [ -d "$d" ] || continue
  r=$(basename "$d")
  rnanneal-ss "$d/$r.fasta" "$out_root/$r" --top-k 500
done
```

Outputs per RNA:

- `out/ribo_benchmark_rnanneal_ss_top500_opt6/<RNA>/predictions.db` (sequence line + dot-brackets)
- `out/ribo_benchmark_rnanneal_ss_top500_opt6/<RNA>/qc.json` (selection diagnostics)

Optional speedup for re-runs (reuse cached CaCoFold outputs from a previous run):

```bash
rnanneal-ss rna_sets/ribo-benchmark/6WJR/6WJR.fasta out/6WJR --top-k 500 \
  --reuse-cacofold-root out/ribo_benchmark_rnanneal_ss_top500_opt6
```

## (Optional) Run the previous implementation for comparison

If you want a side-by-side comparison vs the pre-optimization pipeline, generate a second set of
predictions at commit `21a5569`:

```bash
git checkout 21a5569
out_root="out/ribo_benchmark_rnanneal_ss_top500"
mkdir -p "$out_root"
for d in rna_sets/ribo-benchmark/*; do
  [ -d "$d" ] || continue
  r=$(basename "$d")
  rnanneal-ss "$d/$r.fasta" "$out_root/$r" --top-k 500
done
git checkout main
```

## Score oracle best-of-K (per RNA)

This writes a per-RNA CSV containing F1/MCC for:

- `rn_prev_*`: `out/ribo_benchmark_rnanneal_ss_top500/<RNA>/predictions.db` (optional)
- `rn_opt6_*`: `out/ribo_benchmark_rnanneal_ss_top500_opt6/<RNA>/predictions.db`
- `ef_*`: `rna_sets/ribo-benchmark/<RNA>/eternafold/ef.ss`
- `lf_*`: `rna_sets/ribo-benchmark/<RNA>/linearfold-v/stochastic.ss`

```bash
python3 - <<'PY'
from __future__ import annotations

import csv
import statistics
import string
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path("benchmark_runner/src").resolve()))
sys.path.insert(0, str(Path(".").resolve()))

from src.lib.validate_structure import parse_structure_to_pairs  # noqa: E402
from ssbench.metrics.pair_metrics import compute_pair_metrics  # noqa: E402

KS = [1, 50, 100, 200, 500]
ALLOWED = set(".()[]{}<>") | set(string.ascii_letters)
SEQ_CHARS = set("ACGUTNXacgutnx")


def read_fasta_seq(path: Path) -> str:
    seq = "".join(
        ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith(">")
    )
    return seq.upper().replace("T", "U")


def read_truth_db(path: Path) -> str:
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if s:
            return s.split()[0]
    raise ValueError(f"No truth structure found in {path}")


def read_struct_lines(path: Path, *, length: int, has_sequence_header: bool = False) -> list[str]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        return []

    if has_sequence_header:
        first = lines[0].split()[0]
        if len(first) == length and set(first) <= SEQ_CHARS:
            lines = lines[1:]

    out: list[str] = []
    for raw in lines:
        tok = raw.split()[0]
        if len(tok) != length:
            continue
        if set(tok) <= ALLOWED:
            out.append(tok)
    return out


@dataclass(frozen=True)
class OracleAtK:
    f1: float
    mcc: float
    eff_k: int


def oracle_best_of_k(structs: list[str], *, ref_pairs: set[tuple[int, int]], length: int, k: int) -> OracleAtK:
    eff_k = min(int(k), len(structs))
    best_f1 = -1.0
    best_mcc = 0.0
    for s in structs[:eff_k]:
        pred_pairs = set(parse_structure_to_pairs(s))
        pm = compute_pair_metrics(pred_pairs, ref_pairs, length)
        if pm.f1 > best_f1:
            best_f1 = float(pm.f1)
            best_mcc = float(pm.mcc)
    if best_f1 < 0:
        best_f1 = 0.0
        best_mcc = 0.0
    return OracleAtK(f1=best_f1, mcc=best_mcc, eff_k=eff_k)


root = Path("rna_sets/ribo-benchmark")
rnas = sorted(p.name for p in root.iterdir() if p.is_dir())

rn_prev_root = Path("out/ribo_benchmark_rnanneal_ss_top500")
rn_opt_root = Path("out/ribo_benchmark_rnanneal_ss_top500_opt6")
out_path = Path("out/ribo_benchmark_rnanneal_opt6_metrics.csv")

rows: list[dict[str, object]] = []
for rna in rnas:
    rna_dir = root / rna
    seq = read_fasta_seq(rna_dir / f"{rna}.fasta")
    L = len(seq)
    truth = read_truth_db(rna_dir / "structures" / "references" / "ref.ss")
    if len(truth) != L:
        raise ValueError(f"{rna}: truth length {len(truth)} != seq length {L}")

    ref_pairs = set(parse_structure_to_pairs(truth))

    rn_opt_structs = read_struct_lines(
        rn_opt_root / rna / "predictions.db", length=L, has_sequence_header=True
    )
    rn_prev_structs = (
        read_struct_lines(rn_prev_root / rna / "predictions.db", length=L, has_sequence_header=True)
        if (rn_prev_root / rna / "predictions.db").exists()
        else []
    )
    ef_structs = read_struct_lines(rna_dir / "eternafold" / "ef.ss", length=L)
    lf_structs = read_struct_lines(rna_dir / "linearfold-v" / "stochastic.ss", length=L)

    row: dict[str, object] = {
        "rna": rna,
        "len": L,
        "rn_prev_n": len(rn_prev_structs),
        "rn_opt6_n": len(rn_opt_structs),
        "ef_n": len(ef_structs),
        "lf_n": len(lf_structs),
    }

    for k in KS:
        rn_opt = oracle_best_of_k(rn_opt_structs, ref_pairs=ref_pairs, length=L, k=k)
        ef = oracle_best_of_k(ef_structs, ref_pairs=ref_pairs, length=L, k=k)
        lf = oracle_best_of_k(lf_structs, ref_pairs=ref_pairs, length=L, k=k)
        row |= {
            f"rn_opt6_f1@{k}": rn_opt.f1,
            f"rn_opt6_mcc@{k}": rn_opt.mcc,
            f"rn_opt6_eff@{k}": rn_opt.eff_k,
            f"ef_f1@{k}": ef.f1,
            f"ef_mcc@{k}": ef.mcc,
            f"ef_eff@{k}": ef.eff_k,
            f"lf_f1@{k}": lf.f1,
            f"lf_mcc@{k}": lf.mcc,
            f"lf_eff@{k}": lf.eff_k,
        }
        if rn_prev_structs:
            rn_prev = oracle_best_of_k(rn_prev_structs, ref_pairs=ref_pairs, length=L, k=k)
            row |= {
                f"rn_prev_f1@{k}": rn_prev.f1,
                f"rn_prev_mcc@{k}": rn_prev.mcc,
                f"rn_prev_eff@{k}": rn_prev.eff_k,
            }

    rows.append(row)

out_path.parent.mkdir(parents=True, exist_ok=True)
keys: list[str] = []
for r in rows:
    for k in r.keys():
        if k not in keys:
            keys.append(k)

with out_path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)

print(f"Wrote {out_path} ({len(rows)} RNAs)")

def mean(col: str) -> float:
    vals = [float(r[col]) for r in rows if col in r and r[col] not in ("", None)]
    return statistics.mean(vals) if vals else float("nan")

for k in KS:
    print(
        f"@{k} mean F1/MCC: rn_opt6 {mean(f'rn_opt6_f1@{k}'):.3f}/{mean(f'rn_opt6_mcc@{k}'):.3f}  "
        f"ef {mean(f'ef_f1@{k}'):.3f}/{mean(f'ef_mcc@{k}'):.3f}  "
        f"lf {mean(f'lf_f1@{k}'):.3f}/{mean(f'lf_mcc@{k}'):.3f}"
    )
PY
```

## Reference results (this repo, commit `f95eb98`)

Mean oracle F1/MCC over `rna_sets/ribo-benchmark` (N=16) from `out/ribo_benchmark_rnanneal_opt6_metrics.csv`:

- `@1` RN `0.729/0.733` | EF `0.598/0.596` | LF `0.747/0.750`
- `@50` RN `0.844/0.846` | EF `0.838/0.839` | LF `0.838/0.840`
- `@100` RN `0.848/0.851` | EF `0.848/0.849` | LF `0.847/0.848`
- `@200` RN `0.856/0.858` | EF `0.861/0.863` | LF `0.853/0.854`
- `@500` RN `0.867/0.868` | EF `0.871/0.872` | LF `0.860/0.861`
