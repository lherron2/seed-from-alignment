#!/usr/bin/env bash
set -euo pipefail

# Generic parallel predictor runner for a manifest.
#
# Required env:
#   MANIFEST=...        CSV with at least target_id,pdb_id,chain_id
#   TRUTH_DIR=...       directory of truth JSONs named {target_id}.json
#   PRED_CMD="... {fasta} {outdir}"   predictor command template
#
# Optional env:
#   PRED_DIR=...        predictions output directory
#   JOBS=...            parallel workers (default: 8)
#   TOP_K=...           number of structures requested from predictor (default: 100)
#   PYTHON_BIN=...      python interpreter (default: benchmark_runner/.venv/bin/python)
#
# Notes:
# - Runs `ssbench.cli predict run` per target with concurrency (xargs -P).
# - Relies on the predictor/CLI being idempotent (skips if predictions already exist).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"

MANIFEST="${MANIFEST:-}"
TRUTH_DIR="${TRUTH_DIR:-}"
PRED_CMD="${PRED_CMD:-}"
PRED_DIR="${PRED_DIR:-$ROOT_DIR/data/predictions/parallel_predict_generic}"
JOBS="${JOBS:-8}"
TOP_K="${TOP_K:-100}"

if [[ -z "$MANIFEST" || -z "$TRUTH_DIR" || -z "$PRED_CMD" ]]; then
  echo "Usage: MANIFEST=... TRUTH_DIR=... PRED_CMD='... {fasta} {outdir}' [PRED_DIR=...] [JOBS=...] [TOP_K=...] $0" >&2
  exit 1
fi

mkdir -p "$PRED_DIR"

TARGETS="$("$PYTHON_BIN" - <<'PY'
import csv, os, sys
path = os.environ["MANIFEST"]
with open(path, newline="") as f:
    r = csv.DictReader(f)
    if not r.fieldnames or "target_id" not in r.fieldnames:
        raise SystemExit("manifest missing target_id column")
    for row in r:
        tid = (row.get("target_id") or "").strip()
        if tid:
            print(tid)
PY
)"

export MANIFEST TRUTH_DIR PRED_DIR TOP_K PRED_CMD PYTHON_BIN ROOT_DIR

echo "$TARGETS" | xargs -P "$JOBS" -n 1 bash -lc '
set -euo pipefail
target_id="$1"
safe_id="${target_id//|/_}"
tmp_manifest="$(mktemp)"
cleanup() { rm -f "$tmp_manifest"; }
trap cleanup EXIT

if ! "$PYTHON_BIN" - <<PY
import csv, os
in_path = os.environ["MANIFEST"]
target = "${target_id}"
with open(in_path, newline="") as f:
    r = csv.DictReader(f)
    rows = [row for row in r if (row.get("target_id") or "").strip() == target]
if not rows:
    raise SystemExit(f"Target not found in manifest: {target}")
with open("${tmp_manifest}", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=r.fieldnames or ["target_id","pdb_id","chain_id"])
    w.writeheader()
    for row in rows:
        w.writerow(row)
PY
then
  echo "[WARN] ${target_id} (failed to build per-target manifest)"
  exit 0
fi

log_path="$PRED_DIR/$safe_id/predict.log"
mkdir -p "$(dirname "$log_path")"

if PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli predict run \
  --manifest "$tmp_manifest" \
  --truth "$TRUTH_DIR" \
  --out-dir "$PRED_DIR" \
  --k "$TOP_K" \
  --predictor-cmd "$PRED_CMD" >"$log_path" 2>&1; then
  echo "[OK] ${target_id}"
else
  echo "[WARN] ${target_id} (see $log_path)"
fi
' _

echo "[DONE] Predictions in $PRED_DIR"
