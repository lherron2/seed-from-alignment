#!/usr/bin/env bash
set -euo pipefail

# Resolve paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="${SCRIPT_DIR}/.."

STO="${SCRIPT_DIR}/combined_1.cacofold.sto"
COV="${SCRIPT_DIR}/combined_1.cacofold.cov"
SEQ="6VMY_1_Chain"
OUT_DB="${SCRIPT_DIR}/6VMY_1_Chain_beta1.0_seed0.db"
SUMMARY="${SCRIPT_DIR}/6VMY_summary.txt"


echo "[INFO] Sampling structures for ${SEQ}"
python "${ROOT_DIR}/sample_cli.py" "${STO}" \
  --seq "${SEQ}" \
  --n-samples 100 \
  --thin 10000 \
  --beta 0.5 \
  --cov-file "${COV}" \
  --cov-mode logE \
  --cov-alpha 1 \
  --cov-min-power 0 \
  --out-db "${OUT_DB}" \
  --summary "${SUMMARY}" \
  --cov-forbid-negative

echo "[INFO] Done. Output: ${OUT_DB}"

