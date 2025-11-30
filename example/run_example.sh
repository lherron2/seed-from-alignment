#!/usr/bin/env bash
set -euo pipefail

# Resolve paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="${SCRIPT_DIR}/.."

STO="${SCRIPT_DIR}/combined_1.cacofold.sto"
COV="${SCRIPT_DIR}/combined_1.cacofold.cov"
SEQ="6VMY_1_Chain"
UNREFINED_DB="${SCRIPT_DIR}/6VMY_1_Chain_beta1.0_seed0.db"
REFINED_DB="${SCRIPT_DIR}/6VMY_1_Chain_beta1.0_seed0_refined.db"
ROSETTA_DB="${SCRIPT_DIR}/6VMY_1_rosetta.db"
SUMMARY="${SCRIPT_DIR}/6VMY_summary.txt"


echo "[INFO] Sampling structures for ${SEQ}"
python "${ROOT_DIR}/sample_cli.py" "${STO}" \
  --seq "${SEQ}" \
  --n-samples 100 \
  --thin 1000 \
  --beta 1 \
  --cov-file "${COV}" \
  --cov-mode logE \
  --cov-alpha 2 \
  --cov-min-power 0 \
  --out-db "${UNREFINED_DB}" \
  --summary "${SUMMARY}" \
  --pk-filter-frac 0.2 \
  --cov-forbid-negative

python "${ROOT_DIR}/refine_unpaired_regions.py" \
  --sto "${STO}" \
  --seq-name "${SEQ}" \
  --db-in "${UNREFINED_DB}"\
  --db-out "${REFINED_DB}" \
  --fold-exe "${PROJECT}/repos/RNAstructure/exe/Fold" \
  --duplex-exe "${PROJECT}/repos/RNAstructure/exe/DuplexFold" \
  --min-unpaired 8 \
  --min-terminal-unpaired 1 
   
python "${ROOT_DIR}/filter_db_for_rosetta.py" \
  --db-in "${REFINED_DB}" \
  --db-out "${ROSETTA_DB}" \
  --sto "${STO}" \
  --seq-name "${SEQ}"
