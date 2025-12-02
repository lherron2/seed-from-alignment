#!/usr/bin/env bash
set -euo pipefail

# Example driver for the pre-wired pipelines defined in src/lib/pipeline.py.
#
# This script assumes it lives in the `example/` directory of the repo and
# that the CaCoFold files:
#   - combined_1.cacofold.sto
#   - combined_1.cacofold.cov
# are present alongside it (as they are in this repo).

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "${SCRIPT_DIR}/.." && pwd )"

CONFIG="${SCRIPT_DIR}/pipeline_config.yaml"

# Choose which wiring to run:
#   legacy       : get_consensus_db → sample_pk → refine → rosetta_ready
#   refine_first : get_consensus_db → refine → sample_pk → ensure_valid → rosetta_ready
MODE="refine_first"

echo "[EXAMPLE] Running pipeline in mode: ${MODE}"
echo "[EXAMPLE]   Config: ${CONFIG}"
echo "[EXAMPLE]   Repo root: ${ROOT_DIR}"

python "${ROOT_DIR}/run_pipeline.py" \
    --config "${CONFIG}" \
    --mode "${MODE}"

echo "[EXAMPLE] Done. See pipeline_output/ for generated .db files."

