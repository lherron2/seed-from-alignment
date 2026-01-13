#!/usr/bin/env bash
set -euo pipefail

# Full u400 benchmark: build top-500 predictions for baselines (LF/EF/RS)
# and an RNAnneal-ss + LF + EF ensemble, then (separately) score @K.
#
# This script is resumable: `ssbench.cli predict run` skips targets that already
# have >=K structures in predictions.db.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"
# Default to a conservative parallelism: RNAstructure(AllSub) can be heavy and
# aggressive concurrency has been observed to destabilize WSL on some setups.
JOBS="${JOBS:-4}"

MANIFEST="${MANIFEST:-$ROOT_DIR/data/manifests/bgsu_rep_current_4p0_truth_u400_structured_no_rrna.csv}"
TRUTH_DIR="${TRUTH_DIR:-$ROOT_DIR/data/truth_bgsu_rep_current_4p0_full}"

# Existing RNAnneal-ss run (top-100) used as the pipeline source in the ensemble.
PIPELINE_PREDS="${PIPELINE_PREDS:-$ROOT_DIR/data/predictions/bgsu_rep_current_4p0_u400_structured_no_rrna_lenadap}"

OUT_ROOT="${OUT_ROOT:-$ROOT_DIR/data/predictions/u400_topk500_v1}"
K="${K:-500}"

LF_DIR="${OUT_ROOT}/linearfold_zuker_k${K}"
EF_DIR="${OUT_ROOT}/eternafold_k${K}"
RS_DIR="${OUT_ROOT}/rnastructure_allsub_k${K}"
# NOTE: The original `${RS_DIR}` may contain legacy/broken outputs from an earlier run
# where `RNASTRUCTURE_DOCKER_IMAGE` was accidentally set to an empty string. Write
# fixed results to a new directory to avoid resumable-run skipping.
RS_DIR_FIXED="${OUT_ROOT}/rnastructure_allsub_k${K}_fixed"

# Default to skipping RNAstructure (it is expensive and can destabilize WSL under concurrency).
RUN_RNASTRUCTURE="${RUN_RNASTRUCTURE:-0}"

ENSEMBLE_DIR="${OUT_ROOT}/ensemble_rnanneal_ss_lf_ef_k${K}"
if [[ "$RUN_RNASTRUCTURE" == "1" ]]; then
  ENSEMBLE_DIR="${OUT_ROOT}/ensemble_rnanneal_ss_lf_ef_rs_k${K}"
fi

mkdir -p "$OUT_ROOT"

echo "[1/4] LinearFold-V (+Zuker) top-${K}"
MANIFEST="$MANIFEST" \
TRUTH_DIR="$TRUTH_DIR" \
PRED_DIR="$LF_DIR" \
TOP_K="$K" \
JOBS="$JOBS" \
PRED_CMD="$(printf '%q ' "$PYTHON_BIN" -m ssbench.predict.linearfold_zuker --k "$K" --beamsize 100 --delta 80.0 --max-seconds 120){fasta} {outdir}" \
  bash "$ROOT_DIR/scripts/run_manifest_predict_parallel_generic.sh"

echo "[2/4] EternaFold top-${K}"
MANIFEST="$MANIFEST" \
TRUTH_DIR="$TRUTH_DIR" \
PRED_DIR="$EF_DIR" \
TOP_K="$K" \
JOBS="$JOBS" \
PRED_CMD="$(printf '%q ' "$PYTHON_BIN" -m ssbench.predict.eternafold_subopt --k "$K"){fasta} {outdir}" \
  bash "$ROOT_DIR/scripts/run_manifest_predict_parallel_generic.sh"

if [[ "$RUN_RNASTRUCTURE" == "1" ]]; then
  echo "[3/4] RNAstructure AllSub top-${K}"
  # Use the docker-wrapped RNAstructure tools so we don't depend on a host install.
  ALLSUB_EXE="${ALLSUB_EXE:-$ROOT_DIR/tools/rnastructure_docker/AllSub}"
  FOLD_EXE="${FOLD_EXE:-$ROOT_DIR/tools/rnastructure_docker/Fold}"

  MANIFEST="$MANIFEST" \
  TRUTH_DIR="$TRUTH_DIR" \
  PRED_DIR="$RS_DIR_FIXED" \
  TOP_K="$K" \
  JOBS="${RS_JOBS:-$JOBS}" \
  DATAPATH="${DATAPATH:-}" \
  RNASTRUCTURE_DOCKER_IMAGE="${RNASTRUCTURE_DOCKER_IMAGE:-rnanneal/rnastructure:v6.4.0}" \
  RNASTRUCTURE_DOCKER_CPUS="${RNASTRUCTURE_DOCKER_CPUS:-1}" \
  RNASTRUCTURE_DOCKER_MEMORY="${RNASTRUCTURE_DOCKER_MEMORY:-2g}" \
  RNASTRUCTURE_DOCKER_PIDS_LIMIT="${RNASTRUCTURE_DOCKER_PIDS_LIMIT:-256}" \
  PRED_CMD="$(printf '%q ' "$PYTHON_BIN" -m ssbench.predict.rnastructure_allsub --k "$K" --allsub-exe "$ALLSUB_EXE" --fold-exe "$FOLD_EXE" --abs 80.0 --max-allsub-len "${MAX_ALLSUB_LEN:-250}" --max-seconds 120){fasta} {outdir}" \
    bash "$ROOT_DIR/scripts/run_manifest_predict_parallel_generic.sh"
else
  echo "[3/4] RNAstructure AllSub skipped (RUN_RNASTRUCTURE=0)"
fi

echo "[4/4] Build ensemble top-${K}"
"$PYTHON_BIN" "$ROOT_DIR/scripts/build_fr3d_u400_ensemble_predictions.py" \
  --manifest "$MANIFEST" \
  --out-dir "$ENSEMBLE_DIR" \
  --pipeline-preds "$PIPELINE_PREDS" \
  --linearfold-preds "$LF_DIR" \
  --eternafold-preds "$EF_DIR" \
  --k "$K"

if [[ "$RUN_RNASTRUCTURE" == "1" ]]; then
  "$PYTHON_BIN" "$ROOT_DIR/scripts/build_fr3d_u400_ensemble_predictions.py" \
    --manifest "$MANIFEST" \
    --out-dir "${OUT_ROOT}/ensemble_rnanneal_ss_lf_ef_rs_k${K}" \
    --pipeline-preds "$PIPELINE_PREDS" \
    --linearfold-preds "$LF_DIR" \
    --eternafold-preds "$EF_DIR" \
    --rnastructure-preds "$RS_DIR_FIXED" \
    --k "$K"
fi

echo "[DONE] u400 top-${K} predictions written under $OUT_ROOT"
