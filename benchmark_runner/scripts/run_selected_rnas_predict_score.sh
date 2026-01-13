#!/usr/bin/env bash
set -euo pipefail

# Fast path: run predict+score+report for selected_rnas using existing truth.
# Intended for quick iteration and hyperparameter sweeps without re-downloading data.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"
DEFAULTS_FILE="${DEFAULTS_FILE:-$ROOT_DIR/defaults.yaml}"
export DEFAULTS_FILE

MANIFEST="${MANIFEST:-$ROOT_DIR/data/manifests/selected_rnas_truth.csv}"
TRUTH_DIR="${TRUTH_DIR:-$ROOT_DIR/data/truth_selected}"
PRED_DIR="${PRED_DIR:-$ROOT_DIR/data/predictions/selected_rnas_quick}"
RESULTS_DIR="${RESULTS_DIR:-$ROOT_DIR/data/results}"

mkdir -p "$PRED_DIR" "$RESULTS_DIR"

RUN_TAG="${RUN_TAG:-selected_rnas_$(date +%s)}"
LOG_FILE="${LOG_FILE:-$RESULTS_DIR/run_${RUN_TAG}.log}"
METRICS_FILE="${METRICS_FILE:-$RESULTS_DIR/metrics_${RUN_TAG}.csv}"
SUMMARY_FILE="${SUMMARY_FILE:-$RESULTS_DIR/summary_${RUN_TAG}.md}"

exec > >(tee "$LOG_FILE") 2>&1

CFG_VARS="$("$PYTHON_BIN" - <<'PY'
import os
import shlex
import yaml
from pathlib import Path

cfg_path = Path(os.environ["DEFAULTS_FILE"])
data = yaml.safe_load(cfg_path.read_text()) or {}
paths = data.get("paths", {})
predict = data.get("predict", {})
rnastructure = paths.get("rnastructure", {})

def emit(key, value):
    if value is None:
        return
    print(f"{key}={shlex.quote(str(value))}")

emit("RFAM_CM", paths.get("rfam_cm"))
emit("INFERNAL_BIN", paths.get("infernal_bin"))
emit("RSCAPE", paths.get("rscape"))
emit("ALLSUB_EXE", rnastructure.get("allsub"))
emit("DUPLEX_EXE", rnastructure.get("duplex"))
emit("FOLD_EXE", rnastructure.get("fold"))
emit("PARTITION_EXE", rnastructure.get("partition"))
emit("PROBPLOT_EXE", rnastructure.get("probplot"))
emit("DATAPATH", rnastructure.get("datapath"))
emit("RNASTRUCTURE_DOCKER_IMAGE", rnastructure.get("docker_image"))

emit("TOP_K", predict.get("top_k"))
emit("MIN_LOOP_SEP", predict.get("min_loop_sep"))
emit("PK_ALPHA", predict.get("pk_alpha"))
emit("COV_MODE", predict.get("cov_mode"))
emit("COV_ALPHA", predict.get("cov_alpha"))
emit("COV_MIN_POWER", predict.get("cov_min_power"))
emit("WEIGHT_CAL_METHOD", predict.get("weight_calibration_method"))
emit("WEIGHT_CAL_ZMAX", predict.get("weight_calibration_zmax"))
emit("WEIGHT_ALPHA_CORE", predict.get("weight_alpha_core"))
emit("WEIGHT_ALPHA_ALT", predict.get("weight_alpha_alt"))
emit("WEIGHT_ALPHA_COV", predict.get("weight_alpha_cov"))
emit("WEIGHT_ALPHA_THERMO", predict.get("weight_alpha_thermo"))
emit("THERMO_MODE", predict.get("thermo_mode"))
emit("THERMO_WEIGHT", predict.get("thermo_weight"))
emit("THERMO_MAX_STRUCTURES", predict.get("thermo_max_structures"))
emit("THERMO_MIN_COUNT", predict.get("thermo_min_count"))
emit("THERMO_MIN_PROB", predict.get("thermo_min_prob"))
emit("THERMO_LOG_EPS", predict.get("thermo_log_eps"))
emit("STEM_START_PENALTY_SCALE", predict.get("stem_start_penalty_scale"))
emit("STEM_LEN1_PENALTY_SCALE", predict.get("stem_len1_penalty_scale"))
emit("STEM_LEN2_PENALTY_SCALE", predict.get("stem_len2_penalty_scale"))
emit("STEM_LOG_REWARD_SCALE", predict.get("stem_log_reward_scale"))
emit("STEM_SUPPORT_QUANTILE", predict.get("stem_support_quantile"))
emit("PAIR_PENALTY_MODE", predict.get("pair_penalty_mode"))
emit("PAIR_PENALTY_C0", predict.get("pair_penalty_c0"))
emit("PAIR_PENALTY_C1", predict.get("pair_penalty_c1"))
emit("PAIR_PENALTY", predict.get("pair_penalty"))
emit("PAIR_PENALTY_SCALE", predict.get("pair_penalty_scale"))
emit("REFINE_MAX_SECONDS_DEFAULT", predict.get("refine_max_seconds"))
PY
)"

eval "$CFG_VARS"

: "${DATAPATH:=}"
: "${RNASTRUCTURE_DOCKER_IMAGE:=}"
: "${TOP_K:=100}"
: "${MIN_LOOP_SEP:=0}"
: "${PK_ALPHA:=0.5}"
: "${COV_MODE:=logE}"
: "${COV_ALPHA:=3.0}"
: "${COV_MIN_POWER:=0.1}"
: "${WEIGHT_CAL_METHOD:=robust_z}"
: "${WEIGHT_CAL_ZMAX:=3.0}"
: "${WEIGHT_ALPHA_CORE:=1.0}"
: "${WEIGHT_ALPHA_ALT:=1.0}"
: "${WEIGHT_ALPHA_COV:=1.0}"
: "${WEIGHT_ALPHA_THERMO:=1.0}"
: "${THERMO_MODE:=off}"
: "${THERMO_WEIGHT:=0.0}"
: "${THERMO_MAX_STRUCTURES:=50}"
: "${THERMO_MIN_COUNT:=2}"
: "${THERMO_MIN_PROB:=0.001}"
: "${THERMO_LOG_EPS:=1.0e-6}"
: "${STEM_START_PENALTY_SCALE:=0.05}"
: "${STEM_LEN1_PENALTY_SCALE:=0.10}"
: "${STEM_LEN2_PENALTY_SCALE:=0.02}"
: "${STEM_LOG_REWARD_SCALE:=0.02}"
: "${STEM_SUPPORT_QUANTILE:=0.5}"
: "${PAIR_PENALTY_MODE:=length_aware}"
: "${PAIR_PENALTY_C0:=0.10}"
: "${PAIR_PENALTY_C1:=0.50}"
: "${PAIR_PENALTY:=}"
: "${PAIR_PENALTY_SCALE:=0.25}"

resolve_from_root() {
  local p="${1:-}"
  if [[ -z "$p" ]]; then
    echo ""
  elif [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$ROOT_DIR/$p"
  fi
}

ALLSUB_EXE="$(resolve_from_root "$ALLSUB_EXE")"
DUPLEX_EXE="$(resolve_from_root "$DUPLEX_EXE")"
FOLD_EXE="$(resolve_from_root "$FOLD_EXE")"
PARTITION_EXE="$(resolve_from_root "$PARTITION_EXE")"
PROBPLOT_EXE="$(resolve_from_root "$PROBPLOT_EXE")"

# Defaults tuned to keep wall time low; override via env vars as needed.
: "${N_SAMPLES:=50}"
: "${BURN_IN:=100}"
: "${THIN:=20}"
: "${SEED:=0}"
: "${REFINE_MAX_SECONDS:=${REFINE_MAX_SECONDS_DEFAULT:-5}}"
: "${REFINE_MAX_STRUCTURES:=50}"
: "${REFINE_MIN_UNPAIRED:=10}"
: "${REFINE_END_MASK_STEP:=30}"
: "${REFINE_MAX_END_MASK_LEN:=10}"
: "${REFINE_MAX_HELICES_SEQUENTIAL:=10}"
: "${REFINE_MAX_HELICES_PAIRWISE:=10}"
: "${REFINE_MAX_REGIONS:=2}"
: "${REFINE_MAX_SEEDS:=8}"
: "${REFINE_MAX_SOLUTIONS:=20}"
: "${REFINE_KISSING_CANDIDATES:=0}"

PRED_ARGS=(
  "DATAPATH=${DATAPATH}"
  "RNASTRUCTURE_DOCKER_IMAGE=${RNASTRUCTURE_DOCKER_IMAGE}"
  "PYTHONPATH=$ROOT_DIR/src"
  "$PYTHON_BIN"
  "-m"
  "ssbench.predict.cacofold_mcmc_pipeline"
  "--cm-db" "${RFAM_CM}"
  "--rscape" "${RSCAPE}"
  "--infernal-bin" "${INFERNAL_BIN}"
  "--allsub-exe" "${ALLSUB_EXE}"
  "--duplex-exe" "${DUPLEX_EXE}"
  "--fold-exe" "${FOLD_EXE}"
)

if [[ -n "${PARTITION_EXE:-}" ]]; then
  PRED_ARGS+=("--partition-exe" "${PARTITION_EXE}")
fi
if [[ -n "${PROBPLOT_EXE:-}" ]]; then
  PRED_ARGS+=("--probplot-exe" "${PROBPLOT_EXE}")
fi

PRED_ARGS+=(
  "--top-k" "${TOP_K}"
  "--n-samples" "${N_SAMPLES}"
  "--burn-in" "${BURN_IN}"
  "--thin" "${THIN}"
  "--beta" "1.0"
  "--seed" "${SEED}"
  "--min-loop-sep" "${MIN_LOOP_SEP}"
  "--pk-alpha" "${PK_ALPHA}"
  "--pair-penalty-mode" "${PAIR_PENALTY_MODE}"
  "--pair-penalty-c0" "${PAIR_PENALTY_C0}"
  "--pair-penalty-c1" "${PAIR_PENALTY_C1}"
  "--cov-mode" "${COV_MODE}"
  "--cov-alpha" "${COV_ALPHA}"
  "--cov-min-power" "${COV_MIN_POWER}"
  "--weight-calibration-method" "${WEIGHT_CAL_METHOD}"
  "--weight-calibration-zmax" "${WEIGHT_CAL_ZMAX}"
  "--weight-alpha-core" "${WEIGHT_ALPHA_CORE}"
  "--weight-alpha-alt" "${WEIGHT_ALPHA_ALT}"
  "--weight-alpha-cov" "${WEIGHT_ALPHA_COV}"
  "--weight-alpha-thermo" "${WEIGHT_ALPHA_THERMO}"
  "--thermo-mode" "${THERMO_MODE}"
  "--thermo-weight" "${THERMO_WEIGHT}"
  "--thermo-max-structures" "${THERMO_MAX_STRUCTURES}"
  "--thermo-min-count" "${THERMO_MIN_COUNT}"
  "--thermo-min-prob" "${THERMO_MIN_PROB}"
  "--thermo-log-eps" "${THERMO_LOG_EPS}"
  "--stem-start-penalty-scale" "${STEM_START_PENALTY_SCALE}"
  "--stem-len1-penalty-scale" "${STEM_LEN1_PENALTY_SCALE}"
  "--stem-len2-penalty-scale" "${STEM_LEN2_PENALTY_SCALE}"
  "--stem-log-reward-scale" "${STEM_LOG_REWARD_SCALE}"
  "--stem-support-quantile" "${STEM_SUPPORT_QUANTILE}"
  "--refine-max-seconds" "${REFINE_MAX_SECONDS}"
  "--refine-max-structures" "${REFINE_MAX_STRUCTURES}"
  "--refine-min-unpaired" "${REFINE_MIN_UNPAIRED}"
  "--refine-end-mask-step" "${REFINE_END_MASK_STEP}"
  "--refine-max-end-mask-len" "${REFINE_MAX_END_MASK_LEN}"
  "--refine-max-helices-sequential" "${REFINE_MAX_HELICES_SEQUENTIAL}"
  "--refine-max-helices-pairwise" "${REFINE_MAX_HELICES_PAIRWISE}"
  "--refine-max-regions" "${REFINE_MAX_REGIONS}"
  "--refine-max-seeds" "${REFINE_MAX_SEEDS}"
  "--refine-max-solutions" "${REFINE_MAX_SOLUTIONS}"
  "--refine-kissing-candidates" "${REFINE_KISSING_CANDIDATES}"
)

if [[ -n "${PAIR_PENALTY}" ]]; then
  PRED_ARGS+=("--pair-penalty" "${PAIR_PENALTY}")
else
  PRED_ARGS+=("--pair-penalty-scale" "${PAIR_PENALTY_SCALE}")
fi

PRED_CMD="$(printf '%q ' "${PRED_ARGS[@]}"){fasta} {outdir}"

echo "[RUN] predict run -> $PRED_DIR"
START_TS=$(date +%s)
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli predict run \
  --manifest "$MANIFEST" \
  --truth "$TRUTH_DIR" \
  --out-dir "$PRED_DIR" \
  --k "$TOP_K" \
  --predictor-cmd "$PRED_CMD"

echo "[RUN] score -> $METRICS_FILE"
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli score \
  --manifest "$MANIFEST" \
  --truth "$TRUTH_DIR" \
  --predictions "$PRED_DIR" \
  --out "$METRICS_FILE"

echo "[RUN] report -> $SUMMARY_FILE"
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli report \
  --metrics "$METRICS_FILE" \
  --out "$SUMMARY_FILE"

END_TS=$(date +%s)
echo "[RUN] done in $((END_TS-START_TS))s"

cp -f "$METRICS_FILE" "$RESULTS_DIR/metrics_selected_rnas_quick.csv"
cp -f "$SUMMARY_FILE" "$RESULTS_DIR/summary_selected_rnas_quick.md"
