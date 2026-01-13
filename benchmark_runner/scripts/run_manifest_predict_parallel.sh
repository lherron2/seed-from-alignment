#!/usr/bin/env bash
set -euo pipefail

# Parallel predictor runner for a manifest.
#
# Runs `ssbench.cli predict run` per target with concurrency (xargs -P),
# relying on the predictor/CLI to be idempotent (skips if predictions already exist).
#
# Required env:
#   MANIFEST=...   CSV with at least target_id,pdb_id,chain_id
#   TRUTH_DIR=...  directory of truth JSONs named {target_id}.json
#
# Optional env:
#   PRED_DIR=...   predictions output directory
#   JOBS=...       parallel workers (default: 8)
#   PYTHON_BIN=... python interpreter (default: benchmark_runner/.venv/bin/python)
#   DEFAULTS_FILE=... config file (default: benchmark_runner/defaults.yaml)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"
DEFAULTS_FILE="${DEFAULTS_FILE:-$ROOT_DIR/defaults.yaml}"
export DEFAULTS_FILE

MANIFEST="${MANIFEST:-}"
TRUTH_DIR="${TRUTH_DIR:-}"
PRED_DIR="${PRED_DIR:-$ROOT_DIR/data/predictions/parallel_predict}"
JOBS="${JOBS:-8}"

if [[ -z "$MANIFEST" || -z "$TRUTH_DIR" ]]; then
  echo "Usage: MANIFEST=... TRUTH_DIR=... [PRED_DIR=...] [JOBS=...] $0" >&2
  exit 1
fi

mkdir -p "$PRED_DIR"

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
emit("LENGTH_ADAPTIVE", predict.get("length_adaptive"))
emit("INCLUDE_UNFIXED_SAMPLING", predict.get("include_unfixed_sampling"))
emit("MAX_SCAFFOLDS", predict.get("max_scaffolds"))
emit("MAX_SAMPLES_PER_SCAFFOLD", predict.get("max_samples_per_scaffold"))
emit("N_SAMPLES", predict.get("n_samples"))
emit("BURN_IN", predict.get("burn_in"))
emit("THIN", predict.get("thin"))
emit("BETA", predict.get("beta"))
emit("SEED", predict.get("seed"))
emit("MIN_LOOP_SEP", predict.get("min_loop_sep"))
emit("PK_ALPHA", predict.get("pk_alpha"))
emit("PAIR_PENALTY", predict.get("pair_penalty"))
emit("PAIR_PENALTY_SCALE", predict.get("pair_penalty_scale"))
emit("PAIR_PENALTY_MODE", predict.get("pair_penalty_mode"))
emit("PAIR_PENALTY_C0", predict.get("pair_penalty_c0"))
emit("PAIR_PENALTY_C1", predict.get("pair_penalty_c1"))
emit("PAIR_PENALTY_MIN", predict.get("pair_penalty_min"))
emit("PAIR_PENALTY_MAX", predict.get("pair_penalty_max"))
emit("COV_MODE", predict.get("cov_mode"))
emit("COV_ALPHA", predict.get("cov_alpha"))
emit("COV_MIN_POWER", predict.get("cov_min_power"))
emit("COV_FORBID_NEG", predict.get("cov_forbid_negative"))
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
emit("REFINE_MAX_SECONDS", predict.get("refine_max_seconds"))
emit("REFINE_MAX_STRUCTURES", predict.get("refine_max_structures"))
emit("REFINE_MIN_UNPAIRED", predict.get("refine_min_unpaired"))
emit("REFINE_END_MASK_STEP", predict.get("refine_end_mask_step"))
emit("REFINE_MAX_END_MASK_LEN", predict.get("refine_max_end_mask_len"))
emit("REFINE_MAX_HELICES_SEQUENTIAL", predict.get("refine_max_helices_sequential"))
emit("REFINE_MAX_HELICES_PAIRWISE", predict.get("refine_max_helices_pairwise"))
emit("REFINE_MAX_REGIONS", predict.get("refine_max_regions"))
emit("REFINE_MAX_SEEDS", predict.get("refine_max_seeds"))
emit("REFINE_MAX_SOLUTIONS", predict.get("refine_max_solutions"))
emit("REFINE_KISSING_CANDIDATES", predict.get("refine_kissing_candidates"))
PY
)"

eval "$CFG_VARS"

: "${DATAPATH:=}"
: "${RNASTRUCTURE_DOCKER_IMAGE:=}"
: "${TOP_K:=20}"
: "${LENGTH_ADAPTIVE:=false}"
: "${INCLUDE_UNFIXED_SAMPLING:=}"
: "${MAX_SCAFFOLDS:=20}"
: "${MAX_SAMPLES_PER_SCAFFOLD:=200}"

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
  "--max-scaffolds" "${MAX_SCAFFOLDS}"
  "--max-samples-per-scaffold" "${MAX_SAMPLES_PER_SCAFFOLD}"
  "--n-samples" "${N_SAMPLES}"
  "--burn-in" "${BURN_IN}"
  "--thin" "${THIN}"
  "--beta" "${BETA}"
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

if [[ -n "${PAIR_PENALTY:-}" ]]; then
  PRED_ARGS+=("--pair-penalty" "${PAIR_PENALTY}")
else
  PRED_ARGS+=("--pair-penalty-scale" "${PAIR_PENALTY_SCALE}")
fi
if [[ -n "${PAIR_PENALTY_MIN:-}" ]]; then
  PRED_ARGS+=("--pair-penalty-min" "${PAIR_PENALTY_MIN}")
fi
if [[ -n "${PAIR_PENALTY_MAX:-}" ]]; then
  PRED_ARGS+=("--pair-penalty-max" "${PAIR_PENALTY_MAX}")
fi
if [[ "${COV_FORBID_NEG:-false}" == "True" || "${COV_FORBID_NEG:-false}" == "true" ]]; then
  PRED_ARGS+=("--cov-forbid-negative")
fi
if [[ "${LENGTH_ADAPTIVE}" == "True" || "${LENGTH_ADAPTIVE}" == "true" ]]; then
  PRED_ARGS+=("--length-adaptive")
fi
if [[ -n "${INCLUDE_UNFIXED_SAMPLING:-}" ]]; then
  if [[ "${INCLUDE_UNFIXED_SAMPLING}" == "True" || "${INCLUDE_UNFIXED_SAMPLING}" == "true" ]]; then
    PRED_ARGS+=("--include-unfixed-sampling")
  else
    PRED_ARGS+=("--no-include-unfixed-sampling")
  fi
fi

PRED_CMD="$(printf '%q ' "${PRED_ARGS[@]}"){fasta} {outdir}"
export PRED_CMD

TARGET_IDS="$("$PYTHON_BIN" - <<'PY'
import csv, os
with open(os.environ["MANIFEST"], newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        tid = (row.get("target_id") or "").strip()
        if tid:
            print(tid)
PY
)"

export ROOT_DIR PYTHON_BIN MANIFEST TRUTH_DIR PRED_DIR TOP_K

echo "[INFO] targets=$(echo \"$TARGET_IDS\" | wc -w) jobs=${JOBS} pred_dir=${PRED_DIR}"

set +e
echo "$TARGET_IDS" | xargs -P "$JOBS" -I {} bash -c '
set -uo pipefail
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
' _ {}
XARGS_EC=$?
set -e
if [[ $XARGS_EC -ne 0 ]]; then
  echo "[WARN] xargs exited with code ${XARGS_EC}; some targets may not have run"
fi

echo "[DONE] Predictions in $PRED_DIR"
