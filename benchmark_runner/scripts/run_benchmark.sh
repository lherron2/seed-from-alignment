#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULTS_FILE="${DEFAULTS_FILE:-$ROOT_DIR/defaults.yaml}"
export DEFAULTS_FILE
LOG_FILE="${LOG_FILE:-$ROOT_DIR/data/results/benchmark.log}"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="${PYTHON_BIN:-python}"
fi

if [[ ! -f "$DEFAULTS_FILE" ]]; then
  echo "Defaults file not found: $DEFAULTS_FILE" >&2
  exit 1
fi

mkdir -p "$(dirname "$LOG_FILE")"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[LOG] Writing to $LOG_FILE"

CFG_VARS="$("$PYTHON_BIN" - <<'PY'
import os
import shlex
import yaml
from pathlib import Path

cfg_path = Path(os.environ["DEFAULTS_FILE"])
data = yaml.safe_load(cfg_path.read_text()) or {}
dataset = data.get("dataset", {})
score = data.get("score", {})
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

emit("DATASET_MIN_LEN", dataset.get("min_length"))
emit("DATASET_MAX_LEN", dataset.get("max_length"))
emit("HELIX_OVERLAP", score.get("helix_overlap_threshold"))
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
emit("COV_FORBID_NEGATIVE", predict.get("cov_forbid_negative"))
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

: "${RFAM_CM:=}"
: "${INFERNAL_BIN:=}"
: "${RSCAPE:=}"
: "${ALLSUB_EXE:=}"
: "${DUPLEX_EXE:=}"
: "${FOLD_EXE:=}"
: "${PARTITION_EXE:=}"
: "${PROBPLOT_EXE:=}"
: "${DATAPATH:=}"

: "${DATASET_MIN_LEN:=}"
: "${DATASET_MAX_LEN:=}"
: "${HELIX_OVERLAP:=}"
: "${TOP_K:=20}"
: "${LENGTH_ADAPTIVE:=false}"
: "${INCLUDE_UNFIXED_SAMPLING:=}"
: "${MAX_SCAFFOLDS:=20}"
: "${MAX_SAMPLES_PER_SCAFFOLD:=200}"
: "${N_SAMPLES:=200}"
: "${BURN_IN:=500}"
: "${THIN:=10}"
: "${BETA:=1.0}"
: "${SEED:=0}"
: "${MIN_LOOP_SEP:=0}"
: "${PK_ALPHA:=0.5}"
: "${PAIR_PENALTY:=}"
: "${PAIR_PENALTY_SCALE:=0.25}"
: "${PAIR_PENALTY_MODE:=legacy}"
: "${PAIR_PENALTY_C0:=0.10}"
: "${PAIR_PENALTY_C1:=0.50}"
: "${PAIR_PENALTY_MIN:=}"
: "${PAIR_PENALTY_MAX:=}"
: "${COV_MODE:=logE}"
: "${COV_ALPHA:=3.0}"
: "${COV_MIN_POWER:=0.1}"
: "${COV_FORBID_NEGATIVE:=false}"
: "${WEIGHT_CAL_METHOD:=none}"
: "${WEIGHT_CAL_ZMAX:=3.0}"
: "${WEIGHT_ALPHA_CORE:=1.0}"
: "${WEIGHT_ALPHA_ALT:=1.0}"
: "${WEIGHT_ALPHA_COV:=1.0}"
: "${WEIGHT_ALPHA_THERMO:=1.0}"
: "${THERMO_MODE:=allsub}"
: "${THERMO_WEIGHT:=1.0}"
: "${THERMO_MAX_STRUCTURES:=50}"
: "${THERMO_MIN_COUNT:=2}"
: "${THERMO_MIN_PROB:=0.001}"
: "${THERMO_LOG_EPS:=1.0e-6}"
: "${STEM_START_PENALTY_SCALE:=0.0}"
: "${STEM_LEN1_PENALTY_SCALE:=0.0}"
: "${STEM_LEN2_PENALTY_SCALE:=0.0}"
: "${STEM_LOG_REWARD_SCALE:=0.0}"
: "${STEM_SUPPORT_QUANTILE:=0.5}"
: "${REFINE_MAX_STRUCTURES:=200}"
: "${REFINE_MIN_UNPAIRED:=6}"
: "${REFINE_END_MASK_STEP:=20}"
: "${REFINE_MAX_END_MASK_LEN:=10}"
: "${REFINE_MAX_HELICES_SEQUENTIAL:=20}"
: "${REFINE_MAX_HELICES_PAIRWISE:=20}"
: "${REFINE_MAX_REGIONS:=4}"
: "${REFINE_MAX_SEEDS:=100}"
: "${REFINE_MAX_SOLUTIONS:=100}"
: "${REFINE_KISSING_CANDIDATES:=100}"

FASTA_PATH="${FASTA_PATH:-$ROOT_DIR/../selected_rnas.fasta}"
DATA_DIR="${DATA_DIR:-$ROOT_DIR/data}"
MANIFEST_DIR="${DATA_DIR}/manifests"
PDB_DIR="${DATA_DIR}/pdb"
TRUTH_DIR="${DATA_DIR}/truth_selected"
PRED_DIR="${DATA_DIR}/predictions/selected_rnas_benchmark"
RESULTS_DIR="${DATA_DIR}/results"

mkdir -p "$MANIFEST_DIR" "$PDB_DIR" "$TRUTH_DIR" "$PRED_DIR" "$RESULTS_DIR"

MANIFEST_RAW="${MANIFEST_DIR}/selected_rnas.csv"
MANIFEST_FILTERED="${MANIFEST_DIR}/selected_rnas_filtered.csv"
MANIFEST_FETCHED="${MANIFEST_DIR}/selected_rnas_fetched.csv"
MANIFEST_TRUTH="${MANIFEST_DIR}/selected_rnas_truth.csv"

if [[ ! -f "$FASTA_PATH" ]]; then
  echo "FASTA not found: $FASTA_PATH" >&2
  exit 1
fi

echo "[1/6] Build manifest from FASTA"
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" "$ROOT_DIR/scripts/manifest_from_fasta.py" \
  --fasta "$FASTA_PATH" \
  --out-manifest "$MANIFEST_RAW" \
  --pdb-dir "$PDB_DIR" \
  --download \
  --allow-unmatched \
  --unmatched-out "$MANIFEST_DIR/selected_rnas_unmatched.txt"

echo "[2/6] Filter manifest by length"
"$PYTHON_BIN" - <<PY
import os
import pandas as pd
df = pd.read_csv("$MANIFEST_RAW")
max_len = os.environ.get("DATASET_MAX_LEN")
min_len = os.environ.get("DATASET_MIN_LEN")
if max_len:
    df = df[df["nts_observed"] <= int(max_len)]
if min_len:
    df = df[df["nts_observed"] >= int(min_len)]
df.to_csv("$MANIFEST_FILTERED", index=False)
PY

echo "[3/6] Fetch PDB chains"
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli dataset fetch-pdb \
  --manifest "$MANIFEST_FILTERED" \
  --out-dir "$PDB_DIR" \
  --out-manifest "$MANIFEST_FETCHED"

echo "[4/6] Build truth with Barnaba"
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli truth build \
  --manifest "$MANIFEST_FETCHED" \
  --pdb-dir "$PDB_DIR" \
  --out-dir "$TRUTH_DIR"

echo "[5/6] Filter manifest to truth"
"$PYTHON_BIN" - <<PY
from pathlib import Path
import pandas as pd
df = pd.read_csv("$MANIFEST_FETCHED")
truth_dir = Path("$TRUTH_DIR")
df = df[df["target_id"].apply(lambda t: (truth_dir / f"{t}.json").exists())]
df.to_csv("$MANIFEST_TRUTH", index=False)
PY

PRED_ARGS=(
  "DATAPATH=${DATAPATH}"
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

if [[ -n "${PARTITION_EXE}" ]]; then
  PRED_ARGS+=("--partition-exe" "${PARTITION_EXE}")
fi

if [[ -n "${PROBPLOT_EXE}" ]]; then
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
)

if [[ -n "${PAIR_PENALTY_MIN}" ]]; then
  PRED_ARGS+=("--pair-penalty-min" "${PAIR_PENALTY_MIN}")
fi

if [[ -n "${PAIR_PENALTY_MAX}" ]]; then
  PRED_ARGS+=("--pair-penalty-max" "${PAIR_PENALTY_MAX}")
fi

PRED_ARGS+=(
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
  "{fasta}" "{outdir}"
)

PRED_CMD="$(printf '%q ' "${PRED_ARGS[@]}")"

if [[ "${COV_FORBID_NEGATIVE}" == "True" || "${COV_FORBID_NEGATIVE}" == "true" ]]; then
  PRED_CMD="${PRED_CMD/--refine-kissing-candidates/--cov-forbid-negative --refine-kissing-candidates}"
fi

if [[ -n "${PAIR_PENALTY}" ]]; then
  PRED_CMD="${PRED_CMD/--refine-kissing-candidates/--pair-penalty ${PAIR_PENALTY} --refine-kissing-candidates}"
else
  PRED_CMD="${PRED_CMD/--refine-kissing-candidates/--pair-penalty-scale ${PAIR_PENALTY_SCALE} --refine-kissing-candidates}"
fi

if [[ "${LENGTH_ADAPTIVE}" == "True" || "${LENGTH_ADAPTIVE}" == "true" ]]; then
  PRED_CMD="${PRED_CMD/--refine-kissing-candidates/--length-adaptive --refine-kissing-candidates}"
fi
if [[ -n "${INCLUDE_UNFIXED_SAMPLING:-}" ]]; then
  if [[ "${INCLUDE_UNFIXED_SAMPLING}" == "True" || "${INCLUDE_UNFIXED_SAMPLING}" == "true" ]]; then
    PRED_CMD="${PRED_CMD/--refine-kissing-candidates/--include-unfixed-sampling --refine-kissing-candidates}"
  else
    PRED_CMD="${PRED_CMD/--refine-kissing-candidates/--no-include-unfixed-sampling --refine-kissing-candidates}"
  fi
fi

echo "[6/6] Predict, score, and summarize"
PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli predict run \
  --manifest "$MANIFEST_TRUTH" \
  --truth "$TRUTH_DIR" \
  --out-dir "$PRED_DIR" \
  --k "$TOP_K" \
  --predictor-cmd "$PRED_CMD" \
  --config "$DEFAULTS_FILE"

"$PYTHON_BIN" -m ssbench.cli score \
  --manifest "$MANIFEST_TRUTH" \
  --truth "$TRUTH_DIR" \
  --predictions "$PRED_DIR" \
  --out "$RESULTS_DIR/metrics_selected_rnas.csv" \
  --config "$DEFAULTS_FILE"

"$PYTHON_BIN" -m ssbench.cli report \
  --metrics "$RESULTS_DIR/metrics_selected_rnas.csv" \
  --out "$RESULTS_DIR/summary_selected_rnas.md"

echo "Done. Results in $RESULTS_DIR"
