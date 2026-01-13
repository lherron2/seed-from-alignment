#!/usr/bin/env bash
set -euo pipefail

# Streamed benchmark runner: predict+score per target (report as each finishes).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"
DEFAULTS_FILE="${DEFAULTS_FILE:-$ROOT_DIR/defaults.yaml}"
export DEFAULTS_FILE

MANIFEST="${MANIFEST:-}"
TRUTH_DIR="${TRUTH_DIR:-}"
PRED_DIR="${PRED_DIR:-$ROOT_DIR/data/predictions/stream_benchmark}"
RESULTS_DIR="${RESULTS_DIR:-$ROOT_DIR/data/results}"
RUN_TAG="${RUN_TAG:-stream_$(date +%s)}"

if [[ -z "$MANIFEST" || -z "$TRUTH_DIR" ]]; then
  echo "Usage: MANIFEST=... TRUTH_DIR=... [PRED_DIR=...] [RESULTS_DIR=...] $0" >&2
  exit 1
fi

mkdir -p "$PRED_DIR" "$RESULTS_DIR"

LOG_FILE="${LOG_FILE:-$RESULTS_DIR/run_${RUN_TAG}.log}"
if [[ -n "${NO_TEE:-}" || ! -t 1 ]]; then
  exec >>"$LOG_FILE" 2>&1
elif [[ -t 1 ]]; then
  exec > >(tee "$LOG_FILE") 2>&1
else
  exec >>"$LOG_FILE" 2>&1
fi

trap 'ec=$?; echo "[EXIT] status=$ec run_tag=$RUN_TAG log=$LOG_FILE"; exit $ec' EXIT

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

# Export key hyperparameters so the per-target result snippet can report them.
export TOP_K N_SAMPLES BURN_IN THIN BETA SEED MIN_LOOP_SEP PK_ALPHA

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

SUMMARY_CSV="$RESULTS_DIR/metrics_${RUN_TAG}_per_target.csv"
DIAG_CSV="$RESULTS_DIR/diagnostics_${RUN_TAG}_per_target.csv"
echo "target_id,precision,recall,f1,top1_f1,best_of_k_f1,best_of_k_precision,best_of_k_recall,top1_mcc,best_of_k_mcc,best_of_k_mcc_rank" > "$SUMMARY_CSV"
echo "target_id,best_f1_predictions,best_f1_refined,best_f1_sampled,best_f1_consensus,rank_best_refined_in_predictions,rank_best_sampled_in_predictions" > "$DIAG_CSV"

for target_id in $TARGET_IDS; do
  safe_id="${target_id//|/_}"
  tmp_manifest="$(mktemp)"
  cleanup() { rm -f "$tmp_manifest"; }
  "$PYTHON_BIN" - <<PY
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

  echo "[RUN] ${target_id}"
  if ! PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli predict run \
    --manifest "$tmp_manifest" \
    --truth "$TRUTH_DIR" \
    --out-dir "$PRED_DIR" \
    --k "$TOP_K" \
    --predictor-cmd "$PRED_CMD"; then
    echo "[WARN] predict failed for ${target_id}; continuing"
    cleanup
    continue
  fi

  metrics_path="$RESULTS_DIR/metrics_${RUN_TAG}_${safe_id}.csv"
  if ! PYTHONPATH="$ROOT_DIR/src" "$PYTHON_BIN" -m ssbench.cli score \
    --manifest "$tmp_manifest" \
    --truth "$TRUTH_DIR" \
    --predictions "$PRED_DIR" \
    --out "$metrics_path"; then
    echo "[WARN] score failed for ${target_id}; continuing"
    cleanup
    continue
  fi

  if ! "$PYTHON_BIN" - <<PY
import json
import os
from pathlib import Path

import pandas as pd

from ssbench.metrics.pair_metrics import compute_pair_metrics
from ssbench.predict.parse_dotbracket import parse_dotbracket, pairs_from_dotbracket

df = pd.read_csv("${metrics_path}")
if df.empty:
    raise SystemExit("No metrics rows produced for ${target_id}")
row = df.iloc[0]

truth_path = Path("${TRUTH_DIR}") / f"{row['target_id']}.json"
pred_path = Path("${PRED_DIR}") / "${safe_id}" / "predictions.db"
truth = json.loads(truth_path.read_text())
ref_pairs = {tuple(p) for p in truth["canonical_pairs"]}
seq_len = len(truth["sequence"])
_seq, dotbrs = parse_dotbracket(pred_path.read_text())
pred_pairs_list = [set(pairs_from_dotbracket(s)) for s in dotbrs] or [set()]

best_f1 = None
best_f1_precision = 0.0
best_f1_recall = 0.0

best_mcc = None
best_mcc_rank = 1

for rank, pairs in enumerate(pred_pairs_list, start=1):
    metrics = compute_pair_metrics(pairs, ref_pairs, seq_len)
    if best_f1 is None or metrics.f1 > best_f1:
        best_f1 = metrics.f1
        best_f1_precision = metrics.precision
        best_f1_recall = metrics.recall
    if best_mcc is None or metrics.mcc > best_mcc:
        best_mcc = metrics.mcc
        best_mcc_rank = rank

best_f1 = float(best_f1 or 0.0)
best_mcc = float(best_mcc or 0.0)

def best_f1_in_db(path: Path) -> tuple[float, str | None, int]:
    if not path.exists():
        return 0.0, None, 0
    _s, dotbrs_local = parse_dotbracket(path.read_text())
    best_local = 0.0
    best_struct = None
    for s in dotbrs_local:
        pairs = set(pairs_from_dotbracket(s))
        m = compute_pair_metrics(pairs, ref_pairs, seq_len)
        if m.f1 > best_local:
            best_local = m.f1
            best_struct = s
    return float(best_local), best_struct, len(dotbrs_local)

refined_path = Path("${PRED_DIR}") / "${safe_id}" / "01_refined.db"
sampled_path = Path("${PRED_DIR}") / "${safe_id}" / "02_sampled.db"
consensus_path = Path("${PRED_DIR}") / "${safe_id}" / "00_consensus.db"
best_refined_f1, best_refined_struct, _n_refined = best_f1_in_db(refined_path)
best_sampled_f1, best_sampled_struct, _n_sampled = best_f1_in_db(sampled_path)
best_consensus_f1, _best_consensus_struct, _n_consensus = best_f1_in_db(consensus_path)

rank_best_refined = -1
if best_refined_struct is not None:
    try:
        rank_best_refined = dotbrs.index(best_refined_struct) + 1
    except ValueError:
        rank_best_refined = -1
rank_best_sampled = -1
if best_sampled_struct is not None:
    try:
        rank_best_sampled = dotbrs.index(best_sampled_struct) + 1
    except ValueError:
        rank_best_sampled = -1

gap_refined = best_refined_f1 - best_f1
gap_sampled = best_sampled_f1 - best_f1
if gap_refined > 1e-6 or gap_sampled > 1e-6:
    print(
        f"[DIAG] {row['target_id']} best_f1(predictions)={best_f1:.3f} "
        f"best_f1(refined)={best_refined_f1:.3f} rank_refined_in_top{len(dotbrs)}={rank_best_refined} "
        f"best_f1(sampled)={best_sampled_f1:.3f} rank_sampled_in_top{len(dotbrs)}={rank_best_sampled} "
        f"best_f1(consensus)={best_consensus_f1:.3f}"
    )

print(
    f"[RESULT] {row['target_id']} f1={row['f1']:.3f} p={row['precision']:.3f} r={row['recall']:.3f} "
    f"top1={row['top1_f1']:.3f} best_k={row['best_of_k_f1']:.3f} "
    f"best_k_p={best_f1_precision:.3f} best_k_r={best_f1_recall:.3f} "
    f"top1_mcc={row['mcc']:.3f} best_k_mcc={best_mcc:.3f} best_k_mcc_rank={best_mcc_rank} "
    f"params: top_k={int(os.environ.get('TOP_K', '0') or 0)} n_samples={int(os.environ.get('N_SAMPLES', '0') or 0)} "
    f"burn_in={int(os.environ.get('BURN_IN', '0') or 0)} thin={int(os.environ.get('THIN', '0') or 0)} "
    f"beta={float(os.environ.get('BETA', '0') or 0)} seed={int(os.environ.get('SEED', '0') or 0)} "
    f"min_loop_sep={int(os.environ.get('MIN_LOOP_SEP', '0') or 0)} pk_alpha={float(os.environ.get('PK_ALPHA', '0') or 0)}"
)
with open("${SUMMARY_CSV}", "a") as f:
    f.write(
        "{},{},{},{},{},{},{},{},{},{},{}\n".format(
            row["target_id"],
            row["precision"],
            row["recall"],
            row["f1"],
            row["top1_f1"],
            row["best_of_k_f1"],
            best_f1_precision,
            best_f1_recall,
            row["mcc"],
            best_mcc,
            best_mcc_rank,
        )
    )
with open("${DIAG_CSV}", "a") as f:
    f.write(
        "{},{},{},{},{},{},{}\n".format(
            row["target_id"],
            best_f1,
            best_refined_f1,
            best_sampled_f1,
            best_consensus_f1,
            rank_best_refined,
            rank_best_sampled,
        )
    )
PY
  then
    echo "[WARN] postprocess failed for ${target_id}; continuing"
    cleanup
    continue
  fi

  cleanup
done

echo "[DONE] Per-target metrics: $SUMMARY_CSV"
