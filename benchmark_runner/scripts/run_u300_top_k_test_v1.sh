#!/usr/bin/env bash
set -euo pipefail

# Top-K output-size test for RNAnneal-ss (v3-high settings).
#
# Runs RNAnneal-ss once with --top-k 500, then downstream scoring/reporting can
# evaluate best-of-{1,50,100,200,500} by truncating the ranked output prefix.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"
JOBS="${JOBS:-6}"

MANIFEST="${MANIFEST:-$ROOT_DIR/data/manifests/bgsu_rep_current_4p0_truth_u300_representative50.csv}"
TRUTH_DIR="${TRUTH_DIR:-$ROOT_DIR/data/truth_bgsu_rep_current_4p0_full}"

OUT_DIR="${OUT_DIR:-$ROOT_DIR/data/predictions/top_k_test_v1/rnanneal_ss_ef_scaff_top500}"

# Reuse CaCoFold outputs to avoid cmscan/R-scape reruns.
REUSE_CACOFOLD_ROOT="${REUSE_CACOFOLD_ROOT:-$ROOT_DIR/data/predictions/u300_50_scaffold_validation/rnanneal_ss_rs_scaff}"

RFAM_CM="${RFAM_CM:-$ROOT_DIR/data/rfam/Rfam.cm}"
INFERNAL_BIN="${INFERNAL_BIN:-$REPO_ROOT/infernal-1.1.5/src}"
RSCAPE="${RSCAPE:-$REPO_ROOT/rscape_v2.6.4/src/R-scape}"

ALLSUB_EXE="${ALLSUB_EXE:-$ROOT_DIR/tools/scaffold_backends/eternafold/AllSub}"
FOLD_EXE="${FOLD_EXE:-$ROOT_DIR/tools/scaffold_backends/eternafold/Fold}"
DUPLEX_EXE="${DUPLEX_EXE:-$REPO_ROOT/RNAstructure/exe/DuplexFold}"

mkdir -p "$OUT_DIR"

base_args=(
  "$PYTHON_BIN" -m ssbench.predict.cacofold_mcmc_pipeline
  --cm-db "$RFAM_CM"
  --rscape "$RSCAPE"
  --infernal-bin "$INFERNAL_BIN"
  --reuse-cacofold-root "$REUSE_CACOFOLD_ROOT"
  --allsub-exe "$ALLSUB_EXE"
  --duplex-exe "$DUPLEX_EXE"
  --fold-exe "$FOLD_EXE"
  --top-k 500
  --n-samples 400
  --burn-in 200
  --thin 5
  --beta 1.0
  --seed 0
  --min-loop-sep 0
  --pk-alpha 0.5
  --pair-penalty-mode length_aware
  --pair-penalty-c0 0.10
  --pair-penalty-c1 0.50
  --cov-mode logE_power
  --cov-alpha 3.0
  --cov-min-power 0.1
  --weight-calibration-method robust_z
  --weight-calibration-zmax 3.0
  --weight-alpha-core 1.0
  --weight-alpha-alt 1.0
  --weight-alpha-cov 1.0
  --weight-alpha-thermo 1.0
  --thermo-mode allsub
  --thermo-weight 1.0
  --thermo-max-structures 120
  --thermo-min-count 2
  --thermo-min-prob 0.001
  --thermo-log-eps 1e-6
  --stem-start-penalty-scale 0.05
  --stem-len1-penalty-scale 0.10
  --stem-len2-penalty-scale 0.02
  --stem-log-reward-scale 0.02
  --stem-support-quantile 0.5
  --refine-max-seconds 6
  --refine-max-structures 200
  --refine-min-unpaired 6
  --refine-end-mask-step 10
  --refine-max-end-mask-len 30
  --refine-max-helices-sequential 20
  --refine-max-helices-pairwise 20
  --refine-max-regions 4
  --refine-max-seeds 50
  --refine-max-solutions 50
  --refine-kissing-candidates 50
  --max-scaffolds 15
  --max-samples-per-scaffold 120
  --length-adaptive
  --no-include-unfixed-sampling
  --inject-allsub-scaffolds
  --inject-allsub-scaffolds-max 25
  --inject-allsub-timeout-s 30
  --force-allsub-output 95
)

pred_cmd=""
for a in "${base_args[@]}"; do
  pred_cmd+="$(printf '%q ' "$a")"
done
pred_cmd+='{fasta} {outdir}'

echo "[RUN] top_k=500 -> ${OUT_DIR}"
MANIFEST="$MANIFEST" \
TRUTH_DIR="$TRUTH_DIR" \
PRED_DIR="$OUT_DIR" \
TOP_K=500 \
JOBS="$JOBS" \
PRED_CMD="$pred_cmd" \
  bash "$ROOT_DIR/scripts/run_manifest_predict_parallel_generic.sh"

echo "[DONE] Predictions written to $OUT_DIR"

