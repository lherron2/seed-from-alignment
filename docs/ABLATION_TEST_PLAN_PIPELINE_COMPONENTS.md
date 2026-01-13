# Ablation Plan: Candidate Generation vs Ranking/Scoring

This plan is designed to answer two questions separately:

1) **Generation:** which pipeline components are necessary to *produce* at least one “benchmark-like” structure among the candidates?
2) **Ranking:** which scoring/ranking elements are necessary to *prioritize* those good candidates near the top (top-1 / top-k)?

The core idea is to avoid confounding by separating:
- **oracle (generation) metrics** over an existing candidate pool, and
- **reranking-only (scoring) metrics** on a frozen candidate pool.

## Metrics to record (per target and mean over targets)

From `ssbench.cli score` (already produced in `metrics_*.csv`):
- `top1_f1`, `best_of_k_f1`, and `best_of_k_f1 - top1_f1` (ranking gap)
- `precision`, `recall`, `f1`, `mcc`
- `lonely_pair_rate` (fragmentation proxy)
- `coverage_{1,2,5,10,20,50,100}` (how quickly good candidates appear)

From **oracle candidate evaluation** (new; no ranking involved):
- `oracle_refined_best_f1`: best F1 over structures in `01_refined.db`
- `oracle_sampled_best_f1`: best F1 over structures in `02_sampled.db`
- `oracle_union_best_f1`: best F1 over `refined ∪ sampled`
- candidate counts per stage (to detect “good but absent” vs “good but rare”)

Interpretation:
- If `oracle_union_best_f1` is low: generation is the bottleneck.
- If `oracle_union_best_f1` is high but `top1_f1` is low: ranking is the bottleneck.

## Step 0: Establish a baseline candidate pool (frozen)

Run once to create the candidate DBs you’ll reuse for reranking ablations:

`benchmark_runner/scripts/run_selected_rnas_predict_score.sh` with a dedicated `PRED_DIR`.

Recommended (fast) baseline knobs:
- `TOP_K=100`
- `N_SAMPLES=50`, `BURN_IN=100`, `THIN=20`
- `REFINE_MAX_SECONDS=5` (tight runtime cap)

Example:
`RUN_TAG=ablation_candidates_baseline_$(date +%s) PRED_DIR=benchmark_runner/data/predictions/ablation_candidates_baseline benchmark_runner/scripts/run_selected_rnas_predict_score.sh`

## Step 1: Oracle candidate evaluation (generation attribution)

Compute “best achievable” F1 from each stage’s candidate DBs (no ranking):

`benchmark_runner/scripts/oracle_candidates_selected_rnas.py`
- Inputs: `--manifest`, `--truth-dir`, `--candidates-dir`
- Output: per-target CSV with `oracle_*` columns + summary aggregates.

This tells you whether:
- refinement is adding good structures,
- sampling is adding good structures, and/or
- both are failing to generate.

## Step 2: Reranking-only ablations (scoring attribution)

Freeze candidates and rerank them with different scoring settings.

Use:
- `benchmark_runner/src/ssbench/predict/rerank_existing_candidates.py` to write `predictions.db`
- then `ssbench.cli score` + `ssbench.cli report` as usual.

`benchmark_runner/scripts/ablation_rerank_selected_rnas.py` orchestrates this for a set of named ablations.

### Ranking ablations (recommended initial set)

Each ablation should change **one concept** at a time relative to the baseline scoring.

1) **No helix priors**
   - `--stem-start-penalty-scale 0`
   - `--stem-len1-penalty-scale 0`
   - `--stem-len2-penalty-scale 0`
   - `--stem-log-reward-scale 0`

2) **No pair penalty**
   - `--pair-penalty 0` (or `--pair-penalty-scale 0`)

3) **No calibration**
   - `--weight-calibration-method none`

4) **Component dropouts (reranking only)**
   - `--weight-alpha-cov 0`
   - `--weight-alpha-thermo 0`
   - `--weight-alpha-alt 0`
   - (optionally) `--weight-alpha-core 0` as a sanity check (should degrade severely)

5) **Thermo mode swap (reranking only)**
   - `--thermo-mode off` vs `--thermo-mode pf` vs `--thermo-mode allsub`

6) **Cov mode swap (reranking only)**
   - `--cov-mode off` vs your baseline cov mode

For each ablation, compare:
- mean `top1_f1` vs baseline
- change in ranking gap (`best_of_k_f1 - top1_f1`)
- per-target deltas (some targets will be very sensitive to cov / thermo; others to helix priors)

## Step 3: Generation ablations (optional, more expensive)

Once ranking sensitivity is understood, test “generation-side” components by regenerating candidates:

1) **Refinement off**
   - minimize refinement by setting `REFINE_MAX_SECONDS=0` and strict caps (or add a dedicated `--skip-refine` flag).

2) **Thermo off in generator**
   - `--thermo-mode off` and/or `--thermo-weight 0`

3) **Cov off in generator**
   - `--cov-mode off`

Compare oracle metrics (`oracle_*`) across these runs to see which components are responsible for making good candidates exist.

## Practical notes (speed + fairness)

- Keep `TOP_K` fixed (recommend 100) across ablations so `best_of_k_f1` is comparable.
- Keep MCMC budget fixed for generation ablations (`N_SAMPLES`, `BURN_IN`, `THIN`).
- For scoring-only ablations, do *not* regenerate candidates (rerank the baseline candidate pool).
- Use paired comparisons (same targets) and focus on `Δtop1_f1` and `Δgap`.

