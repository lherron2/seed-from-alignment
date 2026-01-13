# Benchmark Improvement Plan: Ranker + Evidence Calibration + Thermo PF + Stable Pair Penalty

This is an execution-focused plan for implementing improvements **(2)** helix-aware ranking priors (carefully, not overdone), **(3)** evidence calibration, **(4)** replace AllSub-count thermo support with partition-function base-pair probabilities (RNAstructure), and **(5)** a target-stable / length-aware `pair_penalty`.

Context for expected upside:
- Current `selected_rnas` mean `f1≈0.48` and mean `best_of_k_f1≈0.547` (see `benchmark_runner/data/results/metrics_selected_rnas.csv`), so **ranker-only** improvements have bounded mean headroom (~0.067), but can produce large wins on a subset of targets with big gaps.

## Goals (what “success” means)
- Raise **top-1** `f1` (and thus the reported benchmark metric) without needing major increases in sampling budget.
- Reduce overpairing + helix fragmentation in the **selected** structure.
- Keep scoring comparable across targets, especially when MSA/cov evidence is weak or absent.
- Keep changes “RNA-real” (general priors), avoid target-specific hacks.

## Non-goals / guardrails
- Do not tune hyperparameters directly on the 22 scored targets without a held-out split.
- Avoid aggressive helix priors that penalize true “short-helix” motifs.
- Keep defaults conservative; make new terms opt-in or weak-by-default, with clear ablations.

---

# Phase 0: Locate the current scoring path (for the implementing agent)

There are two distinct “energies” in play:
- **Sampler energy** (MCMC): `src/lib/energy.py` and delta updates.
- **Final ranker energy** (selects `top_k` dotbrackets written to `predictions.db`): implemented in `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py` inside `score_struct()`, currently:
  - calls `src/lib/energy.compute_energy(pairs, weights, params)` and
  - adds `pair_penalty * len(pairs)`.

For fastest impact on **top-1**, implement helix priors and improved `pair_penalty` first in the **final ranker** (full recompute per candidate is fine; no delta complexity).

---

# Phase 1 (Point 3): Evidence calibration + weight blending

## Problem
Weights are built from multiple evidence sources (SS_cons/core, alternates, covariation, thermo augmentation). Their scales drift across targets (and “MSA empty” vs “MSA present”), so a fixed `pair_penalty` and ranker objective behave inconsistently.

## Design
Introduce a “components → normalize → blend” pipeline:

1) Build per-pair components:
- `core`: base evidence from `SS_cons` (binary / count)
- `alt`: evidence from `SS_cons_1`, `SS_cons_2`, … (count)
- `cov`: covariation score (float; can be 0 if absent)
- `thermo`: thermodynamic support (later from PF probabilities)

2) Normalize each component **within a target** to a comparable scale:
- Preferred: robust z-score using median / MAD with clipping:
  - `z = clip((x - median)/MAD, -zmax, +zmax)`
  - Handle `MAD==0` by falling back to rank-based percentile or “no-op”.
- Keep a small set of normalization methods behind a config flag.

3) Blend with fixed coefficients:
- `w_total = α_core*z_core + α_alt*z_alt + α_cov*z_cov + α_thermo*z_thermo`
- Start with simple defaults (`α_* = 1.0`) and only tune if you can do a proper holdout.

## Implementation steps
1) Add a small utility module for robust normalization and blending:
   - New file suggestion: `src/lib/weight_calibration.py`
   - Provide:
     - `robust_center_scale(values) -> (median, mad)`
     - `robust_z(x, median, mad, zmax)`
     - `normalize_component(dict[(i,j)]=x, method=...) -> dict[(i,j)]=z`
     - `blend_components(components: dict[name->dict[(i,j)]=val], alphas: ...)`

2) Add a “component extraction” function without breaking existing APIs:
   - Option A (preferred): add a new function in `src/lib/sample_cacofold_structures.py`:
     - `extract_candidate_pair_components(...) -> (L, components_dict)`
     - Keep existing `extract_candidate_pairs()` unchanged.
   - Option B: implement the extraction logic inside `src/lib/pipeline.sample_pk()` and in the benchmark predictor (less reuse, but smaller blast radius).

3) Wire the calibrated weights into both:
   - MCMC sampling: `src/lib/pipeline.sample_pk()` (so the sampler sees stable weights)
   - Final ranking: `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py` (so ranker sees the same calibrated weights)

## Validation / safety checks
- Unit test normalization edge cases:
  - constant component values (MAD=0)
  - sparse components (few non-zeros)
- Ensure calibrated weights preserve ordering when only one source is present.

---

# Phase 2 (Point 2): Helix-aware ranking priors (conservative, evidence-aware)

## Problem
A major failure mode is **helix fragmentation**: many short stems and isolated pairs. The current ranker objective is mostly “sum of pair weights + generic PK/lonely penalties + pair-count penalty”, which doesn’t directly prefer helix-contiguous, RNA-real stem bundles.

## Key caution (don’t overdo it)
Some hard targets have genuine short helices; a strong penalty on len-1/len-2 stems can penalize the truth.

## Design principles
- Apply helix priors only in the **final ranking stage** first (cheap, isolated change).
- Make the len-1 penalty **evidence-conditioned**:
  - only penalize a len-1 stem if its pair has low support under calibrated weights (e.g., below median or below a quantile).
  - do not penalize strongly supported single pairs (which may be real motifs).

## Proposed stem statistics
Given a pair set `P` (with i<j), define stems as maximal runs `(i+k, j-k)` for k=0..(len-1).
Compute:
- `n_stems`
- `n_len1`, `n_len2`
- `sum_log1p_len = Σ log(1+len(stem))` (optional; concave)

## New ranker terms (start weak; scale by weight scale)
Add to the ranker energy:
- `+ stem_start_penalty * n_stems`
- `+ stem_len1_penalty * n_len1_unsupported` (unsupported only)
- `+ stem_len2_penalty * n_len2` (very small / optional)
- `- stem_log_reward * sum_log1p_len` (optional; keep small)

Where “unsupported len1” is defined using calibrated per-pair support:
- `unsupported = (w_total(i,j) < q_support)` where `q_support` might be median or 25th percentile of nonzero weights.

## Implementation steps
1) Add a stem parser/stats helper:
   - New file suggestion: `src/lib/stem_stats.py`
   - Public function: `stem_stats(pairs: set[tuple[int,int]]) -> StemStats`

2) Update final ranker scoring in:
   - `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py`:
     - Extend CLI/config with new hyperparameters:
       - `--stem-start-penalty-scale`
       - `--stem-len1-penalty-scale`
       - `--stem-support-quantile` (or threshold mode)
     - Compute `w_scale` (robust, from calibrated weights) and set penalties as `scale * w_scale`.

3) Keep the sampler energy unchanged initially (no delta complexity).

## Validation
- Add a unit test for `stem_stats()`:
  - simple nested stems, multiple stems, crossing pairs should not crash.
- Add a small ranker regression test (synthetic):
  - two candidate structures with same #pairs but one fragmented: verify the helix-prior score prefers the contiguous helix when evidence is similar.

---

# Phase 3 (Point 4): Replace AllSub counts with PF base-pair probabilities (RNAstructure)

## Why
AllSub “pair frequency” is a rough proxy and varies with parameter choices and structure enumeration details. PF probabilities are smoother, better calibrated, and generally reduce spurious isolated pairs.

## RNAstructure binaries to consider (beyond AllSub)
Within RNAstructure, the typical PF pipeline is:
- `partition`: compute partition function; outputs `.pfs`
- `ProbabilityPlot`: extract base-pair probabilities from `.pfs` (dot plot)
Optional (future, not required for point 4):
- `MaxExpect`: MEA structure from PF (can be used as an additional scaffold seed)
- `stochastic`: sample structures from PF (candidate diversity for sampling-limited targets)
- `ProbKnot`: PK prediction from probabilities (could provide low-weight PK hints)

## Implementation approach (PF → per-pair probabilities)
1) Add wrappers in `src/lib/refine_unpaired_regions.py` or a dedicated module:
   - Suggested new file: `src/lib/rnastructure_pf.py`
   - Functions:
     - `run_partition(partition_exe, fasta_path, out_pfs, env={DATAPATH})`
     - `run_probability_plot(probplot_exe, pfs_path, out_txt)`
     - `parse_probability_plot(out_txt) -> dict[(i,j)] = p_ij`
2) Integrate into thermo augmentation in both places:
   - `src/lib/pipeline.sample_pk()` (replaces current AllSub-count add-on)
   - `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py` (same)

## How to convert probabilities into a weight component
Start simple and stable:
- `thermo_raw(i,j) = log(p_ij + eps)` (or `p_ij` if you want bounded)
- then run the same per-source normalization as other components.

## Feature gating (avoid drowning other evidence)
- Only include canonical pairs.
- Ignore very small probabilities: `p_ij < p_min` (e.g. 1e-3) to keep sparse.
- Keep `α_thermo` small initially (e.g., 0.5–1.0 after normalization).

## Practical notes for the implementing agent
- The exact CLI signatures for `partition` / `ProbabilityPlot` vary by installation; verify via `--help` on the target machine.
- Ensure `DATAPATH` is set in the subprocess environment (benchmarks already do this in predictor commands).
- Keep AllSub as a fallback mode: `thermo_mode = {pf, allsub, off}`.

## Validation
- Unit test `parse_probability_plot()` with a small fixture file (no external binaries).
- Add a smoke-test path that runs the pipeline with `thermo_mode=off` (no dependency on RNAstructure PF during tests).

---

# Phase 4 (Point 5): Make `pair_penalty` length-aware and target-stable

## Problem
Current `pair_penalty` is `pair_penalty_scale * median(weight)` (per target). This drifts with:
- evidence availability (MSA/cov present vs absent)
- weight normalization changes
- sequence length (longer RNAs admit more potential pairings)

## Design
Compute a penalty in “weight units”, but modulated by length:
- Let `w_scale` be a robust scale of the calibrated weights (median or MAD of nonzero weights).
- Define:
  - `pair_penalty = w_scale * (c0 + c1 / sqrt(L))`
  - clamp to `[min_penalty, max_penalty]` in weight units
- Preserve an override:
  - if user sets `--pair-penalty`, use it directly.

## Implementation steps
1) Update `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py`:
   - Keep existing flags for backwards compatibility.
   - Add new flags:
     - `--pair-penalty-c0`, `--pair-penalty-c1`, `--pair-penalty-min`, `--pair-penalty-max`
     - or a single `--pair-penalty-mode {legacy,length_aware}`.
2) Update `benchmark_runner/defaults.yaml` to include the new settings.

## Validation
- Add a small unit test for the penalty computation function (pure math).

---

# Fast ablation plan (to measure what matters)

Run small controlled experiments on the same target set and sampling budget:
1) **Ranker-only helix priors** (no weight calibration, thermo off):
   - Expect: `f1↑`, `best_of_k_f1≈` unchanged.
2) **Weight calibration only** (no helix priors, thermo off):
   - Expect: `f1↑` especially on “MSA-empty-like” targets.
3) **PF thermo only** (calibration on, helix priors off):
   - Expect: fewer isolated pairs; `f1↑` on some hard targets.
4) **Length-aware pair penalty** (with calibration):
   - Expect: less overpairing; modest `f1↑`.

Track per-target changes, not just mean, to ensure you’re not hurting short-helix truth cases.

