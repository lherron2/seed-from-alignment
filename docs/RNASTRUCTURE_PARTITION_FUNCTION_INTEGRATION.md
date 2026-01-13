# RNAstructure PF integration notes (Partition / ProbabilityPlot / MaxExpect / Stochastic)

This file is a tactical companion to `docs/BENCHMARK_IMPROVEMENT_PLAN_RANKER_CALIBRATION.md`, focused on point (4): replacing AllSub-count “thermo support” with partition-function (PF) base-pair probabilities.

## 1) Verify available binaries + CLI

On the target machine, check which of these exist:
- `Partition` (sometimes `partition`)
- `ProbabilityPlot`
- `MaxExpect`
- `Stochastic`
- `ProbKnot`

Record exact `--help` output and required positional args; RNAstructure packaging differs.

### Environment
RNAstructure usually requires `DATAPATH` set to the thermodynamic tables directory.
In this repo, the benchmark predictor already sets `DATAPATH` in the command (see `benchmark_runner/defaults.yaml` and `benchmark_runner/data/results/benchmark.log`).

## 2) Minimal PF pipeline for base-pair probabilities

Goal: compute `p_ij` for the full ungapped sequence.

Typical flow (verify exact arguments):
1) `Partition <input_fasta> <out.pfs> [options]`
2) `ProbabilityPlot <out.pfs> <out.txt> [options]`

Notes:
- Some versions accept FASTA directly; others prefer `.seq`. If `.seq` is required:
  - use RNAstructure’s `Sequence`/conversion tool if available, or write a `.seq` file if format is documented.
- If `ProbabilityPlot` outputs `-log10(p)` rather than `p`, detect and convert:
  - if values are mostly > 1.0, treat them as `-log10(p)` and set `p = 10**(-v)`.

## 3) Parsing ProbabilityPlot output

Implement a parser that is robust to:
- header/comment lines
- 1-based indices
- different column layouts (common patterns include `i j p` or `i j -log10(p)`)

Parsing heuristics:
- Accept lines with at least 3 tokens; parse `i`, `j`, `val`.
- Convert to 0-based `(i-1, j-1)` with `(min,max)` ordering.
- Convert `val` → `p`:
  - if `0 <= val <= 1`: interpret as probability
  - else if `val > 1`: interpret as `-log10(p)`
  - else: skip
- Clamp: `p = min(1.0, max(0.0, p))`

Filtering:
- only canonical pairs (reuse `src/lib/sample_cacofold_structures.is_canonical_pair`)
- optionally skip `p < p_min` to keep the weight map sparse

## 4) Converting PF probabilities into a weight component

Recommended:
- `thermo_raw(i,j) = log(p_ij + eps)` with `eps≈1e-6`

Why:
- makes “very likely” pairs meaningfully stand out
- avoids the tiny dynamic range of raw probabilities

Then:
- normalize this component (robust z-score / percentile) before blending.

## 5) Additional RNAstructure ways to use PF (optional)

These are optional follow-ons (not required to ship point 4), but worth considering for the few sampling-limited targets:

### A) `MaxExpect` (MEA scaffold)
- Use PF to generate a maximum expected accuracy structure.
- Add MEA dot-bracket as an extra scaffold seed (diversity) when cov/MSA is weak.

### B) `Stochastic` (PF structure samples)
- Draw a small number of PF samples (e.g., 10–50) and:
  - use them as extra scaffolds, or
  - use them only to estimate a “pair frequency” feature, similar to AllSub but PF-calibrated.

### C) `ProbKnot` (PK hints)
- Can produce crossing pairs from probabilities.
- Treat as a low-weight component (or scaffold-only), since pk evaluation is tricky and pk_f1 is currently ~0.

## 6) Integration points in this repo

Places currently using AllSub counts:
- `src/lib/pipeline.py` inside `sample_pk()` (thermo augmentation of `candidate_pairs`)
- `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py` inside `build_topk_predictions()` (adds AllSub counts to `weights` before ranking)

Suggested interface:
- Add `thermo_mode: {pf, allsub, off}`
- Add `partition_exe` and `probabilityplot_exe` paths (alongside AllSub/Fold/DuplexFold)
- Keep AllSub as fallback if PF binaries are missing or fail.

## 7) Test strategy (no external binaries)

- Add a fixture text file that mimics `ProbabilityPlot` output (both formats: `p` and `-log10(p)`).
- Unit test the parser and conversion logic only.

