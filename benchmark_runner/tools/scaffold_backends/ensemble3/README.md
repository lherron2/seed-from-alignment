# ensemble3 scaffold backend (RS + LF + EF)

Drop-in `Fold` / `AllSub` replacements that aggregate candidate structures from:

- **RNAstructure** (native `Fold` / `AllSub`)
- **LinearFold-V** (wrapper `linearfold_v/Fold` / `linearfold_v/AllSub`)
- **EternaFold** (wrapper `eternafold/Fold` / `eternafold/AllSub`)

These scripts are intended to be passed to `ssbench.predict.cacofold_mcmc_pipeline` via:

- `--fold-exe benchmark_runner/tools/scaffold_backends/ensemble3/Fold`
- `--allsub-exe benchmark_runner/tools/scaffold_backends/ensemble3/AllSub`

Controls (env vars)

- `RNANNEAL_SS_SUBOPT_MAX`: max number of structures to emit from `AllSub` (default: `200`).
- `ENSEMBLE3_RS_FOLD`: path to RNAstructure `Fold` (default: `RNAstructure/exe/Fold`).
- `ENSEMBLE3_RS_ALLSUB`: path to RNAstructure `AllSub` (default: `RNAstructure/exe/AllSub`).
- `ENSEMBLE3_LF_FOLD`: path to LinearFold-V `Fold` wrapper (default: `benchmark_runner/tools/scaffold_backends/linearfold_v/Fold`).
- `ENSEMBLE3_LF_ALLSUB`: path to LinearFold-V `AllSub` wrapper (default: `benchmark_runner/tools/scaffold_backends/linearfold_v/AllSub`).
- `ENSEMBLE3_EF_FOLD`: path to EternaFold `Fold` wrapper (default: `benchmark_runner/tools/scaffold_backends/eternafold/Fold`).
- `ENSEMBLE3_EF_ALLSUB`: path to EternaFold `AllSub` wrapper (default: `benchmark_runner/tools/scaffold_backends/eternafold/AllSub`).

Notes

- `AllSub` preserves **diversity across sources** by interleaving candidates rather than concatenating.
- Energies come from each source when available; they should be treated as **source-local** (not cross-comparable).

