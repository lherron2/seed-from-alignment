# Scaffold backend wrappers (Fold / AllSub)

These wrappers emulate RNAstructure-style `Fold` and `AllSub` CLIs so RNAnneal-ss can swap in
alternative single-strand secondary structure backends while **still using RNAstructure `DuplexFold`**
for interaction sampling.

Backends:

- `linearfold_v/`: uses LinearFold-V (`--zuker` for suboptimals) and writes a CT file.
- `eternafold/`: uses EternaFold (CONTRAfold) for MEA + samples, then scores energies with
  RNAstructure `efn2` (fallback energy model when the backend lacks an energy).

Common env vars:

- `RNASTRUCTURE_DOCKER_IMAGE` (for `efn2` scoring): default `rnanneal/rnastructure:v6.4.0`
- `LINEARFOLD_DOCKER_IMAGE`: default `rnanneal/linearfold:latest`
- `ETERNAFOLD_DOCKER_IMAGE`: default `rnanneal/eternafold:latest`

Notes:

- These are intended for **benchmarking and ablations**. They do not reproduce the full behavior of
  RNAstructure `AllSub` (e.g., energy-window enumeration); they generate a bounded list of candidates.

