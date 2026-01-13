# RNAnneal-ss

Primary entrypoint (recommended) is the **experiment-proven top-K pipeline**:

- Infernal (`cmscan/cmfetch/cmalign`) + R-scape CaCoFold to build an alignment + covariation tracks
- Ensemble scaffold backends (RNAstructure + LinearFold + EternaFold)
- RNAnneal-ss MCMC sampling
- Writes a **single ranked list** of structures so you can score oracle best-of-`K` by truncating prefixes

## Install (fresh clone)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## Run (single FASTA → `predictions.db`)

```bash
rnanneal-ss rna_sets/ribo-benchmark/6WJR/6WJR.fasta out/6WJR --top-k 500
```

Output: `out/6WJR/predictions.db` (first line is the RNA sequence; remaining lines are dot-brackets).

Notes:
- Uses repo-default tool paths: `benchmark_runner/data/rfam/Rfam.cm`, `infernal-1.1.5/src`, `rscape_v2.6.4/src/R-scape`,
  and scaffold backends in `benchmark_runner/tools/scaffold_backends/ensemble3/`.
- If RNAstructure complains about thermodynamic tables, set `DATAPATH=$(pwd)/RNAstructure/data_tables`.

## Reproducing the benchmark experiments

- u300 representative50 top-K test (v1): `benchmark_runner/scripts/run_u300_top_k_test_v1.sh`
- u400 ensemble vs baselines (v1): `benchmark_runner/scripts/run_u400_topk500_v1.sh`

