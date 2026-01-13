# ssbench

Reproducible benchmark runner for RNA secondary structure prediction using 3D-derived ground truth.

See `benchmark_runner/USAGE.md` for a practical end-to-end guide (single-sequence prediction vs multi-target benchmarking, config editing, and provided runner scripts).

## Quickstart

This example downloads a small representative set, fetches PDBs, builds Barnaba truth, runs a dummy predictor, and scores results.

```bash
cd benchmark_runner
PYTHONPATH=src python -m ssbench.cli --help

PYTHONPATH=src python -m ssbench.cli dataset download --release-id current --resolution 4.0 --limit 25 --out-manifest data/manifests/rep.csv
PYTHONPATH=src python -m ssbench.cli dataset fetch-pdb --manifest data/manifests/rep.csv --out-dir data/pdb

# Build truth (requires barnaba installed)
PYTHONPATH=src python -m ssbench.cli truth build --manifest data/manifests/rep.csv --pdb-dir data/pdb --out-dir data/truth

# Example predictor command template:
# It receives {fasta} and {outdir}. Write K dot-bracket structures to outdir.
PYTHONPATH=src python -m ssbench.cli predict run --manifest data/manifests/rep.csv --truth data/truth \
  --predictor-cmd "python scripts/echo_predictor.py {fasta} {outdir}" --k 1 --out-dir data/predictions

PYTHONPATH=src python -m ssbench.cli score --manifest data/manifests/rep.csv --truth data/truth --predictions data/predictions --out data/results/metrics.csv
PYTHONPATH=src python -m ssbench.cli report --metrics data/results/metrics.csv --out data/results/summary.md
```

## Notes

- Ground truth uses Barnaba (`bb.annotate` + `bb.dot_bracket`) to extract canonical pairs.
- No network calls are made during unit tests; downloads are explicit `ssbench dataset download` and `ssbench dataset fetch-pdb` steps.
- Single-chain IFEs are selected by default (`ife_id` without '+').


## Using the legacy predictor

The legacy predictor in this repo is `sample_cacofold_structures.py`, which requires a CaCoFold `.sto` and (optionally) `.cov` plus the sequence name inside that alignment. For benchmarking against 3D truth, you must provide a mapping from each target to its CaCoFold inputs.

1) Create a mapping CSV with columns:

```
target_id,sto_path,cov_path,seq_name
9E5I|1|A,/path/to/combined_1.cacofold.sto,/path/to/combined_1.cacofold.cov,9E5I_1_Chains
```

2) Use the legacy wrapper as the predictor command:

```
PYTHONPATH=src python -m ssbench.cli predict run \
  --manifest data/manifests/rep_fetched.csv \
  --truth data/truth \
  --predictor-cmd "python3 scripts/legacy_predictor.py --mapping data/legacy_map.csv {fasta} {outdir}" \
  --k 20 --out-dir data/predictions
```

This will generate `.db` files using the legacy sampler and should produce non-zero scores when the CaCoFold inputs are available for the benchmark targets.

## CaCoFold pipeline predictor

To use Infernal + R-scape CaCoFold as the predictor:

```
PYTHONPATH=src python -m ssbench.cli predict run \
  --manifest data/manifests/rep_fetched.csv \
  --truth data/truth \
  --predictor-cmd "python3 -m ssbench.predict.cacofold_pipeline --cm-db /path/to/Rfam.cm --rscape /path/to/rscape {fasta} {outdir}" \
  --k 1 --out-dir data/predictions
```

Requirements:
- Infernal (`cmscan`, `cmfetch`, `cmalign`) in PATH or `--infernal-bin`.
- R-scape binary (`rscape`) in PATH or `--rscape`.
- Local Rfam CM database file.

## CaCoFold + domain masking + MCMC predictor

This pipeline uses Infernal + R-scape CaCoFold, then applies legacy domain masking and MCMC sampling from `src/lib/pipeline.py`, and outputs the top 20 structures by energy score.

```
PYTHONPATH=src python -m ssbench.cli predict run \
  --manifest data/manifests/rep_fetched.csv \
  --truth data/truth \
  --predictor-cmd "python3 -m ssbench.predict.cacofold_mcmc_pipeline \
    --cm-db /path/to/Rfam.cm \
    --rscape /path/to/rscape \
    --infernal-bin /path/to/infernal/bin \
    --allsub-exe /path/to/RNAstructure/AllSub \
    --duplex-exe /path/to/RNAstructure/DuplexFold \
    --fold-exe /path/to/RNAstructure/Fold \
    --top-k 20 \
    {fasta} {outdir}" \
  --k 20 --out-dir data/predictions
```

Notes:
- `AllSub` is required for domain masking (legacy refinement).
- `Fold` is used as a fallback to produce an MFE structure when no coevolution data is found.
- The pipeline scores structures using the legacy energy function and writes the top 20 to `predictions.db`.
