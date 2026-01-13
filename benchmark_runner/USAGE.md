# ssbench usage guide

This directory contains two related things:

1) `ssbench` (a small benchmark CLI) that can build a dataset manifest, fetch PDB chains, build 3D-derived “truth” base pairs, run a predictor, and score results.
2) A predictor implementation (`ssbench.predict.cacofold_mcmc_pipeline`) that takes a **single RNA sequence** and outputs a ranked list of secondary structures.

The glue for most “end-to-end” runs is `benchmark_runner/scripts/run_benchmark.sh`, which wires the CLI + predictor together using `benchmark_runner/defaults.yaml`.

## Installation / invoking the CLI

You can run `ssbench` without installing it by setting `PYTHONPATH` to `benchmark_runner/src`:

```bash
cd benchmark_runner
PYTHONPATH=src python -m ssbench.cli --help
```

If you prefer a `ssbench` executable in your environment, install it from this directory:

```bash
cd benchmark_runner
python -m pip install -e .
ssbench --help
```

## Key concepts / file formats

### Targets and `target_id`

- A “target” is one RNA chain to predict.
- In BGSU/FR3D-style naming, `target_id` often looks like `6XA1|1|Bv`.
- The CLI writes per-target files by replacing `|` with `_` when creating directories.

### Manifests (CSV)

Many commands take `--manifest path/to.csv`. The manifest must have at least:

- `target_id`
- `pdb_id`
- `chain_id`

Other columns (e.g. `rfam`, `nts_observed`, `split`) are optional and used for reporting.

### Truth records (JSON)

Truth JSONs live in a directory like `benchmark_runner/data/truth/...` and are named:

- `{target_id}.json`

Each contains:

- `sequence`: the extracted RNA sequence (canonicalized to A/C/G/U where possible)
- `canonical_pairs`: list of canonical base pairs `(i, j)` (0-based indices)
- `nested_pairs`, `pk_pairs`, `helices`: derived annotations

### Predictions (directory per target)

Predictions live under an output directory like `benchmark_runner/data/predictions/.../{safe_id}/`.

The main expected output is:

- `predictions.db`

`predictions.db` format:

1) First line: sequence (same length as truth)
2) Next lines: dot-bracket structures (same length as sequence)

Pseudoknot bracket alphabets like `[]{}<>` (and `a/A`, `b/B`, …) are supported by the scorer.

## “One sequence” vs “multi sequence”

### One sequence (prediction only; no benchmark/truth)

If you only have an RNA sequence and want structures, run the predictor directly.

Create a single-sequence FASTA:

```bash
cd benchmark_runner
cat > my.fa <<'EOF'
>my_rna
ACGUACGUACGU
EOF
```

Run the predictor (edit paths/flags to match your environment):

```bash
DATAPATH=/path/to/RNAstructure/data_tables \
PYTHONPATH=src .venv/bin/python -m ssbench.predict.cacofold_mcmc_pipeline \
  --cm-db /path/to/Rfam.cm \
  --rscape /path/to/R-scape \
  --infernal-bin /path/to/infernal/bin \
  --allsub-exe /path/to/RNAstructure/AllSub \
  --duplex-exe /path/to/RNAstructure/DuplexFold \
  --fold-exe /path/to/RNAstructure/Fold \
  --top-k 100 \
  my.fa out/my_rna
```

Outputs will appear in `out/my_rna/`, including `out/my_rna/predictions.db`.

Important: `ssbench.predict.cacofold_mcmc_pipeline` expects a **single-sequence FASTA**. If you pass a multi-FASTA file, it will not parse correctly.

### Multi sequence (multiple targets)

There are two common “multi target” workflows:

1) You already have a **manifest CSV** with many targets → run `ssbench` over the manifest.
2) You have a **multi-FASTA** where each entry corresponds to a target with a known PDB → convert the FASTA to a manifest, then run.

#### Multi-target via manifest + truth (recommended)

Given:

- `MANIFEST.csv` with `target_id,pdb_id,chain_id,...`
- `TRUTH_DIR/` with `{target_id}.json`

Run the predictor across all targets:

```bash
cd benchmark_runner
PYTHONPATH=src .venv/bin/python -m ssbench.cli predict run \
  --manifest MANIFEST.csv \
  --truth TRUTH_DIR \
  --predictor-cmd "PYTHONPATH=src .venv/bin/python -m ssbench.predict.cacofold_mcmc_pipeline ... {fasta} {outdir}" \
  --k 100 \
  --out-dir data/predictions/my_run \
  --config defaults.yaml
```

Notes:

- `predict run` writes one `{safe_id}.fa` per target and calls your `--predictor-cmd` once per target.
- The predictor command is a **template**; it must accept `{fasta}` and `{outdir}` placeholders.
- If `predictions.db` already exists with at least `k` structures, `predict run` skips that target (resume-friendly).

#### Multi-FASTA → manifest (for selected hand-picked targets)

If you have a multi-FASTA where headers include a PDB id, you can build a manifest using:

`benchmark_runner/scripts/manifest_from_fasta.py`

Expected header format:

- `>TARGET_ID|PDB_ID|...`

Example:

```bash
cd benchmark_runner
PYTHONPATH=src .venv/bin/python scripts/manifest_from_fasta.py \
  --fasta ../selected_rnas.fasta \
  --out-manifest data/manifests/selected.csv \
  --pdb-dir data/pdb \
  --download \
  --allow-unmatched
```

This matches each FASTA entry to the best PDB chain by sequence identity, and writes a manifest that downstream steps can use.

## Running a full benchmark end-to-end (recommended script)

The easiest “do everything” entry point is:

- `benchmark_runner/scripts/run_benchmark.sh`

It performs:

1) Build manifest from a (possibly multi-)FASTA (`manifest_from_fasta.py`)
2) Fetch PDB chains (`ssbench.cli dataset fetch-pdb`)
3) Build truth from 3D (`ssbench.cli truth build`)
4) Run prediction (`ssbench.cli predict run`)
5) Score and summarize (`ssbench.cli score`, `ssbench.cli report`)

Minimal usage:

```bash
cd benchmark_runner
./scripts/run_benchmark.sh
```

By default it looks for:

- `FASTA_PATH=../selected_rnas.fasta`
- `DEFAULTS_FILE=benchmark_runner/defaults.yaml`

Override those paths:

```bash
cd benchmark_runner
FASTA_PATH=/path/to/targets.fasta \
DEFAULTS_FILE=/path/to/my_defaults.yaml \
LOG_FILE=data/results/my_run.log \
./scripts/run_benchmark.sh
```

Outputs (by default) go under `benchmark_runner/data/`:

- `data/manifests/*`
- `data/pdb/*`
- `data/truth_selected/*`
- `data/predictions/selected_rnas_benchmark/*`
- `data/results/metrics_selected_rnas.csv`
- `data/results/summary_selected_rnas.md`

## Configuration (`defaults.yaml`)

### What to edit

`benchmark_runner/defaults.yaml` has four main sections:

- `dataset`: dataset download + length buckets (used by `ssbench cli dataset download` and some scripts)
- `score`: scoring configuration (e.g. `helix_overlap_threshold`, default `k`)
- `paths`: locations of external tools and databases (Infernal/R-scape/RNAstructure)
- `predict`: predictor hyperparameters (top-k, MCMC budgets, ranking settings, etc)

### How to change config cleanly

1) Copy the defaults:

```bash
cp benchmark_runner/defaults.yaml benchmark_runner/my_defaults.yaml
```

2) Edit `benchmark_runner/my_defaults.yaml`.
3) Point scripts/CLI at it:

- For `run_benchmark.sh` / `run_manifest_predict_parallel.sh`: set `DEFAULTS_FILE=...`
- For `ssbench` subcommands: pass `--config my_defaults.yaml`

### What reads what

- `ssbench.cli` uses `--config` mainly for:
  - dataset download parameters (release, resolution, length buckets, split seed)
  - scoring parameters (helix overlap threshold; default `k`)
- The predictor (`ssbench.predict.cacofold_mcmc_pipeline`) does **not** read `defaults.yaml` directly; it only reads its own CLI flags.
- The provided shell runners (`run_benchmark.sh`, `run_manifest_predict_parallel.sh`) read `DEFAULTS_FILE` and turn YAML keys into predictor CLI flags automatically.

## Common pitfalls / troubleshooting

- **External dependencies**: the CaCoFold-based predictors require Infernal (`cmscan/cmfetch/cmalign`), R-scape, and RNAstructure executables; set them in `defaults.yaml` under `paths:`.
- **RNAstructure thermodynamic tables**: set `DATAPATH` (the runners export it) to avoid PF tools failing.
- **Modified residues / unknown bases**: truth building canonicalizes many modified residues; prediction input is sanitized to a strict `ACGU` alphabet when necessary.
- **rRNA exclusion**: if you use `benchmark_runner/scripts/bgsu_make_structured_manifest.py`, it can exclude ribosomal RNAs via `cmscan` (requires `Rfam.cm`). Disable with `--no-exclude-ribosomal` if you don’t have Rfam available.

## Parallelizing prediction across a manifest

For large manifests, use:

- `benchmark_runner/scripts/run_manifest_predict_parallel.sh`

It runs one target per process (via `xargs -P`) and relies on `predict run` being resume-friendly. Minimal usage:

```bash
cd benchmark_runner
MANIFEST=data/manifests/MANIFEST.csv \
TRUTH_DIR=data/truth/TRUTH_DIR \
PRED_DIR=data/predictions/my_parallel_run \
JOBS=8 \
DEFAULTS_FILE=defaults.yaml \
./scripts/run_manifest_predict_parallel.sh
```

