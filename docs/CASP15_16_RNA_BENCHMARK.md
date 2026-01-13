# CASP15/16 RNA Benchmark (≤300 nt, single-chain)

Goal: benchmark the pipeline on RNA targets from CASP15 and CASP16, restricted to:
- **single-chain RNAs**
- **≤ 300 nucleotides** (as extracted from the truth structure)

This repo’s benchmark runner (`ssbench`) expects:
- a **manifest CSV** with `target_id`, `pdb_id`, `chain_id` columns at minimum
- a **PDB-chain directory** containing files named `{pdb_id}_{chain_id}.pdb`
- a **truth directory** containing `{target_id}.json` built from those chain PDBs

## Directory layout

Recommended:
- `benchmark_runner/data/casp/pdb_raw/` (your downloaded CASP experimental structures)
- `benchmark_runner/data/casp/pdb_chains/` (chain-extracted PDBs named `{pdb_id}_{chain}.pdb`)
- `benchmark_runner/data/manifests/casp15_16_rna.csv` (manifest for ssbench)
- `benchmark_runner/data/truth_casp15_16/` (truth JSONs)
- `benchmark_runner/data/predictions/casp15_16_rna/` (predictions)
- `benchmark_runner/data/results/metrics_casp15_16_rna.csv` (metrics)

## Step 1: Create a target list (CASP15/16 RNAs)

You need a CSV listing the experimental targets you want to include, with at least:
- `target_id`: a unique identifier (e.g. `CASP16_RNA_01`)
- `pdb_id`: an identifier you control (used only for filenames in this benchmark)
- `chain_id`: chain identifier in the PDB (e.g. `A`)
- `pdb_path`: path to the source PDB file (in `pdb_raw/` or anywhere)

Use this helper to inventory chains + rough residue counts from PDB ATOM records:
`benchmark_runner/scripts/casp_list_pdb_chains.py --pdb-dir benchmark_runner/data/casp/pdb_raw --out benchmark_runner/data/casp/chain_inventory.csv`

Then create `benchmark_runner/data/casp/casp15_16_targets.csv` (your curated list).

## Step 2: Extract single-chain PDBs into `{pdb_id}_{chain}.pdb`

`benchmark_runner/scripts/casp_prepare_chain_pdbs.py --targets benchmark_runner/data/casp/casp15_16_targets.csv --out-dir benchmark_runner/data/casp/pdb_chains`

This writes one chain per target as:
- `benchmark_runner/data/casp/pdb_chains/{pdb_id}_{chain_id}.pdb`

## Step 3: Build the ssbench manifest (and filter to ≤300 nt)

Build a manifest suitable for `ssbench truth build`:
`benchmark_runner/scripts/casp_build_manifest.py --targets benchmark_runner/data/casp/casp15_16_targets.csv --out benchmark_runner/data/manifests/casp15_16_rna.csv`

Filtering policy:
- initial filter: keep only targets whose extracted chain has a plausible nucleotide count (from ATOM records)
- strict filter: after truth build, drop any target whose truth sequence length is `>300` or empty

## Step 4: Build truth from 3D (Barnaba)

`benchmark_runner/.venv/bin/python -m ssbench.cli truth build --manifest benchmark_runner/data/manifests/casp15_16_rna.csv --pdb-dir benchmark_runner/data/casp/pdb_chains --out-dir benchmark_runner/data/truth_casp15_16`

Then run:
`benchmark_runner/scripts/casp_filter_truth_len.py --manifest benchmark_runner/data/manifests/casp15_16_rna.csv --truth-dir benchmark_runner/data/truth_casp15_16 --max-nt 300 --out benchmark_runner/data/manifests/casp15_16_rna_len300.csv`

Use the `*_len300.csv` manifest for benchmarking.

## Step 5: Run predictions + scoring

Use the same `ssbench` pipeline as other benchmarks:
- `ssbench.cli predict run`
- `ssbench.cli score`
- `ssbench.cli report`

Or copy the pattern from `benchmark_runner/scripts/run_selected_rnas_predict_score.sh`.

