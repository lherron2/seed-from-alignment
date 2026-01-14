# PKOPT3: targeted pseudoknot exploration in unfixed sampling (ribo-benchmark)

This note documents a small pipeline variant aimed at improving pseudoknot (PK) coverage without
regressing the strong opt6 fixed-scaffold behavior.

## Change summary

PKOPT3 keeps the opt6 fixed-scaffold sampling and top-K selection behavior, but makes the optional
**unfixed** sampling pass (the one with `fix_scaffold_pairs=False`, enabled for `len>=80`) more
PK-friendly:

- `benchmark_runner/src/ssbench/predict/cacofold_mcmc_pipeline.py`: when running the unfixed pass:
  - If CaCoFold `.cov` exists but contains no usable pairs (`cov_empty`), increase
    `thermo_max_structures` to at least `200` and set `thermo_min_count=1` for candidate-pair
    augmentation (sampling-only).
  - Use the CLI scorer settings for sampling weights (`robust_z` normalization + component alphas).
  - Relax PK filters for the unfixed sampler (`pk_filter_frac=0.35`,
    `pk_filter_max_cross_per_pair≈0.20*L` capped at 30, `pk_filter_max_total_cross≈1.50*L` capped at
    500).
  - Reduce the unfixed sampler PK penalty (`pk_alpha *= 0.5`) to increase exploration.
- `src/lib/pipeline.py`: refined-scaffold downselection stays unchanged for the fixed pass, but uses
  a small PK-aware scaffold quota for the unfixed pass (helps seed PK regimes when sampling is
  allowed to modify scaffolds).

## Results (ribo-benchmark oracle best-of-K)

Evaluation target: `rna_sets/ribo-benchmark`, oracle best-of-K over `K ∈ {1,50,100,200,500}`.

Mean oracle F1/MCC (N=16):

- `@500` RN (pkopt3) `0.867299/0.868501` vs RN (opt6) `0.866600/0.868140` (ΔF1 `+0.000699`)
- Other `K` means are unchanged at 3 decimal places vs opt6.

Per-RNA changes at `@500` (RN pkopt3 vs RN opt6):

- `4YAZ` improves from `0.780488` → `0.791667` (Δ `+0.011179`)
- All other RNAs are unchanged within float equality under this scorer script.

The full per-RNA breakdown is written to:

- `out/ribo_benchmark_rnanneal_pkopt3_metrics.csv`

## Reproduce

1) Generate predictions:

```bash
export DATAPATH="$(pwd)/RNAstructure/data_tables"
out_root="out/ribo_benchmark_rnanneal_ss_top500_pkopt3"
mkdir -p "$out_root"
for d in rna_sets/ribo-benchmark/*; do
  [ -d "$d" ] || continue
  r=$(basename "$d")
  rnanneal-ss "$d/$r.fasta" "$out_root/$r" --top-k 500 \
    --reuse-cacofold-root out/ribo_benchmark_rnanneal_ss_top500_opt7
done
```

2) Score oracle best-of-K:

- Use `docs/RIBO_BENCHMARK_ORACLE_TOPK.md` and point the RN root at
  `out/ribo_benchmark_rnanneal_ss_top500_pkopt3`.

## Next PK-focused ideas

The remaining hard cases (e.g. `6WJR`, `3Q3Z`, `6VMY`, `6Q57`) appear to need more than “more
sampling”; candidates worth trying next:

- Short-stem injection: explicitly propose and weight helix segments down to length 2 as candidate
  pairs (e.g. `6WJR` has a 2-bp PK helix), instead of relying on AllSub coverage.
- Dedicated PK moves: add a move that adds/removes an entire helix segment (or a crossing helix
  pair) in one proposal to cross energy barriers.
- Dual-penalty sampling: run an extra chain/pass with `pk_alpha=0` (or very low) and keep a small
  quota of its outputs in the final top-K via the existing PK bucket.
- Thermo PF mode: optionally use RNAstructure `Partition`/`ProbabilityPlot` to include low-probability
  long-range pairs that AllSub may miss.

