# Algorithm And Hyperparameters

This document summarizes the CaCoFold + domain masking + MCMC pipeline and the
current tunable hyperparameters exposed in the CLI.

## Overall Algorithm

1) Input: a target sequence in FASTA.
2) Rfam search (Infernal):
   - `cmscan` against the local Rfam CM database.
   - If a top hit exists, `cmfetch` pulls the model and `cmalign` builds a
     Stockholm alignment.
3) CaCoFold (R-scape):
   - Run `R-scape --cacofold` to produce a CaCoFold Stockholm and optional
     covariation (.cov) file.
   - If CaCoFold fails or drops residues, fall back to MFE-only scaffold.
4) MFE fallback (RNAstructure Fold):
   - Build an MFE dot-bracket and synthesize a minimal Stockholm with SS_cons.
5) Domain masking + refinement (legacy pipeline):
   - Generate masked variants of the scaffold.
   - Refine each mask with RNAstructure AllSub and DuplexFold.
6) MCMC sampling (legacy sampler):
   - Sample alternative structures guided by candidate pairs and covariation
     weights (if available).
7) Ranking and output:
   - Score structures using the legacy energy with covariation weights.
   - Output top-K dot-brackets to `predictions.db`.

## Hyperparameters (CLI Flags)

### Core MCMC
- `--n-samples` (default 200): total MCMC samples.
- `--burn-in` (default 500): burn-in steps.
- `--thin` (default 10): thinning interval.
- `--beta` (default 1.0): inverse temperature for acceptance.
- `--seed` (default 0): RNG seed.
- `--min-loop-sep` (default 3): minimum loop separation for pairing moves.

### Covariation Weighting
- `--cov-mode` (default "logE_power"): covariation weighting mode.
- `--cov-alpha` (default 3.0): covariation weight scale.
- `--cov-min-power` (default 0.1): minimum covariation power for a pair to
  contribute.
- `--cov-forbid-negative` (default False): drop pairs with strong negative
  evidence.

### Refinement Caps (runtime control)
- `--refine-max-structures` (default 200): max refined structures written.
- `--refine-max-regions` (default 10): max regions to refine per structure.
- `--refine-max-seeds` (default 20): max end-to-end seeds per mask.
- `--refine-max-solutions` (default 200): cap solutions per refinement track.
- `--refine-kissing-candidates` (default 200): cap kissing-loop candidates.

### I/O and Tooling
- `--cm-db`: path to `Rfam.cm` (required).
- `--infernal-bin`: path to Infernal bin directory (optional if on PATH).
- `--rscape`: path to R-scape binary (optional if on PATH).
- `--allsub-exe`: path to RNAstructure AllSub (required).
- `--duplex-exe`: path to RNAstructure DuplexFold (optional but recommended).
- `--fold-exe`: path to RNAstructure Fold (required for MFE fallback).
- `--top-k` (default 20): number of structures emitted to `predictions.db`.

## Notes

- If no Rfam hit is found, the pipeline uses MFE as the initial scaffold and
  proceeds with masking + MCMC sampling.
- If CaCoFold output length does not match the input FASTA, the pipeline falls
  back to the MFE scaffold to avoid indexing errors.
