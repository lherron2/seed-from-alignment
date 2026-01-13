# Secondary Structure Prediction Approach

This repository implements a multi-stage pipeline for generating ensembles of RNA secondary structures, including pseudoknots, from comparative sequence alignments. The approach combines consensus structure information from CaCoFold/R-scape with thermodynamic refinement and stochastic sampling.

## Pipeline Overview

```
Stockholm Alignment (.sto) + Covariation Stats (.cov)
                    │
                    ▼
        ┌──────────────────────┐
        │  Consensus Extraction │  (00_consensus.db)
        └──────────────────────┘
                    │
                    ▼
        ┌──────────────────────┐
        │  Thermodynamic        │  (01_refined.db)
        │  Refinement           │
        └──────────────────────┘
                    │
                    ▼
        ┌──────────────────────┐
        │  Pseudoknot Sampling  │  (02_sampled.db)
        │  (Metropolis MCMC)    │
        └──────────────────────┘
                    │
                    ▼
        ┌──────────────────────┐
        │  Rosetta Conversion   │  (04_rosetta.db)
        │  & Filtering          │
        └──────────────────────┘
```

## Stage 1: Consensus Structure Extraction

**Source:** `src/lib/sample_cacofold_structures.py`

The pipeline begins by parsing a Stockholm alignment file containing:
- Multiple aligned RNA sequences
- Consensus secondary structure annotations (`SS_cons`, `SS_cons_1`, etc.) in WUSS notation

Key operations:
1. **Parse Stockholm format** - Extract sequences and all `#=GC SS_cons*` structure tracks
2. **Project to target sequence** - Map alignment-coordinate base pairs to ungapped sequence coordinates
3. **Filter loop separation** - Remove base pairs with hairpin loops shorter than 2 nucleotides
4. **Generate candidate pairs** - Collect all annotated base pairs across SS_cons tracks, weighted by covariation evidence

## Stage 2: Thermodynamic Refinement

**Source:** `src/lib/refine_unpaired_regions.py`

This stage fills in structure for regions left unpaired by the consensus, using RNAstructure's thermodynamic algorithms.

### 2.1 Consensus Masking

Before refinement, the pipeline generates scaffold variants by systematically "masking" (removing) parts of the consensus structure:

- **End Masking Grid**: Incrementally unpair 5' and 3' terminal regions (step size: 5 nt, max: 40 nt)
- **Sequential Helix Masking**: Remove individual helices one at a time (top N by length)
- **Pairwise Helix Masking**: Remove pairs of helices together

This forces exploration of alternative folds in regions that might be incorrectly locked by the consensus.

### 2.2 Terminal Seeding (5'-3' Interactions)

For each masked scaffold, the pipeline explicitly searches for global closing stems:

- Use RNAstructure **DuplexFold** to identify potential 5'-3' terminal interactions
- Dense windowing across unpaired terminal regions
- Seeds are ranked by interaction size and added to scaffolds

### 2.3 Combinatorial Region Refinement

For each unpaired region exceeding a minimum length threshold:

1. **Local Folding** (RNAstructure AllSub):
   - Run suboptimal enumeration on the unpaired subsequence
   - Parameters: temperature (310.15 K), energy window (absolute or percent)
   - Returns multiple alternative local structures with energies

2. **Loop-Loop Interactions** (DuplexFold):
   - Search for inter-loop base pairing (kissing loops, pseudoknots)
   - Compare all pairs of unpaired regions for potential interactions

3. **Combinatorial Assembly**:
   - Enumerate all compatible combinations of local folds
   - Track indices to prevent overlapping base pairs

### 2.4 Topology-Aware Scoring

Final structures are scored using a layer-based weighting scheme:

```
Score = Σ Weight(Layer(p)) × Energy(p)

Where:
  Layer 0 (Nested):     Weight = 1.0   (full stability credit)
  Layer 1 (PK depth 1): Weight = 0.6   (reduced stability)
  Layer 2+ (deep PK):   Weight = -1.0  (penalty for complexity)
```

Layers are assigned via greedy graph coloring - pairs are placed in the lowest layer that has no crossing arcs.

## Stage 3: Pseudoknot Sampling

**Source:** `src/lib/sample_cacofold_structures.py`

This stage uses Metropolis-Hastings MCMC to sample alternative pseudoknotted structures:

### Candidate Pair Generation
- Collect all base pairs from SS_cons tracks
- Filter for canonical pairing (Watson-Crick + GU wobble)
- Optionally weight by covariation statistics (R-scape E-values/power)

### MCMC Sampling
- **State**: Set of active base pairs (valid matching)
- **Moves**: Add/remove/swap pairs
- **Energy**: Sum of pair weights with pseudoknot penalty (pk_alpha parameter)
- **Annealing**: Inverse temperature (beta) controls exploration

### Complexity Filtering
Samples are filtered to remove excessive pseudoknot complexity:
- Maximum crossing fraction (default: 20% of pairs)
- Maximum crossings per pair (default: 2)
- Maximum total crossings (default: 30)
- Optional depth limit

### Scaffold Integration
When a refined scaffold exists, sampled pairs are constrained to avoid conflicting with locked scaffold positions.

## Stage 4: Rosetta Conversion

**Source:** `src/lib/filter_db_for_rosetta.py`

Final processing for compatibility with Rosetta RNA modeling:

1. **Layer Assignment**: Partition pairs into non-crossing layers
2. **Bracket Notation**: Assign bracket types by layer
   - Layer 0: `( )`
   - Layer 1: `[ ]`
   - Layer 2: `{ }`
   - Layer 3+: `a-z / A-Z`

3. **Canonical Filtering**: Remove non-Watson-Crick/wobble pairs
4. **Deduplication**: Merge identical or near-identical structures (Hamming distance)

## Key Configuration Parameters

| Category | Parameter | Default | Description |
|----------|-----------|---------|-------------|
| **Sampling** | n_samples | 1000 | Number of MCMC samples |
| | burn_in | 1000 | Burn-in iterations |
| | beta | 1.0 | Inverse temperature |
| | pk_alpha | 2.0 | Pseudoknot penalty |
| **Refinement** | min_unpaired | 15 | Min region length to refine |
| | temperature | 310.15 | Folding temperature (K) |
| | max_structures | 1000 | Output limit |
| | max_seeds | 50 | 5'-3' seeds per scaffold |
| **Filtering** | max_total_cross | 30 | Max PK crossings |
| | max_cross_per_pair | 2 | Max crossings any pair participates in |

## Output Format

The `.db` (dot-bracket) format:
```
SEQUENCE
(((....)))...[[[[...]]]]
(((....[[[[..))).]]]]...
```

- Line 1: Ungapped RNA sequence
- Lines 2+: One structure per line in pseudoknot-aware bracket notation

## Dependencies

- **RNAstructure**: AllSub (suboptimal enumeration), DuplexFold (intermolecular)
- **CaCoFold/R-scape**: Input Stockholm files with covariation analysis
- **Python 3.10+**: Type hints, dataclasses, pathlib
