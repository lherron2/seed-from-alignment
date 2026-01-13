# Repository Guidelines

## Project Structure & Module Organization

- `src/lib/`: core Python library (`rnanneal`) implementing RNAnneal-ss for sampling/scoring RNA secondary structures (incl. pseudoknots).
- `tests/`: unit + regression tests (pytest).
- `benchmark_runner/`: `ssbench` benchmark runner (dataset/truth extraction, predictor wiring, scoring, reports).
- `docs/`: design notes, benchmark plans, and research writeups.
- `example/`: example configs and pipeline outputs.
- `RNAstructure/`, `infernal-*/`, `rscape_*/`: vendored external tool sources (treat as third-party; avoid local edits unless required).

## Build, Test, and Development Commands

Create an environment and install dev deps:

- `python -m venv .venv && source .venv/bin/activate`
- `python -m pip install -e '.[dev]'`

Run the repository QC suite:

- `./qc.sh` (compile check, `pytest`, `ruff` lint/format check, `mypy`)

Common tasks:

- `python -m pytest -q` (unit tests)
- `python -m ruff check src tests` (lint)
- `python -m ruff format src tests` (format)
- `run-pipeline --help` (main pipeline CLI; see `[project.scripts]` in `pyproject.toml`)

Benchmark runner quickstart (from `benchmark_runner/`):

- `PYTHONPATH=src python -m ssbench.cli --help`

## Coding Style & Naming Conventions

- Python 3.10+, 4-space indentation, max line length 100 (Ruff).
- Prefer `pathlib.Path` over string paths and explicit names like `target_id`, `dotbrackets`.
- Keep external tool paths/config out of code; prefer YAML (e.g., `benchmark_runner/defaults.yaml`, `example/pipeline_config.yaml`).

## Testing Guidelines

- Framework: `pytest` (tests in `tests/test_*.py`).
- Add focused unit tests for new logic; for external binaries, use skip/contract-style tests that only run when the tool is available.

## Commit & Pull Request Guidelines

- Prefer Conventional Commits-style prefixes when possible (`feat:`, `fix:`, `docs:`, `build:`); keep messages imperative and specific.
- PRs: include a short description, the validation command(s) you ran (e.g., `./qc.sh`), and any new environment/config requirements (e.g., `DATAPATH` for RNAstructure tables).

## Security & Configuration Tips

- Do not commit large downloaded datasets, PDBs, or generated benchmark outputs unless explicitly intended.
- If using RNAstructure executables, ensure `DATAPATH` points to `RNAstructure/data_tables` so thermodynamic tables are discoverable.
