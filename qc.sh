#!/bin/bash
# Quality Control script for MCMC sampler refactor
# Run this after EVERY code change

set -e

echo "=== Compile Check ==="
python -m compileall src/ tests/ -q

echo "=== Unit Tests ==="
python -m pytest tests/ -q

echo "=== Lint ==="
ruff check src/ tests/ --quiet || ruff check src/ tests/

echo "=== Format Check ==="
ruff format --check src/ tests/

echo "=== Type Check ==="
mypy src/lib/ --ignore-missing-imports

echo "=== All QC Passed ==="
