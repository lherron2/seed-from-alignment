#!/bin/bash
# Quality Control script for MCMC sampler refactor
# Run this after EVERY code change

set -e

# Use venv if available
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "$SCRIPT_DIR/.venv" ]; then
    source "$SCRIPT_DIR/.venv/bin/activate"
fi

echo "=== Compile Check ==="
python3 -m compileall src/ tests/ -q

echo "=== Unit Tests ==="
python3 -m pytest tests/ -q

echo "=== Lint ==="
python3 -m ruff check src/ tests/ --quiet || python3 -m ruff check src/ tests/

echo "=== Format Check ==="
python3 -m ruff format --check src/ tests/

echo "=== Type Check ==="
python3 -m mypy src/lib/ --ignore-missing-imports

echo "=== All QC Passed ==="
