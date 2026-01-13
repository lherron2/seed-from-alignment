#!/usr/bin/env python3
"""
Batch runner for the RNA pipeline.

Loops over a list of RNA sequence names, generates a specific configuration
file for each, and runs the pipeline, saving outputs to separate directories.

Usage:
    python scripts/run_batch.py \
        --config example/pipeline_config.yaml \
        --rnas RNA1 RNA2 RNA3 \
        --output-base example/batch_results

    OR using a file list:
    python scripts/run_batch.py \
        --config example/pipeline_config.yaml \
        --rnas-file my_targets.txt \
        --output-base example/batch_results
"""

import argparse
import copy
import sys
import subprocess
from pathlib import Path
import yaml  # Requires: pip install pyyaml

def resolve_path_in_config(cfg: dict, key: str, base_dir: Path):
    """
    Helper to convert relative paths in the config to absolute paths.
    This ensures that when we save the new config file in a subdirectory,
    it still points to the correct input files.
    """
    if key in cfg and cfg[key]:
        p = Path(cfg[key])
        if not p.is_absolute():
            # Resolve relative to the original config file's directory
            abs_path = (base_dir / p).resolve()
            cfg[key] = str(abs_path)

def main():
    parser = argparse.ArgumentParser(description="Run RNA pipeline on a batch of sequences.")
    
    # Input selection
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--rnas", nargs="+", help="Space-separated list of RNA names (IDs).")
    group.add_argument("--rnas-file", help="Path to a text file containing one RNA name per line.")
    
    # Configuration
    parser.add_argument("-c", "--config", required=True, help="Path to the base pipeline_config.yaml template.")
    parser.add_argument("-o", "--output-base", required=True, help="Root directory where per-RNA output folders will be created.")
    
    # Optional overrides
    parser.add_argument("--pipeline-script", default="run_pipeline.py", help="Path to the main pipeline script (default: run_pipeline.py).")
    parser.add_argument("--mode", default="refine_first", choices=["refine_first", "legacy"], help="Pipeline mode to run.")
    
    args = parser.parse_args()

    # 1. Parse RNA list
    rnas = []
    if args.rnas:
        rnas = args.rnas
    elif args.rnas_file:
        with open(args.rnas_file, 'r') as f:
            rnas = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not rnas:
        print("No RNAs found to process.")
        sys.exit(1)

    # 2. Load Base Config
    config_path = Path(args.config).resolve()
    with open(config_path, 'r') as f:
        try:
            base_config = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Error parsing YAML config: {exc}")
            sys.exit(1)

    # 3. Resolve Input Paths to Absolute
    # Because we will write new config files into subdirectories, relative paths
    # like "../inputs/seq.sto" in the original config might break. 
    # We resolve them to absolute paths based on the location of the template config.
    base_config_dir = config_path.parent
    
    resolve_path_in_config(base_config, 'sto', base_config_dir)
    resolve_path_in_config(base_config, 'cov', base_config_dir)
    resolve_path_in_config(base_config, 'allsub_exe', base_config_dir)
    resolve_path_in_config(base_config, 'duplex_exe', base_config_dir)

    # 4. Loop and Process
    output_root = Path(args.output_base).resolve()
    pipeline_script = Path(args.pipeline_script).resolve()

    if not pipeline_script.exists():
        print(f"Error: Pipeline script not found at {pipeline_script}")
        sys.exit(1)

    print(f"Found {len(rnas)} RNAs to process.")
    print(f"Output root: {output_root}")
    print("-" * 60)

    failures = []

    for i, rna in enumerate(rnas, 1):
        print(f"[{i}/{len(rnas)}] Processing sequence: {rna}")
        
        # Create directory: output_base/RNA_NAME/
        rna_dir = output_root / rna
        rna_dir.mkdir(parents=True, exist_ok=True)

        # Clone and modify config
        rna_config = copy.deepcopy(base_config)
        rna_config['seq_name'] = rna
        
        # Update IO section to point to specific dir
        if 'io' not in rna_config:
            rna_config['io'] = {}
        
        # We explicitly set work_dir to the new folder
        rna_config['io']['work_dir'] = str(rna_dir)

        # Write the temporary config file for this RNA
        temp_config_path = rna_dir / f"{rna}_config.yaml"
        with open(temp_config_path, 'w') as f:
            yaml.dump(rna_config, f)

        # Construct command
        cmd = [
            sys.executable, 
            str(pipeline_script),
            "--config", str(temp_config_path),
            "--mode", args.mode
        ]

        # Execute
        try:
            # Check=True will raise CalledProcessError if script fails
            subprocess.run(cmd, check=True)
            print(f"  -> Success. Results in {rna_dir}")
        except subprocess.CalledProcessError as e:
            print(f"  -> FAILED. Exit code: {e.returncode}")
            failures.append(rna)
        except Exception as e:
            print(f"  -> FAILED with error: {e}")
            failures.append(rna)
        
        print("-" * 60)

    # 5. Summary
    if failures:
        print(f"\nBatch completed with {len(failures)} failures:")
        for f in failures:
            print(f" - {f}")
        sys.exit(1)
    else:
        print("\nBatch completed successfully.")

if __name__ == "__main__":
    main()
