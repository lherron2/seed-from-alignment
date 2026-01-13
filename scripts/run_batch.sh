#!/bin/bash

python "/home/lukas/research/seed-from-alignment/seed-from-alignment/scripts/run_batch.py" \
--rnas-file "/home/lukas/research/seed-from-alignment/seed-from-alignment/cbl/rnas.txt" \
--config "/home/lukas/research/seed-from-alignment/seed-from-alignment/example/pipeline_config.yaml" \
--output-base "/home/lukas/research/seed-from-alignment/seed-from-alignment/cbl" \
--pipeline-script "/home/lukas/research/seed-from-alignment/seed-from-alignment/run_pipeline.py"
