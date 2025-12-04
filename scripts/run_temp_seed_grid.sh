export DATAPATH='../../RNAstructure/data_tables'

SCRIPT_DIR=$(dirname "$(realpath "$0")")
ROOT_DIR=$(dirname "$SCRIPT_DIR")
EXAMPLE_DIR="${ROOT_DIR}/example"

# Changed .sh to .py
python "${SCRIPT_DIR}/run_temp_seed_grid.py" \
    --config "${EXAMPLE_DIR}/pipeline_config.yaml" \
    --temps 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,4,5,10,20 \
    --seeds 0,1,2 \
    --out 'gridsearch.db'
