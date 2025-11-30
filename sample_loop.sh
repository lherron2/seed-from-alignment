#!/usr/bin/env bash
set -euo pipefail

# Resolve directory of this script so paths work even if called from elsewhere
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# ---- CaCoFold input ----
STO="${SCRIPT_DIR}/cacofold/combined_1.cacofold.sto"
COV="${SCRIPT_DIR}/cacofold/combined_1.cacofold.cov"
SEQ="9E5I_1_Chains"                # CaCoFold sequence name in the .sto

# Derive prefix (9E5I_1) and seq FASTA path (9E5I_1/seq.fa)
SEQ_PREFIX="${SEQ%_Chains}"        # -> 9E5I_1
SEQ_DIR="${SCRIPT_DIR}/${SEQ_PREFIX}"
SEQ_FASTA="${SEQ_DIR}/seq.fa"

# ---- RNAstructure Fold executable ----
FOLD_EXE="${PROJECT}/repos/RNAstructure/exe/Fold"

# ---- Sampling parameters ----
N_SEEDS=5
TEMPS=(0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 2.0 3.0 4.0 5.0)
#TEMPS=(1.0)

# ---- Filter parameters (for filter_db_for_rosetta.py) ----
HAMMING_THRESHOLD=3
#MAX_STRUCTURES=200

# Put all per-seed outputs in 9E5I_1/db
OUT_DIR="${SEQ_DIR}/db"
mkdir -p "${OUT_DIR}"

###############################################
#  CLEAN OLD FILES FOR THIS SEQUENCE/RUN   #
###############################################
# Delete old per-seed raw and refined DBs, plus combined/filtered DBs
rm -f "${OUT_DIR}/${SEQ}"_beta*_seed*.db
rm -f "${OUT_DIR}/${SEQ}"_beta*_seed*_refined.db
rm -f "${OUT_DIR}/${SEQ}_all_refined.db"
rm -f "${OUT_DIR}/${SEQ}_filtered_rosetta.db"

# --------------------------------------------------------------------
# 1) SAMPLE + REFINE PER (beta, seed)
# --------------------------------------------------------------------
for beta in "${TEMPS[@]}"; do
    for (( seed=0; seed<${N_SEEDS}; seed++ )); do
        # Base path for this (beta, seed)
        base="${OUT_DIR}/${SEQ}_beta${beta}_seed${seed}"
        out_db="${base}.db"
        summary="${base}.summary.txt"
        refined_db="${base}_refined.db"

        echo "[INFO] beta=${beta}, seed=${seed} -> ${out_db}"

        # 1) Sample structures from CaCoFold
        python "${SCRIPT_DIR}/sample_cacofold_structures.py" \
            "${STO}" \
            --seq "${SEQ}" \
            --beta "${beta}" \
            --seed "${seed}" \
            --n-samples 100 \
            --thin 1000 \
            --min-loop-sep 1 \
            --out-db "${out_db}" \
	    --summary $summary \
	    --cov-file $COV \
	    --cov-mode 'logE_power' \
	    --cov-alpha 3.0 \
	    --cov-min-power 0.1 \
	    --cov-forbid-negative

        # 2) Refine long unpaired regions (>15 nt) with RNAstructure Fold
        echo "[INFO] Refining long unpaired regions in ${out_db} -> ${refined_db}"
        python "${SCRIPT_DIR}/refine_unpaired_regions.py" \
            --seq-file "${SEQ_FASTA}" \
            --db-in "${out_db}" \
            --db-out "${refined_db}"
    done
done

# --------------------------------------------------------------------
# 2) CONCATENATE ALL REFINED STRUCTURES INTO A SINGLE DB
# --------------------------------------------------------------------
ALL_REFINED_DB="${OUT_DIR}/${SEQ}_all_refined.db"
echo "[INFO] Concatenating all refined .db files into ${ALL_REFINED_DB}"

# This will gather every *_refined.db produced above (fresh each run now)
cat "${OUT_DIR}/${SEQ}"_beta*_seed*_refined.db > "${ALL_REFINED_DB}"

# Optional sanity check:
# echo "[INFO] Total structures in combined DB: $(wc -l < "${ALL_REFINED_DB}")"

# --------------------------------------------------------------------
# 3) RUN FILTER_DB_FOR_ROSETTA ON THE CONCATENATED FILE
# --------------------------------------------------------------------
FILTERED_DB="${OUT_DIR}/${SEQ}_filtered_rosetta.db"

echo "[INFO] Running filter_db_for_rosetta.py on ${ALL_REFINED_DB}"
python "${SCRIPT_DIR}/filter_db_for_rosetta.py" \
    --db-in "${ALL_REFINED_DB}" \
    --db-out "${FILTERED_DB}" \
    --hamming-threshold "${HAMMING_THRESHOLD}" \
    --seq-file "${SEQ_FASTA}"
    # --max-structures "${MAX_STRUCTURES}"

echo "[INFO] Done. Filtered Rosetta-compatible structures:"
echo "       ${FILTERED_DB}"

