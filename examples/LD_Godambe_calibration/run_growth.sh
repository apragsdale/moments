#!/bin/bash
# Growth / asymmetric-migration model (growth_asym_calibration.py): one dataset.
# Arguments: stage, scenario, seed, number of regions, bootstrap refits, number of
# processes. Stages:
#   run      simulate and analyze a new dataset (both methods)
#   analyze  shrink an existing dataset of the second run (compact_data.py), then
#            reanalyze it with both methods
#   compact  only shrink an existing dataset
set -e

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "stage $1, scenario $2, seed $3, $4 regions, $5 bootstrap refits, $6 processes"

if [ "$1" = "run" ]; then
    python3 growth_asym_calibration.py run --scenario "$2" --seed "$3" \
        --regions "$4" --boot-refits "$5" --procs "$6"
else
    # skips files that are already compact
    python3 compact_data.py "data/$2_$3.pkl"
    if [ "$1" = "analyze" ]; then
        for method in ld sfs; do
            python3 growth_asym_calibration.py analyze --scenario "$2" --seed "$3" \
                --method $method --boot-refits "$5" --procs "$6"
        done
    fi
fi

echo "finished stage $1, scenario $2, seed $3"
