#!/bin/bash
# Simulates and analyzes one dataset.
# Arguments: scenario (mig, weak or null), seed, number of regions, bootstrap refits.
set -e

# one core per job: keep numpy and its linear algebra libraries to one thread
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "running scenario $1, seed $2, $3 regions, $4 bootstrap refits"

python3 godambe_calibration.py run --scenario "$1" --seed "$2" --regions "$3" \
    --boot-refits "$4" --procs 1

echo "finished running scenario $1, seed $2"
