#!/bin/bash
# Simulates and analyzes one dataset.
# Arguments: scenario (mig, weak or null), seed, number of regions, bootstrap refits,
# number of processes.
set -e

# one thread per process: the script runs its own processes in parallel
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "running scenario $1, seed $2, $3 regions, $4 bootstrap refits, $5 processes"

python3 godambe_calibration.py run --scenario "$1" --seed "$2" --regions "$3" \
    --boot-refits "$4" --procs "$5"

echo "finished running scenario $1, seed $2"
