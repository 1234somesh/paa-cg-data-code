#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_ibi.sh LAST_STEP
#
# Example:
#   ./run_ibi.sh 5
# will run steps: step_0X starting from the first *incomplete* step
# (e.g. if step_02 is already done, it will start from step_03)
#
# Requirements:
#   - step_01 exists and contains: in.lmp and updated_potentials/
#   - ibi.py is in the parent directory of step_* (i.e. here)
#   - fit_short_range.py is in the same directory as run_ibi.sh (called as ../fit_short_range.py from step_*)
#     when run from inside step_* directories
#   - lmp (LAMMPS) and python3 are in PATH

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 LAST_STEP_INDEX"
    echo "Example: $0 5   # runs from first incomplete step up to step_05"
    exit 1
fi

LAST_STEP="$1"

# Check step_01 exists
if [[ ! -d "step_01" ]]; then
    echo "ERROR: step_01 directory not found. You must have step_01 ready with updated_potentials/."
    exit 1
fi

if [[ ! -d "step_01/updated_potentials" ]]; then
    echo "ERROR: step_01/updated_potentials directory not found."
    echo "Run the first iteration manually before using this script."
    exit 1
fi

# -------------------------------------------------------------
# Determine FIRST_STEP to run:
# - Start from 2
# - If step_02 + updated_potentials/ exists, skip it → try step_03
# - If step_03 + updated_potentials/ exists, skip it → try step_04
# - ...
# -------------------------------------------------------------
FIRST_STEP=2
for ((i=2; i<=LAST_STEP; i++)); do
    s=$(printf "step_%02d" "$i")
    if [[ -d "$s" && -d "$s/updated_potentials" ]]; then
        echo "Found completed ${s} (has updated_potentials/), skipping."
        FIRST_STEP=$((i+1))
    else
        # first step that is either missing or incomplete → start here
        break
    fi
done

if (( FIRST_STEP > LAST_STEP )); then
    echo "All steps up to step_$(printf "%02d" "$LAST_STEP") already completed. Nothing to do."
    exit 0
fi

echo "Will start IBI iterations from step_$(printf "%02d" "$FIRST_STEP") up to step_$(printf "%02d" "$LAST_STEP")."
echo

# Loop from FIRST_STEP to LAST_STEP
for ((i=FIRST_STEP; i<=LAST_STEP; i++)); do
    step=$(printf "step_%02d" "$i")
    prev_step=$(printf "step_%02d" "$((i-1))")

    echo "============================================================"
    echo "  IBI iteration for ${step}  (previous: ${prev_step})"
    echo "============================================================"

    # Make new step directory if it doesn't exist
    if [[ -d "${step}" ]]; then
        echo "WARNING: ${step} already exists, reusing it."
    else
        mkdir "${step}"
    fi

    cd "${step}"

    # Copy in.lmp from previous step
    if [[ ! -f "../${prev_step}/in.lmp" ]]; then
        echo "ERROR: ../${prev_step}/in.lmp not found!"
        exit 1
    fi
    cp "../${prev_step}/in.lmp" ./in.lmp

    echo "Updating pair_coeff lines in ${step}/in.lmp to use ../${prev_step}/updated_potentials/"

    # Update the pair_coeff lines to point to previous step's updated_potentials
    sed -i "s|^pair_coeff[[:space:]]\+1[[:space:]]\+1[[:space:]]\+table.*|pair_coeff      1 1 table ../${prev_step}/updated_potentials/pair_PAA_PAA.table   PAA_PAA|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+1[[:space:]]\+2[[:space:]]\+table.*|pair_coeff      1 2 table ../${prev_step}/updated_potentials/pair_PAA_NA.table  PAA_NA|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+2[[:space:]]\+2[[:space:]]\+table.*|pair_coeff      2 2 table ../${prev_step}/updated_potentials/pair_NA_NA.table NA_NA|" in.lmp

    echo "Running LAMMPS in ${step}..."
    lmp -in in.lmp

    echo "Running IBI update via ibi.py..."
    python3 ../ibi.py   # uses RDFs in current step_* and previous step's updated_potentials

    echo "Smoothing updated potentials in ${step}/updated_potentials/ ..."
    python3 ../fit_short_range.py ./updated_potentials/pair_PAA_PAA.table ./updated_potentials/pair_PAA_PAA.table 
    python3 ../shift_to_zero.py ./updated_potentials/pair_PAA_PAA.table ./updated_potentials/pair_PAA_PAA.table --rcut 20.00
    python3 ../fit_short_range.py ./updated_potentials/pair_PAA_NA.table ./updated_potentials/pair_PAA_NA.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_PAA_NA.table ./updated_potentials/pair_PAA_NA.table --rcut 20.0
    python3 ../fit_short_range.py ./updated_potentials/pair_NA_NA.table ./updated_potentials/pair_NA_NA.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_NA_NA.table ./updated_potentials/pair_NA_NA.table --rcut 20.0

    # Guard point at r = 0.0001 A (continues the exponential wall) so a close
    # contact can't abort LAMMPS with "Pair distance < table inner cutoff".
    # Must run last: it rebuilds the guard after fit_short_range/shift_to_zero.
    echo "Adding inner guard point (r = 0.0001 A) to updated potentials ..."
    for p in PAA_PAA PAA_NA NA_NA; do
        python3 ../add_inner_guard.py ./updated_potentials/pair_${p}.table ./updated_potentials/pair_${p}.table --rmin 0.0001
    done

    echo "Finished ${step}."
    echo

    cd ..
done

echo "All iterations up to step_$(printf "%02d" "$LAST_STEP") completed."

