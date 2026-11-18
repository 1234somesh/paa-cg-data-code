#!/usr/bin/env bash
set -euo pipefail

# IBI driver -- update_6_potential variant.
#
# Pairs updated simultaneously each step (6 of 6):
#   PAA_PAA, PAA_K, PAA_CL, K_CL, K_K, CL_CL
# The remaining pairs are copied through unchanged by ibi.py, but ALL SIX
# tables are written into each step's updated_potentials/, so every
# pair_coeff below can point at the previous step consistently.
#
# Usage:
#   ./run_ibi.sh LAST_STEP
#
# Requirements:
#   - step_01 exists and contains: in.lmp and updated_potentials/
#   - ibi.py and shift_to_zero.py are in the parent dir of step_*
#   - fit_short_range.py is in the same directory as run_ibi.sh (called as ../fit_short_range.py from step_*)
#     (one level deeper than the original Try_* layout)
#   - lmp (LAMMPS) and python3 in PATH

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 LAST_STEP_INDEX"
    exit 1
fi

LAST_STEP="$1"

if [[ ! -d "step_01" ]]; then
    echo "ERROR: step_01 directory not found."
    exit 1
fi
if [[ ! -d "step_01/updated_potentials" ]]; then
    echo "ERROR: step_01/updated_potentials not found."
    exit 1
fi

FIRST_STEP=2
for ((i=2; i<=LAST_STEP; i++)); do
    s=$(printf "step_%02d" "$i")
    if [[ -d "$s" && -d "$s/updated_potentials" ]]; then
        echo "Found completed ${s}, skipping."
        FIRST_STEP=$((i+1))
    else
        break
    fi
done

if (( FIRST_STEP > LAST_STEP )); then
    echo "All steps up to step_$(printf "%02d" "$LAST_STEP") already completed."
    exit 0
fi

echo "Starting from step_$(printf "%02d" "$FIRST_STEP") up to step_$(printf "%02d" "$LAST_STEP")."
echo

for ((i=FIRST_STEP; i<=LAST_STEP; i++)); do
    step=$(printf "step_%02d" "$i")
    prev_step=$(printf "step_%02d" "$((i-1))")

    echo "============================================================"
    echo "  IBI iteration for ${step}  (previous: ${prev_step})"
    echo "============================================================"

    if [[ -d "${step}" ]]; then
        echo "WARNING: ${step} already exists, reusing it."
    else
        mkdir "${step}"
    fi

    cd "${step}"

    if [[ ! -f "../${prev_step}/in.lmp" ]]; then
        echo "ERROR: ../${prev_step}/in.lmp not found!"
        exit 1
    fi
    cp "../${prev_step}/in.lmp" ./in.lmp

    echo "Pointing all 6 pair_coeff lines at ../${prev_step}/updated_potentials/"
    sed -i "s|^pair_coeff[[:space:]]\+1[[:space:]]\+1[[:space:]]\+table.*|pair_coeff      1 1 table ../${prev_step}/updated_potentials/pair_PAA_PAA.table   PAA_PAA|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+1[[:space:]]\+2[[:space:]]\+table.*|pair_coeff      1 2 table ../${prev_step}/updated_potentials/pair_PAA_CS.table  PAA_CS|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+1[[:space:]]\+3[[:space:]]\+table.*|pair_coeff      1 3 table ../${prev_step}/updated_potentials/pair_PAA_CL.table PAA_CL|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+2[[:space:]]\+3[[:space:]]\+table.*|pair_coeff      2 3 table ../${prev_step}/updated_potentials/pair_CS_CL.table CS_CL|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+2[[:space:]]\+2[[:space:]]\+table.*|pair_coeff      2 2 table ../${prev_step}/updated_potentials/pair_CS_CS.table CS_CS|" in.lmp
    sed -i "s|^pair_coeff[[:space:]]\+3[[:space:]]\+3[[:space:]]\+table.*|pair_coeff      3 3 table ../${prev_step}/updated_potentials/pair_CL_CL.table CL_CL|" in.lmp

    echo "Running LAMMPS in ${step}..."
    lmp -in in.lmp

    echo "Running IBI update via ibi.py..."
    python3 ../ibi.py

    echo "Post-processing potentials (monotonic core + U(20A)=0)..."
    # Enforce a smooth, strictly repulsive core on EVERY table LAMMPS will
    # load next step -- including pairs this variant does not update, whose
    # stock tables have a FLAT capped core (zero force inside = overlap risk).
    # fit_short_range.py is idempotent, so re-applying it is a no-op.
    python3 ../fit_short_range.py ./updated_potentials/pair_PAA_PAA.table ./updated_potentials/pair_PAA_PAA.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_PAA_PAA.table ./updated_potentials/pair_PAA_PAA.table --rcut 20.00
    python3 ../fit_short_range.py ./updated_potentials/pair_PAA_CS.table ./updated_potentials/pair_PAA_CS.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_PAA_CS.table ./updated_potentials/pair_PAA_CS.table --rcut 20.00
    python3 ../fit_short_range.py ./updated_potentials/pair_PAA_CL.table ./updated_potentials/pair_PAA_CL.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_PAA_CL.table ./updated_potentials/pair_PAA_CL.table --rcut 20.00
    python3 ../fit_short_range.py ./updated_potentials/pair_CS_CL.table ./updated_potentials/pair_CS_CL.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_CS_CL.table ./updated_potentials/pair_CS_CL.table --rcut 20.00
    python3 ../fit_short_range.py ./updated_potentials/pair_CS_CS.table ./updated_potentials/pair_CS_CS.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_CS_CS.table ./updated_potentials/pair_CS_CS.table --rcut 20.00
    python3 ../fit_short_range.py ./updated_potentials/pair_CL_CL.table ./updated_potentials/pair_CL_CL.table
    python3 ../shift_to_zero.py ./updated_potentials/pair_CL_CL.table ./updated_potentials/pair_CL_CL.table --rcut 20.00

    echo "Finished ${step}."
    echo

    cd ..
done

echo "All iterations up to step_$(printf "%02d" "$LAST_STEP") completed."
