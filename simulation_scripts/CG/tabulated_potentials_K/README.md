# Tabulated CG Potentials (K⁺ counterions)

## Directory Structure
* `/salt_free` -- bonded (`bond_M_M.table`, `angle_MMM.table`, `dihedral_MMMM.table`) and non-bonded (`pair_PAA_PAA.table`, `pair_PAA_K.table`, `pair_K_K.table`) CG potentials for salt-free simulations with K⁺ counterions
  * bonded potentials are not updated during IBI and are taken from `CG_parametrization/box_l_6nm/K/csalt_0.0/update_PAA_PAA_and_PAA_K/initial`
  * non-bonded potentials are the final IBI potentials, `step_79/updated_potentials` of the same parametrization

## Usage
The LAMMPS inputs (`in.lmp`) read the tables from the run directory. Before running a
salt-free K⁺ simulation, copy the tables into that directory, e.g.

    cp tabulated_potentials_K/salt_free/*.table <run_directory>/
