# Tabulated CG Potentials (Cs⁺ counterions)

## Directory Structure
* `/salt_free` -- bonded (`bond_M_M.table`, `angle_MMM.table`, `dihedral_MMMM.table`) and non-bonded (`pair_PAA_PAA.table`, `pair_PAA_CS.table`, `pair_CS_CS.table`) CG potentials for salt-free simulations with Cs⁺ counterions
  * bonded potentials are not updated during IBI and are taken from `CG_parametrization/box_l_6nm/Cs/csalt_0.0/update_PAA_PAA_and_PAA_Cs/initial`
  * non-bonded potentials are the final IBI potentials, `step_83/updated_potentials` of the same parametrization

## Usage
The LAMMPS inputs (`in.lmp`) read the tables from the run directory. Before running a
salt-free Cs⁺ simulation, copy the tables into that directory, e.g.

    cp tabulated_potentials_Cs/salt_free/*.table <run_directory>/
