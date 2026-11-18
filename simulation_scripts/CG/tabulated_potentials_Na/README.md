# Tabulated CG Potentials (Na⁺ counterions)

## Directory Structure
* `/salt_free` -- bonded (`bond_M_M.table`, `angle_MMM.table`, `dihedral_MMMM.table`) and non-bonded (`pair_PAA_PAA.table`, `pair_PAA_NA.table`, `pair_NA_NA.table`) CG potentials shared by the following salt-free simulations:
  * `../free_chain/Na/N_9`, `N_11`, `N_13`, `N_15`, `N_17`, `N_19`, `N_25`
  * `../single_chain_CGM/N=20/LZ_3.2`, `LZ_3.4`, `LZ_3.8`, `LZ_4.4`
  * `../stretching_free_energy`

All other simulation directories contain their own tables.

## Usage
The LAMMPS inputs (`in.lmp`) read the tables from the run directory. Before running one of
the simulations listed above, copy the tables into that directory, e.g.

    cp tabulated_potentials_Na/salt_free/*.table free_chain/Na/N_9/
