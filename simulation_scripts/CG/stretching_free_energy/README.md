# Stretching Free Energy of a Single CG PAA Chain (N = 20, Na⁺ counterions)

Potential of mean force (PMF) along the end-to-end distance $R_\mathrm{ee}$ of a single
coarse-grained PAA chain, obtained by umbrella sampling in LAMMPS and unbiased with WHAM.

## Directory Structure
* `in.lmp` -- LAMMPS input for one umbrella window; a harmonic spring (`fix spring couple`) between the two end beads (ids 1 and 20) restrains $R_\mathrm{ee}$ around `r0`, and $R_\mathrm{ee}$ is written to `ree.dat`
* `paa_cg.data` -- LAMMPS data file of the CG chain with its counterions
* tabulated CG potentials (`pair_*.table`, `bond_M_M.table`, `angle_MMM.table`, `dihedral_MMMM.table`) are not stored here; they are in `../tabulated_potentials_Na/salt_free`
* `launch.sh` -- SLURM array job that runs 266 windows with `r0` = 18.5 ... 45.0 Å (step 0.1 Å) and spring constant `k = 0.5` (LAMMPS `real` units), each in its own `win_<r0>/` directory
* `make_metadat.py` -- converts each window's `ree.dat` from Å to nm (`ree_nm.dat`) and writes the WHAM metadata file
* `metadata.dat` -- WHAM metadata file (time series file, window centre in nm, spring constant) for the 266 windows
* `pmf_grossfield.dat` -- resulting PMF from the WHAM code of A. Grossfield (columns: $R_\mathrm{ee}$ in nm, free energy, error, probability, error), 200 bins from 1.91 to 4.49 nm

## How to Reproduce
1. Copy the tables into this directory (`cp ../tabulated_potentials_Na/salt_free/*.table .`), then run `launch.sh` (LAMMPS version 19 Nov 2024) from here: `sbatch launch.sh`
2. Run `python3 make_metadat.py` in the same directory to generate `win_*/ree_nm.dat` and `metadata.dat`
3. Run WHAM ([A. Grossfield](http://membrane.urmc.rochester.edu/?page_id=126)) with `metadata.dat` at T = 300 K to obtain `pmf_grossfield.dat`

The per-window time series (`win_*/ree.dat`, `win_*/ree_nm.dat`) and LAMMPS logs are not included because of their size.
