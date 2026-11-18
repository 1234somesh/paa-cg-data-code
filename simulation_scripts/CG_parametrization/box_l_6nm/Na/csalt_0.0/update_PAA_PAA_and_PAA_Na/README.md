# IBI parametrization: PAA–PAA and PAA–Na potentials (box_l = 6 nm, c_salt = 0.0 M)

## Directory Structure
* `ibi.py`, `run_ibi.sh` -- iterative Boltzmann inversion (IBI) driver scripts
* `paa_cg.data` -- LAMMPS data file of the coarse-grained system
* `/initial` -- initial pair potentials
* `/targets` -- target RDFs from the atomistic reference simulations
* `/step_01` ... `/step_56` -- IBI iterations; each contains the LAMMPS input (`in.lmp`), log (`log.lammps`), the measured RDFs (`rdf_*.dat`) and the potentials for the next iteration (`/updated_potentials`)
* `merit.py` -- computes the merit function of the PAA–Na and PAA–PAA RDFs for every iteration and plots it (`merit_series.pdf`)
* `plot_rdf_comparison.py`, `shift_to_zero.py` -- helper scripts

Trajectories (`*.lammpstrj`) are not included because of their size.

## Note on step_46 (rerun required)
The LAMMPS run of `step_46` terminated with the error
`Pair distance < table inner cutoff` (`log.lammps`), so its RDF files are empty.
A rerun of the parametrization was therefore required; iterations `step_47` ... `step_56`
come from this rerun. In `merit_series.pdf`, iteration 46 is skipped and the points of
iterations 45 and 47 are connected
