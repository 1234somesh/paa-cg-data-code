# Figure 3: short-range CG pair potentials and RDFs (Na⁺, salt-free)

Plotted with `plotting_scripts/Fig3_nonbonded_potentials_RDFs/plot_vshort.py`.
All files: first column r (nm); potentials V_short in kJ/mol (without the Coulomb term); RDFs dimensionless.

| File | Panel | Content |
|---|---|---|
| `a_PAA_Na_potential_initial.csv`, `a_PAA_Na_potential_final.csv` | (a) | PAA–Na potential, initial and final IBI iteration |
| `b_PAA_PAA_potential_initial.csv`, `b_PAA_PAA_potential_final.csv` | (b) | PAA–PAA potential, initial and final |
| `c_Na_Na_potential.csv` | (c) | Na–Na potential (from the ion–ion PMF, not updated by IBI) |
| `d_PAA_Na_rdf_*.csv` | (d) | PAA–Na RDF: AA target, CG initial, CG final |
| `e_PAA_PAA_rdf_*.csv` | (e) | PAA–PAA RDF: AA target, CG initial, CG final |
| `f_Na_Na_rdf_*.csv` | (f) | Na–Na RDF: AA target, CG initial, CG final |

Sources (in `simulation_scripts/CG_parametrization/box_l_6nm/Na/csalt_0.0/update_PAA_PAA_and_PAA_Na/`):
* initial potentials: `initial/pair_*.table`
* final potentials: `step_55/updated_potentials/pair_*.table` (identical to `simulation_scripts/CG/tabulated_potentials_Na/salt_free`)
* AA target RDFs: `targets/*.dist.dat`
* final CG RDFs: `step_56/rdf_*.dat`
* initial CG RDFs: a CG run with the initial potentials (first iteration of the rerun of the parametrization)

Tables were converted from Å and kcal/mol.
