# Figure 4: stretching free energy and Rg scaling (Na⁺, salt-free)

Plotted with `plotting_scripts/Fig4_stretching_free_energy_Rg_scaling/plot_FE_scaling.py`,
which also fits R_g ∝ N^ν and prints ν (AA: weighted fit, CG: unweighted fit).

| File | Panel | Columns |
|---|---|---|
| `a_PMF_AA.csv` | (a) | R_ee (nm), free energy (kJ/mol), atomistic, PAA 20-mer |
| `a_PMF_CG.csv` | (a) | R_ee (nm), free energy (kJ/mol), CG model |
| `a_PMF_Kremer_Grest.csv` | (a) | R_ee (nm), free energy (kJ/mol), Kremer–Grest bead-spring model |
| `b_Rg_vs_N.csv` | (b) | N, ⟨R_g⟩ AA (nm), error AA (nm), ⟨R_g⟩ CG (nm) |

The free energies are given as in the source files; the plotting script shifts each curve to zero at its minimum.

Sources:
* `a_PMF_AA.csv`: `simulation_scripts/atomistic/stretching_free_energy/profile.xvg`
* `a_PMF_CG.csv`: `simulation_scripts/CG/stretching_free_energy/pmf_grossfield.dat`
* `b_Rg_vs_N.csv`: free-chain simulations in `simulation_scripts/atomistic/free_chain/N_*` and `simulation_scripts/CG/free_chain/Na/N_*`
