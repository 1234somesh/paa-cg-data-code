# Figure 6: ion-specific RDFs and chain size, AA vs CG (Na⁺, K⁺, Cs⁺, salt-free)

Plotted with `plotting_scripts/Fig6_ion_specific_RDFs_chain_size/plot_rdf.py`.

| File | Panel | Columns |
|---|---|---|
| `a_PAA_ion_rdf_<ion>_AA.csv`, `a_PAA_ion_rdf_<ion>_CG.csv` | (a) | r (nm), PAA–ion RDF |
| `b_PAA_PAA_rdf_<ion>_AA.csv`, `b_PAA_PAA_rdf_<ion>_CG.csv` | (b) | r (nm), PAA–PAA RDF |
| `c_ion_ion_rdf_<ion>_AA.csv`, `c_ion_ion_rdf_<ion>_CG.csv` | (c) | r (nm), ion–ion RDF |
| `d_R2_of_n.csv` | (d) | n (number of CG bonds), ⟨R²(n)⟩ (nm²) for each ion, AA and CG; the plot shows ⟨R²(n)⟩^(1/2) |

`<ion>` is `Na`, `K` or `Cs`.

Sources (in `simulation_scripts/CG_parametrization/box_l_6nm/<ion>/csalt_0.0/update_PAA_PAA_and_PAA_<ion>/`):
* AA RDFs: `targets/*.dist.dat`
* CG RDFs and CG ⟨R²(n)⟩: final IBI iteration, `step_56` (Na), `step_80` (K), `step_84` (Cs): `rdf_*.dat`, `R2_of_n_lammps_all.csv`
* AA ⟨R²(n)⟩: atomistic reference simulations at box length 6 nm
