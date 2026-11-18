# Figure 7: radius of gyration vs salt concentration (Na⁺, K⁺, Cs⁺)

Plotted with `plotting_scripts/Fig7_Rg_vs_csalt/plot_combined.py`.

`Rg_vs_csalt.csv` columns:
* `ion` -- counterion (`Na`, `K`, `Cs`)
* `csalt_M` -- added salt concentration (mol/L)
* `Rg_AA_nm`, `Rg_AA_err_nm` -- ⟨R_g⟩ and error, atomistic
* `Rg_CG_nm`, `Rg_CG_err_nm` -- ⟨R_g⟩ and error, CG with the salt-free potentials
* `Rg_CG_reparam_nm`, `Rg_CG_reparam_err_nm` -- ⟨R_g⟩ and error, CG reparametrized at the given salt concentration (empty where not computed)

Simulation inputs: `simulation_scripts/CG/free_chain/<ion>/N_20/csalt_*` (CG with salt) and
`simulation_scripts/CG_parametrization/box_l_6nm/<ion>/csalt_*` (reparametrization).
