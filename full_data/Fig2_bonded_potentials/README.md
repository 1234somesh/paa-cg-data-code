# Figure 2: CG bonded potentials and AA vs CG bonded distributions (Na⁺, salt-free)

Plotted with `plotting_scripts/Fig2_bonded_potentials/plot_with_dist.py`.

| File | Panel | Columns |
|---|---|---|
| `a_bond_potential.csv` | (a) | bond length r (nm), V (kJ/mol) |
| `b_angle_potential.csv` | (b) | angle θ (deg), V (kJ/mol) |
| `c_dihedral_potential.csv` | (c) | dihedral φ (deg), V (kJ/mol) |
| `d_bond_distribution_AA.csv`, `d_bond_distribution_CG.csv` | (d) | r (nm), P(r) |
| `e_angle_distribution.csv` | (e) | θ (deg), P_AA, P_CG |
| `f_dihedral_distribution.csv` | (f) | φ (deg), P_AA, P_CG |

Sources:
* potentials: `simulation_scripts/CG/tabulated_potentials_Na/salt_free/{bond_M_M,angle_MMM,dihedral_MMMM}.table`, converted from Å and kcal/mol
* distributions: `simulation_scripts/CG_parametrization/box_l_6nm/Na/csalt_0.0/update_PAA_PAA_and_PAA_Na/step_56/bond_dist/*.dist.new`; angles and dihedrals converted from rad to deg, probabilities as in the source files
