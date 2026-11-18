# Figure 5: single-chain cell-gel model, PAA 30-mer with Na⁺ (L_z = 6.1 nm), AA vs CG

Plotted with `plotting_scripts/Fig5_single_chain_CGM_N30/plot_6_1.py`
(the schematic in the manuscript figure was added separately).

| File | Panel | Columns |
|---|---|---|
| `a_bond_AA.csv`, `a_bond_CG.csv` | (a) | bond length r (nm), P(r) |
| `b_angle_AA.csv`, `b_angle_CG.csv` | (b) | angle θ (deg), P(θ) |
| `c_dihedral_AA.csv`, `c_dihedral_CG.csv` | (c) | dihedral φ (deg), P(φ) |
| `d_PAA_Na_rdf_AA.csv`, `d_PAA_Na_rdf_CG.csv` | (d) | r (nm), PAA–Na RDF |
| `e_PAA_PAA_rdf_AA.csv`, `e_PAA_PAA_rdf_CG.csv` | (e) | r (nm), PAA–PAA RDF |
| `f_Na_Na_rdf_AA.csv`, `f_Na_Na_rdf_CG.csv` | (f) | r (nm), Na–Na RDF |

Simulation inputs: `simulation_scripts/atomistic/single_chain_CGM/N_30/LZ_6.1` (AA) and
`simulation_scripts/CG/single_chain_CGM/N=30/LZ_6.1` (CG).
