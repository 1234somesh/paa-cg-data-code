# Data of the ESI / ECC (Dang force field with electronic continuum correction) figures

The `ESI_*` folders hold the data series plotted in the corresponding figures, exported from
the plotting scripts in `analysis/ECC_dang/` (values exactly as plotted, i.e. in the plotted units).

| Folder | Figure | Plotting script (`analysis/ECC_dang/...`) |
|---|---|---|
| `ESI_Monomer_Vfull_RDF_ECC` | `Monomer_Vfull_RDF_ECC.pdf` | `potentials/NA/plot_potentials_with_coulomb.py` |
| `ESI_RDF_evolution_box_5_1_6_0` | `RDF_evolution_box_5_1_6_0.pdf` | `potentials/NA/plot_rdf_evolution_box_sizes.py` |
| `ESI_Vfull_PAA_NA_PAA_PAA_box_sizes` | `Vfull_PAA_NA_PAA_PAA_box_sizes.pdf` | `potentials/NA/plot_full_potential_box_sizes.py` |
| `ESI_Rg_vs_box_and_dist_ECC` | `Rg_vs_box_and_dist_ECC.pdf` | `Rgs/plot_rg_vs_box_and_dist.py` |
| `ESI_rdfs_NA_K_CS_1x3_ECC` | `rdfs_NA_K_CS_1x3.pdf` | `RDFs/NA_K_CS/plot_rdf.py` |
| `ESI_qeff_AA_vs_CG_ECC` | `qeff_rc0p72nm_AA_vs_CG_scaled.pdf` | `qeff/free_chain/plot_both.py` |

## File naming
`<panel>_<panel title>_<nn>_<legend label>.csv`
* `<panel>` -- panel letter (a, b, c, ... in reading order)
* `<nn>` -- order in which the series is drawn in that panel
* `n_<i>` in the label -- IBI iteration i; `Target` -- atomistic target
* the column headers are the axis labels of the panel; a third column `*_err` holds the error bars

Series without a legend entry are named by colour and line style:
`color1f77b4` (blue) = AA, `colorff7f0e` (orange) = CG; in `ESI_Rg_vs_box_and_dist_ECC`
the line styles `solid` / `dashed` / `dotted` correspond to box lengths L = 4.3 / 5.1 / 6.0 nm.
