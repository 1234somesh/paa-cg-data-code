# Data and scripts for "Capturing Ion-Specific Effects in Polyacrylic Acid via an Implicit-Solvent, Explicit-Electrostatics Coarse-Grained Model"

Somesh Kurahatti, Svyatoslav Kondrat, David Beyer, Christian Holm
Institute for Computational Physics, University of Stuttgart

This repository contains the simulation inputs, the coarse-grained (CG) parametrization runs,
the analysis scripts and the plotted data needed to reproduce the results of the article. Trajectories are not included because of their size;
they can be regenerated from the inputs provided here.

## Repository layout

```
Kurahatti26a-data/
├── simulation_scripts/   inputs for all simulations (atomistic and CG) and the parametrization runs
├── full_data/            numerical data of every figure panel (CSV), one folder per figure
├── plotting_scripts/     scripts that read full_data/ and produce the figures
├── figures/              figures produced by plotting_scripts/
└── analysis/             analysis scripts used to compute the data from the trajectories
```

### `simulation_scripts/`

| Folder | Content |
|---|---|
| `atomistic/free_chain/N_*` | GROMACS inputs (structure, topology) of a single PAA chain in water, N = 9–25 monomers; force field in `amber03_dang.ff/` |
| `atomistic/single_chain_CGM/` | atomistic single-chain cell-gel model (N = 20, 30) at different box lengths `LZ_*` |
| `atomistic/stretching_free_energy/` | umbrella-sampling inputs and PMF (`profile.xvg`) of the end-to-end stretching free energy |
| `CG_parametrization/box_l_6nm/{Na,K,Cs}/csalt_*/` | iterative Boltzmann inversion (IBI) of the non-bonded CG potentials against atomistic RDFs, for each counterion and salt concentration: driver scripts (`ibi.py`, `run_ibi.sh`), initial potentials, atomistic targets and every IBI iteration (`step_*`) |
| `Implicit_CG_parametrization/` | VOTCA (`csg_inverse`) IBI run of the implicit CG model used for comparison in the ESI: settings, mapping files, targets, and the updated potentials (`*.pot.new`) and RDFs (`*.dist.new`) of every iteration |
| `CG/tabulated_potentials_{Na,K,Cs}/salt_free/` | final salt-free CG potentials (bonded and non-bonded LAMMPS tables) for each counterion |
| `CG/free_chain/{Na,K,Cs}/N_*/csalt_*/` | LAMMPS inputs (`in.lmp`, `paa_cg.data`, tables) of the CG production runs, free chain, with and without added salt |
| `CG/single_chain_CGM/` | CG single-chain cell-gel model runs |
| `CG/stretching_free_energy/` | CG umbrella sampling of the stretching free energy and its PMF (`pmf_grossfield.dat`) |

Folders with additional details contain their own `README.md`.

### `full_data/` and `plotting_scripts/`

| Figure | Data (`full_data/`) | Script (`plotting_scripts/`) | Output (`figures/`) |
|---|---|---|---|
| Fig. 2 – bonded potentials and distributions | `Fig2_bonded_potentials/` | `Fig2_bonded_potentials/plot_with_dist.py` | `bond_ang_dih_pot_dist_combined.pdf` |
| Fig. 3 – non-bonded potentials and RDFs (Na⁺) | `Fig3_nonbonded_potentials_RDFs/` | `Fig3_nonbonded_potentials_RDFs/plot_vshort.py` | `Monomer_Vshort_RDF_2x3.pdf` |
| Fig. 4 – stretching free energy, R_g scaling | `Fig4_stretching_free_energy_Rg_scaling/` | `Fig4_stretching_free_energy_Rg_scaling/plot_FE_scaling.py` | `FE_scaling.pdf` |
| Fig. 5 – single-chain cell-gel model, N = 30 | `Fig5_single_chain_CGM_N30/` | `Fig5_single_chain_CGM_N30/plot_6_1.py` | `bond_ang_dih_and_rdfs_sc_cgm_6_1.pdf` |
| Fig. 6 – ion-specific RDFs and chain size | `Fig6_ion_specific_RDFs_chain_size/` | `Fig6_ion_specific_RDFs_chain_size/plot_rdf.py` | `rdfs_NA_K_CS_2x2.pdf` |
| Fig. 7 – R_g vs salt concentration | `Fig7_Rg_vs_csalt/` | `Fig7_Rg_vs_csalt/plot_combined.py` | `Rg_vs_csalt.pdf` |
| ESI figures (ECC force field) | `ESI_*/` (see `full_data/README_ESI_ECC.md`) | `analysis/ECC_dang/` | – |

Each `full_data/Fig*/README.md` lists the files, their columns and units, and the simulation
folder the data were computed from. All CSV files have a one-line header.

To regenerate a figure:

```bash
cd plotting_scripts/Fig2_bonded_potentials
python3 plot_with_dist.py        # writes ../../figures/bond_ang_dih_pot_dist_combined.pdf
```

The scripts locate the data relative to their own position and can be run from any directory.

### `analysis/`

Scripts used to compute the data from the trajectories, kept as they were used for the
article; file paths inside them refer to the original working directories. Lines that read a
trajectory or topology are marked with a comment `# INPUT: ...` naming the expected file type
(GROMACS `.xtc`/`.trr`/`.tpr`/`.gro` or LAMMPS `.lammpstrj`).

| Folder / file | Analysis |
|---|---|
| `R2_of_n_lammps.py` | mean-squared internal distance ⟨R²(n)⟩ from CG trajectories (Fig. 6d) |
| `cylindrical_analysis/` | cylindrical counterion RDFs and coordination numbers around the chain |
| `spherical_analysis/` | spherical counterion distributions around the chain |
| `qeff/` | effective charge of the chain, free chain and single-chain cell-gel model |
| `compare_with_implicit_CG/` | comparison with the implicit CG model (VOTCA parametrization above) |
| `ECC_dang/` | potentials, RDFs, R_g and effective charge with the ECC-scaled force field (ESI) |

## Software

| Software | Used for |
|---|---|
| GROMACS 2019.6 | atomistic simulations, umbrella sampling |
| LAMMPS (June 2022) | CG simulations, IBI iterations |
| VOTCA 2019.1 | mapping atomistic trajectories to CG, implicit CG parametrization |
| Python 3 with NumPy, pandas, Matplotlib | IBI drivers, analysis and plotting |

## Licence

Data are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0);
scripts are released under the MIT licence.

## Citation

If you use these data or scripts, please cite the article (reference to be added on publication)
and this dataset.
