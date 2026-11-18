"""Figure 5 (main text, data panels): bonded distributions (top row) and RDFs
(bottom row) of the single-chain cell-gel model, PAA 30-mer with Na+, L_z = 6.1 nm;
AA vs CG.

Reads the CSV files in ../../full_data/Fig5_single_chain_CGM_N30 and writes
../../figures/bond_ang_dih_and_rdfs_sc_cgm_6_2.pdf
(the schematic in the manuscript figure was added separately).
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

here = Path(__file__).resolve().parent
data = here.parent.parent / "full_data" / "Fig5_single_chain_CGM_N30"
outdir = here.parent.parent / "figures"


def load(name):
    return np.loadtxt(data / name, delimiter=",", skiprows=1)


fig, axs = plt.subplots(2, 3, figsize=(16, 9), sharey=False)

panels = [
    # (axis, file prefix, title, xlabel, ylabel, xlim, xticks)
    (axs[0, 0], "a_bond", "Bond", "r (nm)", "P(r)", (0.22, 0.5), None),
    (axs[0, 1], "b_angle", "Angle", r"$\theta$ (deg)", r"P($\theta$)", None, np.arange(0, 181, 45)),
    (axs[0, 2], "c_dihedral", "Dihedral", r"$\phi$ (deg)", r"P($\phi$)", None, np.arange(-180, 181, 90)),
    (axs[1, 0], "d_PAA_Na_rdf", "PAA–Na", "r (nm)", "RDF", (0.20, 1.7), None),
    (axs[1, 1], "e_PAA_PAA_rdf", "PAA–PAA", "r (nm)", "RDF", (0.35, 1.7), None),
    (axs[1, 2], "f_Na_Na_rdf", "Na–Na", "r (nm)", "RDF", (0.25, 1.7), None),
]

for ax, prefix, title, xlabel, ylabel, xlim, xticks in panels:
    for tag in ("AA", "CG"):
        d = load(f"{prefix}_{tag}.csv")
        ax.plot(d[:, 0], d[:, 1], label=tag, linewidth=2.0)
    ax.set_title(title, fontsize=16, fontweight="normal")
    ax.set_xlabel(xlabel, fontsize=16, fontweight="normal")
    ax.set_ylabel(ylabel, fontsize=16, fontweight="normal")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.tick_params(labelsize=16, length=8.0, width=2.0)
    ax.legend(fontsize=15)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout()
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "bond_ang_dih_and_rdfs_sc_cgm_6_1.pdf", dpi=300, format="pdf")
