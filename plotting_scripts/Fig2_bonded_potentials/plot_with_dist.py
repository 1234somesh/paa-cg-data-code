"""Figure 2 (main text): CG bonded potentials (top row) and AA vs CG bonded
distributions (bottom row) for PAA with Na+ counterions.

Reads the CSV files in ../../full_data/Fig2_bonded_potentials and writes
../../figures/bond_ang_dih_pot_dist_combined.pdf
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

here = Path(__file__).resolve().parent
data = here.parent.parent / "full_data" / "Fig2_bonded_potentials"
outdir = here.parent.parent / "figures"


def load(name):
    return np.loadtxt(data / name, delimiter=",", skiprows=1)


# --- Potentials (CG tables, converted to nm / deg and kJ/mol) ---
bond_pot = load("a_bond_potential.csv")
ang_pot = load("b_angle_potential.csv")
dih_pot = load("c_dihedral_potential.csv")

# --- Distributions (AA vs CG) ---
bond_aa = load("d_bond_distribution_AA.csv")
bond_cg = load("d_bond_distribution_CG.csv")
ang = load("e_angle_distribution.csv")
dih = load("f_dihedral_distribution.csv")

# --- Figure: 2 rows (potentials, distributions) x 3 columns (bond, angle, dihedral) ---
fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex='col')

# ===================== TOP ROW: POTENTIALS =====================
ax = axs[0, 0]
ax.plot(bond_pot[:, 0], bond_pot[:, 1], linewidth=2.0)
ax.set_title("Bond", fontsize=16, fontweight="normal")
ax.set_xlim([0.22, 0.5])
ax.set_ylim([-1, 20])
ax.set_ylabel("V(r) (kJ/mol)", fontsize=16, fontweight="normal")
ax.tick_params(labelsize=16, length=8.0, width=2.0)

ax = axs[0, 1]
ax.set_ylim([-1.0, 20.0])
ax.plot(ang_pot[:, 0], ang_pot[:, 1], linewidth=2.0)
ax.set_title("Angle", fontsize=16, fontweight="normal")
ax.set_ylabel(r"V($\theta$) (KJ/mol)", fontsize=16, fontweight="normal")
ax.tick_params(labelsize=16, length=8.0, width=2.0)

ax = axs[0, 2]
ax.plot(dih_pot[:, 0], dih_pot[:, 1], linewidth=2.0)
ax.set_title("Dihedral", fontsize=16, fontweight="normal")
ax.set_ylabel(r"V($\phi$) (KJ/mol)", fontsize=16, fontweight="normal")
ax.tick_params(labelsize=16, length=8.0, width=2.0)

# ===================== BOTTOM ROW: DISTRIBUTIONS =====================
ax = axs[1, 0]
ax.plot(bond_aa[:, 0], bond_aa[:, 1], label="AA", linewidth=2.0)
ax.plot(bond_cg[:, 0], bond_cg[:, 1], label="CG", linewidth=2.0)
ax.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
ax.set_ylabel("P(r)", fontsize=16, fontweight="normal")
ax.set_xlim([0.22, 0.5])
ax.tick_params(labelsize=16, length=8.0, width=2.0)
ax.legend(fontsize=16)

ax = axs[1, 1]
ax.plot(ang[:, 0], ang[:, 1], label="AA", linewidth=2.0)
ax.plot(ang[:, 0], ang[:, 2], label="CG", linewidth=2.0)
ax.set_xlabel(r"$\theta$ (Deg)", fontsize=16, fontweight="normal")
ax.set_ylabel(r"P($\theta$)", fontsize=16, fontweight="normal")
ax.tick_params(labelsize=16, length=8.0, width=2.0)
ax.set_xticks(np.arange(0, 181, 45))
ax.legend(fontsize=16)

ax = axs[1, 2]
ax.plot(dih[:, 0], dih[:, 1], label="AA", linewidth=2.0)
ax.plot(dih[:, 0], dih[:, 2], label="CG", linewidth=2.0)
ax.set_xlabel(r"$\phi$ (Deg)", fontsize=16, fontweight="normal")
ax.set_ylabel(r"P($\phi$)", fontsize=16, fontweight="normal")
ax.tick_params(labelsize=16, length=8.0, width=2.0)
ax.set_xticks(np.arange(-180, 181, 90))
ax.legend(fontsize=16)

panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
for k, ax in enumerate(axs.flat):
    ax.text(-0.07, 1.135, panel_labels[k], transform=ax.transAxes,
            fontsize=22, va="top", ha="left", clip_on=False)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout()
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "bond_ang_dih_pot_dist_combined.pdf", dpi=300, format="pdf")
