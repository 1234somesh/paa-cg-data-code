#!/usr/bin/env python3
"""
Root-mean-squared distance between monomers separated by n bonds,
sqrt(<R^2(n)>), comparing CG and atomistic data.
"""
import numpy as np
import matplotlib.pyplot as plt

CG_PATH = "R2_of_n_lammps_all.csv"
AA_PATH = "/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/K/box_l_6_0_last/R2_of_n_atomistic.csv"

n_cg, r2_cg = np.loadtxt(CG_PATH, delimiter=",", skiprows=1, unpack=True)
n_aa, r2_aa = np.loadtxt(AA_PATH, delimiter=",", skiprows=1, unpack=True)

fig, ax = plt.subplots(figsize=(6.5, 5))
ax.plot(n_aa, np.sqrt(r2_aa), "-o", linewidth=2.0, markersize=6, label="AA")
ax.plot(n_cg, np.sqrt(r2_cg), "-o", linewidth=2.0, markersize=6, label="CG")

ax.set_xlabel("n (bonds)", fontsize=15)
ax.set_ylabel(r"$\langle R^2(n) \rangle^{1/2}$ (nm)", fontsize=15)
ax.tick_params(labelsize=13, length=6, width=1.3)
ax.legend(fontsize=13)
for spine in ax.spines.values():
    spine.set_linewidth(1.3)

plt.tight_layout()
plt.savefig("rms_R_vs_n.pdf", dpi=300, format="pdf")
plt.show()
