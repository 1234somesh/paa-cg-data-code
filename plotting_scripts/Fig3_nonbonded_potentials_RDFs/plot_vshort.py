"""Figure 3 (main text): short-range CG pair potentials (top row) and RDFs
(bottom row) for PAA with Na+ counterions, salt-free.

Reads the CSV files in ../../full_data/Fig3_nonbonded_potentials_RDFs and writes
../../figures/Monomer_Vshort_RDF_2x3.pdf
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

here = Path(__file__).resolve().parent
data = here.parent.parent / "full_data" / "Fig3_nonbonded_potentials_RDFs"
outdir = here.parent.parent / "figures"


def load(name):
    return np.loadtxt(data / name, delimiter=",", skiprows=1)


V_MNA0, V_MNA1 = load("a_PAA_Na_potential_initial.csv"), load("a_PAA_Na_potential_final.csv")
V_MM0, V_MM1 = load("b_PAA_PAA_potential_initial.csv"), load("b_PAA_PAA_potential_final.csv")
V_NN = load("c_Na_Na_potential.csv")

def rdfs(tag):
    return (load(f"{tag}_rdf_target_AA.csv"), load(f"{tag}_rdf_CG_initial.csv"),
            load(f"{tag}_rdf_CG_final.csv"))

fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=False, sharey=False)
axV_MNA, axV_MM, axV_NN = axes[0]   # row 1: potentials
axg_MNA, axg_MM, axg_NN = axes[1]   # row 2: RDFs

# ---------- ROW 1: SHORT-RANGE POTENTIALS (no Coulomb) ----------
axV_MNA.plot(V_MNA0[:, 0], V_MNA0[:, 1], label="initial", linewidth=2.0)
axV_MNA.plot(V_MNA1[:, 0], V_MNA1[:, 1], label="final", linewidth=2.0)
axV_MNA.set_xlim(0.265, 1.27)
axV_MNA.set_ylim(-1.0, 5.5)
axV_MNA.set_title("PAA–Na", fontsize=16, fontweight="normal")
axV_MNA.set_ylabel(r"$V_{\mathrm{short}}(r)$ (kJ/mol)", fontsize=16, fontweight="normal")
axV_MNA.tick_params(labelsize=16, length=8.0, width=2.0)
axV_MNA.legend(fontsize=12)

axV_MM.plot(V_MM0[:, 0], V_MM0[:, 1], label="initial", linewidth=2.0)
axV_MM.plot(V_MM1[:, 0], V_MM1[:, 1], label="final", linewidth=2.0)
axV_MM.set_xlim(0.4, 2.1)
axV_MM.set_ylim(-1.0, 4.0)
axV_MM.set_title("PAA–PAA", fontsize=16, fontweight="normal")
axV_MM.tick_params(labelsize=16, length=8.0, width=2.0)
axV_MM.legend(fontsize=12)

axV_NN.plot(V_NN[:, 0], V_NN[:, 1], label="initial", linewidth=2.0)
axV_NN.set_xlim(0.27, 1.3)
axV_NN.set_ylim(-4.0, 5.0)
axV_NN.set_title("Na–Na", fontsize=16, fontweight="normal")
axV_NN.tick_params(labelsize=16, length=8.0, width=2.0)

# ---------- ROW 2: RDFs ----------
target_style = dict(color="black", marker="o", markerfacecolor="none",
                    markersize=6.5, linewidth=0.5, label="target")

for ax, tag in [(axg_MNA, "d_PAA_Na"), (axg_MM, "e_PAA_PAA"), (axg_NN, "f_Na_Na")]:
    t, g0, g1 = rdfs(tag)
    ax.plot(t[:, 0], t[:, 1], **target_style)
    ax.plot(g0[:, 0], g0[:, 1], label="initial", linewidth=1.8)
    ax.plot(g1[:, 0], g1[:, 1], label="final", linewidth=1.8)
    ax.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
    ax.tick_params(labelsize=16, length=8.0, width=2.0)
    ax.legend(fontsize=12)

axg_MNA.set_xlim(0.24, 1.27)
axg_MNA.set_ylim(-1.0, 50.0)
axg_MNA.set_ylabel("RDF", fontsize=16, fontweight="normal")
axg_MM.set_xlim(0.4, 2.1)
axg_MM.set_ylim(-0.4, 12.5)
axg_NN.set_xlim(0.27, 1.3)

# ---------- Panel labels and spines ----------
for k, ax in enumerate(axes.flat):
    ax.text(-0.04, 1.125, f"({'abcdef'[k]})", transform=ax.transAxes,
            fontsize=22, va="top", ha="left", clip_on=False)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout(pad=1.2, h_pad=3.0, w_pad=2.0)
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "Monomer_Vshort_RDF_2x3.pdf", dpi=300, format="pdf", bbox_inches="tight")
