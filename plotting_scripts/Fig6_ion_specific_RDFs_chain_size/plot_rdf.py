"""Figure 6 (main text): ion-specific RDFs (PAA–ion, PAA–PAA, ion–ion) and
internal chain size <R^2(n)>^(1/2), AA vs CG, for Na+, K+ and Cs+ (salt-free).

Reads the CSV files in ../../full_data/Fig6_ion_specific_RDFs_chain_size and writes
../../figures/rdfs_NA_K_CS_2x2.pdf
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

here = Path(__file__).resolve().parent
data = here.parent.parent / "full_data" / "Fig6_ion_specific_RDFs_chain_size"
outdir = here.parent.parent / "figures"


def load(name):
    return np.loadtxt(data / name, delimiter=",", skiprows=1)


ions = ["Na", "K", "Cs"]
colors = {"Na": "tab:blue", "K": "tab:orange", "Cs": "tab:green"}
ls_aa = "-"    # AA solid
ls_cg = "--"   # CG dashed
r2 = pd.read_csv(data / "d_R2_of_n.csv")

# ---- 2x2 layout
fig, axs = plt.subplots(2, 2, figsize=(10, 9), sharey=False)
ax1, ax2 = axs[0]
ax3, ax4 = axs[1]

# --- (a)-(c) RDFs ---
for ax, panel, title, xlim in [(ax1, "a_PAA_ion", "PAA–ion", (0.25, 1.3)),
                               (ax2, "b_PAA_PAA", "PAA–PAA", (0.35, 2.0)),
                               (ax3, "c_ion_ion", "ion–ion", (0.25, 1.5))]:
    for ion in ions:
        c = colors[ion]
        aa = load(f"{panel}_rdf_{ion}_AA.csv")
        cg = load(f"{panel}_rdf_{ion}_CG.csv")
        ax.plot(aa[:, 0], aa[:, 1], color=c, linestyle=ls_aa, linewidth=2.0, label=f"{ion} AA")
        ax.plot(cg[:, 0], cg[:, 1], color=c, linestyle=ls_cg, linewidth=2.0, label=f"{ion} CG")
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("r (nm)", fontsize=16)
    ax.set_ylabel("RDF", fontsize=16)
    ax.set_xlim(*xlim)
    ax.tick_params(labelsize=16, length=8.0, width=2.0)

# --- (d) <R2(n)>^(1/2) ---
axins = ax4.inset_axes([0.5, 0.15, 0.34, 0.32])  # zoom on small n
for ion in ions:
    c = colors[ion]
    for tag, ls in [("AA", ls_aa), ("CG", ls_cg)]:
        y = np.sqrt(r2[f"R2_{ion}_{tag}_nm2"])
        ax4.plot(r2["n"], y, color=c, linestyle=ls, linewidth=2.0, label=f"{ion} {tag}")
        axins.plot(r2["n"], y, color=c, linestyle=ls, linewidth=1.5)

ax4.set_title("Chain size", fontsize=16)
ax4.set_xlabel("# of CG bonds, n", fontsize=16)
ax4.set_ylabel(r"$\langle R^2(n)\rangle^{1/2}$ (nm)", fontsize=16)
ax4.tick_params(labelsize=16, length=8.0, width=2.0)
axins.set_xlim(0, 5)
axins.set_ylim(0.2, 1.1)
axins.tick_params(labelsize=10, length=4, width=1.2)
for spine in axins.spines.values():
    spine.set_linewidth(1.2)

# --- Thicken spines ---
for ax in (ax1, ax2, ax3, ax4):
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

# --- One shared legend for the whole figure ---
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=6, fontsize=15, frameon=False)

for ax, label in zip((ax1, ax2, ax3, ax4), ["(a)", "(b)", "(c)", "(d)"]):
    ax.text(-0.12, 1.02, label, transform=ax.transAxes, fontsize=26, va="bottom", ha="left")

# Leave space for the shared legend
plt.tight_layout(rect=[0, 0, 1, 0.92])
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "rdfs_NA_K_CS_2x2.pdf", dpi=300, format="pdf")
