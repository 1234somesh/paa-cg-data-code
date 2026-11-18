"""Figure 7 (main text): radius of gyration of a PAA 20-mer vs salt concentration
for Na+, K+ and Cs+; AA, CG (salt-free potentials) and CG reparametrized with salt.

Reads ../../full_data/Fig7_Rg_vs_csalt/Rg_vs_csalt.csv and writes
../../figures/Rg_vs_csalt.pdf
"""
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

here = Path(__file__).resolve().parent
data = pd.read_csv(here.parent.parent / "full_data" / "Fig7_Rg_vs_csalt" / "Rg_vs_csalt.csv")
outdir = here.parent.parent / "figures"

# ============================================================
# Shared style
# ============================================================
ions = ["Na", "K", "Cs"]
ion_titles = {"Na": "Na$^+$", "K": "K$^+$", "Cs": "Cs$^+$"}
fmts = {"Na": "o", "K": "s", "Cs": "d"}

capsize = 3
lw = 2.0
ms = 8
title_fs = 18
label_fs = 18
tick_fs = 18
legend_fs = 15

xticks_csalt = [0.0, 0.1, 0.2, 0.5]
xlim_csalt = (-0.02, 0.55)

# ============================================================
# Figure: 1 row × 3 cols
# ============================================================
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

for col, ion in enumerate(ions):
    ax = axs[col]
    d = data[data["ion"] == ion]
    re = d.dropna(subset=["Rg_CG_reparam_nm"])
    fmt = fmts[ion]
    ax.errorbar(d["csalt_M"], d["Rg_AA_nm"], yerr=d["Rg_AA_err_nm"],
                fmt=fmt, mfc='none', mec='C0', ecolor='C0',
                linestyle='-', linewidth=lw, markersize=ms, capsize=capsize, label='AA')
    ax.errorbar(d["csalt_M"], d["Rg_CG_nm"], yerr=d["Rg_CG_err_nm"],
                fmt=fmt, mfc='C1', mec='C1', ecolor='C1',
                linestyle='-', linewidth=lw, markersize=ms, capsize=capsize, label='CG')
    ax.errorbar(re["csalt_M"], re["Rg_CG_reparam_nm"], yerr=re["Rg_CG_reparam_err_nm"],
                fmt=fmt, mfc='none', mec='C1', ecolor='C1', color='C1',
                linestyle='--', linewidth=lw, markersize=ms, capsize=capsize,
                label='CG (reparam)')
    ax.set_title(ion_titles[ion], fontsize=title_fs, fontfamily="serif")
    ax.set_xlabel(r"$c_\mathrm{salt}$ (M)", fontsize=label_fs, fontfamily="serif")
    ax.set_xticks(xticks_csalt)
    ax.set_xlim(xlim_csalt)
    ax.set_ylim(0.82, 1.12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=tick_fs, length=6.0, width=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
    if col == 0:
        ax.legend(fontsize=legend_fs)

axs[0].set_ylabel(r"$R_g$ (nm)", fontsize=label_fs, fontfamily="serif")

for col, label in enumerate(["(a)", "(b)", "(c)"]):
    axs[col].text(-0.12, 1.02, label, transform=axs[col].transAxes,
                  fontsize=26, va="bottom", ha="left")

fig.tight_layout()
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "Rg_vs_csalt.pdf", dpi=300, format="pdf", bbox_inches="tight")
plt.close()
