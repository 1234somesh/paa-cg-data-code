#!/usr/bin/env python3
"""
Single-chain cell-gel model: cylindrical RDF and cumulative coordination
number (CN), AA (solid) vs CG (dashed), one column per ion.

Layout (2x3):
  col 0 = Na, col 1 = K, col 2 = Cs
  top row    = g_cyl(r)  (RDF)
  bottom row = CN(r)     (cumulative)
Vertical dotted line at the cell boundary r_b = 3.105 nm.

Equivalent to, per panel:
  xmgrace <ion>/sc_cgm_5.12/cyl_rdf_aa.dat <ion>/sc_cgm_5.12/cyl_rdf_cg.dat
and the same with "rdf" -> "cn" for the CN column.
Data files are produced by each ion's own cylindrical_analysis/<ion>/sc_cgm_5.12/plot_aa_cg.py.
"""

import numpy as np
import matplotlib.pyplot as plt

BASE = "/work/skurahatti/multivalent/2d-gel/PSS-protocol-1site/amber/figures_Dang/cylindrical_analysis"

ion_dirs = {
    "Na": f"{BASE}/sc_cgm_5.12",
    "K":  f"{BASE}/K/sc_cgm_5.12",
    "Cs": f"{BASE}/Cs/sc_cgm_5.12",
}
colors = {"Na": "tab:blue", "K": "tab:orange", "Cs": "tab:green"}

data = {
    ion: dict(
        rdf_aa=np.loadtxt(f"{d}/cyl_rdf_aa.dat"),
        rdf_cg=np.loadtxt(f"{d}/cyl_rdf_cg.dat"),
        cn_aa=np.loadtxt(f"{d}/cyl_cn_aa.dat"),
        cn_cg=np.loadtxt(f"{d}/cyl_cn_cg.dat"),
    )
    for ion, d in ion_dirs.items()
}

r_b = 3.105  # nm, cell boundary

ions_order = ["Na", "K", "Cs"]
fig, axs = plt.subplots(2, 3, figsize=(16, 9), sharex="col")

for col, ion in enumerate(ions_order):
    d = data[ion]
    c = colors[ion]
    ax_rdf, ax_cn = axs[:, col]

    ax_rdf.plot(d["rdf_aa"][:, 0], d["rdf_aa"][:, 1], color=c, ls="-", lw=2.0, label="AA")
    ax_rdf.plot(d["rdf_cg"][:, 0], d["rdf_cg"][:, 1], color=c, ls="--", lw=2.0, label="CG")
    ax_cn.plot(d["cn_aa"][:, 0], d["cn_aa"][:, 1], color=c, ls="-", lw=2.0, label="AA")
    ax_cn.plot(d["cn_cg"][:, 0], d["cn_cg"][:, 1], color=c, ls="--", lw=2.0, label="CG")

    ax_rdf.axvline(r_b, color="gray", ls=":", lw=1.5, label=rf"$r_b = {r_b}$ nm")
    ax_cn.axvline(r_b, color="gray", ls=":", lw=1.5)

    ax_rdf.set_title(f"{ion}$^+$", fontsize=17, fontfamily="serif", fontweight="semibold")
    ax_cn.set_xlabel("r (nm)", fontsize=16, fontfamily="serif", fontweight="semibold")
    if col == 0:
        ax_rdf.set_ylabel(r"$g_{\mathrm{cyl}}(r)$", fontsize=16, fontfamily="serif", fontweight="semibold")
        ax_cn.set_ylabel("CN(r)", fontsize=16, fontfamily="serif", fontweight="semibold")

    for ax in (ax_rdf, ax_cn):
        ax.tick_params(labelsize=13, length=7.0, width=1.5)
        ax.legend(fontsize=12)
        for spine in ax.spines.values():
            spine.set_linewidth(2.0)

for i, ax in enumerate(axs.flat):
    ax.text(-0.02, 1.02, f"({'abcdef'[i]})", transform=ax.transAxes,
            fontsize=20, va="bottom", ha="right")

plt.tight_layout()
plt.savefig("cyl_rdf_cn_NA_K_CS.pdf", dpi=300, bbox_inches="tight")
plt.show()

print("Saved to cyl_rdf_cn_NA_K_CS.pdf")
