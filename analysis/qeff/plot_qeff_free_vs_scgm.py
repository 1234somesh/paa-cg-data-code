#!/usr/bin/env python3
"""
Combines free-chain/plot.py's qeff_vs_j.pdf (potentials/NA/plot_qeff.py's copy)
and SC-CGM/plot_both.py's qeff_rc0p72nm_AA_vs_CG_scaled.pdf into one 1x2 figure:
q_eff(r_c=0.72 nm, j) for AA vs CG, free chain (a) next to the single-chain
cell-gel model (b).

CG uses scaled charges (e.g. +/-0.117851 in the LAMMPS dump); both panels
convert CG q_eff -> "e" units by the same factor 1/0.117851, as in
SC-CGM/plot_both.py and free-chain's own version of that script.
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q_CG_MONO = 0.117851  # |charge| on a CG PAA bead in the LAMMPS model
SCALE = 1.0 / Q_CG_MONO

PANELS = [
    ("(a)", "Free chain", HERE / "free-chain"),
    ("(b)", "Single-chain cell-gel model", HERE / "SC-CGM"),
]

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)

for ax, (label, title, d) in zip(axes, PANELS):
    aa = np.loadtxt(d / "qeff_rc0p72nm_vs_j.dat")
    cg = np.loadtxt(d / "qeff_rc0p72nm_vs_j_CG_LAMMPS.dat")

    j_aa, qeff_aa = aa[:, 0], aa[:, 1]
    j_cg, qeff_cg = cg[:, 0], cg[:, 1]
    qeff_cg_scaled = qeff_cg * SCALE

    if aa.shape[1] >= 3:
        aa_artist = ax.errorbar(j_aa, qeff_aa, yerr=aa[:, 2], marker="o",
                                linewidth=2.0, capsize=3, label="AA")
    else:
        aa_artist, = ax.plot(j_aa, qeff_aa, marker="o", linewidth=2.0, label="AA")
    cg_artist, = ax.plot(j_cg, qeff_cg_scaled, marker="o", linewidth=2.0, label="CG")
    # legend(): mixing an errorbar container with a plain Line2D reorders
    # handles by artist type (all Line2D before all containers) regardless
    # of creation order -- pass the order explicitly so AA is always first.

    ax.set_title(title, fontsize=15, fontfamily="serif")
    ax.set_xlabel("Bead index $j$", fontsize=16, fontfamily="serif",
                  fontweight="semibold")
    ax.tick_params(labelsize=14, length=6.0, width=2.0)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
    ax.legend([aa_artist, cg_artist], ["AA", "CG"], fontsize=13)
    ax.text(-0.06, 1.12, label, transform=ax.transAxes,
            fontsize=20, va="top", ha="left", clip_on=False)

axes[0].set_ylabel(r"$q_{\mathrm{eff}}(r_c=0.72\ \mathrm{nm},\, j)$  (e)",
                   fontsize=16, fontfamily="serif", fontweight="semibold")

plt.tight_layout(pad=1.2, w_pad=2.5)
plt.savefig(HERE / "qeff_free_vs_scgm.pdf", dpi=300, bbox_inches="tight")
print(f"CG scale factor = {SCALE:.6f}")
print(f"Saved {HERE / 'qeff_free_vs_scgm.pdf'}")
