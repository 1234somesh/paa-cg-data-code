#!/usr/bin/env python3
"""
Figures for the free-chain spherical analysis (run sph_posm_analysis.py first):

  posm_bar_sph_NA_K_CS.pdf  -- CVT osmotic pressure at r_beta = L/2, Na/K/Cs,
                               AA vs CG, with block-averaged error bars
  sph_rdf_cn_NA_K_CS.pdf    -- g_sph(r) and CN(r) from the chain CoM, AA vs CG

Everything is read from the files written by sph_posm_analysis.py; nothing is
recomputed here, so the bars and the error bars describe the same estimator.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
IONS = ["Na", "K", "Cs"]
ION_TITLE = {"Na": "Na$^+$", "K": "K$^+$", "Cs": "Cs$^+$"}

summary = json.loads((HERE / "sph_posm_summary.json").read_text())
missing = [f"{i} {t}" for i in IONS for t in ("AA", "CG")
           if t not in summary.get(i, {})]
if missing:
    raise SystemExit(f"missing results for: {', '.join(missing)} -- "
                     f"run sph_posm_analysis.py first")

# ============================================================
# Figure 1: P_osm bars with error bars
# ============================================================
aa_vals = [summary[i]["AA"]["P_osm"] for i in IONS]
cg_vals = [summary[i]["CG"]["P_osm"] for i in IONS]
aa_errs = [summary[i]["AA"]["P_osm_err"] for i in IONS]
cg_errs = [summary[i]["CG"]["P_osm_err"] for i in IONS]

print("=== free chain, spherical CVT at r_beta = L/2 ===")
for i in IONS:
    for t in ("AA", "CG"):
        e = summary[i][t]
        print(f"  {i} {t}: P_osm = {e['P_osm']:.4f} +/- {e['P_osm_err']:.4f} bar "
              f"({100*e['P_osm_err']/e['P_osm']:.1f}%)  "
              f"[{e['n_frames']} frames, rho(r_beta) = {e['rho_mean']:.4e} ions/A^3]")

x = np.arange(len(IONS))
width = 0.35

fig, ax = plt.subplots(figsize=(7, 5.5))
bars_aa = ax.bar(x - width / 2, aa_vals, width, yerr=aa_errs, capsize=5,
                 error_kw=dict(elinewidth=1.8, capthick=1.8, ecolor="black"),
                 color="black", alpha=0.75, edgecolor="black", linewidth=1.2,
                 label="AA")
bars_cg = ax.bar(x + width / 2, cg_vals, width, yerr=cg_errs, capsize=5,
                 error_kw=dict(elinewidth=1.8, capthick=1.8, ecolor="black"),
                 color="red", alpha=0.75, edgecolor="black", linewidth=1.2,
                 label="CG")

for b, v, e in zip(list(bars_aa) + list(bars_cg), aa_vals + cg_vals, aa_errs + cg_errs):
    ax.text(b.get_x() + b.get_width() / 2, v + e, f"{v:.3f}",
            ha="center", va="bottom", fontsize=12, fontfamily="serif")

r_beta_nm = summary[IONS[0]]["AA"]["r_beta_nm"]
ax.set_xticks(x)
ax.set_xticklabels([ION_TITLE[i] for i in IONS], fontsize=16, fontfamily="serif")
ax.set_ylabel(r"$P_{\mathrm{osm}}$ (bar)", fontsize=16, fontfamily="serif",
              fontweight="semibold")
ax.set_title(rf"Free chain, CVT (spherical), $r_b = 0.5\,L = {r_beta_nm:.3f}$ nm",
             fontsize=14, fontfamily="serif")
ax.legend(fontsize=13)
ax.tick_params(labelsize=14, length=8.0, width=2.0)
for spine in ax.spines.values():
    spine.set_linewidth(2.0)

plt.tight_layout()
plt.savefig(HERE / "posm_bar_sph_NA_K_CS.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("\nSaved posm_bar_sph_NA_K_CS.pdf")

# ============================================================
# Figure 2: g_sph(r) and CN(r), AA vs CG, one column per ion
# ============================================================
# sharey per row: with independently auto-scaled log axes the curve HEIGHTS
# are not comparable between panels, which reads as (for example) K sitting
# above Cs at r_beta when the numbers say the opposite. The g(r_beta) value
# is also printed in each panel so the boundary ordering can be read off
# without measuring against the axis at all.
fig, axs = plt.subplots(2, 3, figsize=(16, 9), sharey="row")

for col, ion in enumerate(IONS):
    r_aa, g_aa = np.loadtxt(HERE / f"sph_rdf_{ion}_aa.dat", unpack=True)
    r_cg, g_cg = np.loadtxt(HERE / f"sph_rdf_{ion}_cg.dat", unpack=True)
    _, cn_aa = np.loadtxt(HERE / f"sph_cn_{ion}_aa.dat", unpack=True)
    _, cn_cg = np.loadtxt(HERE / f"sph_cn_{ion}_cg.dat", unpack=True)
    n_ion = summary[ion]["AA"]["n_ion"]

    ax = axs[0, col]
    ax.plot(r_aa, g_aa, "k-", lw=2.0, label="AA")
    ax.plot(r_cg, g_cg, "r--", lw=2.0, label="CG")
    ax.axvline(r_beta_nm, color="gray", ls=":", lw=1.5,
               label=r"$r_b = L/2$")
    ax.set_yscale("log")
    ax.set_title(ION_TITLE[ion], fontsize=17, fontfamily="serif")
    ax.set_xlabel("r (nm)", fontsize=15, fontfamily="serif")
    if col == 0:
        ax.set_ylabel(r"$g_{\mathrm{sph}}(r)$", fontsize=16, fontfamily="serif")
    g_beta_aa = summary[ion]["AA"]["rho_mean"] / summary[ion]["AA"]["rho_bulk"]
    g_beta_cg = summary[ion]["CG"]["rho_mean"] / summary[ion]["CG"]["rho_bulk"]
    ax.text(0.03, 0.06,
            rf"$g(r_b)$: AA {g_beta_aa:.3f}   CG {g_beta_cg:.3f}",
            transform=ax.transAxes, fontsize=12, fontfamily="serif")
    ax.legend(fontsize=12)

    ax = axs[1, col]
    ax.plot(r_aa, cn_aa, "k-", lw=2.0, label="AA")
    ax.plot(r_cg, cn_cg, "r--", lw=2.0, label="CG")
    ax.axvline(r_beta_nm, color="gray", ls=":", lw=1.5)
    ax.axhline(n_ion, color="gray", ls="-.", lw=1.2, label=f"{n_ion} ions")
    ax.set_xlabel("r (nm)", fontsize=15, fontfamily="serif")
    if col == 0:
        ax.set_ylabel("CN(r)", fontsize=16, fontfamily="serif")
    ax.legend(fontsize=12)

for i, ax in enumerate(axs.flat):
    ax.tick_params(labelsize=13, length=6.0, width=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    ax.text(-0.02, 1.02, f"({'abcdef'[i]})", transform=ax.transAxes,
            fontsize=17, fontfamily="serif", ha="right", va="bottom")

plt.tight_layout()
plt.savefig(HERE / "sph_rdf_cn_NA_K_CS.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved sph_rdf_cn_NA_K_CS.pdf")

# ============================================================
# Figure 3: all three ions on ONE axis around r_beta
# ============================================================
# The per-ion panels above cannot be compared by eye across columns even
# with a shared axis, because each curve sits in its own panel. This puts
# them on a single linear axis over the approach to r_beta, which is the
# comparison the CVT actually rests on: the ion that condenses least leaves
# the most counterions at the boundary and so sits highest here.
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
colors = {"Na": "tab:blue", "K": "tab:orange", "Cs": "tab:green"}

for ax, (tag, low) in zip(axes, (("AA", "aa"), ("CG", "cg"))):
    for ion in IONS:
        r, g = np.loadtxt(HERE / f"sph_rdf_{ion}_{low}.dat", unpack=True)
        g_beta = summary[ion][tag]["rho_mean"] / summary[ion][tag]["rho_bulk"]
        ax.plot(r, g, "-", color=colors[ion], lw=2.0,
                label=rf"{ION_TITLE[ion]}   $g(r_b)$={g_beta:.3f}")
        ax.plot([r_beta_nm], [g_beta], "o", color=colors[ion], ms=8,
                mec="black", mew=1.0, zorder=5)
    ax.axvline(r_beta_nm, color="gray", ls=":", lw=1.5)
    ax.set_xlim(1.5, 3.2)
    ax.set_ylim(0, 2.0)
    ax.set_xlabel("r (nm)", fontsize=15, fontfamily="serif")
    ax.set_title(tag, fontsize=16, fontfamily="serif")
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=13, length=6.0, width=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)

axes[0].set_ylabel(r"$g_{\mathrm{sph}}(r)$", fontsize=16, fontfamily="serif")
plt.tight_layout()
plt.savefig(HERE / "sph_gr_boundary_NA_K_CS.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved sph_gr_boundary_NA_K_CS.pdf")
