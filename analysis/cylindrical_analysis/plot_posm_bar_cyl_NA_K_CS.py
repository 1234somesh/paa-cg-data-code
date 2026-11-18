#!/usr/bin/env python3
"""
Single-chain cell-gel model (cylindrical SC-CGM): Contact-Value-Theorem
osmotic pressure (P_osm) for Na/K/Cs, AA vs CG, as grouped bars.

These RDFs are reproduced by the cell-gel model constrained at Lx=Ly=
62.10120 Ang (6.21012 nm), Lz=60.0 Ang (6.0 nm), with 20 counterions --
see SC_CGM/LZ_freeze_{Na,K,Cs}/paa_cg.data, identical for all three ions
and for both the AA and CG data alike (that's where rho_bulk and r_beta
come from -- no independent per-system trajectory reads).

CVT: rho(r_beta) = g_cyl(r_beta) * rho_bulk, r_beta = 0.5*Lx = 3.105 nm,
read off the already-computed cylindrical g(r) for each system
(cyl_rdf_aa.dat / cyl_rdf_cg.dat), averaged over the last 10% of the
resolved range below r_beta (reduces bin noise; reports a flatness
diagnostic). P_osm = kB*T*rho(r_beta).
"""
import json
import os

import numpy as np
import matplotlib.pyplot as plt

A2NM = 0.1
WINDOW_FRAC = 0.10

# Statistical error bars come from posm_err_timeseries.py, which rebuilds
# rho(r_beta) frame by frame from the same trajectories and block-averages
# the resulting time series (consecutive frames are correlated, so the naive
# per-frame SEM would be far too small). Run that script first; without its
# summary file the bars are drawn without error bars, as before.
ERR_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "posm_err_summary.json")

T_KELVIN = 300.0
ANG3_TO_L = 1e-27
NA_AVO = 6.02214076e23
ANG3_TO_MOLAR = 1.0 / (NA_AVO * ANG3_TO_L)   # 1 ion/Ang^3 -> mol/L
R_GAS_L_BAR = 0.0831446                      # L*bar / (mol*K)

# shared cell-gel box, from SC_CGM/LZ_freeze_{Na,K,Cs}/paa_cg.data
# (identical xlo/xhi, ylo/yhi, zlo/zhi and ion count across all three)
LX = 62.10120
LY = 62.10120
LZ = 60.00000
N_ION = 20
RHO_BULK = N_ION / (LX * LY * LZ)   # ions/Ang^3
R_BETA_NM = 0.5 * LX * A2NM         # = 3.105 nm

CYL = "/work/skurahatti/multivalent/2d-gel/PSS-protocol-1site/amber/figures_Dang/cylindrical_analysis"

rdf_files = {
    "Na": dict(aa=f"{CYL}/sc_cgm_5.12/cyl_rdf_aa.dat",    cg=f"{CYL}/sc_cgm_5.12/cyl_rdf_cg.dat"),
    "K":  dict(aa=f"{CYL}/K/sc_cgm_5.12/cyl_rdf_aa.dat",  cg=f"{CYL}/K/sc_cgm_5.12/cyl_rdf_cg.dat"),
    "Cs": dict(aa=f"{CYL}/Cs/sc_cgm_5.12/cyl_rdf_aa.dat", cg=f"{CYL}/Cs/sc_cgm_5.12/cyl_rdf_cg.dat"),
}


def rho_at_rbeta(r_nm, g, rho_bulk, r_beta_nm, label=""):
    r_nm = np.asarray(r_nm)
    mask = r_nm <= r_beta_nm
    r = r_nm[mask]
    rho_r = g[mask] * rho_bulk

    window_start = r.max() - WINDOW_FRAC * (r.max() - r.min())
    win = r >= window_start
    rho_beta = rho_r[win].mean()
    slope, _ = np.polyfit(r[win], rho_r[win], 1)

    flat_ok = abs(slope) <= 0.05 * abs(rho_beta + 1e-12)
    print(f"  [{label}] r_beta = {r_beta_nm:.3f} nm | rho(r_beta) = {rho_beta:.5e} ions/Ang^3 | "
          f"slope = {slope:.3e} {'(flat)' if flat_ok else '<-- NOT FLAT, provisional'}", flush=True)
    return rho_beta


def posm_bar(rho_ang3, T=T_KELVIN):
    c_molar = rho_ang3 * ANG3_TO_MOLAR
    return c_molar * R_GAS_L_BAR * T, c_molar


# ============================================================
# Compute rho(r_beta) and P_osm for each ion, AA and CG
# ============================================================
print(f"rho_bulk = {RHO_BULK:.6e} ions/Ang^3 (N={N_ION}, box={LX:.4f}x{LY:.4f}x{LZ:.4f} Ang), "
      f"r_beta = {R_BETA_NM:.4f} nm\n", flush=True)

results = {}
for ion, paths in rdf_files.items():
    print(f"=== {ion} ===", flush=True)
    r_aa, g_aa = np.loadtxt(paths["aa"], unpack=True)
    r_cg, g_cg = np.loadtxt(paths["cg"], unpack=True)

    rho_beta_aa = rho_at_rbeta(r_aa, g_aa, RHO_BULK, R_BETA_NM, label=f"{ion} AA")
    rho_beta_cg = rho_at_rbeta(r_cg, g_cg, RHO_BULK, R_BETA_NM, label=f"{ion} CG")

    P_aa, c_aa = posm_bar(rho_beta_aa)
    P_cg, c_cg = posm_bar(rho_beta_cg)
    print(f"  P_osm (AA) = {P_aa:.4f} bar  ({c_aa*1000:.3f} mM equivalent)", flush=True)
    print(f"  P_osm (CG) = {P_cg:.4f} bar  ({c_cg*1000:.3f} mM equivalent)", flush=True)

    results[ion] = dict(P_aa=P_aa, P_cg=P_cg)

# ============================================================
# Bar plot: P_osm, grouped by ion, AA vs CG
# ============================================================
ions_order = ["Na", "K", "Cs"]
aa_vals = [results[i]["P_aa"] for i in ions_order]
cg_vals = [results[i]["P_cg"] for i in ions_order]

# ---- statistical errors from the per-frame time series -------------------
aa_errs = cg_errs = None
if os.path.exists(ERR_JSON):
    with open(ERR_JSON) as f:
        err = json.load(f)
    aa_errs = [err[i]["AA"]["P_osm_err"] for i in ions_order]
    cg_errs = [err[i]["CG"]["P_osm_err"] for i in ions_order]

    print("\n=== Statistical errors (block-averaged, from trajectories) ===", flush=True)
    for i in ions_order:
        for tag, key in (("AA", "P_aa"), ("CG", "P_cg")):
            e = err[i][tag]
            print(f"  {i} {tag}: P_osm = {results[i][key]:.4f} bar "
                  f"+/- {e['P_osm_err']:.4f} "
                  f"({100*e['P_osm_err']/results[i][key]:.1f}%)  "
                  f"[{e['n_frames']} frames; time-series mean "
                  f"{e['P_osm']:.4f} bar, RDF-based {results[i][key]:.4f} bar, "
                  f"diff {100*abs(e['P_osm']-results[i][key])/results[i][key]:.1f}%]",
                  flush=True)
else:
    print(f"\nNOTE: {ERR_JSON} not found -- run posm_err_timeseries.py to get "
          f"error bars; plotting bars without them.", flush=True)

x = np.arange(len(ions_order))
width = 0.35

fig, ax = plt.subplots(figsize=(7, 5.5))
bars_aa = ax.bar(x - width / 2, aa_vals, width, yerr=aa_errs, capsize=5,
                  error_kw=dict(elinewidth=1.8, capthick=1.8, ecolor="black"),
                  color="black", alpha=0.75,
                  edgecolor="black", linewidth=1.2, label="AA")
bars_cg = ax.bar(x + width / 2, cg_vals, width, yerr=cg_errs, capsize=5,
                  error_kw=dict(elinewidth=1.8, capthick=1.8, ecolor="black"),
                  color="red", alpha=0.75,
                  edgecolor="black", linewidth=1.2, label="CG")

all_errs = ((aa_errs + cg_errs) if aa_errs is not None
            else [0.0] * (len(aa_vals) + len(cg_vals)))
for b, v, e in zip(list(bars_aa) + list(bars_cg), aa_vals + cg_vals, all_errs):
    ax.text(b.get_x() + b.get_width() / 2, v + e, f"{v:.3f}",
            ha="center", va="bottom", fontsize=12, fontfamily="serif")

ax.set_xticks(x)
ax.set_xticklabels([f"{i}$^+$" for i in ions_order], fontsize=16, fontfamily="serif")
ax.set_ylabel(r"$P_{\mathrm{osm}}$ (bar)", fontsize=16, fontfamily="serif", fontweight="semibold")
ax.set_title(r"CVT (cylindrical), $r_b = 0.5\,L_x = 3.105$ nm", fontsize=14, fontfamily="serif")
ax.legend(fontsize=13)
ax.tick_params(labelsize=14, length=8.0, width=2.0)
for spine in ax.spines.values():
    spine.set_linewidth(2.0)

plt.tight_layout()
plt.savefig("posm_bar_cyl_NA_K_CS.pdf", dpi=300, bbox_inches="tight")

print("\nSaved to posm_bar_cyl_NA_K_CS.pdf")
