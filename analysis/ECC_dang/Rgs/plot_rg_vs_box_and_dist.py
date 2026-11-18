import math
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# ===================== Error estimator (identical to calc_avg.py) =====================
def autocorrelation(data, normalized=True):
    """
    Compute autocorrelation using FFT
    """
    nobs = len(data)
    corr_data = data - data.mean()
    n = 2**int(math.log(nobs, 2))
    corr_data = corr_data[:n]
    Frf = np.fft.fft(corr_data)
    acf = np.fft.ifft(Frf * np.conjugate(Frf))/corr_data.shape[0]
    if normalized:
        acf /= acf[0]
    acf = np.real(acf)
    return acf[:int(corr_data.shape[0]/2)]

def calc_error(data):
    """
    Error estimation for time series of simulation observables, taking into
    account that these series are correlated (which enhances the estimated
    statistical error).
    """
    acf = autocorrelation(data)
    tau_int = 0.5
    for i in range(len(acf)):
        tau_int += acf[i]
        if i >= 6 * tau_int:
            break
    data_mean = np.mean(data)
    if tau_int > 0.5:
        N_eff = len(data) / (2.0 * tau_int)
        stat_err = np.sqrt(np.var(data) / N_eff)
    else:
        stat_err = np.sqrt(np.var(data) / len(data))
    return data_mean, stat_err

# ===================== Paths =====================
AA_BASE = "/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang_ECC/NA"
AA_SYSTEMS = [
    (4.3, f"{AA_BASE}/box_l_4_3"),
    (5.1, f"{AA_BASE}/box_l_5_1"),
    (6.0, f"{AA_BASE}/box_l_6_0"),
    (7.5, f"{AA_BASE}/box_l_7_5"),
]

CG_SYSTEMS = [
    (4.3, "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_4_3/Try13/step_60", 60, False),
    (5.1, "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_5_1/Try13/step_78", 78, False),
    (6.0, "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/step_79", 79, False),
]

# Rg distribution panel only: skip L=7.5 nm here (still shown in the Rg-vs-L panel).
SKIP_DIST_L = {7.5}

# ===================== Reference (non-ECC) data, from
# figures_Dang/potentials/NA/plot_with_dist_and_msd.py panels (k) and (l) =====================
def load_xvg(path):
    return np.loadtxt(path, comments=["#", "@"], usecols=(0, 1))

# Panel (l): Rg vs box size -- these are the exact hardcoded values used in
# plot_with_dist_and_msd.py, not recomputed here.
REF_box_lengths = np.array([4.3, 5.1, 6.0, 7.5])
REF_aa_rg_box   = np.array([0.9507, 0.9716, 0.9708, 0.9712])
REF_aa_err_box  = np.array([0.006, 0.0088, 0.0078, 0.0064])
REF_cg_rg_box   = np.array([0.933, 0.9475, 0.9729, 0.9761])
REF_cg_err_box  = np.array([0.0025, 0.0021, 0.0024, 0.0021])

# Panel (k): Rg distribution for different box sizes (non-ECC AA/CG).
REF_aa_small  = load_xvg("/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/NA/box_l_4_3/rg_dist.xvg")
REF_cg_small  = load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/box_l_4_3/Try13/step_100/rg_dist.xvg")
REF_aa_medium = load_xvg("/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/NA/box_l_5_1/rg_dist.xvg")
REF_cg_medium = load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/box_l_5_1nm/rg_dist.xvg")
REF_aa_large  = load_xvg("/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/NA/box_l_6_0/rg_dist.xvg")
REF_cg_large  = load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/box_l_6nm/cs_0.0/rg_dist.xvg")

def load_rg_timeseries_aa(dirpath):
    # Rg time series lives in either Rg_super.xvg or Rg.xvg (or both); prefer
    # Rg_super.xvg (the combined/longer series) when present. Both are
    # GROMACS gyrate output, already in nm -- column 1 is the total Rg
    # (column 0 is time, columns 2+ are the Rgx/Rgy/Rgz components).
    for fname in ("Rg_super.xvg", "Rg.xvg"):
        path = os.path.join(dirpath, fname)
        if os.path.isfile(path):
            data = np.loadtxt(path, comments=("#", "@"))
            return data[:, 1]
    raise FileNotFoundError(f"no Rg time series (Rg_super.xvg/Rg.xvg) in {dirpath}")

def load_rg_timeseries_cg(dirpath):
    # LAMMPS fix output: Rg_PAA.dat (and its Rg_PAA.dat.bak backup) are in
    # Angstrom (native LAMMPS "units real"); Rg_PAA.xvg is a pre-converted
    # nm copy of the same data. Returns None if the run hasn't produced any
    # of these yet (still in progress).
    for fname, needs_ang_to_nm in (("Rg_PAA.dat", True), ("Rg_PAA.xvg", False), ("Rg_PAA.dat.bak", True)):
        path = os.path.join(dirpath, fname)
        if os.path.isfile(path):
            data = np.loadtxt(path, comments=("#", "@"))
            rg = data[:, 1]
            return rg * 0.1 if needs_ang_to_nm else rg
    return None

def load_rg_dist(path):
    # rg_dist.xvg (both AA and CG) is already in nm. Some of these files
    # (e.g. box_l_4_3's AA one) hold multiple "&"-separated xmgrace sets --
    # gmx analyze produced one distribution per input column (total Rg,
    # then RgX/RgY/RgZ). Only the first set is the total-Rg distribution
    # we want; the rest are the components we were told to skip.
    with open(path) as f:
        lines = f.readlines()
    first_block = []
    for line in lines:
        s = line.strip()
        if s.startswith(("#", "@")):
            continue
        if s == "&":
            break
        first_block.append(s)
    return np.loadtxt(first_block, usecols=(0, 1))

# ===================== Rg vs box size =====================
aa_L, aa_mean, aa_err = [], [], []
for L, d in AA_SYSTEMS:
    rg = load_rg_timeseries_aa(d)
    mean, err = calc_error(rg)
    aa_L.append(L); aa_mean.append(mean); aa_err.append(err)
    print(f"AA  L={L:4.1f} nm: Rg = {mean:.4f} +- {err:.4f} nm")

cg_L, cg_mean, cg_err, cg_pending = [], [], [], []
for L, d, step, force_placeholder in CG_SYSTEMS:
    if force_placeholder:
        print(f"CG  L={L:4.1f} nm (step_{step}): still running -- leaving as placeholder "
              f"(a partial Rg_PAA.dat already exists but the run isn't finished)")
        cg_pending.append(L)
        continue
    rg = load_rg_timeseries_cg(d)
    if rg is None:
        print(f"CG  L={L:4.1f} nm (step_{step}): still running -- no Rg data yet, leaving as placeholder")
        cg_pending.append(L)
        continue
    mean, err = calc_error(rg)
    cg_L.append(L); cg_mean.append(mean); cg_err.append(err)
    print(f"CG  L={L:4.1f} nm (step_{step}): Rg = {mean:.4f} +- {err:.4f} nm")

# ===================== Rg distributions =====================
aa_dist = {}
for L, d in AA_SYSTEMS:
    if L in SKIP_DIST_L:
        continue
    p = os.path.join(d, "rg_dist.xvg")
    if os.path.isfile(p):
        aa_dist[L] = load_rg_dist(p)

cg_dist = {}
for L, d, step, force_placeholder in CG_SYSTEMS:
    p = os.path.join(d, "rg_dist.xvg")
    if not force_placeholder and os.path.isfile(p):
        cg_dist[L] = load_rg_dist(p)
    else:
        print(f"CG  L={L:4.1f} nm (step_{step}): no rg_dist.xvg yet -- skipping in distribution panel")

# ===================== Plot: 2 rows x 2 columns =====================
# Row 1 (a,b): this ECC analysis. Row 2 (c,d): the corresponding panels
# reproduced from figures_Dang/potentials/NA/plot_with_dist_and_msd.py
# (non-ECC; panels (l) and (k) there). Each column means the same thing
# top-to-bottom (left = Rg vs box size, right = Rg distribution).
fig, axes = plt.subplots(2, 2, figsize=(13, 10.4))
ax_a, ax_b = axes[0]
ax_l, ax_k = axes[1]

ls_box = {4.3: "-", 5.1: "--", 6.0: ":", 7.5: "-."}

# ---------- (a): ECC Rg vs box size (AA and CG) ----------
ax_a.errorbar(aa_L, aa_mean, yerr=aa_err, fmt="s-", color="tab:blue",
              lw=2.0, ms=8, capsize=4, capthick=1.5, label="AA")
ax_a.errorbar(cg_L, cg_mean, yerr=cg_err, fmt="o--", color="tab:orange",
              lw=2.0, ms=8, capsize=4, capthick=1.5, label="CG")
for pL in cg_pending:
    ax_a.axvline(pL, color="tab:orange", linestyle=":", alpha=0.4)
    # Blended transform (x in data coords, y in axes-fraction) so this label
    # can't stretch the y-autoscale the way a data-coordinate annotation would.
    ax_a.text(pL, 0.05, "CG (still running)", transform=ax_a.get_xaxis_transform(),
              color="tab:orange", fontsize=10, ha="center", va="bottom", rotation=90)

ax_a.set_xlabel("Box length $L$ (nm)", fontsize=16, fontweight="normal")
ax_a.set_ylabel(r"$\langle R_g \rangle$ (nm)", fontsize=16, fontweight="normal")
ax_a.set_title("ECC", fontsize=14, fontweight="normal")
ax_a.tick_params(labelsize=16, length=8.0, width=2.0)
ax_a.legend(fontsize=14, frameon=False, loc="lower right")

# ---------- (b): ECC Rg distributions for different box sizes (AA and CG) ----------
for L, dist in aa_dist.items():
    ax_b.plot(dist[:, 0], dist[:, 1], color="tab:blue", ls=ls_box[L], lw=2.0)
for L, dist in cg_dist.items():
    ax_b.plot(dist[:, 0], dist[:, 1], color="tab:orange", ls=ls_box[L], lw=2.0)

ax_b.set_xlabel("$R_g$ (nm)", fontsize=16, fontweight="normal")
ax_b.set_ylabel("P($R_g$)", fontsize=16, fontweight="normal")
ax_b.set_title("ECC", fontsize=14, fontweight="normal")
ax_b.tick_params(labelsize=16, length=8.0, width=2.0)

color_handles = [
    Line2D([0], [0], color="tab:blue", lw=2.0, label="AA"),
    Line2D([0], [0], color="tab:orange", lw=2.0, label="CG"),
]
leg_color_b = ax_b.legend(handles=color_handles, fontsize=12, frameon=False, loc="upper right")
ax_b.add_artist(leg_color_b)

box_handles_b = [
    Line2D([0], [0], color="black", lw=2.0, ls=ls_box[L], label=f"$L={L}$ nm")
    for L in sorted(set(list(aa_dist.keys()) + list(cg_dist.keys())))
]
ax_b.legend(handles=box_handles_b, fontsize=11, frameon=False, loc="center right")

# ---------- (l): non-ECC Rg vs box size (reproduced from plot_with_dist_and_msd.py) ----------
ax_l.errorbar(REF_box_lengths, REF_aa_rg_box, yerr=REF_aa_err_box,
              fmt="s-", color="tab:blue", lw=2.0, ms=8, capsize=4, capthick=1.5, label="AA")
ax_l.errorbar(REF_box_lengths, REF_cg_rg_box, yerr=REF_cg_err_box,
              fmt="o--", color="tab:orange", lw=2.0, ms=8, capsize=4, capthick=1.5, label="CG")
ax_l.set_xlabel("Box length $L$ (nm)", fontsize=16, fontweight="normal")
ax_l.set_ylabel(r"$\langle R_g \rangle$ (nm)", fontsize=16, fontweight="normal")
ax_l.set_title("non-ECC", fontsize=14, fontweight="normal")
ax_l.tick_params(labelsize=16, length=8.0, width=2.0)
ax_l.legend(fontsize=14, frameon=False)

# ---------- (k): non-ECC Rg distributions for different box sizes ----------
ref_ls_box = {"small": "-", "medium": "--", "large": ":"}
ax_k.plot(REF_aa_small[:, 0],  REF_aa_small[:, 1],  color="tab:blue",   ls=ref_ls_box["small"],  lw=2.0)
ax_k.plot(REF_aa_medium[:, 0], REF_aa_medium[:, 1], color="tab:blue",   ls=ref_ls_box["medium"], lw=2.0)
ax_k.plot(REF_aa_large[:, 0],  REF_aa_large[:, 1],  color="tab:blue",   ls=ref_ls_box["large"],  lw=2.0)
ax_k.plot(REF_cg_small[:, 0],  REF_cg_small[:, 1],  color="tab:orange", ls=ref_ls_box["small"],  lw=2.0)
ax_k.plot(REF_cg_medium[:, 0], REF_cg_medium[:, 1], color="tab:orange", ls=ref_ls_box["medium"], lw=2.0)
ax_k.plot(REF_cg_large[:, 0],  REF_cg_large[:, 1],  color="tab:orange", ls=ref_ls_box["large"],  lw=2.0)
ax_k.set_xlabel("$R_g$ (nm)", fontsize=16, fontweight="normal")
ax_k.set_ylabel("P($R_g$)", fontsize=16, fontweight="normal")
ax_k.set_xlim(0.6, 1.7)
ax_k.set_title("non-ECC", fontsize=14, fontweight="normal")
ax_k.tick_params(labelsize=16, length=8.0, width=2.0)

leg_color_k = ax_k.legend(handles=color_handles, fontsize=12, frameon=False, loc="upper right")
ax_k.add_artist(leg_color_k)
ref_box_handles = [
    Line2D([0], [0], color="black", lw=2.0, ls=ref_ls_box["small"],  label="$L=4.30$ nm"),
    Line2D([0], [0], color="black", lw=2.0, ls=ref_ls_box["medium"], label="$L=5.10$ nm"),
    Line2D([0], [0], color="black", lw=2.0, ls=ref_ls_box["large"],  label="$L=6.00$ nm"),
]
ax_k.legend(handles=ref_box_handles, fontsize=11, frameon=False, loc="center right")

# ---------- Panel labels ----------
for ax, lab in zip((ax_a, ax_b, ax_l, ax_k), ["(a)", "(b)", "(c)", "(d)"]):
    ax.text(-0.06, 1.125, lab, transform=ax.transAxes, fontsize=22,
            va="top", ha="left", clip_on=False)

# ---------- Thicken all spines ----------
for ax in (ax_a, ax_b, ax_l, ax_k):
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout(pad=1.2, h_pad=3.5, w_pad=2.5)
plt.savefig("Rg_vs_box_and_dist_ECC.pdf", dpi=300, format="pdf", bbox_inches="tight")
