"""Figure 4 (main text): (a) stretching free energy of a PAA 20-mer (Na+) along the
end-to-end distance, AA vs CG vs Kremer-Grest; (b) Rg vs N with power-law fits.

Reads the CSV files in ../../full_data/Fig4_stretching_free_energy_Rg_scaling and writes
../../figures/FE_scaling.pdf. The fitted exponents nu are printed to stdout.
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pathlib import Path

here = Path(__file__).resolve().parent
data = here.parent.parent / "full_data" / "Fig4_stretching_free_energy_Rg_scaling"
outdir = here.parent.parent / "figures"


def load(name):
    return np.loadtxt(data / name, delimiter=",", skiprows=1)


# ===================== Panel (a) data: stretching free energy (PMF) =====================
data_ref = load("a_PMF_AA.csv")
data2 = load("a_PMF_CG.csv")
data_prim = load("a_PMF_Kremer_Grest.csv")

# ===================== Panel (b) data: Rg vs N scaling =====================
rg = load("b_Rg_vs_N.csv")
N, Rg_aa, Rg_aa_err, Rg_imp = rg[:, 0], rg[:, 1], rg[:, 2], rg[:, 3]

mask = (N >= 9) & (N <= 25)
N_fit = N[mask]
Rg_imp_fit = Rg_imp[mask]
Rg_aa_fit = Rg_aa[mask]
Rg_aa_err_fit = Rg_aa_err[mask]

x = np.log(N_fit)
y_cg = np.log(Rg_imp_fit)
y_aa = np.log(Rg_aa_fit)
sigma_y_aa = Rg_aa_err_fit / Rg_aa_fit


def linfit_unweighted(x, y):
    X = np.vstack([x, np.ones_like(x)]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    m, b = beta
    resid = y - X @ beta
    dof = len(x) - 2
    s2 = np.sum(resid**2) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    return m, b, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])


def linfit_weighted(x, y, sigma):
    X = np.vstack([x, np.ones_like(x)]).T
    w = 1.0 / (sigma**2)
    XtW = X.T * w
    cov0 = np.linalg.inv(XtW @ X)
    beta = cov0 @ (XtW @ y)
    m, b = beta
    resid = y - X @ beta
    chi2 = np.sum((resid / sigma)**2)
    dof = len(x) - 2
    chi2_red = chi2 / dof
    cov = cov0 * chi2_red
    return m, b, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]), chi2_red


nu_cg, lnA_cg, nu_cg_err, lnA_cg_err = linfit_unweighted(x, y_cg)
nu_aa, lnA_aa, nu_aa_err, lnA_aa_err, chi2red_aa = linfit_weighted(x, y_aa, sigma_y_aa)

print(f"nu_aa = {nu_aa:.4f} +/- {nu_aa_err:.4f}")
print(f"nu_cg = {nu_cg:.4f} +/- {nu_cg_err:.4f}")

# ===================== Plot: 1 row x 2 columns =====================
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
ax0, ax3 = axes   # (a): PMF, (b): Rg vs N scaling

# ---------- Panel (a): PMF vs R_ee (each curve shifted to its minimum) ----------
pos1, pmf1 = data_prim[:, 0], data_prim[:, 1] - data_prim[:, 1].min()
pos2, pmf2 = data2[:, 0], data2[:, 1] - data2[:, 1].min()
pos_ref, pmf_ref = data_ref[:, 0], data_ref[:, 1] - data_ref[:, 1].min()
yerr = 0
yerr_ref = 0

ax0.plot(pos_ref, pmf_ref, linewidth=2.0, label="AA", zorder=3)
ax0.fill_between(pos_ref, pmf_ref - yerr_ref, pmf_ref + yerr_ref, alpha=0.25, linewidth=0, zorder=2)
ax0.plot(pos2, pmf2, linewidth=2.0, label="CG", zorder=3)
ax0.fill_between(pos2, pmf2 - yerr, pmf2 + yerr, alpha=0.25, linewidth=0, zorder=2)
ax0.plot(pos1, pmf1, linewidth=2.0, label="Kremer-Grest", zorder=3)
ax0.fill_between(pos1, pmf1 - yerr, pmf1 + yerr, alpha=0.25, linewidth=0, zorder=2)

ax0.set_xlabel(r"$R_{ee}$", fontsize=16, fontweight="normal")
ax0.set_ylabel("Stretch. free energy (kJ/mol)", fontsize=16, fontweight="normal")
ax0.set_xlim([1.9, 4.3])
ax0.set_ylim([-2, 13])
ax0.tick_params(labelsize=16, length=8.0, width=2.0)
ax0.legend(fontsize=14, loc="upper left")

# ---------- Panel (b): Rg vs N with scaling fits ----------
aa_err = ax3.errorbar(N_fit, Rg_aa_fit, yerr=Rg_aa_err_fit, fmt="o", capsize=3, linewidth=2.0,
                      label=f"AA (ν={nu_aa:.3f}±{nu_aa_err:.2f})")
aa_handle = aa_err[0]
cg_line, = ax3.plot(N_fit, Rg_imp_fit, "o", linewidth=2.0,
                    label=f"CG (ν={nu_cg:.3f}±{nu_cg_err:.2f})")

ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlabel("# of monomers, N", fontsize=16, fontweight="normal")
ax3.set_ylabel(r"$\langle R_g \rangle$ (nm)", fontsize=16, fontweight="normal")
ax3.tick_params(labelsize=16, length=8.0, width=2.0)

ax3.xaxis.set_major_locator(mticker.FixedLocator([5, 10, 15, 20, 25]))
ax3.xaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
ax3.yaxis.set_major_locator(mticker.FixedLocator([0.4, 0.8, 1.2]))
ax3.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
ax3.xaxis.set_minor_locator(mticker.NullLocator())
ax3.yaxis.set_minor_locator(mticker.NullLocator())

# fit lines span the full x-axis range
xlim = ax3.get_xlim()
Nfit = np.linspace(xlim[0], xlim[1], 300)
ax3.plot(Nfit, np.exp(lnA_aa) * Nfit**nu_aa, "-", color="tab:blue", linewidth=2.0)
ax3.plot(Nfit, np.exp(lnA_cg) * Nfit**nu_cg, "-", color="tab:orange", linewidth=2.0)
ax3.set_xlim(xlim)

ax3.text(0.1, 0.76, rf"$R_g \sim N^{{{nu_aa:.3f}}}$", transform=ax3.transAxes,
         fontsize=14, color="tab:blue", va="top", ha="left")
ax3.text(0.1, 0.68, rf"$R_g \sim N^{{{nu_cg:.3f}}}$", transform=ax3.transAxes,
         fontsize=14, color="tab:orange", va="top", ha="left")
ax3.legend([aa_handle, cg_line], ["AA", "CG"], fontsize=14, loc="upper left")

# ---------- Panel labels and spines ----------
for ax, label in zip(axes, ['(a)', '(b)']):
    ax.text(-0.04, 1.125, label, transform=ax.transAxes,
            fontsize=22, va="top", ha="left", clip_on=False)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout(pad=1.2, h_pad=3.0, w_pad=2.0)
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "FE_scaling.pdf", dpi=300, format="pdf", bbox_inches="tight")
