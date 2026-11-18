import matplotlib.pyplot as plt
import numpy as np

# ===================== Paths / constants =====================
BASE = "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_4_3/Try13"

# ECC (electronic continuum correction) run: formal unit charges are kept in
# paa_cg.data (q = -1.0 / +1.0), and the Coulomb interaction is rescaled by
# raising the dielectric constant instead of scaling the charges directly:
#   eps_r = eps_water / lambda^2 = 71.9 / 0.75^2 = 127.822222
# (as set by `dielectric 127.822222` in every step_*/in.lmp). Using q=+-1.0
# together with eps_r=127.822222 reproduces exactly the electrostatics used
# in the simulation; using q=+-0.75 together with eps_r=127.822222 would
# double count the charge-scaling correction.
eps_r = 127.822222  # dielectric constant
q_PAA = -1.0  # monomer charge
q_NA  = +1.0  # Na+ charge

def load_table(path):
    # Header length varies between tables in this dataset (some have extra
    # comment lines), so skip everything up through the "N <npoints>" line
    # instead of assuming a fixed number of header rows.
    with open(path) as f:
        lines = f.readlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("N "):
            start = i + 1
            break
    while start < len(lines) and not lines[start].strip():
        start += 1
    return np.loadtxt(lines[start:])

def load_xvg(path):
    # Some inputs (e.g. gmx gyrate with total + per-axis Rg columns) make
    # gmx analyze emit multiple "&"-separated data sets in one file; only
    # the first set (total Rg) is the distribution we want here.
    with open(path) as f:
        lines = f.read().split("\n&", 1)[0].splitlines()
    lines = [ln for ln in lines if ln and not ln.startswith(("#", "@"))]
    return np.loadtxt(lines, usecols=(0, 1))

def coulomb(r_nm, eps_r=eps_r):
    r = np.array(r_nm, dtype=float)
    eps0 = 8.85e-12       # C^2 / (N m^2)
    q    = 1.602e-19      # C
    NA   = 6.022e23       # 1/mol
    V = NA * q**2 * 1e-3 / (4.0 * np.pi * eps_r * eps0 * r * 1e-9)
    return V  # kJ/mol, for q1=q2=+1

def full_potential_from_table(tab, q1, q2):
    # r_nm is taken directly from the tabulated potential, so the Coulomb
    # term is only added up to the same cutoff as the short-range table.
    r_nm    = tab[:, 1] * 0.1      # Å -> nm
    V_short = tab[:, 2] * 4.184    # kcal/mol -> kJ/mol
    V_coul_unit = coulomb(r_nm)    # kJ/mol for q1=q2=+1
    V_full = V_short + (q1 * q2) * V_coul_unit
    return r_nm, V_full

# ===================== Load tabulated SHORT-RANGE potentials =====================
# PAA–Na (this run only has steps up to step_50, unlike the non-ECC Try13
# which went to step_99, so the snapshots are re-picked accordingly)
M_NA_idx0  = load_table(f"{BASE}/initial/pair_PAA_NA.table")
M_NA_idx5  = load_table(f"{BASE}/step_02/updated_potentials/pair_PAA_NA.table")
M_NA_idx11 = load_table(f"{BASE}/step_08/updated_potentials/pair_PAA_NA.table")
M_NA_idx23 = load_table(f"{BASE}/step_59/updated_potentials/pair_PAA_NA.table")

# PAA–PAA
M_M_idx0   = load_table(f"{BASE}/initial/pair_PAA_PAA.table")
M_M_idx5   = load_table(f"{BASE}/step_02/updated_potentials/pair_PAA_PAA.table")
M_M_idx11  = load_table(f"{BASE}/step_08/updated_potentials/pair_PAA_PAA.table")
M_M_idx23  = load_table(f"{BASE}/step_59/updated_potentials/pair_PAA_PAA.table")

# Na–Na
NA_NA_idx0 = load_table(f"{BASE}/step_02/updated_potentials/pair_NA_NA.table")

# ===================== Build FULL potentials (short-range table + Coulomb) =====================
# PAA–Na (q = -1, +1)
r_MNA0,  V_MNA0  = full_potential_from_table(M_NA_idx0,  q_PAA, q_NA)
r_MNA5,  V_MNA5  = full_potential_from_table(M_NA_idx5,  q_PAA, q_NA)
r_MNA11, V_MNA11 = full_potential_from_table(M_NA_idx11, q_PAA, q_NA)
r_MNA23, V_MNA23 = full_potential_from_table(M_NA_idx23, q_PAA, q_NA)

# PAA–PAA (q = -1, -1)
r_MM0,  V_MM0  = full_potential_from_table(M_M_idx0,  q_PAA, q_PAA)
r_MM5,  V_MM5  = full_potential_from_table(M_M_idx5,  q_PAA, q_PAA)
r_MM11, V_MM11 = full_potential_from_table(M_M_idx11, q_PAA, q_PAA)
r_MM23, V_MM23 = full_potential_from_table(M_M_idx23, q_PAA, q_PAA)

# Na–Na (q = +1, +1)
r_NANA0, V_NANA0 = full_potential_from_table(NA_NA_idx0, q_NA, q_NA)

# ===================== Load RDFs =====================
# RDF at step_(n+1) corresponds to the potential trained through step_n
# (that potential is what drove the step_(n+1) simulation), mirroring the
# non-ECC figures_Dang/potentials/NA/plot_with_dist_and_msd.py convention.
# PAA–Na
M_NA_idx0_rdf  = np.loadtxt(f"{BASE}/step_01/rdf_PAA_NA.dat", skiprows=5)
M_NA_idx5_rdf  = np.loadtxt(f"{BASE}/step_03/rdf_PAA_NA.dat", skiprows=5)
M_NA_idx11_rdf = np.loadtxt(f"{BASE}/step_09/rdf_PAA_NA.dat", skiprows=5)
M_NA_idx23_rdf = np.loadtxt(f"{BASE}/step_60/rdf_PAA_NA.dat", skiprows=5)
M_NA_target    = np.loadtxt(f"{BASE}/targets/PAA-NA.dist.dat", usecols=(0, 1))

# PAA–PAA
M_M_idx0_rdf  = np.loadtxt(f"{BASE}/step_01/rdf_PAA_PAA.dat", skiprows=5)
M_M_idx5_rdf  = np.loadtxt(f"{BASE}/step_03/rdf_PAA_PAA.dat", skiprows=5)
M_M_idx11_rdf = np.loadtxt(f"{BASE}/step_09/rdf_PAA_PAA.dat", skiprows=5)
M_M_idx23_rdf = np.loadtxt(f"{BASE}/step_60/rdf_PAA_PAA.dat", skiprows=5)
M_M_target    = np.loadtxt(f"{BASE}/targets/PAA-PAA.dist.dat", usecols=(0, 1))

# Na–Na
N_N_idx0_rdf  = np.loadtxt(f"{BASE}/step_01/rdf_NA_NA.dat", skiprows=5)
N_N_idx5_rdf  = np.loadtxt(f"{BASE}/step_03/rdf_NA_NA.dat", skiprows=5)
N_N_idx11_rdf = np.loadtxt(f"{BASE}/step_09/rdf_NA_NA.dat", skiprows=5)
N_N_idx23_rdf = np.loadtxt(f"{BASE}/step_60/rdf_NA_NA.dat", skiprows=5)
N_N_target    = np.loadtxt(f"{BASE}/targets/NA-NA.dist.dat", usecols=(0, 1))

# ===================== Row 3 data: Rg distributions (3 box sizes) =====================
AA_RG_BASE = "/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang_ECC/NA"
box_sizes = ["4.3", "5.1", "6.0"]
box_colors = {"4.3": "tab:blue", "5.1": "tab:orange", "6.0": "tab:green"}

Rg_aa_dist = {
    "4.3": load_xvg(f"{AA_RG_BASE}/box_l_4_3/rg_dist.xvg"),
    "5.1": load_xvg(f"{AA_RG_BASE}/box_l_5_1/rg_dist.xvg"),
    "6.0": load_xvg(f"{AA_RG_BASE}/box_l_6_0/rg_dist.xvg"),
}
Rg_cg_dist = {
    "4.3": load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_4_3/Try13/step_50/rg_dist.xvg"),
    "5.1": load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_5_1/Try13/step_78/rg_dist.xvg"),
    "6.0": load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/step_60/rg_dist.xvg"),
}

# ===================== Row 3 data: <R^2(n)> comparison (box_l_4_3) =====================
R2_BASE = "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_4_3/Try13/step_50"
n_aa_r2, r2_aa   = np.loadtxt(f"{R2_BASE}/R2_of_n_atomistic.csv",  delimiter=",", skiprows=1, usecols=(0, 1), unpack=True)
n_cg_r2, r2_cg   = np.loadtxt(f"{R2_BASE}/R2_of_n_lammps_all.csv", delimiter=",", skiprows=1, usecols=(0, 1), unpack=True)

# ===================== Row 3 data: q_eff(r_c=0.72nm, j) AA vs CG =====================
QEFF_BASE = "/work/skurahatti/multivalent/2d-gel/PSS-protocol-1site/amber/figures_Dang_ECC/qeff/free_chain"
j_aa_qeff, qeff_aa = np.loadtxt(f"{QEFF_BASE}/qeff_rc0p72nm_vs_j.dat", unpack=True)
j_cg_qeff, qeff_cg = np.loadtxt(f"{QEFF_BASE}/qeff_rc0p72nm_vs_j_CG_LAMMPS.dat", unpack=True)

# ===================== Plot: 3 rows x 3 columns =====================
fig, axes = plt.subplots(3, 3, figsize=(15, 12.6))
axV_MNA, axV_MM, axV_NN = axes[0]     # row 1: full potential (short-range + Coulomb)
axg_MNA, axg_MM, axg_NN = axes[1]     # row 2: RDFs
ax_rg, ax_r2, ax_qeff   = axes[2]     # row 3: Rg dist (box sizes), <R^2(n)>, q_eff(j)

# ---------- ROW 1: FULL POTENTIAL (Coulomb + short-range) ----------
# PAA–Na
axV_MNA.plot(r_MNA0,  V_MNA0,  label="n=0",  linewidth=2.0)
axV_MNA.plot(r_MNA5,  V_MNA5,  label="n=2",  linewidth=2.0)
axV_MNA.plot(r_MNA11, V_MNA11, label="n=8",  linewidth=2.0)
axV_MNA.plot(r_MNA23, V_MNA23, label="n=59", linewidth=2.0)
axV_MNA.set_xlim(0.2, 1.27)
axV_MNA.set_ylim(-6.5, 8)  # matches non-ECC (regular charges) plot for direct comparison
axV_MNA.set_title("PAA–Na", fontsize=16, fontweight="normal")
axV_MNA.set_ylabel(r"$\mathrm{PMF}(r)$ (kJ/mol)", fontsize=14, fontweight="normal")
axV_MNA.tick_params(labelsize=16, length=8.0, width=2.0)
axV_MNA.legend(fontsize=12)

# PAA–PAA
axV_MM.plot(r_MM0,  V_MM0,  label="n=0",  linewidth=2.0)
axV_MM.plot(r_MM5,  V_MM5,  label="n=2",  linewidth=2.0)
axV_MM.plot(r_MM11, V_MM11, label="n=8",  linewidth=2.0)
axV_MM.plot(r_MM23, V_MM23, label="n=59", linewidth=2.0)
axV_MM.set_xlim(0.35, 2.1)
axV_MM.set_ylim(-0.5, 8.0)  # matches non-ECC (regular charges) plot for direct comparison
axV_MM.set_title("PAA–PAA", fontsize=16, fontweight="normal")
axV_MM.tick_params(labelsize=16, length=8.0, width=2.0)
axV_MM.legend(fontsize=12)

# Na–Na
axV_NN.plot(r_NANA0, V_NANA0, linewidth=2.0)
axV_NN.set_xlim(0.2, 1.3)
axV_NN.set_ylim(-2.0, 8.0)  # matches non-ECC (regular charges) plot for direct comparison
axV_NN.set_title("Na–Na", fontsize=16, fontweight="normal")
axV_NN.tick_params(labelsize=16, length=8.0, width=2.0)

# ---------- ROW 2: RDFs ----------
# PAA–Na
axg_MNA.plot(M_NA_target[:, 0]*0.1, M_NA_target[:, 1],
             color="black", marker="o", markerfacecolor="none", markersize=6.5, linewidth=0.5, label="Target")
axg_MNA.plot(M_NA_idx0_rdf[:, 1]*0.1,  M_NA_idx0_rdf[:, 2],  label="n=0",  linewidth=1.8)
axg_MNA.plot(M_NA_idx5_rdf[:, 1]*0.1,  M_NA_idx5_rdf[:, 2],  label="n=2",  linewidth=1.8)
axg_MNA.plot(M_NA_idx11_rdf[:, 1]*0.1, M_NA_idx11_rdf[:, 2], label="n=8",  linewidth=1.8)
axg_MNA.plot(M_NA_idx23_rdf[:, 1]*0.1, M_NA_idx23_rdf[:, 2], label="n=59", linewidth=1.8)
axg_MNA.set_xlim(0.2, 1.27)
axg_MNA.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
axg_MNA.set_ylabel("RDF", fontsize=16, fontweight="normal")
axg_MNA.tick_params(labelsize=16, length=8.0, width=2.0)
axg_MNA.legend(fontsize=12)

# PAA–PAA
axg_MM.plot(M_M_target[:, 0]*0.1, M_M_target[:, 1],
            color="black", marker="o", markerfacecolor="none", markersize=6.5, linewidth=0.5, label="Target")
axg_MM.plot(M_M_idx0_rdf[:, 1]*0.1,  M_M_idx0_rdf[:, 2],  label="n=0",  linewidth=1.8)
axg_MM.plot(M_M_idx5_rdf[:, 1]*0.1,  M_M_idx5_rdf[:, 2],  label="n=2",  linewidth=1.8)
axg_MM.plot(M_M_idx11_rdf[:, 1]*0.1, M_M_idx11_rdf[:, 2], label="n=8",  linewidth=1.8)
axg_MM.plot(M_M_idx23_rdf[:, 1]*0.1, M_M_idx23_rdf[:, 2], label="n=59", linewidth=1.8)
axg_MM.set_xlim(0.35, 2.1)
axg_MM.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
axg_MM.tick_params(labelsize=16, length=8.0, width=2.0)
axg_MM.legend(fontsize=12)

# Na–Na
axg_NN.plot(N_N_target[:, 0]*0.1, N_N_target[:, 1],
            color="black", marker="o", markerfacecolor="none", markersize=6.5, linewidth=0.5, label="Target")
axg_NN.plot(N_N_idx0_rdf[:, 1]*0.1,  N_N_idx0_rdf[:, 2],  label="n=0",  linewidth=1.8)
axg_NN.plot(N_N_idx5_rdf[:, 1]*0.1,  N_N_idx5_rdf[:, 2],  label="n=2",  linewidth=1.8)
axg_NN.plot(N_N_idx11_rdf[:, 1]*0.1, N_N_idx11_rdf[:, 2], label="n=8",  linewidth=1.8)
axg_NN.plot(N_N_idx23_rdf[:, 1]*0.1, N_N_idx23_rdf[:, 2], label="n=59", linewidth=1.8)
axg_NN.set_xlim(0.2, 1.3)
axg_NN.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
axg_NN.tick_params(labelsize=16, length=8.0, width=2.0)
axg_NN.legend(fontsize=12)

# ---------- ROW 3: Rg distributions (box sizes) | <R^2(n)> | q_eff(j) ----------
# Rg distributions, AA vs CG, for 3 box sizes
for box in box_sizes:
    c = box_colors[box]
    aa_dist = Rg_aa_dist[box]
    cg_dist = Rg_cg_dist[box]
    ax_rg.plot(aa_dist[:, 0], aa_dist[:, 1], color=c, linestyle="-",
               linewidth=1.8, label=f"{box} nm AA")
    ax_rg.plot(cg_dist[:, 0], cg_dist[:, 1], color=c, linestyle="--",
               linewidth=1.8, label=f"{box} nm CG")
ax_rg.set_xlabel(r"$R_g$ (nm)", fontsize=16, fontweight="normal")
ax_rg.set_ylabel(r"P($R_g$)", fontsize=16, fontweight="normal")
ax_rg.tick_params(labelsize=16, length=8.0, width=2.0)
ax_rg.legend(fontsize=10)

# <R^2(n)> comparison, atomistic vs CG (box_l_4_3)
ax_r2.plot(n_aa_r2, np.sqrt(r2_aa), "o-", label="Atomistic", lw=1.8, ms=4)
ax_r2.plot(n_cg_r2, np.sqrt(r2_cg), "s-", label="CG", lw=1.8, ms=4)
ax_r2.set_xlabel("n", fontsize=16, fontweight="normal")
ax_r2.set_ylabel(r"$\langle R^2(n) \rangle^{1/2}$ (nm)", fontsize=16, fontweight="normal")
ax_r2.grid(True, alpha=0.3)
ax_r2.tick_params(labelsize=16, length=8.0, width=2.0)
ax_r2.legend(fontsize=12)

# q_eff(r_c=0.72nm, j), AA vs CG
ax_qeff.plot(j_aa_qeff, qeff_aa, marker="o", linewidth=2.0, label="AA")
ax_qeff.plot(j_cg_qeff, qeff_cg, marker="o", linewidth=2.0, label="CG")
ax_qeff.set_xlabel("Bead index $j$", fontsize=16, fontweight="normal")
ax_qeff.set_ylabel(r"$q_{\mathrm{eff}}(r_c=0.72\ \mathrm{nm},\, j)$  (e)",
                    fontsize=14, fontweight="normal")
ax_qeff.tick_params(labelsize=16, length=8.0, width=2.0)
ax_qeff.legend(fontsize=12)

# ---------- Panel labels (a)-(i) ----------
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
k = 0
for i in range(3):
    for j in range(3):
        axes[i, j].text(
            -0.04, 1.125, labels[k],
            transform=axes[i, j].transAxes,
            fontsize=22, va="top", ha="left",
            clip_on=False
        )
        k += 1

# ---------- Thicken all spines ----------
for row in axes:
    for ax in row:
        for spine in ax.spines.values():
            spine.set_linewidth(2.0)

plt.tight_layout(pad=1.2, h_pad=3.0, w_pad=2.0)
plt.savefig("Monomer_Vfull_RDF_ECC.pdf", dpi=300, format="pdf", bbox_inches="tight")
