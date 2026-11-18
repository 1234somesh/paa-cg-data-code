import matplotlib.pyplot as plt
import numpy as np

# ===================== Paths / constants =====================
# ECC (electronic continuum correction): formal unit charges are kept in
# paa_cg.data (q = -1.0 / +1.0), and the Coulomb interaction is rescaled by
# raising the dielectric constant instead of scaling the charges directly:
#   eps_r = eps_water / lambda^2 = 71.9 / 0.75^2 = 127.822222
eps_r = 127.822222  # dielectric constant
q_PAA = -1.0  # monomer charge
q_NA  = +1.0  # Na+ charge

SYSTEMS = [
    ("L=4.3 nm", "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_4_3/Try13", 49),
    ("L=5.1 nm", "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_5_1/Try13", 77),
    ("L=6.0 nm", "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13", 59),
]

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

def coulomb(r_nm, eps_r=eps_r):
    r = np.array(r_nm, dtype=float)
    eps0 = 8.85e-12       # C^2 / (N m^2)
    q    = 1.602e-19      # C
    NA   = 6.022e23       # 1/mol
    V = NA * q**2 * 1e-3 / (4.0 * np.pi * eps_r * eps0 * r * 1e-9)
    return V  # kJ/mol, for q1=q2=+1

def full_potential_from_table(tab, q1, q2):
    r_nm    = tab[:, 1] * 0.1      # Å -> nm
    V_short = tab[:, 2] * 4.184    # kcal/mol -> kJ/mol
    V_coul_unit = coulomb(r_nm)    # kJ/mol for q1=q2=+1
    V_full = V_short + (q1 * q2) * V_coul_unit
    return r_nm, V_full

# ===================== Load + build full potentials =====================
data = []
for label, base, step in SYSTEMS:
    tab_MNA = load_table(f"{base}/step_{step:02d}/updated_potentials/pair_PAA_NA.table")
    tab_MM  = load_table(f"{base}/step_{step:02d}/updated_potentials/pair_PAA_PAA.table")
    r_MNA, V_MNA = full_potential_from_table(tab_MNA, q_PAA, q_NA)
    r_MM,  V_MM  = full_potential_from_table(tab_MM,  q_PAA, q_PAA)
    data.append((label, step, r_MNA, V_MNA, r_MM, V_MM))

# ===================== Plot: 1 row x 2 columns (PAA-Na, PAA-PAA), one curve per L =====================
fig, (ax_MNA, ax_MM) = plt.subplots(1, 2, figsize=(10, 4.6))

for label, step, r_MNA, V_MNA, r_MM, V_MM in data:
    ax_MNA.plot(r_MNA, V_MNA, linewidth=2.0, label=label)
    ax_MM.plot(r_MM, V_MM, linewidth=2.0, label=label)

# The tabulated pair potentials carry an r->0 guard point (added so LAMMPS
# never crashes on a rare close-contact GCMC insertion) whose Coulomb
# contribution diverges as 1/r; it's a numerical safety net, not physical
# data, so it must be cropped out of the plotted range like the other
# figures_Dang_ECC full-potential plots do.
ax_MNA.set_title("PAA–Na", fontsize=16, fontweight="normal")
ax_MNA.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
ax_MNA.set_ylabel(r"$\mathrm{PMF}(r)$ (kJ/mol)", fontsize=14, fontweight="normal")
ax_MNA.set_xlim(0.2, 1.5)
ax_MNA.set_ylim(-6.5, 5.)
ax_MNA.tick_params(labelsize=16, length=8.0, width=2.0)
ax_MNA.legend(fontsize=12)

ax_MM.set_title("PAA–PAA", fontsize=16, fontweight="normal")
ax_MM.set_xlabel("r (nm)", fontsize=16, fontweight="normal")
ax_MM.set_xlim(0.35, 2.1)
ax_MM.set_ylim(-0.2, 4.0)
ax_MM.tick_params(labelsize=16, length=8.0, width=2.0)
ax_MM.legend(fontsize=12)

# ---------- Panel labels (a)-(b) ----------
for ax, lab in zip((ax_MNA, ax_MM), ['(a)', '(b)']):
    ax.text(
        -0.06, 1.125, lab,
        transform=ax.transAxes,
        fontsize=22, va="top", ha="left",
        clip_on=False
    )

# ---------- Thicken all spines ----------
for ax in (ax_MNA, ax_MM):
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout(pad=1.2, w_pad=2.0)
plt.savefig("Vfull_PAA_NA_PAA_PAA_box_sizes.pdf", dpi=300, format="pdf", bbox_inches="tight")
