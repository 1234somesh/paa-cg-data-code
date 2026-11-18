import matplotlib.pyplot as plt
import numpy as np

def load_xvg(path):
    # Some of these files get re-run by gmx analyze as more seeds finish and
    # end up with several "&"-separated blocks appended in one file instead
    # of being overwritten; take the last (most complete) block.
    with open(path) as f:
        text = f.read()
    last_block = text.split("\n&\n")[-1]
    lines = [l for l in last_block.splitlines() if l and not l[0] in "#@&"]
    return np.loadtxt(lines, usecols=(0, 1))

# Na: ECC (electronic continuum correction) run, box_l_6nm
M_NA_cg  = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/step_60/rdf_PAA_NA.dat",
    skiprows=5
)
M_NA_aa  = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/targets/PAA-NA.dist.dat",
    usecols=(0, 1)
)

# K: non-ECC (regular charges), K_last
M_K_cg   = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/K_last/Try13/step_80/rdf_PAA_K.dat",
    skiprows=5
)
M_K_aa   = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/K_last/Try13/targets/PAA-K.dist.dat",
    usecols=(0, 1)
)

# Cs: non-ECC (regular charges), Cs_last
M_CS_cg  = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/Cs_last/Try13/step_84/rdf_PAA_CS.dat",
    skiprows=5
)
M_CS_aa  = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/Cs_last/Try13/targets/PAA-CS.dist.dat",
    usecols=(0, 1)
)

M_M_NA_cg = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/step_60/rdf_PAA_PAA.dat",
    skiprows=5
)
M_M_NA_aa = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/targets/PAA-PAA.dist.dat",
    usecols=(0, 1)
)

M_M_K_cg = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/K_last/Try13/step_80/rdf_PAA_PAA.dat",
    skiprows=5
)
M_M_K_aa = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/K_last/Try13/targets/PAA-PAA.dist.dat",
    usecols=(0, 1)
)

M_M_CS_cg = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/Cs_last/Try13/step_84/rdf_PAA_PAA.dat",
    skiprows=5
)
M_M_CS_aa = np.loadtxt(
    "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/Cs_last/Try13/targets/PAA-PAA.dist.dat",
    usecols=(0, 1)
)

# --- Rg distributions (AA always available; CG only for Cs so far) ---
Rg_NA_aa = load_xvg("/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang_ECC/NA/box_l_6_0/rg_dist.xvg")
Rg_K_aa  = load_xvg("/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/K/box_l_6_0_last/rg_dist.xvg")
Rg_CS_aa = load_xvg("/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/CS/box_l_6_last/rg_dist.xvg")

Rg_NA_cg = load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13/step_60/rg_dist.xvg")
Rg_K_cg  = load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/K_last/Try13/step_80/rg_dist.xvg")
Rg_CS_cg = load_xvg("/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/Cs_last/Try13/step_84/rg_dist.xvg")

# --- Styling
ls_aa = "-"    # AA solid
ls_cg = "--"   # CG dashed

ions = {
    "Na": {"aa": M_NA_aa, "cg": M_NA_cg, "aa_m": M_M_NA_aa, "cg_m": M_M_NA_cg,
           "rg_aa": Rg_NA_aa, "rg_cg": Rg_NA_cg},
    "K":  {"aa": M_K_aa,  "cg": M_K_cg,  "aa_m": M_M_K_aa,  "cg_m": M_M_K_cg,
           "rg_aa": Rg_K_aa,  "rg_cg": Rg_K_cg},
    "Cs": {"aa": M_CS_aa, "cg": M_CS_cg, "aa_m": M_M_CS_aa, "cg_m": M_M_CS_cg,
           "rg_aa": Rg_CS_aa, "rg_cg": Rg_CS_cg},
}

colors = {"Na": "tab:blue", "K": "tab:orange", "Cs": "tab:green"}

# ---- 1x3 layout (PAA-ion, PAA-PAA, Rg distribution)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)
ax1, ax2, ax3 = axs

# --- (1) PAA–ion RDF ---
for ion, d in ions.items():
    c = colors[ion]
    # AA target
    ax1.plot(
        d["aa"][:, 0] * 0.1,
        d["aa"][:, 1],
        color=c, linestyle=ls_aa,
        linewidth=2.0,
        label=f"{ion} AA"
    )
    # CG RDF (VOTCA-style: r in col 1 (Å), g(r) in col 2)
    ax1.plot(
        d["cg"][:, 1] * 0.1,
        d["cg"][:, 2],
        color=c, linestyle=ls_cg,
        linewidth=2.0,
        label=f"{ion} CG"
    )

ax1.set_title("PAA–ion", fontsize=16,
              fontfamily="serif", fontweight="semibold")
ax1.set_xlabel("r (nm)", fontsize=16,
               fontfamily="serif", fontweight="semibold")
ax1.set_ylabel("RDF", fontsize=16,
               fontfamily="serif", fontweight="semibold")
ax1.set_xlim(0.25, 1.3)
ax1.tick_params(labelsize=16, length=8.0, width=2.0)

# --- (2) PAA–PAA RDF ---
for ion, d in ions.items():
    c = colors[ion]
    ax2.plot(
        d["aa_m"][:, 0] * 0.1,
        d["aa_m"][:, 1],
        color=c, linestyle=ls_aa,
        linewidth=2.0,
        label=f"{ion} AA"
    )
    ax2.plot(
        d["cg_m"][:, 1] * 0.1,
        d["cg_m"][:, 2],
        color=c, linestyle=ls_cg,
        linewidth=2.0,
        label=f"{ion} CG"
    )

ax2.set_title("PAA–PAA", fontsize=16,
              fontfamily="serif", fontweight="semibold")
ax2.set_xlabel("r (nm)", fontsize=16,
               fontfamily="serif", fontweight="semibold")
ax2.set_ylabel("RDF", fontsize=16,
               fontfamily="serif", fontweight="semibold")
ax2.set_xlim(0.35, 2.0)
ax2.tick_params(labelsize=16, length=8.0, width=2.0)

# --- (3) Rg distribution ---
for ion, d in ions.items():
    c = colors[ion]
    ax3.plot(
        d["rg_aa"][:, 0],
        d["rg_aa"][:, 1],
        color=c, linestyle=ls_aa,
        linewidth=2.0,
        label=f"{ion} AA"
    )
    ax3.plot(
        d["rg_cg"][:, 0],
        d["rg_cg"][:, 1],
        color=c, linestyle=ls_cg,
        linewidth=2.0,
        label=f"{ion} CG"
    )

ax3.set_title("Rg distribution", fontsize=16,
              fontfamily="serif", fontweight="semibold")
ax3.set_xlabel(r"$R_g$ (nm)", fontsize=16,
               fontfamily="serif", fontweight="semibold")
ax3.set_ylabel(r"P($R_g$)", fontsize=16,
               fontfamily="serif", fontweight="semibold")
ax3.tick_params(labelsize=16, length=8.0, width=2.0)
ax3.set_xlim(0.5, 1.75)
# --- Thicken spines ---
for ax in (ax1, ax2, ax3):
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

# --- One shared legend for the whole figure ---
handles, labels = ax1.get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc="upper center", ncol=6,
    fontsize=15, frameon=False
)

for ax, label in zip(axs, ["(a)", "(b)", "(c)"]):
    ax.text(
        -0.12, 1.02,          # x, y in axes-fraction coords
        label,
        transform=ax.transAxes,
        fontsize=26,
        va="bottom",
        ha="left"
    )

# Leave space for the shared legend
plt.tight_layout(rect=[0, 0, 1, 0.85])
plt.savefig("rdfs_NA_K_CS_1x3.pdf", dpi=300, format="pdf")
