import matplotlib.pyplot as plt
import numpy as np

# ======================= helper =======================

def coord_number(r, g, rho):
    """CN(r) = 4*pi*rho * integral g(r) r^2 dr (cumulative trapezoid)."""
    r = np.asarray(r, float)
    g = np.asarray(g, float)
    idx = np.argsort(r)
    r, g = r[idx], g[idx]
    m = np.isfinite(r) & np.isfinite(g)
    r, g = r[m], g[m]
    if r.size and r[0] > 0:
        r = np.insert(r, 0, 0.0)
        g = np.insert(g, 0, 0.0)
    integrand = g * r**2
    dr = np.diff(r)
    area = 0.5 * (integrand[1:] + integrand[:-1]) * dr
    cn = 4.0 * np.pi * rho * np.concatenate(([0.0], np.cumsum(area)))
    return r, cn

# ======================= data sets =======================
base_aa = "/work/skurahatti/multivalent/2d-gel/PAA_dang/SC_CGM"
base_cg = "/work/skurahatti/multivalent/2d-gel/PAA_dang/SC_CGM/CG2"

N_ions = 20
Lx_nm  = 32.0048 * 0.1
Ly_nm  = 32.0048 * 0.1
Lz_nm  = 32.0    * 0.1
rho    = N_ions / (Lx_nm * Ly_nm * Lz_nm)  # nm^-3

datasets = [
    {
        "label":      r"$L_z$ = 3.2 nm",
        "targets":    f"{base_aa}/LZ_3.2/rdf_com_NA.xvg",
        "implicit":   f"{base_cg}/LZ_3.2/implicit/rdf_M_NA.xvg",
        "explicit":   f"{base_cg}/LZ_3.2/rdf_PAA_NA.dat",
        "xlim":       (0.2, 1.6),
        "xticks":     [0.2, 0.5, 1.0, 1.5],
        "rdf_ylim":   (-0.5, 11),
        "rdf_yticks": [0, 2, 4, 6, 8, 10],
        "cn_ylim":    (-0.5, 20),
        "cn_yticks":  [0, 5, 10, 15, 20],
    },
    {
        "label":      r"$L_z$ = 3.4 nm",
        "targets":    f"{base_aa}/LZ_3.4/M-NA.dist.new",
        "implicit":   f"{base_cg}/LZ_3.4/implicit/rdf_M_NA.xvg",
        "explicit":   f"{base_cg}/LZ_3.4/rdf_PAA_NA.dat",
        "xlim":       (0.2, 1.6),
        "xticks":     [0.2, 0.5, 1.0, 1.5],
        "rdf_ylim":   (-0.5, 12),
        "rdf_yticks": [0, 2, 4, 6, 8, 10, 12],
        "cn_ylim":    (-0.5, 20),
        "cn_yticks":  [0, 5, 10, 15, 20],
    },
    {
        "label":      r"$L_z$ = 3.8 nm",
        "targets":    f"{base_aa}/LZ_3.8/M-NA.dist.new",
        "implicit":   f"{base_cg}/LZ_3.8/implicit/rdf_M_NA.xvg",
        "explicit":   f"{base_cg}/LZ_3.8/rdf_PAA_NA.dat",
        "xlim":       (0.2, 1.7),
        "xticks":     [0.2, 0.5, 1.0, 1.5],
        "rdf_ylim":   (-0.5, 15),
        "rdf_yticks": [0, 3, 6, 9, 12, 15],
        "cn_ylim":    (-0.5, 20),
        "cn_yticks":  [0, 5, 10, 15, 20],
    },
]

# ======================= figure =======================
panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

fig, axes = plt.subplots(3, 2, figsize=(12, 14), sharex=False)

for row, ds in enumerate(datasets):
    targets  = np.loadtxt(ds["targets"],  usecols=(0, 1), comments=('#', '@'))
    implicit = np.loadtxt(ds["implicit"], comments=('#', '@'), usecols=(0, 1))
    explicit = np.loadtxt(ds["explicit"], usecols=(1, 2), skiprows=4)

    # data already in nm (A2NM = 1.0)
    targets_r  = targets[:, 0]
    implicit_r = implicit[:, 0]
    explicit_r = explicit[:, 0]

    rt, cnt = coord_number(targets_r,  targets[:, 1],  rho)
    re, cne = coord_number(explicit_r, explicit[:, 1], rho)
    ri, cni = coord_number(implicit_r, implicit[:, 1], rho)

    ax_rdf = axes[row, 0]
    ax_cn  = axes[row, 1]

    # --- RDF ---
    ax_rdf.plot(targets_r,  targets[:, 1],
                color="black", marker="o", markersize=3, linewidth=1.8, label="AA")
    ax_rdf.plot(explicit_r, explicit[:, 1], linewidth=2.0, label="CG")
    ax_rdf.plot(implicit_r, implicit[:, 1], linewidth=2.0, label="CG(implicit)")

    ax_rdf.set_xlim(ds["xlim"])
    ax_rdf.set_xticks(ds["xticks"])
    ax_rdf.set_ylim(ds["rdf_ylim"])
    ax_rdf.set_yticks(ds["rdf_yticks"])
    ax_rdf.set_ylabel("RDF", fontsize=14, fontfamily="serif", fontweight="semibold")
    ax_rdf.tick_params(labelsize=13, length=6, width=1.5)
    ax_rdf.set_title(ds["label"], fontsize=14, fontfamily="serif", fontweight="semibold")

    # --- CN ---
    ax_cn.plot(rt, cnt, color="black", marker="o", markersize=3, linewidth=1.8, label="AA")
    ax_cn.plot(re, cne, linewidth=2.0, label="CG")
    ax_cn.plot(ri, cni, linewidth=2.0, label="CG(implicit)")

    ax_cn.set_xlim(ds["xlim"])
    ax_cn.set_xticks(ds["xticks"])
    ax_cn.set_ylim(ds["cn_ylim"])
    ax_cn.set_yticks(ds["cn_yticks"])
    ax_cn.set_ylabel("CN(r)", fontsize=14, fontfamily="serif", fontweight="semibold")
    ax_cn.tick_params(labelsize=13, length=6, width=1.5)
    ax_cn.set_title(ds["label"], fontsize=14, fontfamily="serif", fontweight="semibold")

    # legend only on top row
    if row == 0:
        ax_rdf.legend(fontsize=12)
        ax_cn.legend(fontsize=12, loc="upper left")

    # x-labels only on bottom row
    if row == 2:
        ax_rdf.set_xlabel("r (nm)", fontsize=14, fontfamily="serif", fontweight="semibold")
        ax_cn.set_xlabel("r (nm)",  fontsize=14, fontfamily="serif", fontweight="semibold")

    # panel labels
    ax_rdf.text(-0.04, 1.08, panel_labels[row * 2],
                transform=ax_rdf.transAxes, fontsize=18, va="top", ha="left", clip_on=False)
    ax_cn.text(-0.04, 1.08,  panel_labels[row * 2 + 1],
               transform=ax_cn.transAxes,  fontsize=18, va="top", ha="left", clip_on=False)

    # thicken spines
    for ax in (ax_rdf, ax_cn):
        for spine in ax.spines.values():
            spine.set_linewidth(2.0)

fig.suptitle("PAA–Na RDF and Coordination Number",
             fontsize=18, fontfamily="serif", fontweight="semibold", y=0.98)
plt.tight_layout(rect=[0., 0.0, 1, 0.96])
plt.savefig("combined_rdf_cn.pdf", dpi=300, format="pdf")
plt.show()
