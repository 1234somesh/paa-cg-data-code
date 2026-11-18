import matplotlib.pyplot as plt
import numpy as np

# ===================== Paths =====================
# RDF at step_(n+1) corresponds to the potential trained through step_n
# (that potential is what drove the step_(n+1) simulation), same convention
# used in plot_potentials_with_coulomb.py.
BASE_51 = "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_5_1/Try13"
BASE_60 = "/work/skurahatti/multivalent/2d-gel/PAA_dang_ECC/parametrization/box_l_6nm/Try13"

def load_rdf_series(base, final_potential_step):
    final_rdf_step = final_potential_step + 1
    steps = [
        ("n=0",                     f"{base}/step_01"),
        ("n=2",                     f"{base}/step_03"),
        ("n=8",                     f"{base}/step_09"),
        ("n=15",                    f"{base}/step_16"),
        (f"n={final_potential_step}", f"{base}/step_{final_rdf_step:02d}"),
    ]
    rdf_MNA, rdf_MM, rdf_NN = [], [], []
    for label, step_dir in steps:
        rdf_MNA.append((label, np.loadtxt(f"{step_dir}/rdf_PAA_NA.dat", skiprows=5)))
        rdf_MM.append((label,  np.loadtxt(f"{step_dir}/rdf_PAA_PAA.dat", skiprows=5)))
        rdf_NN.append((label,  np.loadtxt(f"{step_dir}/rdf_NA_NA.dat", skiprows=5)))
    target_MNA = np.loadtxt(f"{base}/targets/PAA-NA.dist.dat", usecols=(0, 1))
    target_MM  = np.loadtxt(f"{base}/targets/PAA-PAA.dist.dat", usecols=(0, 1))
    target_NN  = np.loadtxt(f"{base}/targets/NA-NA.dist.dat", usecols=(0, 1))
    return rdf_MNA, rdf_MM, rdf_NN, target_MNA, target_MM, target_NN

# box_l_5_1: last potential snapshot is step_77 -> RDF at step_78
data_51 = load_rdf_series(BASE_51, final_potential_step=77)
# box_l_6nm: last potential snapshot is step_78 -> RDF at step_79
data_60 = load_rdf_series(BASE_60, final_potential_step=78)

def plot_rdf_row(axes_row, data, row_title):
    rdf_MNA, rdf_MM, rdf_NN, target_MNA, target_MM, target_NN = data
    ax_MNA, ax_MM, ax_NN = axes_row

    ax_MNA.plot(target_MNA[:, 0]*0.1, target_MNA[:, 1],
                color="black", marker="o", markerfacecolor="none", markersize=6.5, linewidth=0.5, label="Target")
    for label, rdf in rdf_MNA:
        ax_MNA.plot(rdf[:, 1]*0.1, rdf[:, 2], label=label, linewidth=1.8)
    ax_MNA.set_title(f"PAA–Na ({row_title})", fontsize=16, fontweight="normal")
    ax_MNA.set_ylabel("RDF", fontsize=16, fontweight="normal")
    ax_MNA.set_xlim(0.2, 1.5)
    ax_MNA.tick_params(labelsize=16, length=8.0, width=2.0)
    ax_MNA.legend(fontsize=12)

    ax_MM.plot(target_MM[:, 0]*0.1, target_MM[:, 1],
               color="black", marker="o", markerfacecolor="none", markersize=6.5, linewidth=0.5, label="Target")
    for label, rdf in rdf_MM:
        ax_MM.plot(rdf[:, 1]*0.1, rdf[:, 2], label=label, linewidth=1.8)
    ax_MM.set_title(f"PAA–PAA ({row_title})", fontsize=16, fontweight="normal")
    ax_MM.tick_params(labelsize=16, length=8.0, width=2.0)
    ax_MM.legend(fontsize=12)

    ax_NN.plot(target_NN[:, 0]*0.1, target_NN[:, 1],
               color="black", marker="o", markerfacecolor="none", markersize=6.5, linewidth=0.5, label="Target")
    for label, rdf in rdf_NN:
        ax_NN.plot(rdf[:, 1]*0.1, rdf[:, 2], label=label, linewidth=1.8)
    ax_NN.set_title(f"Na–Na ({row_title})", fontsize=16, fontweight="normal")
    ax_NN.tick_params(labelsize=16, length=8.0, width=2.0)
    ax_NN.legend(fontsize=12)

# ===================== Plot: 2 rows (box sizes) x 3 columns (pairs) =====================
fig, axes = plt.subplots(2, 3, figsize=(15, 8.4))

plot_rdf_row(axes[0], data_51, "L=5.1 nm")
plot_rdf_row(axes[1], data_60, "L=6.0 nm")

for ax in axes[1]:
    ax.set_xlabel("r (nm)", fontsize=16, fontweight="normal")

# ---------- Panel labels (a)-(f) ----------
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
k = 0
for i in range(2):
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
plt.savefig("RDF_evolution_box_5_1_6_0.pdf", dpi=300, format="pdf", bbox_inches="tight")
