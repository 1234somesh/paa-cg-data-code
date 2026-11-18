import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ===================== Paths =====================
BASE = "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/implicit/NA_4_3"

def load_votca(path):
    """Load a VOTCA IBI file (r, value, flag) — returns r and value columns."""
    return np.loadtxt(path, comments='#', usecols=(0, 1))

# ===================== Steps to display =====================
steps  = ["step_004", "step_008", "step_012", "step_016", "step_028"]
labels = ["n=4",      "n=8",      "n=12",     "n=16",     "n=21"]

# ===================== Load potentials (*-*.pot.cur) =====================
pot_MNA  = [load_votca(f"{BASE}/{s}/M-NA.pot.cur")  for s in steps]
pot_MM   = [load_votca(f"{BASE}/{s}/M-M.pot.cur")   for s in steps]
pot_NN   = [load_votca(f"{BASE}/{s}/NA-NA.pot.cur") for s in steps]

# ===================== Load RDFs (*-*.dist.new) + targets =====================
rdf_MNA  = [load_votca(f"{BASE}/{s}/M-NA.dist.new")  for s in steps]
rdf_MM   = [load_votca(f"{BASE}/{s}/M-M.dist.new")   for s in steps]
rdf_NN   = [load_votca(f"{BASE}/{s}/NA-NA.dist.new") for s in steps]

tgt_MNA  = load_votca(f"{BASE}/M-NA.dist.tgt")
tgt_MM   = load_votca(f"{BASE}/M-M.dist.tgt")
tgt_NN   = load_votca(f"{BASE}/NA-NA.dist.tgt")

# ===================== Plot: 2 rows x 3 columns =====================
fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=False, sharey=False)

axV_MNA, axV_MM, axV_NN = axes[0]   # row 1: potentials
axg_MNA, axg_MM, axg_NN = axes[1]   # row 2: RDFs

# ---------- ROW 1: POTENTIALS ----------
for i, (p, lbl) in enumerate(zip(pot_MNA, labels)):
    axV_MNA.plot(p[:, 0], p[:, 1], linewidth=2.0, label=lbl)
axV_MNA.set_title("M–NA",   fontsize=16)
axV_MNA.set_ylabel(r"$V(r)$ (kJ/mol)", fontsize=16)
axV_MNA.set_xlim(0.25, 2.1)
axV_MNA.set_ylim(-5.0, 5.0)
axV_MNA.tick_params(labelsize=14, length=6, width=1.5)
axV_MNA.legend(fontsize=12)

for i, (p, lbl) in enumerate(zip(pot_MM, labels)):
    axV_MM.plot(p[:, 0], p[:, 1], linewidth=2.0, label=lbl)
axV_MM.set_title("M–M",     fontsize=16)
axV_MM.set_xlim(0.3, 2.1)
axV_MM.set_ylim(-7.0, 12.0)
axV_MM.tick_params(labelsize=14, length=6, width=1.5)
axV_MM.legend(fontsize=12)

for i, (p, lbl) in enumerate(zip(pot_NN, labels)):
    axV_NN.plot(p[:, 0], p[:, 1], linewidth=2.0, label=lbl)
axV_NN.set_title("NA–NA",   fontsize=16)
axV_NN.set_xlim(0.25, 2.1)
axV_NN.set_ylim(-2.7, 6.0)
axV_NN.tick_params(labelsize=14, length=6, width=1.5)
axV_NN.legend(fontsize=12)

# ---------- ROW 2: RDFs — set None to use matplotlib autoscale ----------
rdf_ylim_MNA = (-0.5, 25.0)
rdf_ylim_MM  = (-0.5,  8.0)
rdf_ylim_NN  = (-0.5,  20.0)

kw_tgt = dict(color="black", marker="o", markerfacecolor="none",
              markersize=5, linewidth=0.8, label="Target")

axg_MNA.plot(tgt_MNA[:, 0], tgt_MNA[:, 1], **kw_tgt)
for g, lbl in zip(rdf_MNA, labels):
    axg_MNA.plot(g[:, 0], g[:, 1], linewidth=1.8, label=lbl)
axg_MNA.set_xlim(0.25, 2.1)
if rdf_ylim_MNA: axg_MNA.set_ylim(rdf_ylim_MNA)
axg_MNA.set_xlabel("r (nm)", fontsize=16)
axg_MNA.set_ylabel("RDF",    fontsize=16)
axg_MNA.tick_params(labelsize=14, length=6, width=1.5)
axg_MNA.legend(fontsize=12)

axg_MM.plot(tgt_MM[:, 0], tgt_MM[:, 1], **kw_tgt)
for g, lbl in zip(rdf_MM, labels):
    axg_MM.plot(g[:, 0], g[:, 1], linewidth=1.8, label=lbl)
axg_MM.set_xlim(0.3, 2.1)
if rdf_ylim_MM: axg_MM.set_ylim(rdf_ylim_MM)
axg_MM.set_xlabel("r (nm)", fontsize=16)
axg_MM.tick_params(labelsize=14, length=6, width=1.5)
axg_MM.legend(fontsize=12)

axg_NN.plot(tgt_NN[:, 0], tgt_NN[:, 1], **kw_tgt)
for g, lbl in zip(rdf_NN, labels):
    axg_NN.plot(g[:, 0], g[:, 1], linewidth=1.8, label=lbl)
axg_NN.set_xlim(0.25, 2.1)
if rdf_ylim_NN: axg_NN.set_ylim(rdf_ylim_NN)
axg_NN.set_xlabel("r (nm)", fontsize=16)
axg_NN.tick_params(labelsize=14, length=6, width=1.5)
axg_NN.legend(fontsize=12)

# ---------- Panel labels ----------
panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
for k, ax in enumerate(axes.flat):
    ax.text(-0.04, 1.12, panel_labels[k],
            transform=ax.transAxes, fontsize=20,
            va="top", ha="left", clip_on=False)

# ---------- Thicken spines ----------
for ax in axes.flat:
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout(pad=1.2, h_pad=3.0, w_pad=2.0)
plt.savefig("implicit_V_RDF.pdf", dpi=300, format="pdf", bbox_inches="tight")
plt.show()
