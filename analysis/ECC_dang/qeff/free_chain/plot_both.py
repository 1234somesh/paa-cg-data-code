#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Files ----------------
aa_file = "qeff_rc0p72nm_vs_j.dat"
cg_file = "qeff_rc0p72nm_vs_j_CG_LAMMPS.dat"

# ---------------- Load ----------------
aa = np.loadtxt(aa_file)
cg = np.loadtxt(cg_file)

j_aa, qeff_aa = aa[:, 0], aa[:, 1]
j_cg, qeff_cg = cg[:, 0], cg[:, 1]

# ---------------- Put on same charge scale ----------------
# This ECC run represents the q=0.75 correction via a scaled dielectric
# (127.822222 = 71.9/0.75**2 in ../step_50/in.lmp) while keeping the
# literal per-atom charges at +/-1.0 -- q_eff is a raw charge sum (never
# touches the Coulomb interaction/dielectric), so the dumped charges are
# already genuine "e" units and need no rescaling here.
q_cg_mono = 1.0  # |charge| on a CG PAA bead in your LAMMPS model
scale = 1.0 / q_cg_mono
qeff_cg_scaled = qeff_cg * scale

# ---------------- Plot ----------------
plt.figure(figsize=(7, 4.5))

plt.plot(j_aa, qeff_aa, marker="o", linewidth=2.0, label="AA")
plt.plot(j_cg, qeff_cg_scaled, marker="o", linewidth=2.0,
         label=f"CG")

plt.xlabel("Bead index $j$", fontsize=16, fontfamily="serif", fontweight="semibold")
plt.ylabel(r"$q_{\mathrm{eff}}(r_c=0.72\ \mathrm{nm},\, j)$  (e)",
           fontsize=16, fontfamily="serif", fontweight="semibold")

plt.tick_params(labelsize=14, length=6.0, width=2.0)
for spine in plt.gca().spines.values():
    spine.set_linewidth(2.0)

plt.legend(fontsize=13)
plt.tight_layout()
plt.savefig("qeff_rc0p72nm_AA_vs_CG_scaled.pdf", dpi=300)
plt.show()

# ---------------- Quick sanity output ----------------
print(f"CG scale factor = {scale:.6f}")
print(f"AA mean q_eff = {qeff_aa.mean():.4f} e")
print(f"CG mean q_eff (raw) = {qeff_cg.mean():.4f} (model charge units)")
print(f"CG mean q_eff (scaled) = {qeff_cg_scaled.mean():.4f} e")

