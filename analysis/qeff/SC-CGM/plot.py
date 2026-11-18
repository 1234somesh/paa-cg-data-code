import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt
from util import error_binning

# ------------ Load topology + trajectory ------------
u = mda.Universe(
    "/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/SC_CGM/LZ_4.4/nvt.gro",  # INPUT: atomistic GROMACS structure (.gro) used as topology -- replace with your own
    "/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang/SC_CGM/LZ_4.4/NVT.xtc"  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
)

# ------------ Selections ------------
polymer = u.select_atoms("resname PAA")
ions = u.select_atoms("resname NA")
if ions.n_atoms == 0:
    ions = u.select_atoms("name NA")

print("Number of PAA residues:", len(polymer.residues))
print("Number of NA ions:     ", ions.n_atoms)

if len(polymer.residues) == 0 or ions.n_atoms == 0:
    raise SystemExit("Selection failed: check resname/name for PAA or NA.")

# ------------ Parameters (nm) ------------
rc_nm = 0.72
ANG_TO_NM = 0.1  # MDAnalysis positions are Å

n_monomers = len(polymer.residues)

# charges (in units of e)
monomer_charge = -1.0
ion_charge = +1.0

# ------------ Main loop: store per-frame counts ------------
counts_per_frame = []

for ts in u.trajectory:
    # COM of each residue (Å)
    monomer_coms_A = polymer.center_of_mass(compound="residues")  # (n_monomers, 3)
    ion_pos_A = ions.positions                                    # (n_ions, 3)

    # distances with PBC (Å): shape (n_ions, n_monomers)
    d_A = distance_array(ion_pos_A, monomer_coms_A, box=u.dimensions)

    # closest residue index for each ion
    j_min = np.argmin(d_A, axis=1)
    dmin_nm = d_A[np.arange(d_A.shape[0]), j_min] * ANG_TO_NM

    # ions within rc contribute to that closest bead
    mask = dmin_nm <= rc_nm
    frame_counts = np.zeros(n_monomers, dtype=np.int64)
    np.add.at(frame_counts, j_min[mask], 1)
    counts_per_frame.append(frame_counts)

counts_per_frame = np.array(counts_per_frame)  # (n_frames, n_monomers)
n_frames = counts_per_frame.shape[0]
print("Number of frames:", n_frames)

# ------------ q_eff(rc, j) with error via binning analysis ------------
q_eff = np.zeros(n_monomers)
q_eff_err = np.zeros(n_monomers)

for i in range(n_monomers):
    # per-frame q_eff time series for monomer i
    q_series = monomer_charge + ion_charge * counts_per_frame[:, i].astype(float)
    mean, tau, N_eff, err = error_binning(q_series)
    q_eff[i] = mean
    q_eff_err[i] = err

j = np.arange(1, n_monomers + 1, dtype=int)
out = np.column_stack((j, q_eff, q_eff_err))
np.savetxt("qeff_rc0p72nm_vs_j.dat", out, header="j  q_eff_rc(0.72nm)  error")

print("Saved: qeff_rc0p72nm_vs_j.dat")

# ------------ Plot ------------
plt.figure(figsize=(7, 4.5))
plt.errorbar(j, q_eff, yerr=q_eff_err, marker="o", linewidth=2.0, capsize=3)
plt.xlabel("Bead index j", fontsize=16, fontfamily="serif", fontweight="semibold")
plt.ylabel(r"$q_{\mathrm{eff}}(r_c=0.72\ \mathrm{nm},\, j)$", fontsize=16, fontfamily="serif", fontweight="semibold")
plt.tick_params(labelsize=14, length=6.0, width=2.0)
for spine in plt.gca().spines.values():
    spine.set_linewidth(2.0)
plt.tight_layout()
plt.savefig("qeff_rc0p72nm_vs_j.pdf", dpi=300)
plt.show()

