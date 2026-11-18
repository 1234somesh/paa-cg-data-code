#!/usr/bin/env python3
"""
Combined plot: AA and CG cylindrical RDF and coordination number on the same axes,
plus Contact-Value-Theorem osmotic pressure (P_osm) at r_beta = 0.5*Lx.
Centroid method only -- the nearest-monomer/bead diagnostic (Option A) and its
Monte Carlo normalization were explored and dropped: they answer a different
question (nearest-neighbor distance to a point cloud) than what the contact
value theorem needs (distance from a well-defined cylindrical axis), so they
were never a valid basis for P_osm. Kept here only in spirit as the reason
r_beta is read from the centroid-based profile below.
"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

A2NM = 0.1
nbins = 200

# ---- CVT / osmotic-pressure constants -----------------------------
T_KELVIN = 300.0          # <-- set to your actual simulation temperature
KB_J     = 1.380649e-23   # Boltzmann constant, J/K
ANG3_TO_L = 1e-27         # 1 Angstrom^3 in liters
NA_AVO    = 6.02214076e23
ANG3_TO_MOLAR = 1.0 / (NA_AVO * ANG3_TO_L)   # 1 ion/Ang^3 -> mol/L (~1660.5)
R_GAS_L_BAR = 0.0831446   # L*bar / (mol*K)

def accessible_area_circle_in_rectangle(r, a, b, n=4000):
    """
    Area of a circle of radius r (centered on the chain) that actually
    lies inside the periodic box, whose transverse cross-section is a
    rectangle of half-widths a=Lx/2, b=Ly/2 centered on the same point.

    Why this is needed: a plain circular annulus area, pi*r^2, is only
    correct while the circle is small enough to fit entirely inside the
    box. Once r grows past the nearest box edge (r > a, or r > b), part
    of that circle would stick out through the periodic boundary -- it
    doesn't physically extend there, so pi*r^2 overcounts the available
    area from that point on. This function computes the TRUE area of
    "circle intersected with rectangle" instead of the naive full-circle
    formula, so the density normalization stays correct at every r, not
    just for r < min(a, b).

    How it works, concretely: by the rectangle's mirror symmetry, the
    intersection area is 4x whatever lies in just one quadrant (x,y >= 0).
    Walking across that quadrant in thin vertical strips of width dx
    (from x=0 to x=a), the circle's edge in that strip would reach height
    y = sqrt(r^2 - x^2) -- but the rectangle only extends up to y = b, so
    the ACCESSIBLE height of that strip is capped at b:
        accessible_height(x) = min(sqrt(r^2 - x^2), b)
    Summing (integrating) accessible_height(x) over the strips gives the
    quadrant's area; np.trapz does that numerical integration, and the
    result is multiplied by 4 for all four quadrants. For small r (circle
    fits entirely inside the box), this reduces to the ordinary pi*r^2;
    for large r (circle bigger than the box diagonal), it saturates at
    the full box area 4*a*b -- exactly the two limits you'd expect.
    """
    x = np.linspace(0.0, a, n)
    y = np.sqrt(np.clip(r*r - x*x, 0.0, None))
    y = np.minimum(y, b)
    return 4.0 * np.trapezoid(y, x)

def wrap_to_box(x, lo, L):
    return (x - lo) % L + lo

def pbc_circular_mean(x, lo, L):
    """
    PBC-safe mean position on a periodic axis of length L. A plain
    arithmetic mean of x-coordinates breaks down if the chain straddles
    the periodic boundary (some monomers near lo, others near lo+L): the
    naive mean lands near the middle of the box, nowhere near where the
    (torn-looking but actually contiguous) chain really sits. Treating
    each x as an angle around the periodic box and averaging on the unit
    circle instead sidesteps that -- the result is the correct wrapped
    position regardless of where the chain sits relative to the box edges.
    """
    ang = 2.0 * np.pi * (x - lo) / L
    mean_ang = np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang)))
    if mean_ang < 0:
        mean_ang += 2.0 * np.pi
    return lo + (mean_ang / (2.0 * np.pi)) * L

# ============================================================
# AA: atomistic cylindrical RDF (from rdf.py logic)
# ============================================================
print("=== Processing AA ===", flush=True)
u_aa = mda.Universe(
    "/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_f/periodic_5.12/NVT_nonwater.gro",  # INPUT: atomistic GROMACS structure (.gro) used as topology -- replace with your own
    "/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_freeze_Na/NVT.xtc"  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
)

paa_atoms = u_aa.select_atoms("name CH CO2M CH2 OE1 OE2 H1 H2A H2B")
na_atoms  = u_aa.select_atoms("name NA")
monomers  = paa_atoms.residues

print(f"  AA: {u_aa.atoms.n_atoms} atoms, {len(u_aa.trajectory)} frames", flush=True)
print(f"  PAA residues: {len(monomers)}, Na ions: {len(na_atoms)}", flush=True)

u_aa.trajectory[0]
Lx_aa, Ly_aa, Lz_aa = u_aa.dimensions[:3]
a_aa, b_aa = Lx_aa / 2.0, Ly_aa / 2.0
r_max_aa = np.sqrt(a_aa**2 + b_aa**2)   # same convention as CG: box-diagonal reach

bin_edges_aa   = np.linspace(0.0, r_max_aa, nbins + 1)
bin_centers_aa = 0.5 * (bin_edges_aa[:-1] + bin_edges_aa[1:])
rdf_counts_aa  = np.zeros(nbins)
n_frames_aa = 0

for ts in u_aa.trajectory:
    Lx, Ly, Lz = ts.dimensions[:3]
    mon_com = np.array([res.atoms.center_of_mass() for res in monomers])
    # PBC-safe chain center (matches the CG half below): a plain
    # mon_com[:, 0].mean() would be wrong whenever the chain drifts
    # near a box edge and its monomers end up wrapped across the
    # periodic boundary (some near x=0, some near x=Lx) -- the naive
    # mean would then land near the box middle, not where the chain
    # actually is, corrupting r_xy (and the RDF/CN/P_osm) for that frame.
    xw = wrap_to_box(mon_com[:, 0], 0.0, Lx)
    yw = wrap_to_box(mon_com[:, 1], 0.0, Ly)
    xc = pbc_circular_mean(xw, 0.0, Lx)
    yc = pbc_circular_mean(yw, 0.0, Ly)

    dx = na_atoms.positions[:, 0] - xc
    dy = na_atoms.positions[:, 1] - yc
    dx -= np.rint(dx / Lx) * Lx
    dy -= np.rint(dy / Ly) * Ly
    r_xy = np.sqrt(dx*dx + dy*dy)

    hist, _ = np.histogram(r_xy, bins=bin_edges_aa)
    rdf_counts_aa += hist
    n_frames_aa += 1
    if n_frames_aa % 500 == 0:
        print(f"  AA frame {n_frames_aa}", flush=True)

rho_aa = na_atoms.n_atoms / (Lx_aa * Ly_aa * Lz_aa)
A_edges_aa = np.array([accessible_area_circle_in_rectangle(r, a_aa, b_aa) for r in bin_edges_aa])
dA_aa = A_edges_aa[1:] - A_edges_aa[:-1]
dV_aa = dA_aa * Lz_aa   # rectangle-clipped normalization, correct at every r
g_cyl_aa = rdf_counts_aa / (n_frames_aa * rho_aa * dV_aa)
CN_aa = np.cumsum(rdf_counts_aa) / n_frames_aa
r_nm_aa = bin_centers_aa * A2NM

print(f"  AA done: {n_frames_aa} frames", flush=True)

np.savetxt("cyl_rdf_aa.dat", np.column_stack([r_nm_aa, g_cyl_aa]),
           header="r(nm) g_cyl(r)  [AA]")
np.savetxt("cyl_cn_aa.dat", np.column_stack([r_nm_aa, CN_aa]),
           header="r(nm) CN(r)  [AA]")

# ============================================================
# CG: implicit solvent cylindrical RDF (from rdf_imp.py logic)
# ============================================================
print("=== Processing CG ===", flush=True)

traj_cg = "/work/skurahatti/multivalent/2d-gel/PAA_dang/SC_CGM/LZ_freeze_Na/traj_paa_cg.lammpstrj"  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own

# Read first frame to set up bins
with open(traj_cg) as f:
    f.readline()  # ITEM: TIMESTEP
    f.readline()  # timestep
    f.readline()  # ITEM: NUMBER OF ATOMS
    natoms_cg = int(f.readline().strip())
    f.readline()  # ITEM: BOX BOUNDS
    xlo_cg, xhi_cg = map(float, f.readline().split())
    ylo_cg, yhi_cg = map(float, f.readline().split())
    zlo_cg, zhi_cg = map(float, f.readline().split())

Lx_cg = xhi_cg - xlo_cg
Ly_cg = yhi_cg - ylo_cg
Lz_cg = zhi_cg - zlo_cg
a_cg, b_cg = Lx_cg/2.0, Ly_cg/2.0
r_max_cg = np.sqrt(a_cg**2 + b_cg**2)

bin_edges_cg   = np.linspace(0.0, r_max_cg, nbins + 1)
bin_centers_cg = 0.5 * (bin_edges_cg[:-1] + bin_edges_cg[1:])
A_edges_cg = np.array([accessible_area_circle_in_rectangle(r, a_cg, b_cg) for r in bin_edges_cg])
dA_cg = A_edges_cg[1:] - A_edges_cg[:-1]
dV_cg = dA_cg * Lz_cg

g_sum_cg  = np.zeros(nbins)
cn_sum_cg = np.zeros(nbins)
n_frames_cg = 0
rho_cg_list = []   # track per-frame bulk Na density for a clean post-hoc average

# Full read
with open(traj_cg) as f:
    while True:
        line = f.readline()
        if not line:
            break
        if "ITEM: TIMESTEP" not in line:
            continue
        f.readline()  # timestep
        f.readline()  # ITEM: NUMBER OF ATOMS
        natoms = int(f.readline().strip())
        f.readline()  # ITEM: BOX BOUNDS
        xlo, xhi = map(float, f.readline().split())
        ylo, yhi = map(float, f.readline().split())
        zlo, zhi = map(float, f.readline().split())
        Lx, Ly, Lz = xhi-xlo, yhi-ylo, zhi-zlo
        f.readline()  # ITEM: ATOMS header

        typ = np.empty(natoms, dtype=int)
        x = np.empty(natoms)
        y = np.empty(natoms)

        for i in range(natoms):
            parts = f.readline().split()
            typ[i] = int(parts[1])
            x[i] = float(parts[3])
            y[i] = float(parts[4])

        xw = wrap_to_box(x, xlo, Lx)
        yw = wrap_to_box(y, ylo, Ly)

        paa_mask = (typ == 1)
        na_mask  = (typ == 2)

        xc = pbc_circular_mean(xw[paa_mask], xlo, Lx)
        yc = pbc_circular_mean(yw[paa_mask], ylo, Ly)

        dx = xw[na_mask] - xc
        dy = yw[na_mask] - yc
        dx -= np.rint(dx / Lx) * Lx
        dy -= np.rint(dy / Ly) * Ly
        r_xy = np.sqrt(dx*dx + dy*dy)

        hist, _ = np.histogram(r_xy, bins=bin_edges_cg)
        rho = np.count_nonzero(na_mask) / (Lx * Ly * Lz)
        g_sum_cg += hist / (rho * dV_cg)
        cn_sum_cg += np.cumsum(hist).astype(float)
        rho_cg_list.append(rho)
        n_frames_cg += 1

        if n_frames_cg % 20000 == 0:
            print(f"  CG frame {n_frames_cg}", flush=True)

g_cyl_cg = g_sum_cg / n_frames_cg
CN_cg = cn_sum_cg / n_frames_cg
r_nm_cg = bin_centers_cg * A2NM

print(f"  CG done: {n_frames_cg} frames", flush=True)

np.savetxt("cyl_rdf_cg.dat", np.column_stack([r_nm_cg, g_cyl_cg]),
           header="r(nm) g_cyl(r)  [CG]")
np.savetxt("cyl_cn_cg.dat", np.column_stack([r_nm_cg, CN_cg]),
           header="r(nm) CN(r)  [CG]")

# ============================================================
# Contact Value Theorem: P_osm from rho(r_beta), r_beta = 0.5*Lx
# ============================================================
# r_beta = 0.5*Lx is the distance to the nearest periodic-image face --
# the field-free symmetry plane where the CVT boundary condition applies.
# This requires "distance from a well-defined cylindrical axis," which the
# centroid method approximates (one averaged transverse point per frame,
# standing in for the axis) -- the only method used throughout this script.

def rho_at_rbeta(bin_centers, g_cyl, rho_bulk, r_beta, window_frac=0.10, label=""):
    """
    Local ion density at the cell boundary: rho(r) = g_cyl(r) * rho_bulk,
    averaged over a small window just inside r_beta (reduces bin noise
    and stays clear of r_beta itself). Also reports the local slope as
    a flatness diagnostic -- if it's not small, there isn't a real
    plateau yet at this box size and the resulting P_osm should be
    treated as provisional (rerun with a larger transverse box).
    """
    bin_centers = np.asarray(bin_centers)
    mask = bin_centers <= r_beta
    if mask.sum() < 3:
        raise ValueError(f"[{label}] r_beta={r_beta:.3f} Ang falls below the "
                          f"resolved bin range -- check units.")
    r = bin_centers[mask]
    rho_r = g_cyl[mask] * rho_bulk

    window_start = r.max() - window_frac * (r.max() - r.min())
    win = r >= window_start

    rho_beta = rho_r[win].mean()
    rho_beta_std = rho_r[win].std()
    slope, _ = np.polyfit(r[win], rho_r[win], 1)

    flat_ok = abs(slope) <= 0.05 * abs(rho_beta + 1e-12)
    print(f"  [{label}] r_beta = {r_beta*A2NM:.3f} nm | "
          f"rho(r_beta) = {rho_beta:.6f} +/- {rho_beta_std:.6f} ions/Ang^3 | "
          f"slope = {slope:.3e} "
          f"{'(flat)' if flat_ok else '<-- NOT FLAT, treat P_osm as provisional'}",
          flush=True)
    return rho_beta

r_beta_aa = 0.5 * Lx_aa
r_beta_cg = 0.5 * Lx_cg
rho_cg_bulk = np.mean(rho_cg_list)
print(f"  CG bulk density check: mean={rho_cg_bulk:.6e}, "
      f"std={np.std(rho_cg_list):.2e} ions/Ang^3 (should be ~constant under NVT)",
      flush=True)

print("=== Contact Value Theorem: rho(r_beta) and P_osm ===", flush=True)
rho_beta_aa = rho_at_rbeta(bin_centers_aa, g_cyl_aa, rho_aa,      r_beta_aa, label="AA")
rho_beta_cg = rho_at_rbeta(bin_centers_cg, g_cyl_cg, rho_cg_bulk, r_beta_cg, label="CG")

# P_osm = kB*T * sum_species rho_species(r_beta) -- single mobile species
# (Na+, no added salt) here; if you add a co-ion later, compute its own
# rho_at_rbeta(...) and add it into rho_beta_total before converting.
def posm_bar(rho_beta_ang3, T=T_KELVIN):
    c_molar = rho_beta_ang3 * ANG3_TO_MOLAR   # mol/L
    return c_molar * R_GAS_L_BAR * T, c_molar

P_osm_aa, c_aa = posm_bar(rho_beta_aa)
P_osm_cg, c_cg = posm_bar(rho_beta_cg)

print(f"  P_osm (AA) = {P_osm_aa:.4f} bar  ({c_aa*1000:.3f} mM equivalent)", flush=True)
print(f"  P_osm (CG) = {P_osm_cg:.4f} bar  ({c_cg*1000:.3f} mM equivalent)", flush=True)

# ============================================================
# Plot: AA vs CG
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# RDF
ax1.plot(r_nm_aa, g_cyl_aa, 'k-', lw=2.0, label="AA")
ax1.plot(r_nm_cg, g_cyl_cg, 'r--', lw=2.0, label="CG")
ax1.set_xlabel("r (nm)", fontsize=16, fontfamily="serif", fontweight="semibold")
ax1.set_ylabel(r"$g_{\mathrm{cyl}}(r)$", fontsize=16, fontfamily="serif", fontweight="semibold")
ax1.legend(fontsize=14)
ax1.tick_params(labelsize=14, length=8.0, width=2.0)
ax1.text(-0.04, 1.10, "(a)", transform=ax1.transAxes, fontsize=22, va="top", ha="left")

# CN
ax2.plot(r_nm_aa, CN_aa, 'k-', lw=2.0, label="AA")
ax2.plot(r_nm_cg, CN_cg, 'r--', lw=2.0, label="CG")
ax2.axhline(20.0, linestyle=":", color="gray", lw=1.5, label="20 ions")
ax2.set_xlabel("r (nm)", fontsize=16, fontfamily="serif", fontweight="semibold")
ax2.set_ylabel("CN(r)", fontsize=16, fontfamily="serif", fontweight="semibold")
ax2.legend(fontsize=14)
ax2.tick_params(labelsize=14, length=8.0, width=2.0)
ax2.text(-0.04, 1.10, "(b)", transform=ax2.transAxes, fontsize=22, va="top", ha="left")

for ax in (ax1, ax2):
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.tight_layout()
plt.savefig("cyl_rdf_aa_vs_cg.pdf", dpi=300, bbox_inches="tight")
#plt.savefig("cyl_rdf_aa_vs_cg.png", dpi=300, bbox_inches="tight")
plt.show()

print("Saved to cyl_rdf_aa_vs_cg.pdf and cyl_rdf_aa_vs_cg.png")

# ============================================================
# Bar plot: CVT osmotic pressure, AA vs CG
# ============================================================
fig2, ax3 = plt.subplots(figsize=(4.5, 5.0))

labels = ["AA", "CG"]
values = [P_osm_aa, P_osm_cg]
colors = ["black", "red"]

bars = ax3.bar(labels, values, color=colors, alpha=0.75, width=0.55,
                edgecolor="black", linewidth=1.2)
for b, v in zip(bars, values):
    ax3.text(b.get_x() + b.get_width() / 2, v, f"{v:.3f}",
              ha="center", va="bottom", fontsize=13, fontfamily="serif")

ax3.set_ylabel(r"$P_{\mathrm{osm}}$ (bar)", fontsize=16, fontfamily="serif", fontweight="semibold")
ax3.set_title(r"CVT, $r = 0.5\,L_x$", fontsize=14, fontfamily="serif")
ax3.tick_params(labelsize=14, length=8.0, width=2.0)
for spine in ax3.spines.values():
    spine.set_linewidth(2.0)

plt.tight_layout()
plt.savefig("posm_cvt_bar.pdf", dpi=300, bbox_inches="tight")
plt.show()

print("Saved to posm_cvt_bar.pdf")
