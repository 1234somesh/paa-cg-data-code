import numpy as np
import matplotlib.pyplot as plt

traj = "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization/box_l_6nm/Try14/step_11/traj_paa_cg.lammpstrj"  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own

# --- Parameters ---
rc_nm = 0.72
rc_A  = rc_nm * 10.0     # LAMMPS dump here looks like Å (box ~ 43.7 Å = 4.37 nm)

# --- Helpers ---
def wrap_to_box(x, lo, L):
    """Wrap coordinate x into [lo, lo+L). Works even for huge unwrapped values."""
    return (x - lo) % L + lo

def min_image(d, L):
    """Minimum-image displacement."""
    return d - np.rint(d / L) * L

def read_lammpstrj_frames(fname):
    """Stream frames from an orthorhombic LAMMPS dump (no tilt factors)."""
    with open(fname, "r") as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            _ = f.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(f.readline().strip())

            _ = f.readline()  # ITEM: BOX BOUNDS ...
            xlo, xhi = map(float, f.readline().split()[:2])
            ylo, yhi = map(float, f.readline().split()[:2])
            zlo, zhi = map(float, f.readline().split()[:2])
            Lx, Ly, Lz = (xhi - xlo), (yhi - ylo), (zhi - zlo)

            header = f.readline().strip().split()
            cols = header[2:]  # after "ITEM: ATOMS"
            col_index = {c: i for i, c in enumerate(cols)}

            # Support xu/yu/zu or x/y/z
            xkey = "xu" if "xu" in col_index else "x"
            ykey = "yu" if "yu" in col_index else "y"
            zkey = "zu" if "zu" in col_index else "z"

            needed = ["id", "type", "q", xkey, ykey, zkey]
            for k in needed:
                if k not in col_index:
                    raise RuntimeError(f"Missing column '{k}' in dump header: {cols}")

            id_i   = col_index["id"]
            type_i = col_index["type"]
            q_i    = col_index["q"]
            x_i, y_i, z_i = col_index[xkey], col_index[ykey], col_index[zkey]

            ids = np.empty(natoms, dtype=int)
            typ = np.empty(natoms, dtype=int)
            q   = np.empty(natoms, dtype=float)
            x   = np.empty(natoms, dtype=float)
            y   = np.empty(natoms, dtype=float)
            z   = np.empty(natoms, dtype=float)

            for i in range(natoms):
                parts = f.readline().split()
                ids[i] = int(parts[id_i])
                typ[i] = int(parts[type_i])
                q[i]   = float(parts[q_i])
                x[i]   = float(parts[x_i])
                y[i]   = float(parts[y_i])
                z[i]   = float(parts[z_i])

            yield timestep, (xlo, ylo, zlo), (Lx, Ly, Lz), ids, typ, q, x, y, z


# --- Accumulators (set after first frame once we know bead count/order) ---
q_beads = None
sum_assigned_ion_charge = None
n_frames = 0

for timestep, (xlo, ylo, zlo), (Lx, Ly, Lz), ids, typ, q, x, y, z in read_lammpstrj_frames(traj):
    # sort by atom id to keep bead order consistent (ids 1..20 are PAA in your example)
    order = np.argsort(ids)
    ids, typ, q = ids[order], typ[order], q[order]
    x, y, z     = x[order], y[order], z[order]

    # wrap coords into primary box (important for xu/yu/zu)
    xw = wrap_to_box(x, xlo, Lx)
    yw = wrap_to_box(y, ylo, Ly)
    zw = wrap_to_box(z, zlo, Lz)

    paa = (typ == 1)   # first 20 beads in your dump
    na  = (typ == 2)   # last 20 ions in your dump

    if q_beads is None:
        bead_ids = ids[paa]
        # bead index j will follow ascending bead id
        q_beads = q[paa].copy()  # could be per-bead (usually all equal)
        sum_assigned_ion_charge = np.zeros_like(q_beads, dtype=float)
        print(f"Detected {paa.sum()} PAA beads and {na.sum()} Na ions")
        print(f"Bead ids (j=1..N): {bead_ids}")

    # positions
    bead_pos = np.column_stack((xw[paa], yw[paa], zw[paa]))     # (Nb, 3)
    ion_pos  = np.column_stack((xw[na],  yw[na],  zw[na]))      # (Ni, 3)
    ion_q    = q[na]                                            # (Ni,)

    # distances ion->bead with minimum image (Å)
    dx = ion_pos[:, None, 0] - bead_pos[None, :, 0]
    dy = ion_pos[:, None, 1] - bead_pos[None, :, 1]
    dz = ion_pos[:, None, 2] - bead_pos[None, :, 2]
    dx = min_image(dx, Lx)
    dy = min_image(dy, Ly)
    dz = min_image(dz, Lz)
    d  = np.sqrt(dx*dx + dy*dy + dz*dz)  # (Ni, Nb)

    j_min = np.argmin(d, axis=1)
    dmin  = d[np.arange(d.shape[0]), j_min]

    # ions within rc contribute their charge to their closest bead
    mask = (dmin <= rc_A)
    np.add.at(sum_assigned_ion_charge, j_min[mask], ion_q[mask])

    n_frames += 1

if n_frames == 0:
    raise RuntimeError("No frames read. Check the trajectory path and format.")

# Average assigned ion charge within rc for each bead
P_rc = sum_assigned_ion_charge / n_frames

# Effective charge within rc
q_eff = q_beads + P_rc

# Save + plot
j = np.arange(1, len(q_eff) + 1, dtype=int)
np.savetxt("qeff_rc0p72nm_vs_j_CG_LAMMPS.dat",
           np.column_stack((j, q_eff)),
           header="j  q_eff_rc(0.72nm)")

plt.figure(figsize=(7, 4.5))
plt.plot(j, q_eff, marker="o", linewidth=2.0)
plt.xlabel("Bead index j", fontsize=16, fontfamily="serif", fontweight="semibold")
plt.ylabel(r"$q_{\mathrm{eff}}(r_c=0.72\ \mathrm{nm},\, j)$",
           fontsize=16, fontfamily="serif", fontweight="semibold")
plt.tick_params(labelsize=14, length=6.0, width=2.0)
for spine in plt.gca().spines.values():
    spine.set_linewidth(2.0)
plt.tight_layout()
plt.savefig("qeff_rc0p72nm_vs_j_CG_LAMMPS.pdf", dpi=300)
plt.show()

print(f"Frames used: {n_frames}")
print(f"Mean bead charge (type 1): {np.mean(q_beads):.6f}")
print(f"Mean ion  charge (type 2): {np.mean(ion_q):.6f}")

