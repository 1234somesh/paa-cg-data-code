#!/usr/bin/env python3
import numpy as np
import math

# ---------------- USER CONFIG ----------------
LAMMPS_DUMP   = "traj_paa_cg.lammpstrj"   # path to your LAMMPS trajectory
N_MONOMERS    = 20                        # total beads in the chain (1..20)
ATOMTYPE_PAA  = 1                         # polymer beads are type 1
# --------------------------------------------

ANG_TO_NM = 0.1  # 1 Å = 0.1 nm

def compute_R2_from_lammps(filename,
                           n_monomers=N_MONOMERS,
                           atomtype_poly=ATOMTYPE_PAA):
    """
    Compute <R^2(n)> for a single 20-mer chain from a LAMMPS dump file.

    All monomers 1..20 (1-based), i.e. 0-based indices 0..19, are included.
    Max separation n = 19.

    Uses:
      - wrapped coordinates obtained from unwrapped xu,yu,zu via modulo L.
      - minimum image convention for bead–bead vectors.
    Units: output R^2 in nm^2, contour length s in nm.
    """
    # All monomers: 1..20 (1-based) => 0..19 (0-based)
    interior_start = 0              # 0-based index for resid 1
    interior_end   = n_monomers - 1  # 19 for a 20-mer (resid 20)
    n_interior     = interior_end - interior_start + 1  # 20
    max_n          = n_interior - 1                      # 19

    sum_r2 = np.zeros(max_n, dtype=float)
    cnt    = np.zeros(max_n, dtype=int)

    sum_r1 = 0.0   # for mean n=1 bond length (interior only)
    cnt1   = 0
    frames = 0

    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            # Look for the start of a frame
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            # TIMESTEP line
            ts_line = f.readline()
            # ts = int(ts_line.strip())  # not used but could be

            # NUMBER OF ATOMS
            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Expected 'ITEM: NUMBER OF ATOMS', got: " + line)
            natoms = int(f.readline().strip())

            # BOX BOUNDS
            line = f.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Expected 'ITEM: BOX BOUNDS', got: " + line)

            xlo, xhi = map(float, f.readline().split()[:2])
            ylo, yhi = map(float, f.readline().split()[:2])
            zlo, zhi = map(float, f.readline().split()[:2])
            Lx, Ly, Lz = xhi - xlo, yhi - ylo, zhi - zlo

            # ATOMS header
            line = f.readline()
            if not line.startswith("ITEM: ATOMS"):
                raise ValueError("Expected 'ITEM: ATOMS', got: " + line)

            # Positions of the 20 monomers (in Å, wrapped)
            positions = np.zeros((n_monomers, 3), dtype=float)
            present   = np.zeros(n_monomers, dtype=bool)

            for _ in range(natoms):
                parts = f.readline().split()
                if not parts:
                    continue

                atom_id = int(parts[0])
                atype   = int(parts[1])

                # monomer beads: IDs 1..20 of type 1
                if 1 <= atom_id <= n_monomers and atype == atomtype_poly:
                    xu = float(parts[3])
                    yu = float(parts[4])
                    zu = float(parts[5])

                    # Wrap unwrapped coords into primary box [xlo,xhi)
                    x = (xu - xlo) % Lx + xlo
                    y = (yu - ylo) % Ly + ylo
                    z = (zu - zlo) % Lz + zlo

                    positions[atom_id - 1] = (x, y, z)
                    present[atom_id - 1]   = True

            # If any bead missing, skip this frame
            if not present.all():
                continue

            frames += 1

            # --- Compute R^2(n) for this frame, all resid 1..20 ---
            for n in range(1, max_n + 1):         # n = 1..17
                # i0 runs over interior_start..interior_end-n
                for i0 in range(interior_start, interior_end + 1 - n):
                    j0 = i0 + n
                    # j0 is automatically <= interior_end by loop bounds

                    dx, dy, dz = positions[j0] - positions[i0]

                    # Minimum image convention
                    if dx > 0.5 * Lx:
                        dx -= Lx
                    elif dx < -0.5 * Lx:
                        dx += Lx

                    if dy > 0.5 * Ly:
                        dy -= Ly
                    elif dy < -0.5 * Ly:
                        dy += Ly

                    if dz > 0.5 * Lz:
                        dz -= Lz
                    elif dz < -0.5 * Lz:
                        dz += Lz

                    # Convert to nm
                    dx_nm = dx * ANG_TO_NM
                    dy_nm = dy * ANG_TO_NM
                    dz_nm = dz * ANG_TO_NM

                    r2 = dx_nm*dx_nm + dy_nm*dy_nm + dz_nm*dz_nm
                    sum_r2[n - 1] += r2
                    cnt[n - 1]    += 1

                    # n = 1 bond stats: all bonds 1-2 .. 19-20
                    if n == 1:
                        r = math.sqrt(r2)
                        sum_r1 += r
                        cnt1   += 1

    if frames == 0:
        raise RuntimeError("No frames parsed. Check file name and format.")

    R2_mean = sum_r2 / cnt
    ns      = np.arange(1, max_n + 1, dtype=int)

    mean_bond_len_nm = sum_r1 / cnt1
    s_axis_nm        = ns * mean_bond_len_nm

    return ns, R2_mean, mean_bond_len_nm, s_axis_nm, frames


def main():
    ns, R2_mean, L, s_axis_nm, frames = compute_R2_from_lammps(LAMMPS_DUMP)

    print(f"Frames processed: {frames}")
    print(f"Mean n=1 bond length L (all resid 1..20) = {L:.4f} nm")

    # Save CSVs
    np.savetxt(
        "R2_of_n_lammps_all.csv",
        np.c_[ns, R2_mean],
        delimiter=",",
        header="n,<R2(n)>_nm2   (all resid 1..20)",
        comments=""
    )
    np.savetxt(
        "R2_of_s_lammps_all.csv",
        np.c_[s_axis_nm, R2_mean],
        delimiter=",",
        header="s_nm,<R2>_nm2   (all resid 1..20)",
        comments=""
    )
    print("Wrote R2_of_n_lammps_all.csv and R2_of_s_lammps_all.csv")

if __name__ == "__main__":
    main()

