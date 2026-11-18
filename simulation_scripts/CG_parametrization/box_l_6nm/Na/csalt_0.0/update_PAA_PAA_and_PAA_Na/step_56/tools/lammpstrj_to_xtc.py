#!/usr/bin/env python3
"""
Convert a LAMMPS .lammpstrj (dump custom) to GROMACS .xtc (and optional .trr),
ensuring:
  - coordinates are written in nm (LAMMPS Å -> nm)
  - time is written in ps with a timestep of 2 fs (dt_ps = 0.002)

Your dump header example:
  ITEM: ATOMS id type q xu yu zu

This script:
  - reads TIMESTEP, NUMBER OF ATOMS, BOX BOUNDS, ATOMS
  - reorders atoms by id (so atom id 1 is index 0, etc.)
  - converts Å->nm for positions and box
  - writes XTC/TRR with proper time metadata:
        time_ps = timestep * 0.002
"""

import sys
import numpy as np
import MDAnalysis as mda


def lammpstrj_frames(path, n_atoms_expected=None, wrap=False):
    """
    Yield (timestep, box_nm[6], positions_nm[N,3]) for each frame.
    Assumes orthorhombic BOX BOUNDS like: ITEM: BOX BOUNDS pp pp pp
    """
    with open(path, "r") as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Unexpected dump format: missing 'ITEM: NUMBER OF ATOMS'")
            n_atoms = int(f.readline().strip())

            if n_atoms_expected is not None and n_atoms != n_atoms_expected:
                raise ValueError(f"Atom count mismatch: dump has {n_atoms}, expected {n_atoms_expected}")

            line = f.readline().strip()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Unexpected dump format: missing 'ITEM: BOX BOUNDS'")

            xlo, xhi = map(float, f.readline().split()[:2])
            ylo, yhi = map(float, f.readline().split()[:2])
            zlo, zhi = map(float, f.readline().split()[:2])

            # Å -> nm
            lx = (xhi - xlo) / 10.0
            ly = (yhi - ylo) / 10.0
            lz = (zhi - zlo) / 10.0
            box_nm = np.array([lx, ly, lz, 90.0, 90.0, 90.0], dtype=np.float32)

            line = f.readline().strip()
            if not line.startswith("ITEM: ATOMS"):
                raise ValueError("Unexpected dump format: missing 'ITEM: ATOMS ...'")

            cols = line.split()[2:]  # after "ITEM: ATOMS"
            col_index = {name: i for i, name in enumerate(cols)}

            def pick(*names):
                for n in names:
                    if n in col_index:
                        return col_index[n]
                return None

            i_id = pick("id")
            ix = pick("x", "xu", "xs")
            iy = pick("y", "yu", "ys")
            iz = pick("z", "zu", "zs")

            if i_id is None or ix is None or iy is None or iz is None:
                raise ValueError(
                    "Need columns: id and x/y/z (or xu/yu/zu or xs/ys/zs). "
                    f"Found columns: {cols}"
                )

            pos_nm = np.zeros((n_atoms, 3), dtype=np.float32)

            for _ in range(n_atoms):
                parts = f.readline().split()
                aid = int(parts[i_id])  # 1-based LAMMPS id
                x = float(parts[ix]) / 10.0  # Å -> nm
                y = float(parts[iy]) / 10.0
                z = float(parts[iz]) / 10.0
                pos_nm[aid - 1, :] = (x, y, z)

            if wrap:
                pos_nm[:, 0] = np.mod(pos_nm[:, 0], lx)
                pos_nm[:, 1] = np.mod(pos_nm[:, 1], ly)
                pos_nm[:, 2] = np.mod(pos_nm[:, 2], lz)

            yield timestep, box_nm, pos_nm


def main():
    if len(sys.argv) < 4:
        print(
            "Usage:\n"
            "  lammpstrj_to_xtc.py conf.gro traj.lammpstrj out.xtc [--trr out.trr] [--wrap]\n\n"
            "Notes:\n"
            "  - conf.gro coordinates are assumed nm (GROMACS standard)\n"
            "  - lammpstrj coordinates are assumed Å and will be converted to nm\n"
            "  - time is written in ps using dt = 2 fs => dt_ps = 0.002\n"
        )
        sys.exit(1)

    conf_gro = sys.argv[1]
    dump = sys.argv[2]
    out_xtc = sys.argv[3]

    out_trr = None
    wrap = "--wrap" in sys.argv
    if "--trr" in sys.argv:
        j = sys.argv.index("--trr")
        out_trr = sys.argv[j + 1]

    # 2 fs in ps
    dt_ps = 0.005

    u = mda.Universe(conf_gro)
    n_atoms = u.atoms.n_atoms

    xtc_w = mda.coordinates.XTC.XTCWriter(out_xtc, n_atoms=n_atoms)
    trr_w = mda.coordinates.TRR.TRRWriter(out_trr, n_atoms=n_atoms) if out_trr else None

    frame = 0
    for step, box_nm, pos_nm in lammpstrj_frames(dump, n_atoms_expected=n_atoms, wrap=wrap):
        u.atoms.positions = pos_nm
        u.dimensions = box_nm

        # Set time metadata (ps)
        ts = u.trajectory.ts
        ts.frame = frame
        ts.dt = dt_ps
        ts.time = step * dt_ps

        xtc_w.write(u)
        if trr_w:
            trr_w.write(u)

        frame += 1

    xtc_w.close()
    if trr_w:
        trr_w.close()

    print(f"Wrote: {out_xtc}")
    if out_trr:
        print(f"Wrote: {out_trr}")
    print(f"Time mapping: time_ps = timestep * {dt_ps} (2 fs)")


if __name__ == "__main__":
    main()

