#!/usr/bin/env python3
"""
Histogram bonded (adjacent-monomer) distances from a LAMMPS dump (.lammpstrj).

Your case:
- Monomers are the negatively charged beads (q ~ -0.117851) and you said
  atom IDs 81-161 are NOT monomers -> so monomers are IDs 1-80.
- Dump uses: ITEM: ATOMS id type q xu yu zu  (positions in Angstrom, unwrapped)
- Output distances are reported in nm (GROMACS units).

How "adjacent distance" is defined here:
- Uses the bond list from the LAMMPS data file (recommended, robust).
  We take all bonds where BOTH atoms are monomers (atom id <= 80).
  This corresponds to adjacent beads along each polymer chain.

PBC handling:
- Applies minimum-image convention using the box lengths from BOX BOUNDS.

Outputs:
- <prefix>_lengths_nm.csv : one distance per bond per frame
- <prefix>_hist.png       : histogram plot
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def parse_bonds_from_data(data_path: Path):
    """Return list of (i, j) bonded atom-id pairs from a LAMMPS data file."""
    lines = data_path.read_text(errors="ignore").splitlines()

    bonds_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Bonds":
            bonds_idx = i
            break
    if bonds_idx is None:
        raise RuntimeError("No 'Bonds' section found in the data file.")

    bonds = []
    i = bonds_idx + 1
    while i < len(lines) and lines[i].strip() == "":
        i += 1

    while i < len(lines):
        s = lines[i].strip()
        if s == "" or re.match(r"^[A-Za-z]", s):
            break
        parts = s.split()
        if len(parts) >= 4:
            ai = int(parts[2])
            aj = int(parts[3])
            if ai != aj:
                bonds.append((min(ai, aj), max(ai, aj)))
        i += 1

    return sorted(set(bonds))


def lammpstrj_bond_lengths(
    traj_path: Path,
    bonds,
    monomer_max_id: int = 80,
    require_charge=None,
    charge_tol: float = 1e-4,
):
    """
    Stream a LAMMPS dump and yield bond lengths (nm) for each frame.

    bonds: list of (i, j) atom-id pairs (1-based).
    If require_charge is not None, only atoms with q close to that value are used.
    """
    bonds_m = [(i, j) for (i, j) in bonds if i <= monomer_max_id and j <= monomer_max_id]
    if not bonds_m:
        raise RuntimeError("No monomer-monomer bonds found after filtering. Check monomer_max_id.")

    needed_ids = set()
    for i, j in bonds_m:
        needed_ids.add(i)
        needed_ids.add(j)
    needed_ids = sorted(needed_ids)

    with traj_path.open("r") as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            step = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError("Unexpected format: missing 'ITEM: NUMBER OF ATOMS'")
            n_atoms = int(f.readline().strip())

            line = f.readline().strip()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("Unexpected format: missing 'ITEM: BOX BOUNDS'")

            xlo, xhi = map(float, f.readline().split()[:2])
            ylo, yhi = map(float, f.readline().split()[:2])
            zlo, zhi = map(float, f.readline().split()[:2])

            Lx = xhi - xlo
            Ly = yhi - ylo
            Lz = zhi - zlo

            line = f.readline().strip()
            if not line.startswith("ITEM: ATOMS"):
                raise RuntimeError("Unexpected format: missing 'ITEM: ATOMS ...'")

            cols = line.split()[2:]
            col_index = {c: k for k, c in enumerate(cols)}

            def pick(*names):
                for n in names:
                    if n in col_index:
                        return col_index[n]
                return None

            i_id = pick("id")
            i_q = pick("q")
            ix = pick("x", "xu", "xs")
            iy = pick("y", "yu", "ys")
            iz = pick("z", "zu", "zs")
            if None in (i_id, ix, iy, iz):
                raise RuntimeError(f"Need id and x/y/z columns. Found: {cols}")

            pos = {}  # atom-id -> (x,y,z) in Angstrom
            for _ in range(n_atoms):
                parts = f.readline().split()
                aid = int(parts[i_id])

                if aid > monomer_max_id:
                    continue

                if require_charge is not None and i_q is not None:
                    q = float(parts[i_q])
                    if abs(q - float(require_charge)) > charge_tol:
                        continue

                if aid in needed_ids:
                    x = float(parts[ix])
                    y = float(parts[iy])
                    z = float(parts[iz])
                    pos[aid] = (x, y, z)

            lengths_nm = []
            for i, j in bonds_m:
                if i not in pos or j not in pos:
                    continue
                xi, yi, zi = pos[i]
                xj, yj, zj = pos[j]

                dx = xi - xj
                dy = yi - yj
                dz = zi - zj

                dx -= Lx * round(dx / Lx)
                dy -= Ly * round(dy / Ly)
                dz -= Lz * round(dz / Lz)

                rA = math.sqrt(dx * dx + dy * dy + dz * dz)  # Angstrom
                lengths_nm.append(rA / 10.0)                 # nm

            yield step, np.array(lengths_nm, dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data", type=Path, help="LAMMPS data file (to read Bonds section)")
    ap.add_argument("traj", type=Path, help="LAMMPS dump trajectory (.lammpstrj)")
    ap.add_argument("--monomer-max-id", type=int, default=80, help="Monomer atom IDs are <= this value (default: 80)")
    ap.add_argument("--require-charge", type=str, default="-0.117851",
                    help="Safety filter on monomer charge (default: -0.117851). Set to 'none' to disable.")
    ap.add_argument("--bins", type=int, default=80, help="Histogram bins (default: 80)")
    ap.add_argument("--max-frames", type=int, default=0, help="If >0, stop after this many frames (debug)")
    ap.add_argument("--out-prefix", type=str, default="bond", help="Output prefix (default: bond)")
    args = ap.parse_args()

    require_charge = None if args.require_charge.lower() == "none" else float(args.require_charge)

    bonds = parse_bonds_from_data(args.data)

    all_lengths = []
    nframes = 0
    for step, lengths_nm in lammpstrj_bond_lengths(
        args.traj,
        bonds,
        monomer_max_id=args.monomer_max_id,
        require_charge=require_charge,
    ):
        if lengths_nm.size:
            all_lengths.append(lengths_nm)
        nframes += 1
        if args.max_frames and nframes >= args.max_frames:
            break

    if not all_lengths:
        raise RuntimeError("No bond lengths collected. Check monomer selection and dump format.")

    all_lengths = np.concatenate(all_lengths)

    csv_path = Path(f"{args.out_prefix}_lengths_nm.csv")
    np.savetxt(csv_path, all_lengths, delimiter=",", header="bond_length_nm", comments="")
    print(f"Wrote {csv_path} ({all_lengths.size} bond lengths)")

        # Normalized histogram (area = 1)
    hist, bin_edges = np.histogram(all_lengths, bins=args.bins, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    hist_dat = Path(f"{args.out_prefix}_hist_density.dat")
    np.savetxt(hist_dat, np.c_[bin_centers, hist], fmt="%.8f",
               header="bin_center_nm  probability_density_1_per_nm", comments="# ")
    print(f"Wrote {hist_dat}")

    fig_path = Path(f"{args.out_prefix}_hist.png")
    plt.figure()
    plt.hist(all_lengths, bins=args.bins, density=True)
    plt.xlabel("Bond length (nm)")
    plt.ylabel("Probability density (1/nm)")
    plt.title("Bond-length distribution (monomer–monomer bonds)")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=200)
    print(f"Wrote {fig_path}")

    print(f"Frames read: {nframes}")
    print(f"Mean bond length (nm): {all_lengths.mean():.4f}")
    print(f"Std  bond length (nm): {all_lengths.std(ddof=1):.4f}")
    print(f"Min/Max (nm): {all_lengths.min():.4f} / {all_lengths.max():.4f}")


if __name__ == "__main__":
    main()
