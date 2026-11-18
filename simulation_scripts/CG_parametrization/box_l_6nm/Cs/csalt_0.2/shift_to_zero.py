import argparse
import numpy as np


def shift_potential_to_zero_at_rcut(
    infile,
    outfile,
    r_cut=20.0,
    header_lines=5,
    mode="interp",
    snap_first_tail_point=True,
    zero_tail=True,
):
    """
    Shift a LAMMPS pair_style table by a constant so that U(r_cut) = 0.

    This applies:
        U_new(r) = U(r) - U(r_cut)

    In this workflow, the tables typically extend beyond the physical cutoff.
    A pure constant shift makes any already-flat zero tail become a nonzero
    plateau. To avoid that, this helper can also:
      1. snap the first point at/above r_cut to exactly r_cut
      2. force U(r) = 0 and F(r) = 0 for all tabulated points at/above r_cut

    Parameters
    ----------
    infile : str
        Input LAMMPS table file (5-line header like pair_*.table).
    outfile : str
        Output table file.
    r_cut : float
        Cutoff radius where U should become exactly zero (same units as r in table, Å in your workflow).
    header_lines : int
        Number of header lines to copy verbatim.
    mode : {"interp","nearest"}
        How to get U(r_cut) if r_cut is not exactly on the grid:
        - "interp": linear interpolation (recommended)
        - "nearest": use value at nearest tabulated r
    snap_first_tail_point : bool
        If True, move the first point at/above r_cut to exactly r_cut.
    zero_tail : bool
        If True, force U=0 and F=0 for all points at/above r_cut.
    """

    # --- read header and data ---
    with open(infile, "r") as f:
        header = [f.readline() for _ in range(header_lines)]
        data = np.loadtxt(f)

    idx = data[:, 0].astype(int)
    r   = data[:, 1]
    U   = data[:, 2]

    # --- get U(rcut) ---
    if mode == "interp":
        if r_cut < r.min() or r_cut > r.max():
            raise ValueError(f"r_cut={r_cut} is outside table range [{r.min()}, {r.max()}].")
        U_rcut = np.interp(r_cut, r, U)
    elif mode == "nearest":
        j = int(np.argmin(np.abs(r - r_cut)))
        U_rcut = U[j]
    else:
        raise ValueError("mode must be 'interp' or 'nearest'")

    # --- shift ---
    U_new = U - U_rcut

    if zero_tail:
        first_tail = int(np.searchsorted(r, r_cut, side="left"))
        if first_tail >= len(r):
            raise ValueError(
                f"r_cut={r_cut} is not covered by any tabulated point at/above cutoff."
            )

        if snap_first_tail_point and not np.isclose(r[first_tail], r_cut):
            if first_tail > 0 and r[first_tail - 1] >= r_cut:
                raise ValueError(
                    "Table is not strictly increasing around the cutoff; cannot snap safely."
                )
            r[first_tail] = r_cut

        U_new[first_tail:] = 0.0

    # --- recompute forces from the modified potential ---
    F_new = np.zeros_like(U_new)
    n = len(r)
    for i in range(1, n - 1):
        F_new[i] = -(U_new[i + 1] - U_new[i - 1]) / (r[i + 1] - r[i - 1])
    F_new[0] = -(U_new[1] - U_new[0]) / (r[1] - r[0])
    F_new[-1] = -(U_new[-1] - U_new[-2]) / (r[-1] - r[-2])

    if zero_tail:
        F_new[first_tail:] = 0.0

    # --- write out ---
    with open(outfile, "w") as f:
        f.writelines(header)
        for i in range(n):
            f.write(f"{idx[i]:6d} {r[i]:10.6f} {U_new[i]:15.8e} {F_new[i]:15.8e}\n")


def main():
    p = argparse.ArgumentParser(
        description="Shift a LAMMPS table so U(r_cut)=0, with an optional flat zero tail."
    )
    p.add_argument("infile", help="input LAMMPS table file")
    p.add_argument("outfile", help="output shifted table file")
    p.add_argument("--rcut", type=float, default=20.0, help="cutoff radius where U=0 (default 20.0)")
    p.add_argument("--header-lines", type=int, default=5, help="header lines to copy (default 5)")
    p.add_argument("--mode", choices=["interp", "nearest"], default="interp",
                   help="how to evaluate U(rcut) if rcut not on grid")
    p.add_argument(
        "--no-snap-first-tail-point",
        action="store_true",
        help="do not move the first tabulated point at/above rcut onto exactly rcut",
    )
    p.add_argument(
        "--no-zero-tail",
        action="store_true",
        help="do not force the tail at/above rcut to zero after shifting",
    )
    args = p.parse_args()

    shift_potential_to_zero_at_rcut(
        args.infile, args.outfile,
        r_cut=args.rcut,
        header_lines=args.header_lines,
        mode=args.mode,
        snap_first_tail_point=not args.no_snap_first_tail_point,
        zero_tail=not args.no_zero_tail,
    )

if __name__ == "__main__":
    main()
