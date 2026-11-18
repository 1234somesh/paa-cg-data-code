#!/usr/bin/env python3
"""
Extend a LAMMPS pair table down to a tiny inner radius (default r = 0.0001 A)
so that a close contact never trips

    ERROR: Pair distance < table inner cutoff (src/pair_table.cpp)

LAMMPS aborts as soon as any pair in the neighbor list comes closer than the
first tabulated r. The guard point continues the repulsive wall that
fit_short_range.py already put in the core,

    U(r) = U1 * exp(b (r1 - r)),   b = ln(U1/U2) / (r2 - r1),   F = b * U,

fitted on the first two table points, so the wall stays smooth and monotonic.
If the first two points carry no repulsive decay (b <= 0) a flat cap at U1
with F = 0 would be useless, so the script refuses instead.

The table must have the usual 5-line header (comment, blank, keyword,
"N <npts> ...", blank); the N count is bumped and the index column renumbered.

If the table already starts at r_min (a guarded table coming back through
ibi.py -> fit_short_range.py -> shift_to_zero.py), the guard point is not
added again but re-derived from the next two points. fit_short_range.py's
one-sided force across the wide r_min -> r_2 gap otherwise leaves the guard
with almost no force (F ~ 0), i.e. a flat cap.

Usage
-----
    python3 add_inner_guard.py IN_TABLE OUT_TABLE [--rmin 0.0001]
"""
import argparse
import re

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--rmin", type=float, default=1e-4,
                    help="inner guard radius in A (default 0.0001)")
    args = ap.parse_args()

    with open(args.infile) as f:
        header = [f.readline() for _ in range(5)]
        data = np.loadtxt(f)

    # drop an existing guard point so it is rebuilt from the real wall
    rebuilt = data[0, 1] <= args.rmin
    if rebuilt:
        data = data[1:]
        if data[0, 1] <= args.rmin:
            raise RuntimeError(f"{args.infile}: more than one point at r <= {args.rmin} A")

    r, U = data[:, 1], data[:, 2]
    b = np.log(U[0] / U[1]) / (r[1] - r[0]) if U[0] > 0 and U[1] > 0 else -1.0
    if b <= 0.0:
        raise RuntimeError(
            f"{args.infile}: first two points (U={U[0]:.4g}, {U[1]:.4g}) are not a "
            f"repulsive wall; run fit_short_range.py on it first")
    U_g = U[0] * np.exp(b * (r[0] - args.rmin))
    guard = np.array([[0, args.rmin, U_g, b * U_g]])
    new_data = np.vstack([guard, data])
    new_data[:, 0] = np.arange(1, len(new_data) + 1)
    print(f"{args.infile}: {'rebuilt' if rebuilt else 'added'} guard r = {args.rmin} A  "
          f"U = {U_g:.4g} kcal/mol  F = {b * U_g:.4g} kcal/mol/A  (b = {b:.4g} /A "
          f"from r = {r[0]:.6f}, {r[1]:.6f} A)")

    # bump "N <npts>" in the header, keeping any trailing keywords
    for i, line in enumerate(header):
        if re.match(r"^\s*N\s+\d+", line):
            header[i] = re.sub(r"^(\s*N\s+)\d+", rf"\g<1>{len(new_data)}", line)
            break
    else:
        raise RuntimeError(f"{args.infile}: no 'N <npts>' line in the 5-line header")

    with open(args.outfile, "w") as f:
        f.writelines(header)
        np.savetxt(f, new_data, fmt="%8d %10.6f %15.8e %15.8e")


if __name__ == "__main__":
    main()
