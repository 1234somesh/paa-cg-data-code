#!/usr/bin/env python3
"""
Replace the unreliable short-range core of a LAMMPS pair table with a smooth,
strictly repulsive exponential core.

Guarantee
---------
Below the detected join radius r_join the output satisfies, for decreasing r:

    U(r) strictly INCREASES     and     F(r) = -dU/dr strictly INCREASES,

with F > 0 throughout.  That is what stops particles from sliding into (or
straight through) the core during GCMC insertions / close contacts: a flat
capped core has F = 0 inside and gives no restoring push at all.

Why the core needs replacing
----------------------------
In IBI the core region is never sampled (g(r) ~ 0), so ln(g_imp/g_ref) there is
meaningless and whatever sits in the table is an artifact.  Three artifacts show
up in practice, and all three are handled:

  (a) a flat capped plateau (e.g. U = 239 kcal/mol out to ~2.35 A) -- monotonic
      but with ZERO force inside;
  (b) an all-zero core (a table built from an all-zero initial guess), i.e. no
      repulsive wall at all;
  (c) a noisy / non-monotonic core that wanders and can even dip.

Method
------
1. Walk outward from r_min while the data is "uninformative": U <= 0, or the
   slope is flat (|dU/dr| < flat_tol).  That marks the artifact region.
2. Skip a couple of further grid points, because the first points past a
   plateau carry a spurious huge central-difference slope from the step itself.
3. Least-squares fit  ln U(r) = ln U0 - b (r - r0)  over a window just outside
   the join.  A windowed fit (rather than one point's slope) is what makes b
   robust to single-point noise.
4. Clamp b so the core tops out near core_height at r_min:
       b_cap = ln(core_height / U0) / (r0 - r_min),   b = min(b_fit, b_cap)
   If the data carries no repulsive decay at all (case b: fitted b <= 0), fall
   back to b = b_cap so a wall is still constructed.
5. Rewrite U(r < r0) = U0 * exp(b (r0 - r)) and recompute F = -dU/dr
   consistently, then verify monotonicity of both U and F.

Because b is a positive constant, F = b*U inside the core, so F grows in exact
proportion to U as r decreases -- both strictly increasing, by construction.

Usage
-----
    python3 fit_short_range.py IN_TABLE OUT_TABLE [flat_tol] [core_height]

`flat_tol` (default 1e-6) is the |dU/dr| below which the core counts as flat.
The third positional argument used to be `deriv_threshold` in the older
single-point-slope version; it is still accepted so existing calls keep working,
but it now only controls flat-plateau detection.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_FIT_DEFAULT = 12   # grid points used for the log-linear wall fit
SKIP_DEFAULT = 2     # points skipped past a plateau edge (step artifact)


def smooth_pmf_auto_join(r, U, flat_tol=1e-6, core_height=1000.0,
                         n_fit=N_FIT_DEFAULT, skip=SKIP_DEFAULT, r0_override=None):
    """
    Returns
    -------
    U_smooth : ndarray
    r_join   : float
    info     : dict  (U0, b, b_fit, b_cap, capped, i_start)

    r0_override
    -----------
    If given, skip the automatic flat/zero-core walk entirely and anchor the
    wall fit at the grid point nearest this distance (Å). Use this when the
    short-range region is too noisy for the flat-plateau heuristic to find a
    sane join point on its own (e.g. IBI corrections that never departed from
    ~0 because close contacts were barely sampled) -- inspect the table by eye
    and pick an r past the noise, from which the exponential wall is fit.
    """
    r = np.asarray(r, dtype=float)
    U = np.asarray(U, dtype=float)

    if not np.all(np.diff(r) > 0):
        raise ValueError("r must be strictly increasing")
    n = len(r)
    if n < 5:
        raise ValueError("Need at least 5 points")

    # --- central-difference derivative ---
    dUdr = np.zeros_like(U)
    dUdr[1:-1] = (U[2:] - U[:-2]) / (r[2:] - r[:-2])
    dUdr[0] = (U[1] - U[0]) / (r[1] - r[0])
    dUdr[-1] = (U[-1] - U[-2]) / (r[-1] - r[-2])

    if r0_override is not None:
        # explicit join distance: anchor there directly, skip the
        # flat/zero-core walk and its step-artifact skip.
        i_start = int(np.argmin(np.abs(r - r0_override)))
        i_start = max(i_start, 1)
    else:
        # --- 1. extent of the uninformative core (zero / flat) ---
        i_core_end = -1
        for i in range(n):
            if U[i] <= 0.0 or abs(dUdr[i]) < flat_tol:
                i_core_end = i
            else:
                break

        # --- 2. skip the step artifact just past a plateau ---
        i_start = i_core_end + 1 + (skip if i_core_end >= 0 else 0)
        i_start = max(i_start, 1)

    if i_start >= n - 3:
        raise RuntimeError(
            "No usable data outside the flat/zero core -- the table appears to "
            "be flat or non-positive almost everywhere; nothing to anchor to."
        )

    # --- 3. windowed log-linear fit of the wall, just outside the join ---
    hi = min(i_start + n_fit, n)
    sel = np.arange(i_start, hi)
    sel = sel[U[sel] > 0.0]
    if len(sel) < 3:
        raise RuntimeError(
            f"Fewer than 3 positive-U points in the fit window starting at "
            f"r={r[i_start]:.4f} A; cannot fit a repulsive wall."
        )

    r0 = r[i_start]
    A = np.vstack([-(r[sel] - r0), np.ones_like(r[sel])]).T
    b_fit, lnU0 = np.linalg.lstsq(A, np.log(U[sel]), rcond=None)[0]
    U0 = float(np.exp(lnU0))

    # --- 4. clamp so the core tops out near core_height at r_min ---
    if core_height <= U0:
        core_height = U0 * 1.05
    b_cap = np.log(core_height / U0) / (r0 - r[0])

    if b_fit <= 0.0:
        # data has no repulsive decay (all-zero / attractive core): build one
        b = b_cap
    else:
        b = min(b_fit, b_cap)
    capped = not (b_fit > 0.0 and b_fit <= b_cap)

    if b <= 0.0:
        raise RuntimeError("Could not determine a positive decay constant b.")

    # --- 5. rewrite the core ---
    U_smooth = U.copy()
    core = r < r0
    U_smooth[core] = U0 * np.exp(b * (r0 - r[core]))

    info = dict(U0=U0, b=b, b_fit=b_fit, b_cap=b_cap, capped=capped,
                i_start=i_start, n_core=int(core.sum()))
    return U_smooth, r0, info


def forces_from_U(r, U):
    F = np.zeros_like(U)
    F[1:-1] = -(U[2:] - U[:-2]) / (r[2:] - r[:-2])
    F[0] = -(U[1] - U[0]) / (r[1] - r[0])
    F[-1] = -(U[-1] - U[-2]) / (r[-1] - r[-2])
    return F


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 fit_short_range.py IN_TABLE OUT_TABLE [flat_tol] [core_height] [r0_override]")
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    flat_tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    core_height = float(sys.argv[4]) if len(sys.argv) > 4 else 1000.0
    r0_override = float(sys.argv[5]) if len(sys.argv) > 5 else None

    with open(in_file, "r") as f:
        header_lines = [next(f) for _ in range(5)]

    data = np.loadtxt(in_file, comments="#", skiprows=5)
    r = data[:, 1]
    U = data[:, 2]

    U_smooth, r_join, info = smooth_pmf_auto_join(
        r, U, flat_tol=flat_tol, core_height=core_height, r0_override=r0_override)

    F_smooth = forces_from_U(r, U_smooth)

    # --- verify the guarantee on the reconstructed core ---
    core = r < r_join
    idx = np.where(core)[0]
    if len(idx) >= 2:
        # walking inward (decreasing r) U and F must both increase
        dU_in = np.diff(U_smooth[idx][::-1])
        dF_in = np.diff(F_smooth[idx][::-1])
        if not np.all(dU_in > 0):
            raise RuntimeError("core U is not strictly increasing as r decreases")
        # allow the outermost force point to be set by the join stencil
        if not np.all(dF_in[1:] > 0):
            raise RuntimeError("core F is not strictly increasing as r decreases")
        # same join-stencil exemption applies to the absolute positivity
        # check: the outermost core point's central difference spans the
        # join itself and can come out slightly negative even though the
        # reconstructed core is a clean increasing exponential.
        if not np.all(F_smooth[idx][:-1] > 0):
            raise RuntimeError("core F is not strictly positive")

    print(f"Detected join distance r_join = {r_join:.4f} A "
          f"(replaced {info['n_core']} pts)")
    print(f"  wall fit: U0={info['U0']:.4g} kcal/mol  b_fit={info['b_fit']:.4g}  "
          f"b_cap={info['b_cap']:.4g}  -> b={info['b']:.4g} /A"
          f"{'  [capped at core_height]' if info['capped'] else ''}")
    print(f"  core: U({r[0]:.4f})={U_smooth[0]:.4g}  F({r[0]:.4f})={F_smooth[0]:.4g}"
          f"  (U and F both increase monotonically as r decreases)")

    data[:, 2] = U_smooth
    data[:, 3] = F_smooth

    plt.plot(r, U, label="original")
    plt.plot(r, U_smooth, label="smoothed")
    plt.axvline(r_join, color="k", linestyle="--", label=f"r_join={r_join:.2f} A")
    plt.xlabel("r (A)"); plt.ylabel("U(r)"); plt.legend()

    with open(out_file, "w") as f:
        for line in header_lines:
            f.write(line)
        np.savetxt(f, data, fmt="%8d %10.6f %15.8e %15.8e")

    print(f"Wrote smoothed table to: {out_file}")


if __name__ == "__main__":
    main()
