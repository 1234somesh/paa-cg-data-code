#!/usr/bin/env python3
"""
Statistical error bars for the CVT osmotic pressure in posm_bar_cyl_NA_K_CS.pdf.

The bar script computes a single number per system from the time-AVERAGED
cylindrical RDF:

    rho(r_beta) = <g_cyl(r_beta)> * rho_bulk,   P_osm = kB*T*rho(r_beta)

which carries no uncertainty. Here the same quantity is rebuilt frame by
frame straight from the trajectories, so rho(r_beta) becomes a time series
whose fluctuations give the error bar.

Per frame, exactly mirroring the averaged calculation in plot_aa_cg.py:
  * chain axis   = PBC-safe circular mean of the monomer/bead transverse
                   positions (one point per frame standing in for the axis)
  * r_xy         = minimum-image transverse distance of each counterion
  * histogram    = same 200 bins over [0, sqrt(a^2+b^2)]
  * rho_frame    = mean over the window bins of counts/dV, with dV from the
                   rectangle-clipped annulus area (so the frame average of
                   this series reproduces the published rho(r_beta) exactly)

The window is the same "last WINDOW_FRAC of the resolved range below
r_beta" used by plot_posm_bar_cyl_NA_K_CS.py.

Error bar: block averaging. The raw per-frame SEM would be far too small --
consecutive frames are strongly correlated -- so the series is blocked over
increasing block lengths until the block-mean SEM stops growing; that
plateau is the reported statistical error, and it is propagated linearly
into P_osm (P_osm is linear in rho).

Run with the MDAnalysis venv, e.g.
    .../radial_density/mda-env/bin/python posm_err_timeseries.py
Writes posm_err_summary.json (+ per-system .npy time series) next to this file.
"""
import itertools
import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

A2NM = 0.1
NBINS = 200
WINDOW_FRAC = 0.10

T_KELVIN = 300.0
ANG3_TO_L = 1e-27
NA_AVO = 6.02214076e23
ANG3_TO_MOLAR = 1.0 / (NA_AVO * ANG3_TO_L)
R_GAS_L_BAR = 0.0831446

MIN_BLOCKS = 20          # smallest number of blocks kept in the block scan
AA_PAA_SEL = "name CH CO2M CH2 OE1 OE2 H1 H2A H2B"

SYSTEMS = {
    "Na": dict(
        aa_top="/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_f/periodic_5.12/NVT_nonwater.gro",  # INPUT: atomistic GROMACS structure (.gro) used as topology -- replace with your own
        aa_trj="/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_freeze_Na/NVT.xtc",  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
        aa_ion="NA",
        cg_trj="/work/skurahatti/multivalent/2d-gel/PAA_dang/SC_CGM/LZ_freeze_Na/traj_paa_cg.lammpstrj",  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own
        aa_dat="sc_cgm_5.12/cyl_rdf_aa.dat",
        cg_dat="sc_cgm_5.12/cyl_rdf_cg.dat",
    ),
    "K": dict(
        aa_top="/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_f/periodic_5.12_K/test/conf.gro",  # INPUT: atomistic GROMACS structure (.gro) used as topology -- replace with your own
        aa_trj="/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_freeze_K/NVT.xtc",  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
        aa_ion="K",
        cg_trj="/work/skurahatti/multivalent/2d-gel/PAA_dang/SC_CGM/LZ_freeze_K/traj_paa_cg.lammpstrj",  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own
        aa_dat="K/sc_cgm_5.12/cyl_rdf_aa.dat",
        cg_dat="K/sc_cgm_5.12/cyl_rdf_cg.dat",
    ),
    "Cs": dict(
        aa_top="/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_f/periodic_5.12_Cs/test/conf.gro",  # INPUT: atomistic GROMACS structure (.gro) used as topology -- replace with your own
        aa_trj="/work/skurahatti/multivalent/SC-CGM/PAA20mer/P-V_protocol/Compress/PULL_freeze_Cs/NVT.xtc",  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
        aa_ion="Cs",
        cg_trj="/work/skurahatti/multivalent/2d-gel/PAA_dang/SC_CGM/LZ_freeze_Cs/traj_paa_cg.lammpstrj",  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own
        aa_dat="Cs/sc_cgm_5.12/cyl_rdf_aa.dat",
        cg_dat="Cs/sc_cgm_5.12/cyl_rdf_cg.dat",
    ),
}


# ---------- geometry helpers (identical to plot_aa_cg.py) ----------

def accessible_area_circle_in_rectangle(r, a, b, n=4000):
    """Area of the circle of radius r that actually lies inside the box."""
    x = np.linspace(0.0, a, n)
    y = np.sqrt(np.clip(r * r - x * x, 0.0, None))
    y = np.minimum(y, b)
    return 4.0 * np.trapezoid(y, x)


def wrap_to_box(x, lo, L):
    return (x - lo) % L + lo


def pbc_circular_mean(x, lo, L):
    """PBC-safe mean position on a periodic axis (chain may straddle the edge)."""
    ang = 2.0 * np.pi * (x - lo) / L
    mean_ang = np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang)))
    if mean_ang < 0:
        mean_ang += 2.0 * np.pi
    return lo + (mean_ang / (2.0 * np.pi)) * L


def build_grid(Lx, Ly, Lz, dat_path=None):
    """Bins, per-bin accessible volume, and the r_beta window mask.

    The window mask is taken from the very .dat grid that
    plot_posm_bar_cyl_NA_K_CS.py averages, so the error bar describes exactly
    the same estimator as the plotted bar. That matters here: the 10% window
    start lands *exactly* on a bin center, so whether that boundary bin is
    included flips on rounding -- the .dat stores it 4e-07 nm low and drops
    it (14 bins), a freshly rebuilt grid keeps it (15 bins). Harmless for CG
    (flat there) but a 0.5% shift for AA, whose profile is ~3x steeper at
    r_beta. Rebuilding the mask independently would silently put the error
    bar on a slightly different number than the bar it annotates.
    """
    a, b = Lx / 2.0, Ly / 2.0
    r_max = np.sqrt(a * a + b * b)
    edges = np.linspace(0.0, r_max, NBINS + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    A_edges = np.array([accessible_area_circle_in_rectangle(r, a, b) for r in edges])
    dV = (A_edges[1:] - A_edges[:-1]) * Lz

    r_beta = 0.5 * Lx

    if dat_path is not None:
        r_dat_nm = np.loadtxt(dat_path, usecols=0)
        if r_dat_nm.size != centers.size:
            raise ValueError(f"{dat_path}: {r_dat_nm.size} bins, expected {centers.size}")
        misfit = np.abs(r_dat_nm - centers * A2NM).max()
        if misfit > 1e-4:
            raise ValueError(f"{dat_path}: grid mismatch {misfit:.2e} nm -- "
                             f"different box or nbins than this trajectory")
        r_grid = r_dat_nm            # nm, as the bar script sees it
        r_beta_grid = r_beta * A2NM
    else:
        r_grid = centers
        r_beta_grid = r_beta

    below = r_grid <= r_beta_grid
    r_below = r_grid[below]
    window_start = r_below.max() - WINDOW_FRAC * (r_below.max() - r_below.min())
    win = below & (r_grid >= window_start)

    if dat_path is not None:
        window_start = window_start / A2NM   # back to Angstrom for reporting

    return edges, centers, dV, win, r_beta, window_start


# ---------- error estimation ----------

def block_scan(x, min_blocks=MIN_BLOCKS):
    """SEM of the mean vs block length; the plateau is the correlated-data SEM."""
    x = np.asarray(x, dtype=float)
    n = x.size
    out = []
    b = 1
    while n // b >= min_blocks:
        m = n // b
        block_means = x[:m * b].reshape(m, b).mean(axis=1)
        out.append((b, m, block_means.std(ddof=1) / np.sqrt(m)))
        b *= 2
    return out


def sem_from_scan(scan):
    """Plateau estimate: the largest SEM in the scan (conservative)."""
    return max(s for _, _, s in scan)


# ---------- per-frame rho(r_beta) ----------

def rho_series_cg(traj, dat_path=None, max_frames=None, report_every=50000):
    """Per-frame rho(r_beta) from a LAMMPS text dump (id type q xu yu zu)."""
    n_header = 9
    with open(traj) as f:
        head = [f.readline() for _ in range(n_header)]
    natoms = int(head[3])
    xlo, xhi = map(float, head[5].split())
    ylo, yhi = map(float, head[6].split())
    zlo, zhi = map(float, head[7].split())
    ncols = len(head[8].split()) - 2          # "ITEM: ATOMS" + column names
    Lx, Ly, Lz = xhi - xlo, yhi - ylo, zhi - zlo

    edges, centers, dV, win, r_beta, window_start = build_grid(Lx, Ly, Lz, dat_path)
    dV_win = dV[win]
    lines_per_frame = n_header + natoms

    rho, nwin = [], []
    t0 = time.time()
    with open(traj) as f:
        while True:
            block = list(itertools.islice(f, lines_per_frame))
            if len(block) < lines_per_frame:
                break
            arr = np.array("".join(block[n_header:]).split(),
                           dtype=np.float64).reshape(natoms, ncols)
            typ = arr[:, 1].astype(int)
            xw = wrap_to_box(arr[:, 3], xlo, Lx)
            yw = wrap_to_box(arr[:, 4], ylo, Ly)

            paa = typ == 1
            ion = typ == 2
            xc = pbc_circular_mean(xw[paa], xlo, Lx)
            yc = pbc_circular_mean(yw[paa], ylo, Ly)

            dx = xw[ion] - xc
            dy = yw[ion] - yc
            dx -= np.rint(dx / Lx) * Lx
            dy -= np.rint(dy / Ly) * Ly
            r_xy = np.sqrt(dx * dx + dy * dy)

            hist, _ = np.histogram(r_xy, bins=edges)
            rho.append(np.mean(hist[win] / dV_win))
            nwin.append(hist[win].sum())

            if len(rho) % report_every == 0:
                print(f"    CG frame {len(rho)} ({time.time()-t0:.0f}s)", flush=True)
            if max_frames is not None and len(rho) >= max_frames:
                break

    return np.array(rho), np.array(nwin, dtype=int), dict(
        Lx=Lx, Ly=Ly, Lz=Lz, r_beta=r_beta, window_start=window_start,
        n_win_bins=int(win.sum()))


def rho_series_aa(top, trj, ion_name, dat_path=None, max_frames=None, report_every=20000):
    """Per-frame rho(r_beta) from the atomistic trajectory."""
    import MDAnalysis as mda

    u = mda.Universe(top, trj)
    paa_atoms = u.select_atoms(AA_PAA_SEL)
    ions = u.select_atoms(f"name {ion_name}")
    if ions.n_atoms == 0:
        raise RuntimeError(f"no atoms matched 'name {ion_name}' in {top}")
    residues = paa_atoms.residues

    # Vectorised replacement for [res.atoms.center_of_mass() for res in residues]
    idx = np.concatenate([r.atoms.ix for r in residues])
    masses = np.concatenate([r.atoms.masses for r in residues])
    counts = np.array([r.atoms.n_atoms for r in residues])
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    mass_sum = np.add.reduceat(masses, starts)

    u.trajectory[0]
    Lx, Ly, Lz = u.dimensions[:3]
    edges, centers, dV, win, r_beta, window_start = build_grid(Lx, Ly, Lz, dat_path)
    dV_win = dV[win]

    ion_ix = ions.ix
    rho, nwin = [], []
    t0 = time.time()
    for ts in u.trajectory:
        Lx_t, Ly_t = ts.dimensions[:2]
        pos = ts.positions[idx]
        com = np.add.reduceat(pos * masses[:, None], starts, axis=0) / mass_sum[:, None]

        xw = wrap_to_box(com[:, 0], 0.0, Lx_t)
        yw = wrap_to_box(com[:, 1], 0.0, Ly_t)
        xc = pbc_circular_mean(xw, 0.0, Lx_t)
        yc = pbc_circular_mean(yw, 0.0, Ly_t)

        ipos = ts.positions[ion_ix]
        dx = ipos[:, 0] - xc
        dy = ipos[:, 1] - yc
        dx -= np.rint(dx / Lx_t) * Lx_t
        dy -= np.rint(dy / Ly_t) * Ly_t
        r_xy = np.sqrt(dx * dx + dy * dy)

        hist, _ = np.histogram(r_xy, bins=edges)
        rho.append(np.mean(hist[win] / dV_win))
        nwin.append(hist[win].sum())

        if len(rho) % report_every == 0:
            print(f"    AA frame {len(rho)} ({time.time()-t0:.0f}s)", flush=True)
        if max_frames is not None and len(rho) >= max_frames:
            break

    return np.array(rho), np.array(nwin, dtype=int), dict(
        Lx=float(Lx), Ly=float(Ly), Lz=float(Lz), r_beta=float(r_beta),
        window_start=float(window_start), n_win_bins=int(win.sum()))


def posm_bar(rho_ang3):
    """rho (ions/Ang^3) -> P_osm (bar); linear, so errors scale by the same factor."""
    return rho_ang3 * ANG3_TO_MOLAR * R_GAS_L_BAR * T_KELVIN


def summarise(tag, rho, nwin, geom):
    scan = block_scan(rho)
    sem = sem_from_scan(scan)
    mean = rho.mean()
    P, dP = posm_bar(mean), posm_bar(sem)

    print(f"  [{tag}] frames={rho.size}  rho={mean:.6e} +/- {sem:.3e} ions/Ang^3 "
          f"(raw std={rho.std(ddof=1):.3e}, naive SEM={rho.std(ddof=1)/np.sqrt(rho.size):.3e})",
          flush=True)
    print(f"        window {geom['window_start']*A2NM:.3f}-{geom['r_beta']*A2NM:.3f} nm "
          f"({geom['n_win_bins']} bins), <ions in window> = {nwin.mean():.3f}", flush=True)
    print(f"        P_osm = {P:.4f} +/- {dP:.4f} bar   "
          f"(block scan: " + ", ".join(f"b={b}:{s:.2e}" for b, _, s in scan[::3]) + ")",
          flush=True)

    return dict(
        n_frames=int(rho.size),
        rho_mean=float(mean),
        rho_sem=float(sem),
        rho_std=float(rho.std(ddof=1)),
        rho_sem_naive=float(rho.std(ddof=1) / np.sqrt(rho.size)),
        P_osm=float(P),
        P_osm_err=float(dP),
        n_ions_window_mean=float(nwin.mean()),
        block_scan=[[int(b), int(m), float(s)] for b, m, s in scan],
        **{k: float(v) if not isinstance(v, int) else v for k, v in geom.items()},
    )


def main():
    results = {}
    for ion, spec in SYSTEMS.items():
        print(f"=== {ion} ===", flush=True)
        results[ion] = {}

        print("  reading AA trajectory ...", flush=True)
        rho_aa, nwin_aa, geom_aa = rho_series_aa(spec["aa_top"], spec["aa_trj"], spec["aa_ion"],
                                                 dat_path=HERE / spec["aa_dat"])
        np.save(HERE / f"posm_rho_series_{ion}_aa.npy", rho_aa)
        results[ion]["AA"] = summarise(f"{ion} AA", rho_aa, nwin_aa, geom_aa)

        print("  reading CG trajectory ...", flush=True)
        rho_cg, nwin_cg, geom_cg = rho_series_cg(spec["cg_trj"], dat_path=HERE / spec["cg_dat"])
        np.save(HERE / f"posm_rho_series_{ion}_cg.npy", rho_cg)
        results[ion]["CG"] = summarise(f"{ion} CG", rho_cg, nwin_cg, geom_cg)

    out = HERE / "posm_err_summary.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out}", flush=True)


if __name__ == "__main__":
    main()
