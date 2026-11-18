#!/usr/bin/env python3
"""
Free chain (no cell-gel constraint): spherical RDF between the PAA 20mer
centre of mass and its 20 counterions, and the Contact-Value-Theorem
osmotic pressure at the cell boundary r_beta = L/2, with error bars.

This is the spherical counterpart of ../cylindrical_analysis. The free chain
has no imposed axis, so the natural reference is a single point -- the chain
CoM -- and the natural symmetry is spherical:

    rho(r_beta) = g_sph(r_beta) * rho_bulk,   P_osm = kB*T*rho(r_beta)

Per frame:
  1. the chain is made whole across PBC and its mass-weighted CoM computed;
  2. every counterion is taken at its nearest periodic image of that CoM;
  3. the CoM of the whole system (chain + counterions) is shifted to the box
     centre, as specified for this analysis;
  4. ion distances from the chain CoM are histogrammed into the spherical
     shells.
Step 3 does not change any distance -- minimum-image separations are
translation invariant -- so it cannot move the RDF or P_osm. It is done
because it defines the frame the "boundary at half the box length" refers
to, and the per-frame chain-CoM offset from the box centre is reported as a
diagnostic.

Shell volumes are the volume of (sphere of radius r) INTERSECTED WITH the
minimum-image cube of side L centred on the reference point -- not the bare
4/3*pi*r^3. Below r = L/2 the two agree exactly (the sphere still fits), so
the CVT value at r_beta = L/2 is unaffected; beyond it the bare formula would
overcount the space that actually exists, which would corrupt the tail of
g(r) that gets plotted next to it.

Error bar: block averaging of the per-frame rho(r_beta) series, exactly as
in ../cylindrical_analysis/posm_err_timeseries.py -- consecutive frames are
correlated, so the naive per-frame SEM is far too small.

Run with the MDAnalysis venv:
    .../radial_density/mda-env/bin/python sph_posm_analysis.py [ion ...]
"""
import itertools
import json
import os
import sys
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

MIN_BLOCKS = 20

# Set SPH_MAX_FRAMES=N to process only the first N frames of each
# trajectory (quick trial run); unset/0 processes everything.
MAX_FRAMES = int(os.environ.get("SPH_MAX_FRAMES", "0")) or None

AA_BASE = "/beegfs/ant_work/skurahatti/work/multivalent/2d-gel/PAA_dang"
CG_BASE = "/work/skurahatti/multivalent/2d-gel/PAA_dang/parametrization"

# The stripped 182-atom topology labels every counterion "NA", whichever salt
# it really is, so ion masses are set explicitly here. They only enter the
# (distance-invariant) system-CoM shift, but leaving them wrong would make
# that diagnostic meaningless.
ION_MASS = {"Na": 22.98977, "K": 39.0983, "Cs": 132.90545}
CG_BEAD_MASS = 71.0

SYSTEMS = {
    "Na": dict(
        aa_top=f"{AA_BASE}/NA/box_l_6_0/NVT_hbonds.tpr",  # INPUT: atomistic GROMACS run input (.tpr) used as topology -- replace with your own
        aa_trj=f"{AA_BASE}/NA/box_l_6_0/combined_NVT.xtc",  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
        cg_trj=f"{CG_BASE}/box_l_6nm/Try14/step_11/traj_paa_cg.lammpstrj",  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own
    ),
    "K": dict(
        aa_top=f"{AA_BASE}/K/box_l_6_0_last/NVT_hbonds.tpr",  # INPUT: atomistic GROMACS run input (.tpr) used as topology -- replace with your own
        aa_trj=f"{AA_BASE}/K/box_l_6_0_last/combined_NVT.xtc",  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
        # Try14/step_10 supersedes Try13/step_80: same frame count (400001),
        # but P_osm = 1.361 +/- 0.006 bar vs 1.508 +/- 0.007 -- roughly halves
        # the gap to the AA reference (1.266 +/- 0.009), consistent with
        # Try14 reflecting a later, better-converged/fixed K parametrization.
        cg_trj=f"{CG_BASE}/K_last/Try14/step_10/traj_paa_cg.lammpstrj",  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own
    ),
    "Cs": dict(
        aa_top=f"{AA_BASE}/CS/box_l_6_last/NVT_hbonds.tpr",  # INPUT: atomistic GROMACS run input (.tpr) used as topology -- replace with your own
        aa_trj=f"{AA_BASE}/CS/box_l_6_last/combined_NVT.xtc",  # INPUT: atomistic GROMACS trajectory (.xtc/.trr) -- replace with the path to your own AA trajectory
        cg_trj=f"{CG_BASE}/Cs_last/Try13/step_84/traj_paa_cg.lammpstrj",  # INPUT: CG LAMMPS dump trajectory (.lammpstrj; columns id type q x y z) -- replace with your own
    ),
}


# ---------------- geometry ----------------

def sphere_cube_volume(r, a, n=400):
    """Volume of {|x| < r} intersected with the cube of half-width a.

    The minimum-image convention puts every counterion inside the cube of
    side 2a centred on the reference point, so this -- not 4/3*pi*r^3 -- is
    the space a shell of radius r actually has available. Equal to the full
    sphere for r <= a, and saturating at the cube volume for r >= a*sqrt(3).
    """
    if r <= 0.0:
        return 0.0
    if r <= a:
        # sphere still fits inside the cube: exact, no quadrature error
        return 4.0 / 3.0 * np.pi * r ** 3
    if r >= a * np.sqrt(3.0):
        return 8.0 * a ** 3
    xs = np.linspace(0.0, min(r, a), n)
    inner = np.empty(n)
    for i, x in enumerate(xs):
        rho2 = r * r - x * x
        if rho2 <= 0.0:
            inner[i] = 0.0
            continue
        ys = np.linspace(0.0, min(np.sqrt(rho2), a), n)
        h = np.minimum(np.sqrt(np.clip(rho2 - ys * ys, 0.0, None)), a)
        inner[i] = np.trapezoid(h, ys)
    return 8.0 * np.trapezoid(inner, xs)


def build_grid(L):
    """Bins out to the box-corner reach, shell volumes, and the r_beta window."""
    a = L / 2.0
    r_max = a * np.sqrt(3.0)
    edges = np.linspace(0.0, r_max, NBINS + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    V_edges = np.array([sphere_cube_volume(r, a) for r in edges])
    dV = V_edges[1:] - V_edges[:-1]

    r_beta = a
    below = centers <= r_beta
    r_below = centers[below]
    window_start = r_below.max() - WINDOW_FRAC * (r_below.max() - r_below.min())
    win = below & (centers >= window_start)

    return edges, centers, dV, win, r_beta, window_start


def chain_com(pos, masses, L, already_unwrapped):
    """Mass-weighted chain CoM, PBC-safe.

    Wrapped coordinates are made whole by walking along the atom ordering and
    accumulating minimum-image steps between CONSECUTIVE atoms. Consecutive
    atoms are always within a bond length of each other, so each step is
    unambiguous and the result is exact however far the chain is stretched.

    The obvious cheaper alternative -- taking every atom at its nearest image
    of atom 0 -- is only valid while the chain spans less than L/2, which a
    free chain in a 6 nm box does violate (~0.3% of frames here). On exactly
    those frames it was found to misplace the CoM by up to 7.9 A, whereas
    this walk reproduces MDAnalysis's bond-graph unwrap() to ~5e-06 A.

    Returns (com, max_consecutive_step). That second value is the validity
    check: it must stay far below L/2, otherwise a step was ambiguous and the
    unwrap cannot be trusted.
    """
    if already_unwrapped:
        return (pos * masses[:, None]).sum(0) / masses.sum(), 0.0
    d = np.diff(pos, axis=0)
    d -= L * np.rint(d / L)
    whole = np.empty_like(pos)
    whole[0] = pos[0]
    whole[1:] = pos[0] + np.cumsum(d, axis=0)
    return ((whole * masses[:, None]).sum(0) / masses.sum(),
            float(np.abs(d).max()))


def frame_distances(pos_paa, m_paa, pos_ion, m_ion, L, already_unwrapped):
    """Ion distances from the chain CoM, after centring the system CoM.

    Returns (r_ion, chain_offset_from_box_centre, max_consecutive_step).
    """
    com, step = chain_com(pos_paa, m_paa, L, already_unwrapped)

    # counterions at their nearest image of the chain CoM
    d_ion = pos_ion - com
    d_ion -= L * np.rint(d_ion / L)

    # system CoM relative to the chain CoM, then shift it onto the box centre.
    # Distances below are computed from d_ion, which this shift leaves
    # untouched; it only fixes where the configuration sits in the box.
    m_tot = m_paa.sum() + m_ion.sum()
    com_sys_rel = (d_ion * m_ion[:, None]).sum(0) / m_tot
    chain_offset = -com_sys_rel        # chain CoM relative to box centre after the shift

    return np.sqrt((d_ion ** 2).sum(1)), np.linalg.norm(chain_offset), step


# ---------------- error estimation ----------------

def block_scan(x, min_blocks=MIN_BLOCKS):
    x = np.asarray(x, dtype=float)
    n = x.size
    out = []
    b = 1
    while n // b >= min_blocks:
        m = n // b
        bm = x[:m * b].reshape(m, b).mean(axis=1)
        out.append((b, m, bm.std(ddof=1) / np.sqrt(m)))
        b *= 2
    return out


def sem_from_scan(scan):
    return max(s for _, _, s in scan)


def posm_bar(rho_ang3):
    return rho_ang3 * ANG3_TO_MOLAR * R_GAS_L_BAR * T_KELVIN


# ---------------- per-system drivers ----------------

def run_cg(traj, ion, report_every=100000):
    """LAMMPS text dump, columns: id type q xu yu zu (unwrapped)."""
    n_header = 9
    with open(traj) as f:
        head = [f.readline() for _ in range(n_header)]
    natoms = int(head[3])
    xlo, xhi = map(float, head[5].split())
    ylo, yhi = map(float, head[6].split())
    zlo, zhi = map(float, head[7].split())
    ncols = len(head[8].split()) - 2
    L = np.array([xhi - xlo, yhi - ylo, zhi - zlo])
    if not np.allclose(L, L[0]):
        raise ValueError(f"{traj}: box is not cubic ({L}) -- r_beta = L/2 is ambiguous")

    edges, centers, dV, win, r_beta, window_start = build_grid(L[0])
    dV_win = dV[win]
    lines_per_frame = n_header + natoms

    counts = np.zeros(NBINS)
    rho_series, offsets = [], []
    n_frames = 0
    t0 = time.time()
    with open(traj) as f:
        while True:
            block = list(itertools.islice(f, lines_per_frame))
            if len(block) < lines_per_frame:
                break
            arr = np.array("".join(block[n_header:]).split(),
                           dtype=np.float64).reshape(natoms, ncols)
            typ = arr[:, 1].astype(int)
            pos = arr[:, 3:6]
            paa = typ == 1
            ionm = typ == 2

            m_paa = np.full(paa.sum(), CG_BEAD_MASS)
            m_ion = np.full(ionm.sum(), ION_MASS[ion])
            r_ion, off, _ = frame_distances(pos[paa], m_paa, pos[ionm], m_ion,
                                            L, already_unwrapped=True)

            h, _ = np.histogram(r_ion, bins=edges)
            counts += h
            rho_series.append(np.mean(h[win] / dV_win))
            offsets.append(off)
            n_frames += 1
            if n_frames % report_every == 0:
                print(f"    CG frame {n_frames} ({time.time()-t0:.0f}s)", flush=True)
            if MAX_FRAMES and n_frames >= MAX_FRAMES:
                break

    n_ion = int(ionm.sum())
    return _finish(counts, np.array(rho_series), np.array(offsets), n_frames,
                   n_ion, L[0], edges, centers, dV, win, r_beta, window_start,
                   n_bad_span=0)


def run_aa(top, trj, ion, report_every=200000):
    import MDAnalysis as mda

    u = mda.Universe(top, trj)
    paa = u.select_atoms("resname PAA")
    ions = u.atoms - paa
    if paa.n_atoms == 0 or ions.n_atoms == 0:
        raise RuntimeError(f"{top}: PAA={paa.n_atoms} ions={ions.n_atoms}")

    u.trajectory[0]
    L = np.array(u.dimensions[:3], dtype=float)
    if not np.allclose(L, L[0], rtol=1e-4):
        raise ValueError(f"{trj}: box is not cubic ({L})")

    edges, centers, dV, win, r_beta, window_start = build_grid(L[0])
    dV_win = dV[win]

    ix_paa, ix_ion = paa.ix, ions.ix
    m_paa = paa.masses.copy()
    m_ion = np.full(ions.n_atoms, ION_MASS[ion])
    half = L[0] / 2.0

    counts = np.zeros(NBINS)
    rho_series, offsets = [], []
    n_frames = 0
    worst_step = 0.0
    n_bad_span = 0
    t0 = time.time()
    for ts in u.trajectory:
        p = ts.positions
        r_ion, off, step = frame_distances(p[ix_paa], m_paa, p[ix_ion], m_ion,
                                           L, already_unwrapped=False)
        worst_step = max(worst_step, step)
        if step >= half:
            n_bad_span += 1
        h, _ = np.histogram(r_ion, bins=edges)
        counts += h
        rho_series.append(np.mean(h[win] / dV_win))
        offsets.append(off)
        n_frames += 1
        if n_frames % report_every == 0:
            print(f"    AA frame {n_frames} ({time.time()-t0:.0f}s)", flush=True)
        if MAX_FRAMES and n_frames >= MAX_FRAMES:
            break

    res = _finish(counts, np.array(rho_series), np.array(offsets), n_frames,
                  ions.n_atoms, L[0], edges, centers, dV, win, r_beta,
                  window_start, n_bad_span)
    res["worst_unwrap_step"] = worst_step
    return res


def _finish(counts, rho_series, offsets, n_frames, n_ion, L,
            edges, centers, dV, win, r_beta, window_start, n_bad_span):
    rho_bulk = n_ion / L ** 3
    g = counts / (n_frames * rho_bulk * dV)
    cn = np.cumsum(counts) / n_frames

    scan = block_scan(rho_series)
    sem = sem_from_scan(scan)
    mean = rho_series.mean()

    return dict(
        r_nm=centers * A2NM, g=g, cn=cn, rho_series=rho_series,
        n_frames=n_frames, n_ion=n_ion, L=L, rho_bulk=rho_bulk,
        r_beta_nm=r_beta * A2NM, window_start_nm=window_start * A2NM,
        n_win_bins=int(win.sum()),
        rho_mean=mean, rho_sem=sem, rho_std=rho_series.std(ddof=1),
        rho_sem_naive=rho_series.std(ddof=1) / np.sqrt(rho_series.size),
        P_osm=posm_bar(mean), P_osm_err=posm_bar(sem),
        chain_offset_mean_nm=float(offsets.mean() * A2NM),
        chain_offset_max_nm=float(offsets.max() * A2NM),
        n_bad_span=n_bad_span,
        block_scan=[[int(b), int(m), float(s)] for b, m, s in scan],
    )


def report(tag, res):
    print(f"  [{tag}] frames={res['n_frames']}  rho(r_beta)={res['rho_mean']:.6e} "
          f"+/- {res['rho_sem']:.3e} ions/Ang^3 "
          f"(naive SEM {res['rho_sem_naive']:.2e})", flush=True)
    print(f"        r_beta={res['r_beta_nm']:.3f} nm, window "
          f"{res['window_start_nm']:.3f}-{res['r_beta_nm']:.3f} nm "
          f"({res['n_win_bins']} bins), rho_bulk={res['rho_bulk']:.4e}", flush=True)
    print(f"        chain CoM offset from box centre: mean "
          f"{res['chain_offset_mean_nm']:.3f} nm, max {res['chain_offset_max_nm']:.3f} nm"
          + (f"  | WARNING: {res['n_bad_span']} frames with an ambiguous "
             f"unwrap step >= L/2" if res['n_bad_span'] else ""), flush=True)
    if "worst_unwrap_step" in res:
        print(f"        unwrap check: largest consecutive-atom step "
              f"{res['worst_unwrap_step']:.2f} A (must be << L/2 = {res['L']/2:.1f} A)",
              flush=True)
    print(f"        P_osm = {res['P_osm']:.4f} +/- {res['P_osm_err']:.4f} bar", flush=True)


def main():
    wanted = sys.argv[1:] or list(SYSTEMS)
    summary_path = HERE / "sph_posm_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}

    for ion in wanted:
        spec = SYSTEMS[ion]
        print(f"=== {ion} ===", flush=True)
        summary.setdefault(ion, {})

        for tag, fn in (("AA", lambda: run_aa(spec["aa_top"], spec["aa_trj"], ion)),
                        ("CG", lambda: run_cg(spec["cg_trj"], ion))):
            print(f"  reading {tag} trajectory ...", flush=True)
            res = fn()
            report(f"{ion} {tag}", res)

            low = tag.lower()
            np.savetxt(HERE / f"sph_rdf_{ion}_{low}.dat",
                       np.column_stack([res["r_nm"], res["g"]]),
                       header=f"r(nm) g_sph(r)  [{ion} {tag}, CoM-of-chain to counterions]")
            np.savetxt(HERE / f"sph_cn_{ion}_{low}.dat",
                       np.column_stack([res["r_nm"], res["cn"]]),
                       header=f"r(nm) CN(r)  [{ion} {tag}]")
            np.save(HERE / f"sph_rho_series_{ion}_{low}.npy", res["rho_series"])

            summary[ion][tag] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                                 for k, v in res.items()
                                 if k not in ("r_nm", "g", "cn", "rho_series")}
            summary_path.write_text(json.dumps(summary, indent=2))

    print(f"\nWrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
