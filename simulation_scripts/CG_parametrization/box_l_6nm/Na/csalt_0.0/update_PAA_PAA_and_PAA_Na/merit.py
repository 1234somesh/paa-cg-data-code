#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Use np.trapezoid if present (no deprecation warning), else np.trapz
_trapz = getattr(np, "trapezoid", np.trapz)

# ---- your existing knobs (unchanged) ----
# IBI iterations step_01..step_56 (step_46 crashed and has empty RDFs; it
# is skipped and the lines connect step_45 directly to step_47)
step_sources = list(range(1, 57))
base = Path(__file__).parent
skiprows_rdf = 4
usecols_rdf = (1, 2)
usecols_target = (0, 1)
weight_sigma = None  # set a float to enable w(r)=exp(-r/sigma)

def merit_from_curves(r, g_cg, r_t, g_t, sigma=None):
    g_ti = np.interp(r, r_t, g_t)
    w = 1.0 if sigma is None else np.exp(-r / sigma)
    num = _trapz((g_cg - g_ti)**2 * w, x=r)
    den = _trapz((g_ti**2) * w, x=r)
    return num / den if den != 0 else np.nan

def load_step_rdf(step_dir, rdf_name):
    # LAMMPS appends a new block header each time the run restarts, producing
    # multiple "timestep nrows" lines with only 2 columns inside the file.
    # Reading the last 200 lines always gives the final averaged block.
    with open(step_dir / rdf_name) as f:
        lines = f.readlines()
    if not lines:  # failed run (e.g. step_46) -> no data
        return None, None
    r, g = np.loadtxt(lines[-200:], usecols=usecols_rdf, unpack=True)
    return r, g

# Targets
rt_na, gt_na = np.loadtxt(base / "targets" / "PAA-NA.dist.dat",
                          usecols=usecols_target, unpack=True)
rt_mm, gt_mm = np.loadtxt(base / "targets" / "PAA-PAA.dist.dat",
                          usecols=usecols_target, unpack=True)

iters = np.arange(1, len(step_sources) + 1)
merit_na, merit_mm = [], []

for step_id in step_sources:
    step_dir_na = base / f"step_{step_id:02d}"
    r_na, g_na = load_step_rdf(step_dir_na, "rdf_PAA_NA.dat")
    merit_na.append(np.nan if r_na is None else
                    merit_from_curves(r_na, g_na, rt_na, gt_na, sigma=weight_sigma))

    step_dir_mm = base / f"step_{step_id:02d}"
    r_mm, g_mm = load_step_rdf(step_dir_mm, "rdf_PAA_PAA.dat")
    merit_mm.append(np.nan if r_mm is None else
                    merit_from_curves(r_mm, g_mm, rt_mm, gt_mm, sigma=weight_sigma))

merit_na = np.asarray(merit_na)
merit_mm = np.asarray(merit_mm)

# Drop failed iterations so the lines connect across them without a gap
ok = np.isfinite(merit_na) & np.isfinite(merit_mm)
iters, merit_na, merit_mm = iters[ok], merit_na[ok], merit_mm[ok]

#np.savetxt("merit_series.dat",
#           np.column_stack([iters, merit_na, merit_mm]),
#           header="iter  merit_PAA-NA  merit_PAA-PAA", comments="")

plt.figure(figsize=(5.2, 3.8), dpi=140)
plt.plot(iters, merit_na, marker="o", lw=1.5, label="PAA–NA")
plt.plot(iters, merit_mm, marker="s", lw=1.5, label="PAA–PAA")
plt.xlabel("Iteration", fontsize=16, fontfamily="serif", fontweight="semibold")
plt.ylabel(r"$f_\mathrm{merit}$", fontsize=15, fontfamily="serif", fontweight="semibold")
plt.yscale("log")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("merit_series.pdf", dpi=300, format="pdf")

