#!/usr/bin/env python3
"""
Compare CG RDFs from a given step against the target (AA) RDFs, for every
pair present in the system.

Usage:
    python3 plot_rdf_comparison.py step_02
    python3 plot_rdf_comparison.py step_05 --save

Which pairs get plotted is auto-detected from the rdf_*.dat files actually
present in the step directory, so the same script works for the salt-free
runs (PAA-PAA, PAA-K, K-K) and for the added-salt runs where all six pairs
(+ PAA-CL, K-CL, CL-CL) are computed.

Each panel is annotated with that pair's merit value

    f_merit = int (g_CG - g_target)^2 dr / int g_target^2 dr

using the same definition as merit.py, so you can see at a glance which pairs
are converged (< 0.003) and which are not.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# numpy >= 2.0 removed np.trapz (renamed np.trapezoid); numpy < 2.0 has only
# np.trapz. Prefer np.trapz, fall back to np.trapezoid, so this runs unchanged
# under both -- same idea as the shim already used in merit.py.
try:
    _trapz = np.trapz
except AttributeError:  # numpy >= 2.0
    _trapz = np.trapezoid

if len(sys.argv) < 2:
    print("Usage: python3 plot_rdf_comparison.py <step_dir> [--save]")
    print("Example: python3 plot_rdf_comparison.py step_02")
    sys.exit(1)

step = sys.argv[1].rstrip("/")
save = "--save" in sys.argv
step_dir = Path(step)
targets_dir = Path("targets")

if not step_dir.is_dir():
    sys.exit(f"ERROR: {step_dir} not found")


# --- read functions ---

def read_lammps_rdf(path):
    """Read LAMMPS fix ave/time RDF output, return (r, g_avg) in Angstrom.
    Handles both multi-block format (uncommented headers) and
    single-block format (headers commented out with #)."""
    rs = None
    g_blocks = []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if (not line) or line.startswith("#") or line.startswith("Time-averaged"):
                continue
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                nrows = int(parts[1])
                r_block, g_block = [], []
                for _ in range(nrows):
                    row = f.readline()
                    if not row:
                        break
                    cols = row.split()
                    if len(cols) < 3:
                        continue
                    r_block.append(float(cols[1]))
                    g_block.append(float(cols[2]))
                r_block = np.array(r_block)
                g_block = np.array(g_block)
                if rs is None:
                    rs = r_block
                g_blocks.append(g_block)
    if g_blocks:
        return rs, np.mean(np.stack(g_blocks), axis=0)
    # Fallback: plain columnar data (block headers commented out)
    arr = np.loadtxt(path, comments="#")
    return arr[:, 1], arr[:, 2]


def read_target_rdf(path):
    """Read target RDF (cols 0,1), return r in Angstrom, g(r)."""
    arr = np.loadtxt(path, usecols=(0, 1))
    return arr[:, 0], arr[:, 1]


def find_target(candidates):
    """First existing target file from a list of candidate names.
    .dist.dat is preferred; .dist.new is an accepted fallback because in these
    trees the two are identical (convert_dist.py's x10 was undone), and some
    pairs were never converted."""
    for name in candidates:
        p = targets_dir / name
        if p.exists():
            return p
    return None


def merit(r, g_cg, r_t, g_t):
    """Same definition as merit.py."""
    g_ti = np.interp(r, r_t, g_t)
    num = _trapz((g_cg - g_ti) ** 2, x=r)
    den = _trapz(g_ti ** 2, x=r)
    return num / den if den else np.nan


# --- pair definitions: (label, rdf file, candidate target files) ---
ALL_PAIRS = [
    ("PAA-PAA", "rdf_PAA_PAA.dat", ["PAA-PAA.dist.dat", "PAA-PAA.dist.new", "M-M.dist.new"]),
    ("PAA-K",   "rdf_PAA_K.dat",   ["PAA-K.dist.dat",   "PAA-K.dist.new",   "M-K.dist.new"]),
    ("PAA-Cl",  "rdf_PAA_CL.dat",  ["PAA-CL.dist.dat",  "PAA-CL.dist.new",  "M-CL.dist.new"]),
    ("K-Cl",    "rdf_K_CL.dat",    ["K-CL.dist.dat",    "K-CL.dist.new"]),
    ("K-K",     "rdf_K_K.dat",     ["K-K.dist.dat",     "K-K.dist.new"]),
    ("Cl-Cl",   "rdf_CL_CL.dat",   ["CL-CL.dist.dat",   "CL-CL.dist.new"]),
]

# keep only pairs this system actually computed
pairs = [p for p in ALL_PAIRS if (step_dir / p[1]).exists()]
if not pairs:
    sys.exit(f"ERROR: no rdf_*.dat files found in {step_dir}")

ncol = 3
nrow = int(np.ceil(len(pairs) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 4.3 * nrow), squeeze=False)
flat = axes.flatten()

print(f"{step}:  {len(pairs)} pair(s) found")
for ax, (label, rdf_name, tgt_names) in zip(flat, pairs):
    tgt_path = find_target(tgt_names)
    r_tgt = g_tgt = None
    if tgt_path is not None:
        r_tgt, g_tgt = read_target_rdf(tgt_path)
        ax.plot(r_tgt * 0.1, g_tgt, "k-", linewidth=2.0, label="Target (AA)")

    r_cg, g_cg = read_lammps_rdf(step_dir / rdf_name)
    ax.plot(r_cg * 0.1, g_cg, "r-", linewidth=1.8, alpha=0.85, label=f"CG ({step})")

    if r_tgt is not None:
        f = merit(r_cg, g_cg, r_tgt, g_tgt)
        col = "green" if f < 0.003 else "darkred"
        ax.text(0.97, 0.93, rf"$f_\mathrm{{merit}}$ = {f:.4f}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=12, color=col)
        print(f"   {label:8s} f_merit = {f:.5f}"
              f"{'   (converged)' if f < 0.003 else ''}")
    else:
        ax.text(0.97, 0.93, "no target", transform=ax.transAxes,
                ha="right", va="top", fontsize=11, color="gray")
        print(f"   {label:8s} (no target file)")

    ax.set_title(label, fontsize=15)
    ax.set_xlabel("r (nm)", fontsize=13)
    ax.set_xlim(0, 2.1)
    ax.tick_params(labelsize=12, length=6, width=1.5)
    ax.legend(fontsize=10, loc="center right")
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

for j in range(len(pairs), len(flat)):
    flat[j].axis("off")
for i in range(nrow):
    axes[i, 0].set_ylabel("g(r)", fontsize=13)

plt.tight_layout()
plt.show()
