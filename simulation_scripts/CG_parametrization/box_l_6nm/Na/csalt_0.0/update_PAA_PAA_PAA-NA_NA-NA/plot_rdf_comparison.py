#!/usr/bin/env python3
"""
Compare CG RDFs from a given step with target (AA) RDFs.

Usage:
    python3 plot_rdf_comparison.py step_02
    python3 plot_rdf_comparison.py step_05
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: python3 plot_rdf_comparison.py <step_dir>")
    print("Example: python3 plot_rdf_comparison.py step_02")
    sys.exit(1)

step = sys.argv[1]
step_dir = Path(step)
targets_dir = Path("targets")

if not step_dir.is_dir():
    sys.exit(f"ERROR: {step_dir} not found")

# --- read functions ---

def read_lammps_rdf(path):
    """Read LAMMPS fix ave/time RDF output, return (r, g_avg) in Å.
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
    """Read target RDF (cols 0,1), return r in Å, g(r)."""
    arr = np.loadtxt(path, usecols=(0, 1))
    return arr[:, 0], arr[:, 1]

# --- pair definitions ---
pairs = [
    ("PAA–Na",  "rdf_PAA_NA.dat",  "PAA-NA.dist.dat"),
    ("PAA–PAA", "rdf_PAA_PAA.dat", "PAA-PAA.dist.dat"),
    ("Na–Na",   "rdf_NA_NA.dat",   "NA-NA.dist.dat"),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

for ax, (label, rdf_name, tgt_name) in zip(axes, pairs):
    rdf_path = step_dir / rdf_name
    tgt_path = targets_dir / tgt_name

    if tgt_path.exists():
        r_tgt, g_tgt = read_target_rdf(tgt_path)
        ax.plot(r_tgt * 0.1, g_tgt, "k-", linewidth=2.0, label="Target (AA)")

    if rdf_path.exists():
        r_cg, g_cg = read_lammps_rdf(rdf_path)
        ax.plot(r_cg * 0.1, g_cg, "r-", linewidth=1.8, alpha=0.85, label=f"CG ({step})")
    else:
        ax.text(0.5, 0.5, f"No {rdf_name}", transform=ax.transAxes,
                ha="center", fontsize=12, color="gray")

    ax.set_title(label, fontsize=16)
    ax.set_xlabel("r (nm)", fontsize=14)
    ax.set_xlim(0, 2.1)
    ax.tick_params(labelsize=13, length=6, width=1.5)
    ax.legend(fontsize=12)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

axes[0].set_ylabel("g(r)", fontsize=14)

plt.tight_layout()
outname = f"rdf_comparison_{step}.png"
plt.show()
print(f"Saved {outname}")
