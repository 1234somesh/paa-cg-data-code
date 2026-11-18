#!/usr/bin/env python3
"""Convert .dist.new files (r in nm) to .dist.dat files (r in Angstrom)."""

import os

conversions = [
    ("Cs-CL.dist.new",  "CS-CL.dist.dat"),
    ("CL-CL.dist.new",  "CL-CL.dist.dat"),
    ("Cs-Cs.dist.new",  "CS-CS.dist.dat"),
    ("M-Cs.dist.new",   "PAA-CS.dist.dat"),
    ("M-M.dist.new",    "PAA-PAA.dist.dat"),
    ("M-CL.dist.new",    "PAA-CL.dist.dat")
]

base = os.path.dirname(os.path.abspath(__file__))

for src, dst in conversions:
    src_path = os.path.join(base, src)
    dst_path = os.path.join(base, dst)
    with open(src_path) as f_in, open(dst_path, 'w') as f_out:
        for line in f_in:
            parts = line.split()
            r_ang = float(parts[0]) * 1.0
            f_out.write(f"{r_ang:.4f} {parts[1]} {parts[2]}\n")
    print(f"{src} -> {dst}")

print("Done.")
