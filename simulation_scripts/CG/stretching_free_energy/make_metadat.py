#!/usr/bin/env python3
import re
from pathlib import Path

K = 209.2                 # as you requested (units consistent with your WHAM run)
INFILE = "ree.dat"
OUT_TS = "ree_nm.dat"     # converted timeseries: time  Ree(nm)
OUT_META = "metadata.dat"

# win_019.7  -> 19.7 Å  -> 1.97 nm
pat = re.compile(r"^win_(\d+\.\d+)$")

num_re = re.compile(r'^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$')
def is_num(s: str) -> bool:
    return bool(num_re.match(s))

rows = []
skipped = []

for d in sorted(Path(".").glob("win_*")):
    if not d.is_dir():
        continue
    m = pat.match(d.name)
    if not m:
        continue

    r0_A = float(m.group(1))
    r0_nm = r0_A / 10.0

    inp = d / INFILE
    if not inp.exists():
        skipped.append((d.name, "missing ree.dat"))
        continue

    out = d / OUT_TS
    n_keep = 0

    with inp.open("r", errors="ignore") as fin, out.open("w") as fout:
        for line in fin:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("@"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            if not (is_num(parts[0]) and is_num(parts[1])):
                continue

            t = float(parts[0])
            ree_A = float(parts[1])
            ree_nm = ree_A / 10.0

            fout.write(f"{t:.6f} {ree_nm:.8f}\n")
            n_keep += 1

    if n_keep == 0:
        skipped.append((d.name, "no numeric data"))
        continue

    rows.append((r0_nm, out))

# sort by R0
rows.sort(key=lambda x: x[0])

with open(OUT_META, "w") as f:
    for r0_nm, ts in rows:
        # Grossfield/Bauer-style 1D metadata: file  R0  k
        f.write(f"{ts.as_posix()} {r0_nm:.6f} {K:.6f}\n")

print(f"Wrote {OUT_META} with {len(rows)} windows.")
print(f"Wrote converted timeseries as win_*/{OUT_TS} (nm).")
if skipped:
    print(f"Skipped {len(skipped)} windows (first 10):")
    for w, why in skipped[:10]:
        print(f"  {w}: {why}")

