#!/usr/bin/env python3
"""
Iterative Boltzmann inversion update for LAMMPS tabulated nonbonded potentials.

Usage
-----
From inside a step directory (e.g. step_01) that contains RDFs:

    rdf_PAA_PAA.dat
    rdf_PAA_NA.dat
    rdf_NA_NA.dat

and with ibi.py, initial/, and targets/ one level above:

    python3 ../ibi.py [T_K]

Arguments
---------
T_K : optional temperature in K (default 300.0)

Behaviour
---------
The script will:
  * read RDFs from the current directory
  * compare to targets in ../targets/
  * choose base tables from ../initial/ or previous step's updated_potentials
  * update **PAA-PAA and PAA-NA in every step**
  * never update NA-NA (always copied)
  * write tables (updated or copied) to ./updated_potentials/pair_*.table
"""

from pathlib import Path
import sys
import numpy as np

# --- constants ---
kB_kcal_per_mol_K = 0.0019872041  # kcal/mol/K


# ---------- Utilities ----------

def read_lammps_rdf(path):
    """
    Read a LAMMPS 'fix ave/time ... compute rdf' output in mode vector
    and return (r, g_avg) where:
      - r is the distance grid (Å)
      - g_avg is the time-averaged RDF

    Assumes format like:
        Time-averaged data for fix ...
        # TimeStep Number-of-rows
        # Row ...
        200 200
        1  r  g  ...
        ...
        400 200
        ...

    We use column 2 (index 1) as r, column 3 (index 2) as g(r).
    """
    path = Path(path)
    rs = None
    g_blocks = []

    with path.open() as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if (not line) or line.startswith("#") or line.startswith("Time-averaged"):
                continue

            parts = line.split()
            # header of a block: "<timestep> <nrows>"
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                nrows = int(parts[1])
                r_block = []
                g_block = []
                for _ in range(nrows):
                    row = f.readline()
                    if not row:
                        break
                    cols = row.split()
                    if len(cols) < 3:
                        continue
                    r_block.append(float(cols[1]))
                    g_block.append(float(cols[2]))
                r_block = np.asarray(r_block, dtype=float)
                g_block = np.asarray(g_block, dtype=float)

                if rs is None:
                    rs = r_block
                else:
                    if not np.allclose(rs, r_block, rtol=1e-6, atol=1e-6):
                        raise ValueError(f"Inconsistent r-grid in RDF file {path}")

                g_blocks.append(g_block)

    if rs is None or len(g_blocks) == 0:
        raise RuntimeError(f"No RDF blocks found in {path}")

    g_avg = np.mean(np.stack(g_blocks, axis=0), axis=0)
    return rs, g_avg


def read_target_rdf(path):
    """
    Read target RDF: assume 2+ columns, r and g(r) in the first two.
    r is assumed to be in nm if max(r) < 5, else Å.
    We return r in Å, g as-is.
    """
    arr = np.loadtxt(path, usecols=(0, 1))
    if arr.ndim == 1:
        arr = arr.reshape(-1, len(arr))
    r = arr[:, 0].astype(float)
    g = arr[:, 1].astype(float)

    # crude unit detection: nm -> Å
    if r.max() < 5.0:
        r = r * 10.0  # nm -> Å

    return r, g


def read_lammps_table(path):
    """
    Read a LAMMPS pair_style table file and split into header + numeric data.

    Returns:
        header_lines : list of strings (kept as-is)
        data         : ndarray of shape (N, 4): [idx, r, U, F]
    """
    path = Path(path)
    with path.open() as f:
        lines = f.readlines()

    header_lines = []
    data_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped == "":
            header_lines.append(line)
            continue
        # data lines start with a digit (index)
        if stripped[0].isdigit() or (stripped[0] == '-' and len(stripped) > 1 and stripped[1].isdigit()):
            data_lines.append(line)
        else:
            header_lines.append(line)

    if not data_lines:
        raise RuntimeError(f"No numeric data lines found in table {path}")

    data = np.loadtxt(data_lines)
    return header_lines, data


def write_lammps_table(path, header_lines, data):
    """
    Write a LAMMPS pair table back out with the original header and
    updated numeric data (idx, r, U, F).
    """
    path = Path(path)
    with path.open("w") as f:
        for line in header_lines:
            f.write(line)
        np.savetxt(f, data, fmt="%8d %10.6f %15.8e %15.8e")


def compute_force_from_potential(r, U):
    """
    Compute F(r) = -dU/dr via central finite differences.
    r, U are 1D arrays of same length.
    """
    r = np.asarray(r, dtype=float)
    U = np.asarray(U, dtype=float)

    dUdr = np.zeros_like(U)
    # central differences
    dUdr[1:-1] = (U[2:] - U[:-2]) / (r[2:] - r[:-2])
    # forward/backward for edges
    dUdr[0] = (U[1] - U[0]) / (r[1] - r[0])
    dUdr[-1] = (U[-1] - U[-2]) / (r[-1] - r[-2])

    F = -dUdr
    return F


def ibi_update_on_table(r_table, U_table, r_rdf, g_imp, r_ref, g_ref,
                        T_K=300.0, r_cut_corr=20.0, g_floor=1e-6, prefactor=0.04):
    """
    Apply one IBI update step to U_table(r) on the same r grid as the table.

    U_new(r) = U_table(r) + prefactor * kT * ln(g_imp/g_ref) for r <= r_cut_corr
    U_new(r) = U_table(r) for r >  r_cut_corr

    r_table : array (Å) – r-grid of the LAMMPS table
    U_table : array (kcal/mol) – short-range potential
    r_rdf   : array (Å) – r-grid from LAMMPS RDF
    g_imp   : array – RDF from current CG simulation
    r_ref   : array (Å) – r-grid of target RDF
    g_ref   : array – target RDF
    """
    r_table = np.asarray(r_table, dtype=float)
    U_table = np.asarray(U_table, dtype=float)

    # interpolate RDFs onto the table grid
    g_imp_interp = np.interp(r_table, r_rdf, g_imp,
                             left=g_imp[0], right=g_imp[-1])
    g_ref_interp = np.interp(r_table, r_ref, g_ref,
                             left=g_ref[0], right=g_ref[-1])

    # avoid log(0)
    g_imp_clip = np.clip(g_imp_interp, g_floor, None)
    g_ref_clip = np.clip(g_ref_interp, g_floor, None)

    kT = kB_kcal_per_mol_K * T_K
    dU = prefactor * kT * np.log(g_imp_clip / g_ref_clip)

    # only apply correction up to r_cut_corr
    mask_corr = (r_table <= r_cut_corr)
    dU_eff = np.zeros_like(U_table)
    dU_eff[mask_corr] = dU[mask_corr]

    U_new = U_table + dU_eff
    return U_new


# ---------- Main driver ----------
def main():
    # Run from inside step_* directory:
    #   python3 ../ibi.py [T_K]
    #
    # Current working dir = step_X
    step_dir = Path.cwd()
    parent = step_dir.parent  # where step_* live

    # root directory where ibi.py resides (and where initial/ and targets/ live)
    root = Path(__file__).resolve().parent

    # optional temperature argument
    T_K = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0

    # --- determine where to read base potentials from ---
    # If we are in step_01 -> base ../initial
    # If we are in step_0N (N>1) -> base ../step_0(N-1)/updated_potentials if it exists
    import re
    step_name = step_dir.name
    base_pot_dir = None

    m = re.match(r"step_(\d+)", step_name)
    if m:
        idx = int(m.group(1))
        if idx > 1:
            prev_name = f"step_{idx-1:02d}"
            prev_updated = parent / prev_name / "updated_potentials"
            if prev_updated.is_dir():
                base_pot_dir = prev_updated
            else:
                # fallback if previous updated_potentials is missing
                base_pot_dir = parent / "initial"
        else:
            # step_01
            idx = 1
            base_pot_dir = parent / "initial"
    else:
        # not in a step_* folder: just default to ../initial, treat as idx=1
        idx = 1
        base_pot_dir = parent / "initial"

    # --- update policy per step ---
    # Now: always update BOTH PAA_PAA and PAA_NA in each step.
    # NA_NA is never updated – only copied.
    pairs_to_update = {"PAA_PAA", "PAA_NA"}

    initial_dir = base_pot_dir.resolve()
    targets_dir = (root / "targets").resolve()
    updated_dir = (step_dir / "updated_potentials").resolve()
    updated_dir.mkdir(parents=True, exist_ok=True)

    print(f"Step directory          : {step_dir}")
    print(f"Step index              : {idx}")
    print(f"Base potentials (input) : {initial_dir}")
    print(f"Target RDFs             : {targets_dir}")
    print(f"Output updated pots     : {updated_dir}")
    print(f"Temperature             : {T_K} K")
    print(f"Pairs updated this step : {sorted(pairs_to_update)}")
    print("Corrections applied only up to r_cut_corr for each pair; tails unchanged.\n")

    # pair definitions: (label, rdf_file, target_rdf, table_name)
    pairs = [
        ("PAA_PAA", "rdf_PAA_PAA.dat", "PAA-PAA.dist.dat", "pair_PAA_PAA.table"),
        ("PAA_NA",  "rdf_PAA_NA.dat",  "PAA-NA.dist.dat",  "pair_PAA_NA.table"),
        ("NA_NA",   "rdf_NA_NA.dat",   "NA-NA.dist.dat",   "pair_NA_NA.table"),
    ]

    for label, rdf_name, tgt_name, table_name in pairs:
        print(f"--- Handling {label} ---")
        rdf_path = step_dir / rdf_name
        tgt_path = targets_dir / tgt_name
        table_path = initial_dir / table_name
        out_table_path = updated_dir / table_name

        # Always need a base table to either update or copy
        if not table_path.exists():
            print(f"  WARNING: base table {table_path} not found. Skipping {label}.")
            continue

        # Read base table
        header_lines, data = read_lammps_table(table_path)
        idx_col = data[:, 0].astype(int)
        r_table = data[:, 1]
        U_table = data[:, 2]

        # Decide whether to actually update this pair
        if label not in pairs_to_update:
            # Just copy the table unchanged
            print(f"  Not updating {label}; copying base table to {out_table_path}")
            write_lammps_table(out_table_path, header_lines, data)
            continue

        # For labels we DO want to update (PAA_PAA, PAA_NA):
        if not rdf_path.exists():
            print(f"  WARNING: RDF file {rdf_path} not found. Copying base table for {label}.")
            write_lammps_table(out_table_path, header_lines, data)
            continue
        if not tgt_path.exists():
            print(f"  WARNING: target RDF {tgt_path} not found. Copying base table for {label}.")
            write_lammps_table(out_table_path, header_lines, data)
            continue

        print(f"  Updating {label} via IBI")

        # 1) Read RDF from LAMMPS
        r_rdf, g_imp = read_lammps_rdf(rdf_path)

        # 2) Read target RDF
        r_ref, g_ref = read_target_rdf(tgt_path)

        # 3) Do IBI update (pair-specific correction cutoff)
        if label == "PAA_NA":
            r_cut_corr = 20.00  # Å
        elif label == "PAA_PAA":
            r_cut_corr = 20.0  # Å (adjust if you want 12 Å here as well)
        else:
            r_cut_corr = 20.0  # fallback, not used for NA_NA currently

        U_new = ibi_update_on_table(
            r_table, U_table,
            r_rdf, g_imp,
            r_ref, g_ref,
            T_K=T_K,
            r_cut_corr=r_cut_corr,
            g_floor=1e-6,
        )

        # 4) Recompute force to be consistent with updated potential
        F_new = compute_force_from_potential(r_table, U_new)

        new_data = np.column_stack([idx_col, r_table, U_new, F_new])

        write_lammps_table(out_table_path, header_lines, new_data)
        print(f"  -> wrote updated table to {out_table_path}")

    print("\nDone. For the next iteration, point LAMMPS to the tables in:")
    print(f"   {updated_dir}")


if __name__ == "__main__":
    main()

