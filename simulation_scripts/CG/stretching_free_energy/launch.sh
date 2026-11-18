#!/usr/bin/env bash
#SBATCH --job-name=UmbRee
#SBATCH --time=09:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-265%96
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

set -euo pipefail

module load lammps/20241119
LMPEXE="${LMPEXE:-lmp}"
K=0.5

# --- map task id -> r0 (Å) exactly: 18.5..45.0 step 0.1 ---
tenth=$((185 + SLURM_ARRAY_TASK_ID))          # 18.5 Å = 185 tenths
r0=$(awk -v t="$tenth" 'BEGIN{printf "%.1f", t/10.0}')
win=$(printf "win_%05.1f" "$r0")

# --- required files (match your uploaded in.lmp) ---
req=(in.lmp paa_cg.data
     pair_PAA_PAA.table pair_PAA_NA.table pair_NA_NA.table
     bond_M_M.table angle_MMM.table dihedral_MMMM.table)

for f in "${req[@]}"; do
  [[ -f "$f" ]] || { echo "Missing $f in submit dir: $(pwd)"; exit 2; }
done

mkdir -p "$win"
cp -f in.lmp paa_cg.data "$win/"

# symlink tables into window dir (so relative filenames work)
for f in pair_PAA_PAA.table pair_PAA_NA.table pair_NA_NA.table bond_M_M.table angle_MMM.table dihedral_MMMM.table; do
  ln -sf "../$f" "$win/$f"
done

cd "$win"

echo "JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID r0=$r0 k=$K host=$(hostname)"
echo "Running: $LMPEXE -in in.lmp -var r0 $r0 -var k $K"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# Use srun (safe on most clusters even for 1 task)
srun -n 1 "$LMPEXE" -in in.lmp -var r0 "$r0" -var k "$K" \
  -log lammps.log -screen screen.out -echo both
