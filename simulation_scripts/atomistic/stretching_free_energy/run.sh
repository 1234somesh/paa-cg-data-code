#!/bin/bash -l
#SBATCH --partition=gpu-single
#SBATCH --job-name=PULL    # Job name
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00           # Time limit (hh:mm:ss)
#SBATCH --mem=2gb
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

crdnt=$crdnt
wd=$wd

# Load GROMACS module
#source /etc/profile.d/lmod.sh
module load devel/cuda/12.1
module load compiler/gnu/12.1
module load chem/gromacs/2024.3-cuda-12.6

# Change to the correct working directory
echo "Running gromacs for coordinate $crdnt in $(pwd)"

unset OMP_NUM_THREADS

gmx grompp -f ../pro.mdp -c ../InitConfig/conf${crdnt}.gro -p ../topol.top -r ../InitConfig/conf${crdnt}.gro -o OUT$crdnt -n ../index.ndx -maxwarn 1
gmx mdrun -v -deffnm OUT$crdnt


