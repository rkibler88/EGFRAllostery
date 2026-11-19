#!/bin/bash
#SBATCH --job-name=EGFR_100ns
#SBATCH --output=md_%j.out
#SBATCH --error=md_%j.err
#SBATCH --partition=hpg-turin
#SBATCH --account=raquel.dias     
#SBATCH --gres=gpu:1               
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30gb               
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rkibler@ufl.edu

# Load the modules
module load cuda/12.8.1 gcc/14.2.0 openmpi/5.0.7 gromacs/2025.2

# 1. Pre-process
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_production.tpr

# 2. Run with GPU offloading
gmx mdrun -v -deffnm md_production \
    -nb gpu \
    -pme gpu \
    -bonded gpu \
    -update gpu \
    -cpt 15
