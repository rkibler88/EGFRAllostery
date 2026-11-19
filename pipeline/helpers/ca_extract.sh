#!/bin/bash
#SBATCH --job-name=ca_extract
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=00:30:00
#SBATCH --output=ca_extract_%j.out

module purge
module load cuda/12.8.1 gcc/14.2.0 openmpi/5.0.7 gromacs/2025.2
cd /blue/raquel.dias/rkibler/EGFR

echo C-alpha | gmx trjconv -s md_production.tpr -f md_clean.xtc -o md_ca.xtc
echo C-alpha | gmx trjconv -s md_production.tpr -f md_clean.xtc -o ca_ref.pdb -dump 0
gmx check -f md_clean.xtc
