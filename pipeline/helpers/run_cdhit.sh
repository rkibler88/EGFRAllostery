#!/bin/bash
#SBATCH --job-name=EGFR_CDHIT
#SBATCH --output=cdhit_%j.out
#SBATCH --error=cdhit_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rkibler@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=04:00:00

module load cdhit/4.6.8

# -i: input (your new fasta)
# -o: output (the filtered file)
# -c: 0.9 (90% identity threshold)
# -n: 5 (word size for 90%)
# -d: 0 (keep long headers)
# -T: 8 (matches cpus-per-task)
# -M: 64000 (matches memory in MB)

cd-hit -i kinase_ty_raw.fasta -o kinase_ty_62.fasta -c 0.62 -n 4 -M 32000 -T 8
