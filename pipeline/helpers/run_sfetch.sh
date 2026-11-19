#!/bin/bash
#SBATCH --job-name=EGFR_sfetch
#SBATCH --output=sfetch_%j.out
#SBATCH --mem=16gb
#SBATCH --time=01:00:00

module load hmmer

esl-sfetch --index PF00069.alignment.full

esl-sfetch -f PF00069.alignment.full winner_ids.txt > kinase_final_aligned.fasta
