#!/bin/bash
#SBATCH --job-name=EGFR_PLMC
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#
# Stage 2 - direct coupling analysis.
#
# Output: egfr_couplings.txt  (31,878 pairs = C(253,2))
#
# plmc reports 32,672 valid sequences of 32,912 and 253 scored sites of 3,839.
# Effective sample size Neff = 27,612.8, i.e. Neff/L = 109 for L = 253.
# Because the alignment is already 62%-identity representatives, the -t 0.2 reweighting
# is inert: plmc reports a mean 80%-identity neighbourhood of 1.000 sequences.
#
# NOTE ON THE SCORE COLUMN: column 5 (Frobenius norm) is zero throughout this output.
# Column 6, the average-product-corrected score, is what downstream analysis uses.
# It is signed and approximately mean-zero, so it must not be summed over arbitrary
# residue-pair sets as though it were a magnitude (see docs/CORRECTIONS.md Section 5).
set -euo pipefail
module load conda && conda activate egfr_env

plmc -o egfr.plmc \
     -c egfr_couplings.txt \
     -f EGFR_HUMAN/712-968 \
     -t 0.2 -m 500 -le 1.0 -lh 0.01 -n 8 \
     kinase_final_aligned.fasta
