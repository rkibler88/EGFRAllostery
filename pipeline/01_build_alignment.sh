#!/bin/bash
# Stage 1 - redundancy-filtered alignment from Pfam PF07714.
#
# Input : PF07714.alignment.full   (Stockholm, 400,321 sequences x 3,839 columns)
# Output: kinase_final_aligned.fasta (32,912 x 3,839)
#
# One CD-HIT pass at 62% identity. EGFR_HUMAN/712-968 is not selected as a cluster
# representative, so its aligned row is recovered from the Stockholm file and appended:
# 32,911 representatives + 1 rescued EGFR = 32,912.
set -euo pipefail
module load cdhit hmmer

python convert_safe.py PF07714.alignment.full kinase_ty_raw.fasta   # Stockholm -> FASTA, gaps preserved
cd-hit -i kinase_ty_raw.fasta -o kinase_ty_62.fasta -c 0.62 -n 4 -M 32000 -T 8
grep "^>" kinase_ty_62.fasta | sed "s/^>//" | cut -d" " -f1 > winner_ids.txt

esl-sfetch --index PF07714.alignment.full
esl-sfetch -f PF07714.alignment.full winner_ids.txt > kinase_aligned.fasta
python rescue_egfr.py kinase_aligned.fasta kinase_final_aligned.fasta   # append EGFR_HUMAN/712-968
