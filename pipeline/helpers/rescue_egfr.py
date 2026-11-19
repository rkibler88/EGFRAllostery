"""Append the aligned EGFR_HUMAN/712-968 row to the representative alignment (stage 01).

EGFR is not selected as a CD-HIT cluster representative, so its aligned row is
pulled from the full Stockholm alignment and appended: 32,911 representatives
+ 1 rescued EGFR = 32,912 sequences. The same row, written on its own, is the
column -> residue map used by every downstream script
(data/processed/egfr_reference_aligned.fasta).

Original: rescue_egfr.py, which only wrote the reference row and took no
arguments. Behaviour is otherwise unchanged.

Usage: python rescue_egfr.py <representatives.fasta> <out.fasta> [stockholm]
"""
import sys
from Bio import AlignIO

reps, out_file = sys.argv[1], sys.argv[2]
stockholm = sys.argv[3] if len(sys.argv) > 3 else "PF07714.alignment.full"
TARGET = "EGFR_HUMAN/712-968"

row = None
for alignment in AlignIO.parse(stockholm, "stockholm"):
    for record in alignment:
        if record.id == TARGET:
            row = str(record.seq)
            break
    if row is not None:
        break
if row is None:
    sys.exit(f"{TARGET} not found in {stockholm}")

with open("egfr_reference_aligned.fasta", "w") as fh:
    fh.write(f">{TARGET}\n{row}\n")

with open(out_file, "w") as out, open(reps) as src:
    out.write(src.read())
    out.write(f">{TARGET}\n{row}\n")

print(f"Rescued {TARGET} ({len(row)} columns) into {out_file} "
      f"and egfr_reference_aligned.fasta")
