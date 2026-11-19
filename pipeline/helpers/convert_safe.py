"""Stockholm -> FASTA with gap columns preserved (stage 01).

Original: convert_safe.py in the working directory, which hardcoded
PF00069.alignment.full as input. The alignment actually used for the
published couplings is PF07714 (see docs/SCRIPT_PROVENANCE.md), so the
family is now a command-line argument rather than a constant.

Usage: python convert_safe.py <alignment.stockholm> <out.fasta>
"""
import sys
from Bio import AlignIO

in_file, out_file = sys.argv[1], sys.argv[2]
print(f"Converting {in_file} -> {out_file} ...")
n = AlignIO.convert(in_file, "stockholm", out_file, "fasta")
print(f"Converted {n} alignment(s); gap columns preserved.")
