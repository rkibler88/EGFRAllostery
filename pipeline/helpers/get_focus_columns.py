from Bio import SeqIO

# Load the reference sequence we just rescued
record = SeqIO.read("egfr_reference_aligned.fasta", "fasta")
sequence = str(record.seq)

# Find positions (1-indexed) that are letters (A-Z), not dashes (-)
# These represent the 257 residues of the EGFR kinase domain
focus_cols = [i + 1 for i, char in enumerate(sequence) if char.isalpha()]

print(f"EGFR has {len(focus_cols)} amino acids in this alignment.")

# Save as a comma-separated list for the plmc tool
with open("focus_columns.txt", "w") as f:
    f.write(",".join(map(str, focus_cols)))

print("focus_columns.txt has been created.")
