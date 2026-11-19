from Bio import AlignIO
import sys

# 1. Load your winner IDs into a 'set' (very fast for searching)
print("Loading winner IDs...")
winners = set()
with open("winner_ids.txt", "r") as f:
    for line in f:
        winners.add(line.strip())

print(f"Searching for {len(winners)} sequences...")

# 2. Open the input (Stockholm) and output (FASTA)
input_file = "PF07714.alignment.full"  # New input
output_file = "kinase_final_aligned.fasta" # Same output name is fine

# 3. Stream through the alignment
# We use AlignIO.parse to handle the Stockholm file one record at a time
# This saves memory!
count = 0
found = 0

with open(output_file, "w") as out_handle:
    # Stockholm files are parsed as a list of alignments; Pfam has 1 big one
    for alignment in AlignIO.parse(input_file, "stockholm"):
        for record in alignment:
            count += 1
            if record.id in winners:
                # This writes the sequence with gaps intact
                out_handle.write(f">{record.id} {record.description}\n")
                out_handle.write(f"{str(record.seq)}\n")
                found += 1
            
            # Print progress every 100,000 sequences
            if count % 100000 == 0:
                print(f"Processed {count} sequences... found {found}")

print(f"Done! Found {found} out of {len(winners)} winners.")
