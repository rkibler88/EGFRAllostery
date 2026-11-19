import pandas as pd
from Bio import SeqIO

# 1. Load the Map (Seat # -> Residue #)
record = SeqIO.read("egfr_reference_aligned.fasta", "fasta")
mapping = {}
current_res = 712
for i, char in enumerate(str(record.seq)):
    if char.isalpha():
        mapping[i + 1] = current_res
        current_res += 1

# 2. Load the Raw PLMC Data (all 31,878 rows)
# We use Column 5 for the score
df = pd.read_csv("egfr_couplings.txt", sep=" ", header=None, 
                 usecols=[0, 2, 5], names=["col_i", "col_j", "score"])

# 3. Translate and Save the FULL list
df['res_i'] = df['col_i'].map(mapping)
df['res_j'] = df['col_j'].map(mapping)
df = df.dropna(subset=['res_i', 'res_j'])
df['res_i'] = df['res_i'].astype(int)
df['res_j'] = df['res_j'].astype(int)

df.to_csv("egfr_couplings_translated.csv", index=False)
print("Master List Created: egfr_couplings_translated.csv")
