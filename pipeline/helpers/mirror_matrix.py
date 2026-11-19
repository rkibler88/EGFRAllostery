# Save as export_mirror_matrix.py
import pandas as pd
import numpy as np

# Load MD matrix
dccm = pd.read_csv("egfr_dccm_matrix.csv", index_col=0)
res_list = dccm.index.astype(int).tolist()

# Create a blank DCA matrix for the same residues
dca_matrix = np.zeros((len(res_list), len(res_list)))
couplings = pd.read_csv("egfr_couplings_translated.csv")

# Fill matrix
res_to_idx = {res: i for i, res in enumerate(res_list)}
for _, row in couplings.iterrows():
    ri, rj = int(row['res_i']), int(row['res_j'])
    if ri in res_to_idx and rj in res_to_idx:
        dca_matrix[res_to_idx[ri], res_to_idx[rj]] = row['score']
        dca_matrix[res_to_idx[rj], res_to_idx[ri]] = row['score']

np.save("fig7_dca_aligned.npy", dca_matrix)
print("Aligned DCA matrix saved for Figure 7.")
