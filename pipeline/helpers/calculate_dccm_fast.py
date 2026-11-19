import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load your 253 evolutionary residues (The 'Truth' labels)
comm_df = pd.read_csv("egfr_communities_PRODUCTION_1.8L.csv")
evolution_residues = sorted(comm_df['residue'].unique())
print(f"Evolutionary analysis residues: {len(evolution_residues)}")

# 2. Load the GROMACS ASCII matrix
print("Loading covariance data (7.1MB)...")
# Flatten immediately to treat it as a 1D stream of numbers
raw_data = np.loadtxt("covar.ascii").flatten()
total_elements = raw_data.size
print(f"Total values loaded: {total_elements}")

# 3. Reshape based on total size
# D should be 777 (which is 259 residues * 3 dimensions)
D = int(np.sqrt(total_elements))
raw_covar = raw_data.reshape(D, D)
n_md_residues = D // 3
print(f"Matrix successfully reshaped to {D}x{D} ({n_md_residues} residues in MD)")

# 4. Define the MD Residue Range (712 to 970 = 259 residues)
md_residues = list(range(712, 712 + n_md_residues))

# 5. Vectorized Compression (The Stride Trick)
# This converts Cartesian (3N x 3N) to Residue (N x N)
dccm_full = raw_covar[0::3, 0::3] + raw_covar[1::3, 1::3] + raw_covar[2::3, 2::3]

# 6. SUBSETTING: Match the 259 physical residues to the 253 evolutionary residues
indices = [i for i, res in enumerate(md_residues) if res in evolution_residues]
dccm_subset = dccm_full[np.ix_(indices, indices)]
print(f"Subsetted DCCM to {dccm_subset.shape} to match evolutionary data.")

# 7. Normalize to Correlation Matrix [-1, 1]
diag = np.diag(dccm_subset).copy()
diag[diag == 0] = 1e-9 
norm_factor = np.sqrt(np.outer(diag, diag))
dccm_normalized = dccm_subset / norm_factor

# 8. Save and Plot
dccm_df = pd.DataFrame(dccm_normalized, index=evolution_residues, columns=evolution_residues)
dccm_df.to_csv("egfr_dccm_matrix.csv")

plt.figure(figsize=(12, 10))
sns.heatmap(dccm_df, cmap='RdBu_r', vmin=-1, vmax=1, center=0, square=True,
            cbar_kws={"label": "Correlation Coefficient"})

# Formatting for a scientific paper
plt.title("EGFR Kinase Domain: Dynamic Cross-Correlation Matrix (DCCM)")
plt.xlabel("Residue Number")
plt.ylabel("Residue Number")

plt.savefig("EGFR_DCCM_Publication.png", dpi=300)
print("Success! Final Analysis complete. Image: EGFR_DCCM_Publication.png")
