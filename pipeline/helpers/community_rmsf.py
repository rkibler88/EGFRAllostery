import pandas as pd
import numpy as np

# 1. Load your communities (The 'Teams')
comm_df = pd.read_csv("egfr_communities_PRODUCTION_1.8L.csv")

# 2. Load the RMSF data (The 'Vibrations')
# We skip the header lines starting with @ or #
rmsf_data = []
with open("rmsf.xvg", "r") as f:
    for line in f:
        if not (line.startswith("@") or line.startswith("#")):
            parts = line.split()
            rmsf_data.append([int(float(parts[0])), float(parts[1])])

rmsf_df = pd.DataFrame(rmsf_data, columns=['residue', 'rmsf'])

# 3. Merge the two worlds
merged = pd.merge(comm_df, rmsf_df, on='residue')

# 4. Calculate Average Fluctuation per Community
# This is the most objective way to see which modules are 'engines' vs 'scaffolds'
summary = merged.groupby('community').agg({
    'rmsf': ['mean', 'max', 'count']
}).reset_index()

# Flatten columns for easier reading
summary.columns = ['community', 'mean_rmsf', 'max_rmsf', 'residue_count']

# 5. Define our Community Names again for the report
comm_names = {
    0:  "N-lobe Beta Framework",
    1:  "P-loop/Phosphate Binding",
    2:  "Hinge Region (ATP Binding)",
    3:  "C-lobe Core/Tail Base",
    4:  "alphaC-helix (Distal Segment)",
    5:  "Asymmetric Dimer Interface",
    6:  "Central Linker Hub",
    7:  "Regulatory alphaC-Helix (Clutch)",
    8:  "Scaffold-Switch Interface",
    9:  "C-lobe Core (alphaE)",
    10: "Tail Base Anchor A",
    11: "Tail Base Anchor B",
    12: "Activation Loop (Distal)",
    13: "Catalytic Switch (DFG)",
    14: "Main Engine Block (Scaffold)"
}

summary['name'] = summary['community'].map(comm_names)

# 6. Print the Objective Results
print("\n" + "="*80)
print(f"{'COMMUNITY NAME':<40} | {'MEAN RMSF (nm)':<15} | {'SIZE'}")
print("="*80)
sorted_summary = summary.sort_values(by='mean_rmsf', ascending=False)
for _, row in sorted_summary.iterrows():
    print(f"{row['name']:<40} | {row['mean_rmsf']:<15.4f} | {int(row['residue_count'])}")

# Save for your paper's supplemental data
summary.to_csv("community_dynamics_summary.csv", index=False)
