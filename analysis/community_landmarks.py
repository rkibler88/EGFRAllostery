import pandas as pd

# 1. Load results
df = pd.read_csv("egfr_communities_PRODUCTION_1.8L.csv")

# 2. Define Landmarks (UniProt numbers)
landmarks = {
    "P-loop (Phosphate Binding)": range(718, 725),
    "N-lobe Beta Core": range(712, 752),
    "alphaC-helix (The Clutch)": range(752, 769),
    "Hinge Region (ATP Binding)": range(790, 796),
    "alphaD-helix (Structural Link)": range(797, 805),
    "alphaE-helix (C-lobe Core)": range(811, 827),
    "Catalytic Loop": range(830, 840),
    "Activation Loop (The Switch)": range(840, 866),
    "alphaF-helix (Engine Block/Scaffold)": range(875, 896),
    "alphaG-helix (Dimerization Interface)": range(905, 921),
    "alphaH-alphaI (Tail Base)": range(930, 969)
}

print(f"{'Comm':<5} | {'Size':<5} | {'Dominant Structural Landmarks'}")
print("-" * 70)

for comm_id in sorted(df['community'].unique()):
    # FIX: Get the actual values from the 'residue' column
    res_list = df[df['community'] == comm_id]['residue'].tolist()
    
    matches = []
    for name, r_range in landmarks.items():
        # Check how many residues of this community fall into this landmark range
        overlap = set(res_list).intersection(set(r_range))
        if len(overlap) > 1: # If 2 or more residues match
            matches.append(f"{name} ({len(overlap)} res)")
    
    landmark_str = ", ".join(matches) if matches else "Distal/Inter-loop regions"
    print(f"{comm_id:<5} | {len(res_list):<5} | {landmark_str}")
