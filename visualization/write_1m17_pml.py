import pandas as pd

# 1. Load your Production Data (Uses UniProt 712-968 scale)
df = pd.read_csv("egfr_communities_PRODUCTION_1.8L.csv")

# 2. Define the Hubs and Bridges for highlighting
top_hubs = [965, 717, 799, 851, 964]    # Weighted Degree
top_bridges = [792, 838, 793, 853, 852] # Betweenness

# 15 High-contrast colors
colors = ["red", "blue", "green", "yellow", "magenta", "cyan", 
          "orange", "pink", "slate", "purple", "wheat", "brightorange", 
          "limon", "forest", "deepblue"]

with open("viz_final_comparison.pml", "w") as f:
    f.write("reinitialize\n")
    f.write("load 1M17.pdb, active_state\n")
    f.write("load egfr_alphafold.pdb, inactive_state\n")
    
    # 3. FIX NUMBERING FOR 1M17 (+24 offset)
    f.write("alter active_state, resi = str(int(resi) + 24)\n")
    f.write("sort\n")
    
    # 4. ALIGN AND ZOOM
    # Align 1M17 to the kinase domain of AlphaFold
    f.write("cealign inactive_state and resi 712-968, active_state\n")
    f.write("hide everything\n")
    f.write("show cartoon, resi 712-968\n")
    f.write("zoom resi 712-968\n")
    
    # 5. DUAL COLORING
    f.write("color gray80, all\n") # Background
    
    for i, row in df.iterrows():
        res = int(row['residue'])
        comm = int(row['community'])
        color = colors[comm % len(colors)]
        # Color both structures using the same UniProt residue number
        f.write(f"color {color}, resi {res}\n")
    
    # 6. HIGHLIGHT THE "CLUTCH" (alphaC-helix) to see the 25-degree shift
    # We will make this stand out
    f.write("set cartoon_transparency, 0.4, resi 752-768\n")
    
    # 7. SHOW HUBS (Black) and BRIDGES (White)
    for h in top_hubs:
        f.write(f"show spheres, resi {h}\n")
        f.write(f"color black, resi {h}\n")
    for b in top_bridges:
        f.write(f"show spheres, resi {b}\n")
        f.write(f"color white, resi {b}\n")

    f.write("set sphere_scale, 0.8\n")
    f.write("bg_color white\n")
    f.write("set ray_opaque_background, on\n")
    
    # Labeling landmark residues for the paper
    f.write("label (active_state and resi 853 and n. CA), 'DFG-In (Active)'\n")
    f.write("label (inactive_state and resi 853 and n. CA), 'DFG-Out (Inactive)'\n")
    f.write("set label_size, 20\n")
    f.write("set label_color, black\n")
