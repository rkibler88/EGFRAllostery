import pandas as pd

# Load the production results
df = pd.read_csv("egfr_communities_PRODUCTION_1.8L.csv")

# 1. Get Top 5 Hubs (Weighted Degree)
top_hubs = df.sort_values(by='weighted_degree', ascending=False).head(5).index.tolist()

# 2. Get Top 5 Bridges (Betweenness)
top_bridges = df.sort_values(by='betweenness', ascending=False).head(5).index.tolist()

# 15 distinct colors for the communities
colors = ["red", "blue", "green", "yellow", "magenta", "cyan", 
          "orange", "pink", "slate", "purple", "wheat", "brightorange", 
          "limon", "forest", "deepblue"]

with open("viz_hubs_and_bridges.pml", "w") as f:
    f.write("load 1M17.pdb\n")
    f.write("hide all\n")
    f.write("show cartoon\n")
    f.write("set ray_opaque_background, on\n")
    f.write("bg_color white\n")
    f.write("color gray80, all\n") # Neutral background
    
    # 3. Color by Community
    for i, row in df.iterrows():
        res = int(row['residue'])
        comm = int(row['community'])
        color = colors[comm % len(colors)]
        f.write(f"color {color}, resi {res}\n")
    
    # 4. Show Hubs (The Anchors) as Black Spheres
    for h in top_hubs:
        f.write(f"show spheres, resi {h}\n")
        f.write(f"color black, resi {h}\n")
    
    # 5. Show Bridges (The Switches) as White Spheres
    for b in top_bridges:
        f.write(f"show spheres, resi {b}\n")
        f.write(f"color white, resi {b}\n")
    
    f.write("set sphere_scale, 1.0\n")
    # Add labels to see which is which in PyMOL
    f.write("label resi 793, 'Hinge Bridge (793)'\n")
    f.write("label resi 853, 'DFG Switch (853)'\n")
    f.write("label resi 965, 'Activator Hub (965)'\n")
    f.write("set label_size, 20\n")
    f.write("set label_color, black\n")
