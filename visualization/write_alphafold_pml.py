import pandas as pd

# 1. Load your Production Data
df = pd.read_csv("egfr_communities_PRODUCTION_1.8L.csv")

# 2. Get the specific Hubs and Bridges
top_hubs = df.sort_values(by='weighted_degree', ascending=False).head(5).index.tolist()
top_bridges = df.sort_values(by='betweenness', ascending=False).head(5).index.tolist()

colors = ["red", "blue", "green", "yellow", "magenta", "cyan", 
          "orange", "pink", "slate", "purple", "wheat", "brightorange", 
          "limon", "forest", "deepblue"]

# 3. Create the PML file
with open("viz_alphafold_FINAL.pml", "w") as f:
    f.write("reinitialize\n") # Clear everything
    f.write("load egfr_alphafold.pdb\n")
    f.write("hide all\n")
    
    # Show ONLY the Kinase Domain (where our data is)
    f.write("select kinase, resi 712-970\n")
    f.write("show cartoon, kinase\n")
    f.write("color gray90, all\n") # Make the rest of the protein faint gray
    
    # 4. Color by Community
    for i, row in df.iterrows():
        res = int(row['residue'])
        comm = int(row['community'])
        color = colors[comm % len(colors)]
        # We specify 'egfr_alphafold' to make sure it colors the right object
        f.write(f"color {color}, egfr_alphafold and resi {res}\n")
    
    # 5. Show Hubs (The Anchors) as BIG Black Spheres
    for h in top_hubs:
        f.write(f"show spheres, egfr_alphafold and resi {h}\n")
        f.write(f"color black, egfr_alphafold and resi {h}\n")
    
    # 6. Show Bridges (The Switches) as White Spheres
    for b in top_bridges:
        f.write(f"show spheres, egfr_alphafold and resi {b}\n")
        f.write(f"color white, egfr_alphafold and resi {b}\n")

    f.write("set sphere_scale, 1.2\n") # Make spheres larger to find them
    f.write("zoom kinase\n")
    f.write("bg_color white\n")
