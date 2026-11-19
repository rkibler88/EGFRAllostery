from pymol import cmd

def build_fig6():
    # 1. SETUP THE SCENE
    cmd.reinitialize()
    # Using 1M17 (Active) to show the 'Handshake' in its functional state
    cmd.load("1M17.pdb", "egfr")
    # Apply the verified +24 offset for 1M17
    cmd.alter("egfr", "resi = str(int(resi) + 24)")
    cmd.sort()
    
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ambient", 0.5)
    
    # 2. THE TARGET: COMMUNITY 13 (DFG SWITCH)
    # This is the 'Allosteric Integration Hub'
    cmd.select("comm13", "resi 853-855 or resi 793-795 or resi 816-820")
    cmd.show("surface", "comm13")
    cmd.color("brightorange", "comm13")
    cmd.set("transparency", 0.3, "comm13")
    
    # Close-up of the DFG Motif atoms
    cmd.show("sticks", "resi 853-855")
    cmd.set("stick_radius", 0.3)
    
    # 3. THE INFRASTRUCTURE: PROTEIN GHOST
    cmd.show("cartoon", "resi 712-968")
    cmd.color("gray90", "egfr")
    cmd.set("cartoon_transparency", 0.7)
    
    # 4. PATHWAY 1: THE ONCOGENIC SHORT-CIRCUIT (L858R -> DFG)
    # Data from cluster: 858 -> 857 -> 767 -> 763 -> 853
    path_mut = [858, 857, 767, 763, 853]
    for i in range(len(path_mut)-1):
        name = f"path_mut_{i}"
        cmd.distance(name, f"resi {path_mut[i]} and n. CA", f"resi {path_mut[i+1]} and n. CA")
        cmd.set("dash_color", "red", name)
        cmd.set("dash_width", 5, name)
        cmd.set("dash_gap", 0, name) # Continuous tube look
    
    # Highlight the mutation site
    cmd.show("spheres", "resi 858")
    cmd.color("black", "resi 858")
    
    # 5. PATHWAY 2: THE NATURAL CLUTCH (P-loop -> alphaC -> DFG)
    # Representative 'Validation' highway
    path_nat = [722, 755, 762, 853]
    for i in range(len(path_nat)-1):
        name = f"path_nat_{i}"
        cmd.distance(name, f"resi {path_nat[i]} and n. CA", f"resi {path_nat[i+1]} and n. CA")
        cmd.set("dash_color", "marine", name)
        cmd.set("dash_width", 5, name)
        cmd.set("dash_gap", 0.2, name) # Dotted line look
        
    # 6. LABELS AND FOCUS
    cmd.hide("labels")
    cmd.center("resi 853")
    cmd.zoom("resi 712-968", 5)
    
    # Orientation: Looking into the pocket
    cmd.set_view((
     0.36,   -0.54,    0.75,
     0.88,    0.46,   -0.09,
    -0.29,    0.70,    0.64,
     0.00,    0.00, -150.00,
    12.00,   25.00,   35.00,
    110.00,  190.00,    0.00 ))

    print("Figure 6 Scene Prepared. Use 'ray' to render.")

build_fig6()
