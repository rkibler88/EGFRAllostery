import pandas as pd

# 1. Parse the RMSF data
def parse_rmsf(filename):
    res_rmsf = {}
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(("@", "#")):
                parts = line.split()
                # Store as residue_number: rmsf_value
                res_rmsf[int(float(parts[0]))] = float(parts[1])
    return res_rmsf

rmsf_map = parse_rmsf("rmsf.xvg")

# 2. Open your AlphaFold PDB and replace B-factors
with open("egfr_alphafold.pdb", "r") as f_in, open("egfr_rmsf_mapped.pdb", "w") as f_out:
    for line in f_in:
        if line.startswith("ATOM"):
            res_num = int(line[22:26].strip())
            # Get the RMSF for this residue, default to 0 if not found
            val = rmsf_map.get(res_num, 0.0)
            
            # PDB B-factor column is characters 61-66
            # We rewrite the line with the new value formatted to 6.2f
            new_line = line[:60] + f"{val:6.2f}" + line[66:]
            f_out.write(new_line)
        else:
            f_out.write(line)

print("Success! Created 'egfr_rmsf_mapped.pdb' with real MD data.")
