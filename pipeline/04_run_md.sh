#!/bin/bash
# Stage 4 - molecular dynamics and dynamic cross-correlation.
#
# GROMACS 2025.2, AMBER99SB-ILDN, TIP3P, PME, Parrinello-Rahman.
# Starting structure: AlphaFold2 model of P00533, kinase domain, active state.
#
# Output: egfr_dccm_matrix.csv (253 x 253, normalised to unit diagonal)
#         rmsf.xvg
#
# SCOPE, which bounds every dynamic claim in the manuscript:
#   - ONE replicate (n = 1), 100 ns, 50,001 frames
#   - monomer; no dimer, no mutant, no inactive-state simulation
#   - construct spans residues 712-970 (259 C-alpha), so residue 965 lies five
#     residues from a free terminus and the C-terminal tail 971-1210 is absent
#   - no conformational transition was sampled
set -euo pipefail
module load gromacs

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md_production

gmx trjconv -s md.tpr -f md_production.xtc -o md_clean.xtc -pbc mol -center
gmx covar  -s md.tpr -f md_clean.xtc -o eigenval.xvg -ascii covar.ascii -nomwa
gmx rmsf   -s md.tpr -f md_clean.xtc -o rmsf.xvg -res

python calculate_dccm_fast.py   # 3Nx3N covariance -> NxN correlation, subset to the 253 scored residues
