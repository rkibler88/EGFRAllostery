# Upstream pipeline

These stages produced the files in `data/processed/`. They are **not runnable from a clone**:
they need inputs of 2.5 GB (the Pfam alignment) to 8.2 GB (the MD trajectory), which are not
committed. They are recorded here so the provenance is auditable and the work can be repeated
from the original sources. `analysis/` consumes their outputs and needs none of them.

| stage | what it does | key output |
|---|---|---|
| `01_build_alignment.sh` | PF07714 full alignment -> CD-HIT 62% representatives -> EGFR row rescued | `kinase_final_aligned.fasta` (32,912 x 3,839) |
| `02_run_plmc.sh` | pseudo-likelihood maximisation DCA | `egfr_couplings.txt` (31,878 pairs) |
| `03_build_network.py` | 462-edge network, weighted Louvain, 100 seeds | `egfr_communities_PRODUCTION_1.8L.csv` |
| `04_run_md.sh` | 100 ns GROMACS production run and covariance analysis | `egfr_dccm_matrix.csv`, `rmsf.xvg` |
| `05_md_replicates/` | four further 100 ns replicates and their agreement statistics | `replicate_*.csv`, `egfr_dccm_matrix_rep{2..5}.csv` |
| `06_sca_sectors.md` | statistical coupling analysis and sector decomposition (pySCA) | `sca_sectors.csv`, `pair_table.csv.gz` |

`helpers/` holds the author's scripts that the stage drivers call, and `mdp/` the GROMACS
parameter files for minimisation, ion placement, NVT, NPT and production.

## Corrections applied to the original scripts

| stage | original | here |
|---|---|---|
| 01 | `run_sfetch.sh` indexed **PF00069** (protein kinase domain), not the family actually used | PF07714 |
| 01 | `run_cdhit.sh` header comments described a 90% / `-n 5` run; the command executed was `-c 0.62 -n 4` | comments match the command |
| 01 | `convert_safe.py` and `rescue_egfr.py` hardcoded their filenames | both take arguments |
| 03 | `L_ref = 257` (focus columns requested) while the network has 253 nodes (match-state positions scored) | both named explicitly; edge count kept at 462 for continuity with the deposited partition |
| 04 | DCCM was first computed with an O(N^2) Python loop | the vectorised `calculate_dccm_fast.py`, which is what produced the deposited matrix |

Full disposition of every script in the original working directories:
`docs/SCRIPT_PROVENANCE.md`.

## What the manuscript used

A **single** CD-HIT pass at 62%. An 80% -> 60% branch also exists in the original working
directory, terminating in `kinase_60_winners.fasta` (93,773 sequences); it was never aligned
and did not feed the published couplings. Earlier drafts described that branch in Methods.
