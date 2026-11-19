# Data

## Committed — `processed/`

Everything `analysis/` needs. SHA-256 for each file is in `processed/SHA256SUMS`
(`cd data/processed && shasum -c SHA256SUMS`).

| file | size | what it is |
|---|---|---|
| `community_dynamics_summary.csv` | 1 KB | mean, max RMSF and size per community (replicate 1) |
| `contact_qc.json` | 8 KB | the 328 EGFR kinase-domain chains used for the conformational-ensemble contact test, and the chains rejected |
| `egfr_clean.pdb` | 167 KB | AlphaFold2 model of P00533, kinase domain — source of Ca-Ca distances |
| `egfr_communities_PRODUCTION_1.8L.csv` | 9 KB | Louvain partition, weighted degree and betweenness at 462 edges |
| `egfr_couplings.txt` | 799 KB | raw plmc output, 31,878 pairs; residue identities in columns 2 and 4; APC score in column 6 |
| `egfr_couplings_translated.csv` | 952 KB | the same couplings mapped to UniProt P00533 numbering |
| `egfr_dccm_matrix.csv` | 1.3 MB | dynamic cross-correlation matrix, 253 x 253, unit diagonal (replicate 1 — the published trajectory) |
| `egfr_dccm_matrix_mean5.csv` | 867 KB | element-wise mean DCCM over the five replicates |
| `egfr_dccm_matrix_rep2.csv` | 868 KB | replicate 2 DCCM |
| `egfr_dccm_matrix_rep3.csv` | 867 KB | replicate 3 DCCM |
| `egfr_dccm_matrix_rep4.csv` | 866 KB | replicate 4 DCCM |
| `egfr_dccm_matrix_rep5.csv` | 867 KB | replicate 5 DCCM |
| `egfr_reference_aligned.fasta` | 4 KB | the aligned EGFR_HUMAN/712-968 row; uppercase = match-state, lowercase = insert-state. The column -> residue map every script uses |
| `fig7_dca_aligned.npy` | 512 KB | DCA scores arranged in the DCCM residue order, for the mirror heatmap |
| `pair_table.csv.gz` | 2.0 MB | one row per residue pair (31,878): APC coupling, SCA coupling, dynamic correlation, minimum heavy-atom distance over 328 structures, inter-chain distance, Ca-Ca distance, community and sector labels |
| `pair_table_ss.csv.gz` | 2.2 MB | the same table with DSSP secondary-structure class, backbone H-bond occupancy and the cross-strand register flag added |
| `provenance_rep2.txt` | 1 KB | replicate 2: seed, node, job id, GROMACS build, starting coordinates |
| `provenance_rep3.txt` | 1 KB | replicate 3 provenance |
| `provenance_rep4.txt` | 1 KB | replicate 4 provenance |
| `provenance_rep5.txt` | 1 KB | replicate 5 provenance |
| `rmsf.xvg` | 4 KB | per-residue RMSF, replicate 1 |
| `rmsf_rep2.xvg` | 5 KB | per-residue RMSF, replicate 2 |
| `rmsf_rep3.xvg` | 5 KB | per-residue RMSF, replicate 3 |
| `rmsf_rep4.xvg` | 4 KB | per-residue RMSF, replicate 4 |
| `rmsf_rep5.xvg` | 4 KB | per-residue RMSF, replicate 5 |
| `secondary_structure_per_residue.csv` | 5 KB | DSSP class per residue over 501 trajectory frames, with the fraction of frames beta / helix / loop |
| `site_definitions.json` | 0 KB | residue sets for the ATP site, the allosteric (EAI045 / JBJ-04-125-02) pocket and the dimer interface |

Numbering is UniProt **P00533** throughout. The analysed domain is **714–966**: of the 257
residues spanning 712–968, four (712, 713, 967, 968) are insert-state columns and were not
scored by plmc, so the network has 253 nodes.

Five independent 100 ns trajectories exist. **Replicate 1 is the one the main text reports**;
replicates 2–5 are the reproducibility check added in revision (Section "Replicate
validation", `docs/CORRECTIONS.md` Section 9.1).

## Not committed

| file | size | how to obtain |
|---|---|---|
| `PF07714.alignment.full` | 2.5 GB | InterPro/Pfam, family PF07714 (`PK_Tyr_Ser-Thr`; previous ID `Pkinase_Tyr`), full alignment, Stockholm |
| `kinase_final_aligned.fasta` | 122 MB | regenerate with `pipeline/01_build_alignment.sh` |
| `egfr.plmc` | 107 MB | regenerate with `pipeline/02_run_plmc.sh` |
| `md_clean.xtc` | 8.2 GB per replicate | regenerate with `pipeline/04_run_md.sh` (replicate 1) or `pipeline/05_md_replicates/` (2–5) |
| `Vsca.npy`, `Csca.npy`, `Di` | — | regenerate with the pySCA stage, `pipeline/06_sca_sectors.md` |
