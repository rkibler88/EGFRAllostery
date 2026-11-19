# Script provenance

Every script found in the two working directories this repository was assembled from
(`Downloads/EGFR`, 55 files, and `Downloads/EGFR_Final_Results`, 20 files) is listed below:
where it went, or why it did not come along. 37 of the 75 are carried into the repository.

Nothing here is deleted history — the working directories are unchanged. This table exists so
that a reader who finds, say, `fig7spoof.py` in the author's folder can see that it is the
submitted version of a figure whose published form is `analysis/make_figure6.py`.

Status values:

| status | meaning |
|---|---|
| kept | carried over unchanged (path may be renamed) |
| adapted | carried over with the change described |
| superseded | a later script in the repository does the same job, usually with a correction |
| withdrawn | the analysis it implements was retracted in revision |
| excluded | scratch, duplicate or unrelated |
| duplicate | byte-identical to another copy that is kept |

| original | status | in this repository | note |
|---|---|---|---|
| `EGFR/1m17fixed.py` | kept | `visualization/write_1m17_pml.py` | Writes the 1M17 comparison script, applying the +24 offset between 1M17 and UniProt numbering. |
| `EGFR/alphafoldpymol.py` | kept | `visualization/write_alphafold_pml.py` | Writes the AlphaFold community-colouring PyMOL script. |
| `EGFR/ca_extract.sh` | kept | `pipeline/helpers/ca_extract.sh` | C-alpha trajectory extraction before the covariance step. |
| `EGFR/calculate_dccm_fast.py` | kept | `pipeline/helpers/calculate_dccm_fast.py` | 3Nx3N Cartesian covariance to NxN correlation, subset to the 253 scored residues. Produces `egfr_dccm_matrix.csv`. |
| `EGFR/communitydecoderfixed.py` | kept | `analysis/community_landmarks.py` | Assigns each community its dominant structural landmarks. |
| `EGFR/extract_winners.py` | kept | `pipeline/helpers/extract_winners.py` | Biopython alternative to `esl-sfetch` for pulling cluster representatives out of the Stockholm alignment. |
| `EGFR/get_focus_columns.py` | kept | `pipeline/helpers/get_focus_columns.py` | Builds the plmc focus-column list (257 columns requested; 253 scored). |
| `EGFR/ions.mdp` | kept | `pipeline/mdp/ions.mdp` | Ion placement. |
| `EGFR/make_master.py` | kept | `pipeline/helpers/make_master.py` | Maps alignment columns to UniProt P00533 numbering; produces `egfr_couplings_translated.csv`. |
| `EGFR/md.mdp` | kept | `pipeline/mdp/md.mdp` | Production MD parameters. |
| `EGFR/md_replicates/README_replicates.md` | kept | `pipeline/05_md_replicates/README.md` | Cluster configuration and run notes. |
| `EGFR/md_replicates/aggregate_replicates.py` | kept | `pipeline/05_md_replicates/aggregate_replicates.py` | Between-replicate agreement; produces the `replicate_*.csv` tables. |
| `EGFR/md_replicates/dccm_from_covar.py` | kept | `pipeline/05_md_replicates/dccm_from_covar.py` | Replicate DCCM from `gmx covar -ascii` output. |
| `EGFR/md_replicates/diagnose_replicates.sh` | kept | `pipeline/05_md_replicates/diagnose_replicates.sh` | Post-mortem for a failed array. |
| `EGFR/md_replicates/eq_seed.mdp` | kept | `pipeline/05_md_replicates/eq_seed.mdp` | Per-replicate equilibration with an independent velocity seed. |
| `EGFR/md_replicates/md_rep.mdp` | kept | `pipeline/05_md_replicates/md_rep.mdp` | Replicate production parameters. |
| `EGFR/md_replicates/md_replicates.sbatch` | kept | `pipeline/05_md_replicates/md_replicates.sbatch` | Five-replicate array job (CORRECTIONS.md Section 9.1). |
| `EGFR/md_replicates/monitor_replicates.sh` | kept | `pipeline/05_md_replicates/monitor_replicates.sh` | Watches a running array over ssh. |
| `EGFR/md_replicates/rerun_analysis.sbatch` | kept | `pipeline/05_md_replicates/rerun_analysis.sbatch` | Covariance and RMSF re-analysis for the replicate array. |
| `EGFR/minim.mdp` | kept | `pipeline/mdp/minim.mdp` | Energy minimisation. |
| `EGFR/mirror_matrix.py` | kept | `pipeline/helpers/mirror_matrix.py` | Aligns the DCA scores to the DCCM residue order; produces `fig7_dca_aligned.npy`. |
| `EGFR/npt.mdp` | kept | `pipeline/mdp/npt.mdp` | NPT equilibration. |
| `EGFR/nvt.mdp` | kept | `pipeline/mdp/nvt.mdp` | NVT equilibration. |
| `EGFR/pymolvisualizer.py` | kept | `visualization/write_hubs_bridges_pml.py` | Writes the hub/bridge PyMOL script. |
| `EGFR/rmsfcheck.py` | kept | `pipeline/helpers/community_rmsf.py` | Mean RMSF per community; produces `community_dynamics_summary.csv`. |
| `EGFR/run_cdhit.sh` | kept | `pipeline/helpers/run_cdhit.sh` | Stage 01 redundancy filter. Header comments described a 90% / `-n 5` run; the command executed, and used, is `-c 0.62 -n 4`. |
| `EGFR/run_md.sh` | kept | `pipeline/helpers/run_md_sbatch.sh` | The cluster job actually submitted for the 100 ns production run (GROMACS 2025.2, GPU offload). |
| `EGFR/run_plmc.sh` | kept | `pipeline/02_run_plmc.sh` | Stage 02, annotated with the plmc run statistics and the column-5/column-6 score caveat. |
| `EGFR/run_sfetch.sh` | kept | `pipeline/helpers/run_sfetch.sh` | Stage 01 sequence retrieval. The copy here still names PF00069; `pipeline/01_build_alignment.sh` runs it against PF07714, the family the published couplings come from. |
| `EGFR/viz_alphafold_FINAL.pml` | kept | `visualization/viz_alphafold_communities.pml` | Generated PyMOL session for the AlphaFold structural panels. |
| `EGFR/viz_final_comparison.pml` | kept | `visualization/viz_1m17_comparison.pml` | Generated PyMOL session, 1M17 vs AlphaFold. |
| `EGFR/viz_hubs_and_bridges.pml` | kept | `visualization/viz_hubs_and_bridges.pml` | Generated PyMOL session, hubs and bridges. |
| `EGFR_Final_Results/fig2pdbgen.py` | kept | `visualization/map_rmsf_to_bfactor.py` | Writes RMSF into the B-factor column for the rigidity-map panel. |
| `EGFR_Final_Results/makefig6.py` | kept | `visualization/pymol_switch_handshake.py` | PyMOL render of the switch/handshake structural panel on 1M17. |
| `EGFR/convert_safe.py` | adapted | `pipeline/helpers/convert_safe.py` | Stockholm to FASTA. Input family is now an argument; the original hardcoded PF00069. |
| `EGFR/rescue_egfr.py` | adapted | `pipeline/helpers/rescue_egfr.py` | Recovers the aligned EGFR_HUMAN/712-968 row. Now also appends it to the representative alignment, which is what stage 01 needs; the original only wrote the reference row. |
| `EGFR/run_perturbation.py` | adapted | `analysis/network_knockout.py` | Global-efficiency knockout. Rewritten to build the network the same way as the rest of the repository and to write `results/tables/network_knockout.csv`; values are unchanged (793: 2.26%, 858: 1.44%, 965: 1.00%). |
| `EGFR/1_8_egfrcommunity.py` | superseded | `pipeline/03_build_network.py` | The production partition. Rewritten so the two meanings of L (257 focus columns, 253 scored positions) are explicit; edge count unchanged at 462. |
| `EGFR/Qtitration.py` | superseded | — | Earlier, coarser version of the same titration. |
| `EGFR/build_network.py` | superseded | — | Exploratory multiplier sweep that wrote `egfr_communities_{m}L.csv`; the deposited partition is the 1.80L run. |
| `EGFR/calculate_dccm.py` | superseded | — | Same mathematics in an O(N^2) Python loop; the vectorised version replaced it. |
| `EGFR/communitydecoder.py` | superseded | — | Indexing bug; fixed in the next version. |
| `EGFR/connectivityQ.py` | superseded | — | Giant-component / Q table at seven thresholds; replaced by the percolation scan in `analysis/make_figure2.py`. |
| `EGFR/correlation_test.py` | superseded | — | Coupling-vs-DCCM enrichment against a random-pair null. The distance-matched null in `analysis/make_figure6.py` is the test the manuscript now reports. |
| `EGFR/correlation_test_withneighbor.py` | superseded | — | Same test with sequence neighbours included. |
| `EGFR/generate_topology_data.py` | superseded | — | Produced `fig4_percolation_data.csv` with unweighted modularity while the manuscript quoted the weighted value. |
| `EGFR/get_allosteric_paths.py` | superseded | — | Two hardcoded Dijkstra paths; `analysis/compute_pathways.py` is the version the manuscript uses. |
| `EGFR/optimize_q.py` | superseded | — | Q titration over 1.0-4.0L. `analysis/make_figure2.py` recomputes the sweep with weighted Louvain, which is what the manuscript quotes. |
| `EGFR/quantify_allostery.py` | superseded | — | Inter-module coupling by summing signed APC scores. |
| `EGFR/quantifyallostery_normalizationfactor.py` | superseded | — | Adds a density normalisation to the same sum. `analysis/make_figure5.py` reports all three measures and treats the signed sum as signed (CORRECTIONS.md Section 5). |
| `EGFR/rank_couplings.py` | superseded | — | Same column-to-residue mapping as `make_master.py` plus a top-500 slice; the full translated list is what every downstream script reads. |
| `EGFR_Final_Results/1bschematic.py` | superseded | — | Figure 1B schematic. `analysis/make_figure1.py` corrects the DFG motif (D855-F856-G857) and the alphaG span. |
| `EGFR_Final_Results/1c.py` | superseded | — | Figure 1C flowchart. The sequence count (93,773) and the convergence label did not survive revision. |
| `EGFR_Final_Results/backboneremove.py` | superseded | — | Path search with backbone edges forbidden; folded into `analysis/compute_pathways.py`. |
| `EGFR_Final_Results/dijkstra.py` | superseded | — | Shortest paths at 1.80L; replaced by `analysis/compute_pathways.py`. |
| `EGFR_Final_Results/distaldijkstra.py` | superseded | — | Same with a long-range edge bonus. |
| `EGFR_Final_Results/fig7spoof.py` | superseded | — | The submitted version of the mirror figure. Panel D means were hardcoded; `analysis/make_figure6.py` computes them and adds the distance-matched control. |
| `EGFR_Final_Results/handshakefind.py` | superseded | — | Interactive search for coupled pairs with high correlation; folded into `analysis/make_figure6.py`. |
| `EGFR_Final_Results/makefig2.py` | superseded | — | Replaced by `analysis/make_figure4.py` (RMSF). |
| `EGFR_Final_Results/makefig3.py` | superseded | — | Earlier version of the flux figure. |
| `EGFR_Final_Results/makefig3beige.py` | superseded | — | Last author version of the flux figure; replaced by `analysis/make_figure5.py`. |
| `EGFR_Final_Results/makefig4.py` | superseded | — | Replaced by `analysis/make_figure2.py` (percolation, weighted). |
| `EGFR_Final_Results/makefig5.py` | superseded | — | Replaced by `analysis/make_figure3.py` (centrality). |
| `EGFR_Final_Results/makefig7.py` | superseded | — | DCA/DCCM mirror figure; replaced by `analysis/make_figure6.py`. |
| `EGFR/allostericcircuits.py` | withdrawn | — | Built `FINAL_ALLOSTERIC_CENSUS_PAIRS.csv` / `FINAL_ALLOSTERIC_ROBUST_CIRCUITS.csv`. The circuit census rested on thresholded sums of signed couplings; the claim was withdrawn in revision. |
| `EGFR/moduleintegrator.py` | withdrawn | — | Ranked "highways" from the same census. |
| `EGFR/dccmread.py` | excluded | — | Five-line interactive check of one residue row. |
| `EGFR/qscore.py` | excluded | — | Twelve-line interactive modularity check. |
| `EGFR/run_extract.sh` | excluded | — | Ten-line sbatch wrapper around `extract_winners.py`; the module loads are in `pipeline/01_build_alignment.sh`. |
| `EGFR/run_hhfilter.sh` | excluded | — | 80% identity / 50% coverage branch. It terminates in `kinase_60_winners.fasta`, was never aligned, and did not feed the published couplings. |
| `EGFR_Final_Results/auditranking.py` | excluded | — | Prints the hub and bridge rankings; reproduced by `analysis/make_figure3.py`. |
| `EGFR_Final_Results/mmlhw2.py` | excluded | — | Unrelated linear-algebra coursework that happened to sit in the folder. |
| `EGFR_Final_Results/viz_alphafold_FINAL.pml` | duplicate | — | Byte-identical to the copy in `EGFR/`. |
| `EGFR_Final_Results/viz_final_comparison.pml` | duplicate | — | Byte-identical to the copy in `EGFR/`. |
| `EGFR_Final_Results/viz_hubs_and_bridges.pml` | duplicate | — | Byte-identical to the copy in `EGFR/`. |

## Scripts in this repository with no counterpart in the working directories

| path | why it exists |
|---|---|
| `analysis/make_figure1.py` … `make_figure6.py` | corrected replacements for the author's figure scripts, written during the error audit (`docs/CORRECTIONS.md`) |
| `analysis/compute_pathways.py` | shortest-path analysis with the backbone handling the manuscript describes |
| `analysis/verify_numbers.py` | recomputes every manuscript number from `data/processed` and exits non-zero on disagreement |
| `analysis/make_figure7_sectors.py`, `make_figure8_beta_sheet.py`, `make_figureS1_replicates.py` | the three revision figures; see the header of each for what is recomputed and what is read from a committed table |
| `pipeline/01_build_alignment.sh`, `02_run_plmc.sh`, `03_build_network.py`, `04_run_md.sh` | stage drivers that tie the author's helper scripts together in the order they were run |
| `pipeline/06_sca_sectors.md` | the pySCA stage, documented rather than scripted — see that file |
