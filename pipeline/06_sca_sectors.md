# Stage 06 — statistical coupling analysis and sector decomposition

This stage is **documented, not scripted**. It ran interactively during revision; its outputs
are committed, and what follows is enough to repeat it. Figure 7 and the sector rows of
`results/tables/` come from here.

## Inputs

* `kinase_final_aligned.fasta` — the stage 01 alignment (32,912 sequences, not committed)
* `data/processed/egfr_clean.pdb` and the 328 EGFR kinase-domain chains listed in
  `data/processed/contact_qc.json` (RCSB; chain IDs given there)

## Method

pySCA 6.1 (`pysca.scaTools`), on the 253 match-state positions 714–966:

1. `scaMat` on the filtered alignment gives the SCA coupling matrix `Csca`, the positional
   conservation vector `Di`, and the eigendecomposition `Vsca` / `Lsca`; a randomised
   alignment gives `Lrand` for the significance cutoff.
2. `kpos = 6` significant eigenmodes, rotated with `st.rotICA(Vsca, kmax=6)`.
3. `st.icList(Vpica, kpos, Csca, p_cut=0.95)` assigns positions to independent components;
   components 1–6 are reported as sectors in `results/tables/sca_sectors.csv` with their
   cutoffs. Random seed 0.

## Spatial tests (`results/tables/sector_spatial_tests.csv`)

For each sector: contiguity in the 8 Å contact graph, mean pairwise Cα–Cα distance, and two
null models over 2,000 draws each — positions resampled within deciles of `Di`
(conservation-matched), and runs of the same sequence lengths placed at random
(run-length-matched). Sector 2 is the only sector compact beyond both nulls
(10.94 Å observed, p < 0.0005).

## What is committed

`sca_sectors.csv`, `sector_spatial_tests.csv`, `sector2_null_comparison.csv`,
`residue_level_bootstrap.csv`, `results/egfr_sca_sectors.pdb` (sectors written into the
B-factor column), and the per-pair table `data/processed/pair_table.csv.gz`, which carries the
SCA coupling, the APC coupling, the dynamic correlation, the minimum heavy-atom distance over
the 328 structures and the Cα–Cα distance for all 31,878 pairs.

The per-position conservation vector `Di` is **not** committed, so the conservation-matched
null cannot be recomputed from a clone; `analysis/make_figure7_sectors.py` recomputes
everything else and reads that null's published values from `sector_spatial_tests.csv`.
