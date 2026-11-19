# EGFR kinase domain: co-evolutionary architecture and allosteric network

Code, data, figures and manuscript source for a co-evolutionary network analysis of the EGFR
kinase domain (UniProt **P00533**, residues 714–966), compared against molecular dynamics and
against the conformational ensemble of solved EGFR structures.

Direct coupling analysis over Pfam **PF07714** yields a residue interaction network that
partitions into 15 communities (weighted *Q* = 0.8036 at 462 edges). Centrality separates
structural anchors from inter-module bridges. Five 100 ns simulations of the monomeric domain
ask whether co-evolutionary coupling corresponds to correlated motion — and, under a spatial
control, it does not. Statistical coupling analysis of the same alignment recovers one sector
that is spatially compact beyond conservation and sequence contiguity.

## Quickstart

```bash
conda env create -f environment.yml
conda activate egfr-allostery

make verify     # recompute every manuscript number from data/processed
make figures    # regenerate all data-driven figure panels into results/figures
```

`make verify` is the reproducibility test. It recomputes the partition, the modularity, the
percolation threshold, the centrality rankings, the RMSF ratios, the inter-module coupling and
the residue identities straight from the committed data, and exits non-zero if any value
disagrees with what the manuscript states. It takes about a minute.

## Layout

```
data/processed/    every input the analysis needs (committed, with SHA256SUMS)
pipeline/          upstream stages: alignment -> DCA -> network -> MD -> replicates -> SCA
analysis/          figure generation, pathway and knockout analysis, the verification test
visualization/     PyMOL sessions and the scripts that write them (structural panels)
results/           figures and tables produced by analysis/
manuscript/        LaTeX source, bibliography, and the submitted figures
docs/              corrections, revision changelog, script provenance, method selection
```

`pipeline/` is recorded for provenance and is **not** runnable from a clone — it needs the
2.5 GB Pfam alignment and the 8.2 GB-per-replicate trajectories. Everything in `analysis/`
runs from the committed data.

## What the analysis supports, and what it does not

This repository is the corrected version of an earlier preprint. Several headline claims did
not survive recomputation and have been withdrawn rather than patched. The full record is in
[`docs/CORRECTIONS.md`](docs/CORRECTIONS.md); the short version:

- **The co-evolution/dynamics correspondence is spatial.** The reported 1.83-fold enrichment
  in correlated motion among co-evolving distal pairs reproduces exactly — but it is measured
  against a null matched only on *sequence* separation. The top-coupled pairs sit at a median
  11.4 Å Cα–Cα against 24.3 Å for the rest. Against a null matched on Cα–Cα distance the
  enrichment is 1.02-fold (bootstrap 95% CI 0.94–1.11), in every one of the five replicates.
- **There are almost no strong long-range couplings to test.** Of 55 pairs scoring above 2.0,
  one survives a 10 Å cutoff and none survive 12 Å; 92% of the 462 network edges are spatial
  contacts with a median length of 3.9 Å. The network is close to a contact map, which is the
  expected behaviour of DCA. Across the 328 solved EGFR kinase-domain chains, 43% of the
  top-500 "long-range" couplings are in contact in at least one state (4.3% of distal pairs).
- **The alignment is superfamily-wide, not RTK-specific.** PF07714 spans tyrosine and
  serine/threonine kinases. After 62% clustering the retained set is 34.8% Viridiplantae,
  26.2% Metazoa, 24.0% Fungi, and contains **one** ErbB-family sequence — EGFR itself, which
  had to be manually re-inserted. Claims of receptor-tyrosine-kinase specificity were removed.
- **The dynamics are five 100 ns replicates of a monomer** (500 ns aggregate). No conformational
  transition, mutant, or dimer was simulated. Replicate agreement is precision, not accuracy:
  all five share a force field and a starting model.
- **Numbering was wrong in places.** Residue 853 is isoleucine; the DFG motif is D855-F856-G857.
  Earlier figures and text labelled 853 as DFG.

One analysis that would settle the remaining open question — alignments restricted to receptor
tyrosine kinases and to the ErbB family — has not been run. It is specified in
[`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) Section 9.2.

## Scripts

| script | produces | notes |
|---|---|---|
| `analysis/make_figure1.py` | Figure 1 panels B, C | replaces `1bschematic.py`, `1c.py` |
| `analysis/make_figure2.py` | Figure 2 (percolation) | weighted Louvain; replaces `makefig4.py` |
| `analysis/make_figure3.py` | Figure 3 panels A, B (centrality) | replaces `makefig5.py` |
| `analysis/make_figure4.py` | Figure 4 panels A, C (RMSF) | replaces `makefig2.py` |
| `analysis/make_figure5.py` | Figure 5 panels B, C (inter-module coupling) | replaces `makefig3.py` |
| `analysis/make_figure6.py` | Figure 6 (DCA vs DCCM) | computes panel D; replaces `fig7spoof.py` |
| `analysis/make_figure7_sectors.py` | Figure 7 (ensemble contacts, SCA sectors) | panel B's conservation null is read from a table — see the docstring |
| `analysis/make_figure8_beta_sheet.py` | Figure 8 (β-sheet pairing) | fully recomputed from `pair_table_ss.csv.gz` |
| `analysis/make_figureS1_replicates.py` | Figure S1 (replicate agreement) | fully recomputed from the replicate tables |
| `analysis/compute_pathways.py` | shortest-path results | |
| `analysis/network_knockout.py` | `results/tables/network_knockout.csv` | global-efficiency knockout |
| `analysis/community_landmarks.py` | community → structural landmark table | |
| `analysis/verify_numbers.py` | the reproducibility test | |

The 3D panels (1A, 3C, 4B, 5A) are PyMOL renders from `visualization/`. They carry no
numerical annotations and are not regenerated by `make figures`.

Every script in the two original working directories is accounted for in
[`docs/SCRIPT_PROVENANCE.md`](docs/SCRIPT_PROVENANCE.md) — kept, adapted, superseded or
excluded, with the reason.

## Notes for anyone extending this

- The coupling score is **average-product-corrected** and signed, with mean ≈ 0 and 55.5% of
  values negative. Do not sum it over arbitrary residue-pair sets as though it were a
  magnitude — doing so accumulates background and scales with set size (the original
  inter-module metric did this, and 75 of 105 module pairs came out negative). Restrict sums
  to network edges, or use a rank statistic.
- Proximity filtering throughout the original pipeline is by **sequence index** (`|i-j| > 5`),
  never by 3D distance. A sequence filter does not control for spatial proximity.
- `egfr_couplings.txt` column 5 (Frobenius norm) is zero throughout; column 6 is the score.
- Two lengths are in play: 257 focus columns requested from plmc, 253 match-state positions
  actually scored. The network has 253 nodes; the 462-edge threshold is `int(257 × 1.80)`.

## Citing

See `CITATION.cff`. Code is MIT; manuscript text and figures are CC BY 4.0.
