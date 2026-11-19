# Pipeline provenance audit — `~/Downloads/EGFR`

Purpose: establish which files and which script versions actually produced the published results, so corrections are applied to the real pipeline rather than the one the Methods section describes.

## 1. The pipeline that produced the published results

Every downstream artefact traces to this chain. Each link is confirmed by at least two independent pieces of evidence (file content + log).

```
PF07714.alignment.full                 Stockholm, 400,321 sequences, 3,839 columns
  │                                    #=GF ID   PK_Tyr_Ser-Thr
  │                                    #=GF AC   PF07714.24
  │                                    #=GF DE   Protein tyrosine and serine/threonine kinase
  │                                    #=GF PI   Pkinase_Tyr;
  │
  ├─[gaps stripped]──► kinase_ty_raw.fasta        (not retained; 400,321 per CD-HIT log)
  │
  ├─ run_cdhit.sh:  cd-hit -c 0.62 -n 4 -M 32000 -T 8
  │     cdhit_26371612.out: "total seq: 400321" ... "400321 finished 32911 clusters"
  │                              ▼
  │                    kinase_ty_62.fasta         32,911 sequences (unaligned)
  │                              ▼
  │                    winner_ids.txt             cluster representative IDs
  │
  ├─ extract_winners.py  (streams PF07714.alignment.full, keeps winners, gaps intact)
  │  + rescue_egfr.py    (appends EGFR_HUMAN/712-968, not a 62% representative)
  │                              ▼
  │                    kinase_final_aligned.fasta 32,912 sequences × 3,839 columns
  │                              ▼
  └─ run_plmc.sh:  plmc -f EGFR_HUMAN/712-968 -t 0.2 -m 500 -le 1.0 -lh 0.01 -n 8
        plmc_26375071.err:
            "Found focus EGFR_HUMAN/712-968 as sequence 32912"
            "32672 valid sequences out of 32912"
            "253 sites out of 3839"
            "Effective number of samples: 27612.8  (80% identical neighborhood = 1.000 samples)"
                               ▼
                    egfr_couplings.txt            31,878 lines = C(253,2)
                               ▼  make_master.py / rank_couplings.py
                    egfr_couplings_translated.csv 31,878 rows, residues 714–966
                               ▼  1_8_egfrcommunity.py (462 edges, weighted Louvain, 100 seeds)
                    egfr_communities_PRODUCTION_1.8L.csv   253 residues, 15 communities
```

Everything reported in the manuscript — communities, centrality, inter-module flux, mutation mapping, the DCA side of the DCCM comparison — descends from `egfr_couplings_translated.csv`, and therefore from a **32,912-sequence alignment**.

## 2. Abandoned branches present in the folder

These files exist, look authoritative, and did **not** feed the published results. Three of the manuscript's Methods numbers come from them.

| File / script | Size | Status | Evidence |
|---|---|---|---|
| `kinase_80_winners.fasta` | 215,428 seqs | abandoned | no downstream consumer; not the plmc input |
| `kinase_60_winners.fasta` | **93,773 seqs** | abandoned | **this is the manuscript's reported final depth**; unaligned (26 MB ÷ 93,773 ≈ 277 B/seq), so it could not have been a plmc input |
| `convert_safe.py` | — | abandoned | reads `PF00069.alignment.full` (Ser/Thr kinase family), not PF07714 |
| `run_sfetch.sh` | — | abandoned | `esl-sfetch` from `PF00069.alignment.full`; would have built `kinase_final_aligned.fasta` from the wrong family |
| `run_hhfilter.sh` | — | abandoned | `hhfilter -id 80 -cov 50`; the likely origin of the "80% identity" step in Methods 3.1, which was hhfilter, not CD-HIT, and did not feed the final alignment |
| `calculate_dccm.py` | — | superseded | assumes `n_atoms` from `covar.ascii` matches the 253 community residues; it does not (259 vs 253). Superseded by `calculate_dccm_fast.py`, which subsets correctly |
| `communitydecoder.py` | — | superseded | indexes `df.index` instead of `df['residue']`; superseded by `communitydecoderfixed.py` |
| `quantify_allostery.py` | — | superseded | raw-sum flux only; `quantifyallostery_normalizationfactor.py` adds the density metric but still ranks by raw sum |

**The "3,667 columns" figure** is a stale code comment in `rescue_egfr.py` (`# Save it with all 3,667 columns intact`). The file it writes is 3,839 columns wide. **The "1,051,876 sequences" figure** does not correspond to any file in the folder; the PF07714 full alignment contains 400,321 records.

## 3. Verified bugs in the production scripts

These affect published numbers and are fixed in the re-analysis.

**3.1 `generate_topology_data.py` — unweighted Louvain in the percolation sweep.**
```python
q = community_louvain.modularity(community_louvain.best_partition(G), G)   # no weight=
```
The production partition (`1_8_egfrcommunity.py`) uses `weight='weight'` in both `best_partition` and `modularity`. Consequence: `fig4_percolation_data.csv` reports Q = 0.7986 at 1.80L while the manuscript text quotes Q = 0.8036 from the weighted production run. Figure 2 and the Results text are computed differently.

**3.2 `quantifyallostery_normalizationfactor.py` — unordered community pairs split across two rows.**
```python
highway_stats = allosteric_links.groupby(['comm_i','comm_j'])['score'].agg(['sum','count'])
```
`comm_i`/`comm_j` inherit residue order from the couplings file, so a community pair appears as both (A,B) and (B,A). The table has 126 rows where there are only C(15,2) = 105 community pairs. The dimer-interface/tail-anchor pathway is split: (5,10) sum = 54.291 over 294 pairs, (10,5) sum = 3.378 over 14 pairs; 294 + 14 = 308 = |C5| × |C10| = 14 × 22. **The published Φ = 54.29 is a partial sum.**

**3.3 `L = 257` hard-coded where the network has 253 nodes.**
Present in `1_8_egfrcommunity.py`, `Qtitration.py`, `optimize_q.py`, `connectivityQ.py`, `build_network.py`, `generate_topology_data.py`. The edge count at "1.80L" is `int(257 × 1.80) = 462`; with the true L it would be `int(253 × 1.80) = 455`. `connectivityQ.py` additionally divides the giant-component size by 257 rather than 253, under-reporting connectivity — though that script did not produce the published figure (`generate_topology_data.py` correctly uses 253).

**3.4 `correlation_test.py` — no RNG seed.** The 500-pair null is redrawn on every execution, so the published baseline (0.156) and p-value are not reproducible.

**3.5 Sequence-separation filter.** Every script that filters uses `|i−j| > 5` (`correlation_test.py`, `correlation_test_withneighbor.py`, `allostericcircuits.py`, `rank_couplings.py`). The manuscript states `|i−j| > 10` in §2.7 and in the Figure 6D caption. No code in the folder uses 10.

## 4. The 257 → 253 discrepancy, resolved

`EGFR_HUMAN/712-968` occupies 257 alphabetic positions in the alignment: **253 uppercase (match state) and 4 lowercase (insert state)**. `plmc`'s focus mode scores match-state columns only, which is why it reports "253 sites out of 3839" and emits C(253,2) = 31,878 pairs.

The four unscored residues are **712, 713, 967, 968** — the first two and last two residues of the domain definition. `get_focus_columns.py` and `make_master.py` use `char.isalpha()`, which counts both cases and yields 257; the residue-number mapping remains correct because each column maps to its own residue, and the four insert-state columns simply never appear in the couplings file.

This also means the kinase domain as analysed is **714–966**, not 712–968.

## 5. Confirmed-correct items

Not everything is wrong. These check out against the data and need no change:

- Top five by node strength: 965 (24.74), 717, 799, 851, 964 — ranks 1–5 confirmed.
- Top five by betweenness: 792 (0.3118), 838, 793 (0.2651), 853, 852 — ranks 1–5 confirmed.
- Residue 793 dual role: node strength 18.32 (rank 7) and betweenness 0.265 (rank 3) — both in the top 10, as claimed.
- 1.80L is genuinely the first threshold reaching 100% connectivity (1.70L and 1.75L both give 99.6%, i.e. 252/253).
- Per-community mean RMSF values (C7 = 0.2135, C3 = 0.0545, C14 = 0.0636) reproduce exactly.
- The RMSF merge is sound: `rmsf.xvg` is numbered 712–970, matching the community file's residue numbering; 253 residues merge correctly.
- MD construct and length: 259 Cα, residues 712–970, 50,001 frames over 100 ns (`covar.log`).
- `egfr_couplings.txt` → `egfr_couplings_translated.csv` column-to-residue mapping is correct.

## 6. Analyses present in the data but absent from the manuscript

`fig8_perturbation_results.csv` (from `run_perturbation.py`) contains a network knockout analysis never reported:

| Residue | Global-efficiency drop on removal |
|---|---|
| 793 | +2.26% |
| 858 | +1.44% |
| 965 | +1.00% |
| 884 | −0.01% |
| 733 | −0.02% |

The two negative values mean removing residues 884 and 733 marginally *increases* global network efficiency. This matters because the manuscript names 733 and 884 as "the primary high-centrality residues" (§2.3 robustness check) and describes 883/884 as ranking "among the primary infrastructure hubs" (§2.4) — in the production network 733 ranks 34th and 884 ranks 98th of 253 by node strength. Either report this analysis or withdraw the hub claims for 733/883/884.
