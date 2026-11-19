# CORRECTIONS — *Evolutionary Hubs Coordinate the Conformational Switch of the Epidermal Growth Factor Receptor*

Every value below was recomputed from the source files in `~/Downloads/EGFR`. Supporting tables and figures are listed in §9. Where a claim could not be verified or was contradicted, the source file and the computation are named so the result can be checked independently.

**Reproducibility note.** The published community partition reproduces *exactly* from `egfr_couplings_translated.csv` — weighted Louvain, 462 edges, best of 100 seeds gives Q = 0.8036, 15 communities, and 100% co-membership agreement with `egfr_communities_PRODUCTION_1.8L.csv`. The published enrichment reproduces too (1.83× vs 1.84×). **The computations are sound; the problems are in how they are described and what is inferred from them.**

---

## 1. Severity summary

| # | Issue | Class | Effect |
|---|---|---|---|
| 1 | Alignment statistics describe an abandoned pipeline branch | **Factual error** | 4 numbers in Methods 3.1 and Figure 1C wrong |
| 2 | Alignment is 74% non-metazoan; 1 ErbB sequence | **Framing invalid** | Abstract + Introduction RTK claims unsupported |
| 3 | T790M community claim contradicted by the data | **Factual error** | §2.9 and abstract claim false |
| 4 | Φ built on a signed, mean-zero score | **Result does not hold** | §2.6 headline pathways not reproducible under any well-posed metric |
| 5 | DCCM enrichment is a spatial-proximity artefact | **Result does not hold** | §2.7–2.8 and Fig 6D central claim collapses |
| 6 | Residue 853 is Ile, not the DFG motif | **Factual error** | §2.9 shortest-path target, Fig 3A caption, C13 label |
| 7 | Community 5 is not the αG helix; dimer residues elsewhere | **Factual error** | §2.6 pathway annotation |
| 8 | "Dimer–tail" pathway is backbone-adjacent | **Claim inverted** | §2.6 explicitly denies this |
| 9 | Local-coupling-bias model unsupported | **Claim withdrawn** | §2.9 |
| 10 | Two references wrong, one fabricated | **Integrity risk** | Bibliography |
| 11 | Robustness claim false (communities 14–16; betweenness unstable) | **Factual error** | §2.3 |
| 12 | Hub claims for 733/883/884 unsupported | **Factual error** | §2.3, §2.4 |
| 13 | Arithmetic, filter and definition inconsistencies | **Minor** | Several sections |

---

## 2. Alignment and DCA statistics (Methods 3.1–3.2, §2.1–2.2, Figure 1C)

The manuscript describes the **80% → 60% CD-HIT branch**, which terminates in `kinase_60_winners.fasta` (93,773 sequences, unaligned). The published couplings come from a **single 62% pass** terminating in `kinase_final_aligned.fasta`. Evidence: `run_cdhit.sh`, `cdhit_26371612.out`, `run_plmc.sh`, `plmc_26375071.err`.

| Manuscript | Verified | Source |
|---|---|---|
| "over 1,051,876 sequences" | **400,321** | `PF07714.alignment.full` record count |
| "3,667 alignment columns" | **3,839** | plmc: "253 sites out of 3839"; stale comment in `rescue_egfr.py` is the origin of 3,667 |
| "clustered at 80% identity and again at 62%" | **one pass: `cd-hit -c 0.62 -n 4`** | `run_cdhit.sh`; log shows 400,321 → 32,911 clusters |
| "retained 93,773 sequences" | **32,912** (32,672 valid) | plmc: "32672 valid sequences out of 32912" |
| "257 focus columns" | 257 residues, **253 match-state scored** | 253 uppercase + 4 lowercase in `egfr_reference_aligned.fasta` |
| N/L ≈ 370 (Figure 1C) | **N/L = 130.1; Neff/L = 109.1** | Neff = 27,612.8 from plmc |
| Neff | **not reported in the manuscript — 27,612.8** | `plmc_26375071.err` line 5 |

**Replacement text for Methods 3.1 (final paragraph):**

> The PF07714 (PK_Tyr_Ser-Thr) full alignment was downloaded in Stockholm format and contained 400,321 sequences of 3,839 columns. Gap characters were preserved during conversion to FASTA using Biopython v1.81. Phylogenetic redundancy was reduced in a single pass with CD-HIT v4.7 at 62% identity (`-c 0.62 -n 4`), yielding 32,911 cluster representatives. Because EGFR_HUMAN/712-968 was not selected as a representative, its aligned row was recovered from the original Stockholm file and appended, giving a final alignment of 32,912 sequences. Of the 257 EGFR residues spanning 712–968, 253 occupy match-state (uppercase) columns and were carried into the coupling analysis; the four insert-state residues (712, 713, 967 and 968) were not scored. The analysed domain is therefore 714–966.

**Replacement text for Methods 3.2 (add after the parameter list):**

> plmc reported 32,672 valid sequences and an effective sample size of Neff = 27,612.8 after similarity reweighting, i.e. Neff/L = 109 for L = 253 scored positions. Because the alignment had already been reduced to 62%-identity cluster representatives, the `-t 0.2` reweighting threshold had almost no effect: plmc reports a mean 80%-identity neighbourhood of 1.000 sequences.

**Also correct:** the reported coupling score is the **average-product-corrected (APC) score** (column 6 of `egfr_couplings.txt`), not the raw Frobenius norm. Column 5, where plmc writes the Frobenius norm, is zero throughout this file. Methods 3.2 currently describes only the Frobenius norm and must state that APC correction was applied. This matters for §2.6 (see §5 below).

**Figure 1C** must read "32,912 sequences (Neff/L ≈ 109)" and must drop "specialized".

---

## 3. Family scope and alignment composition (Abstract, Introduction, Figure 1C)

Two independent problems.

**3.1 PF07714 is not a tyrosine-kinase family.** The downloaded file's own header reads:

```
#=GF ID   PK_Tyr_Ser-Thr
#=GF AC   PF07714.24
#=GF DE   Protein tyrosine and serine/threonine kinase
#=GF PI   Pkinase_Tyr;
```

`Pkinase_Tyr` is the **previous** identifier. Methods 3.1 cites the retired name.

**3.2 The retained alignment is predominantly non-metazoan.** Organism mnemonics were mapped to NCBI taxonomy for all 32,912 retained and 400,321 raw sequences (`alignment_composition.csv`):

| Group | Raw PF07714 | Retained for DCA | Change |
|---|---|---|---|
| Viridiplantae | 43.9% | **34.8%** | 0.79× |
| Metazoa | 45.4% | **26.2%** | 0.58× |
| Fungi | 7.5% | **24.0%** | 3.21× |
| Other eukaryota | 2.5% | **11.7%** | 4.76× |
| Bacteria | 0.3% | 1.6% | 5.64× |

The 62% clustering did exactly what aggressive clustering does: it collapsed the mutually similar vertebrate sequences and retained the divergent ones.

| | Raw | Retained |
|---|---|---|
| Vertebrate sequences | 143,257 (35.8%) | **1,158 (3.5%)** |
| Human sequences | 403 | **2** |
| ErbB-family sequences | 17 | **1 — EGFR itself, manually rescued** |

**Consequence.** These statements are not supported and must be removed or rewritten:

- Abstract: "a redundancy-filtered **tyrosine kinase** alignment of 93,773 sequences"
- Abstract: "evolutionary analyses have often aggregated diverse kinase families, potentially obscuring signals specific to receptor tyrosine kinases" — the design aggregates diverse kinase families, including serine/threonine kinases of plants and fungi
- Introduction: "we move beyond universal kinase motifs to identify the community structures **specific to the receptor tyrosine kinase signaling mechanism**"
- Figure 1C: "a **specialized** alignment of 93,773 **tyrosine kinase** sequences"

**Replacement text for the abstract's data sentence:**

> Using a redundancy-filtered alignment of 32,912 sequences from the Pfam protein tyrosine and serine/threonine kinase family (PF07714; Neff/L = 109), we applied Direct Coupling Analysis to infer co-evolutionary constraints across the EGFR kinase domain.

**Replacement text for the Introduction's scope claim:**

> Because PF07714 spans tyrosine and serine/threonine kinases across eukaryotes, the couplings recovered here reflect constraints conserved across the kinase superfamily rather than features specific to receptor tyrosine kinases or to the ErbB family. Resolving ErbB-specific architecture would require a subfamily-restricted alignment, which we address in the Limitations section.

**The honest framing is stronger than the current one**: a superfamily-wide analysis that recovers kinase-universal architecture is a legitimate and interpretable result. Claiming RTK specificity from an alignment containing one ErbB sequence is not.

---

## 4. T790M (Abstract, Introduction, §2.9)

**4.1 The community claim is false.** §2.9 states that T790M "maps to the Hinge Region community (Community 2), which harbors the two highest-betweenness residues in the network: residue 792 ... and residue 793". In `egfr_communities_PRODUCTION_1.8L.csv`:

| Residue | Community | Node strength (rank) | Betweenness (rank) |
|---|---|---|---|
| T790 | **2** | 11.678 (61/253) | 0.1166 (**20**/253) |
| L792 | **0** | 14.735 (26/253) | 0.3118 (1/253) |
| M793 | **13** | 18.319 (7/253) | 0.2651 (3/253) |

Neither 792 nor 793 is in Community 2. Across 100 Louvain seeds, T790 shares a community with 792 in **8%** of runs and with 793 in **8%**; 792 and 793 share a community with each other in **8%**. Across thresholds from 404 to 506 edges, 790 is consistently in one community and 792 in another. The claim is robustly false, not marginally so.

The abstract's "placed the T790M resistance mutation at the network's primary information bottleneck" must go.

**4.2 "Distal" is wrong.** T790 is the gatekeeper, lining the ATP cleft — confirmed as Thr at position 790 in the alignment. Introduction lines 50–53 group it with mutations that "act at a distance from the ATP-binding cleft", contradicting §2.9's own correct use of "gatekeeper".

**4.3 The conventional mechanism is mischaracterized.** §2.9 contrasts the network model against "steric occlusion of inhibitor binding, the conventionally cited mechanism". Since Yun et al. (2008, PNAS 105:2070–2075) the accepted account is thermodynamic: T790M restores ATP affinity, shifting the ATP/inhibitor competition.

**Replacement text for §2.9, first paragraph:**

> The T790M gatekeeper resistance mutation lies in Community 2, a 17-residue module spanning the β4–β5 region and the hinge (775–791). T790 itself ranks 20th of 253 by betweenness centrality (0.117) and 61st by node strength (11.68); the network's two highest-betweenness residues, L792 (0.312) and M793 (0.265), fall in Communities 0 and 13 respectively and are not co-assigned with T790 (8% of 100 Louvain initialisations). We therefore do not claim that T790M occupies the network's principal information bottleneck. The established resistance mechanism is a restoration of ATP affinity rather than simple steric occlusion (Yun et al., 2008); the present analysis adds only that T790 sits adjacent to, but not within, the hinge residues that carry the highest inter-module betweenness in this network.

---

## 5. Inter-module pathways (§2.6, Methods 3.4, Figure 5)

**5.1 The metric is ill-posed.** Φ_AB = Σ J_ij sums the APC-corrected score, which is signed with mean ≈ 0 by construction (across all 31,878 pairs: mean 0.0004, min −1.775, max 6.384, **55.5% negative**). Summing it over |A|×|B| pairs accumulates background: **75 of 105 community pairs have a negative Φ**, and Φ correlates *negatively* with the number of pairs (Pearson r = −0.33). Φ cannot be described as "the total evolutionary pressure maintaining coordinated behavior between modules".

**5.2 A grouping bug split pathways.** `quantifyallostery_normalizationfactor.py` groups on the ordered pair `['comm_i','comm_j']`, so each unordered community pair can appear twice. The published table has 126 rows for 105 pairs. The dimer–tail pathway was split as (5,10) = 54.291 over 294 pairs and (10,5) = 3.378 over 14 pairs; 294 + 14 = 308 = 14 × 22. **The published Φ = 54.29 is a partial sum; the correct total is 57.67.**

**5.3 Methods and Results define different statistics.** Methods 3.4 describes normalization "by the number of possible residue pairs between modules"; §2.6 defines a raw sum. Both were computed; only the raw sum was used for ranking.

**5.4 Under no well-posed metric do the two reported pathways dominate.** (`intermodule_flux_corrected.csv`)

| Pathway | Published Φ | Corrected Φ (all pairs) | Rank | Density rank | Network edges | Φ over those edges | Rank |
|---|---|---|---|---|---|---|---|
| C5–C10 "dimer–tail" | 54.29 | **57.67** | 1 | 4 | 3 | 7.58 | 5 |
| C0–C7 (unreported) | — | **56.43** | 2 | 6 | **0** | 0.00 | 28 |
| C1–C7 "P-loop–αC" | 54.66 | 54.66 | 3 | 2 | **1** | 2.05 | 22 |
| C2–C4 (unreported) | — | 37.28 | 5 | **1** | 2 | 3.57 | 13 |
| C3–C6 | — | −26.70 | 74 | — | **8** | **18.01** | **1** |

Only 63 of 462 edges are inter-community, spread over 27 of 105 community pairs. The claim that the P-loop/αC-helix pathway "serves as a positive methodological validation" and "supports the interpretive fidelity of the co-evolutionary network" rests on **one edge out of 462**.

**5.5 The novelty claim is inverted.** §2.6 states the dimer–tail pathway "is not readily explained by sequence adjacency or obvious local structural proximity". Community 5 spans residues **909–922** (contiguous) and Community 10 spans **908 and 923–943** — Community 10 immediately flanks Community 5 on both sides. The three network edges joining them are:

| Edge | \|i−j\| | J |
|---|---|---|
| 922–923 | **1** | 2.698 |
| 908–909 | **1** | 2.537 |
| 918–926 | 8 | 2.343 |

Two of three are immediately sequence-adjacent backbone pairs. The pathway *is* sequence adjacency. Note also that the |i−j| > 5 filter applied in §2.8 "to exclude trivially correlated backbone-adjacent contacts" was never applied to the pathway analysis, where it would have removed this result.

**5.6 The structural annotation is wrong.** Against experimentally observed secondary structure for P00533 (UniProt, PDB-derived):

- **Community 5 (909–922) is not a helix.** It comprises turns (909–911, 914–917) plus a short helix at 920–922 — the αG–αH loop. The long observed helix in this region is **893–908**, which is split across Communities 3, 6, 10 and 14.
- The canonical activator-side asymmetric-dimer residues — V948, M952, I953 in P00533 numbering (= V924, M928, I929 in the mature numbering used by Zhang et al. 2006) — are in Communities **3, 11 and 11**. None is in Community 5 or Community 10. V948 ranks 182/253 by node strength and 201/253 by betweenness.
- The decoder scripts use hand-typed landmark ranges (`range(905, 921)` for "alphaG-helix (Dimer Int)") that do not match the observed secondary structure. Both `quantify_allostery.py` and `quantifyallostery_normalizationfactor.py` carry *different* label sets for the same 15 communities.

**5.7 "Top 0.1%" is arithmetically impossible.** With 15 communities there are C(15,2) = 105 pairs. Two of 105 is **1.9%**; the top 0.1% of 105 is 0.105 pathways.

**Replacement text for §2.6:**

> For each pair of communities we computed the cumulative inter-module coupling Φ_AB = Σ J_ij. Because the coupling scores are average-product-corrected and therefore signed with mean approximately zero, Φ computed over all |A|×|B| residue pairs is dominated by background accumulation: 75 of the 105 community pairs return a negative value and Φ is negatively correlated with the number of contributing pairs (r = −0.33). We therefore report inter-module coupling restricted to the 462 edges that define the network partition. Under this measure the strongest inter-module channels connect Community 3 to Community 6 (8 edges, Σ J = 18.01) and Community 3 to Community 14 (7 edges, Σ J = 12.95). Only 63 of the 462 edges are inter-community, distributed across 27 of the 105 possible module pairs, consistent with the high modularity of the partition (Q = 0.8036). We note that the αG–αH loop (Community 5) and the αH region (Community 10) are joined by three edges of which two, 922–923 and 908–909, are immediately sequence-adjacent; this channel reflects backbone continuity between consecutive C-lobe segments rather than long-range allosteric coupling, and we do not interpret it as a dimerization-to-tail conduit.

Delete the dimerization interpretation (the "external dimerization machinery … at the cell surface" sentence) entirely: the αG–αH region is intracellular, the simulation is a monomer, and a superfamily-wide alignment containing one ErbB sequence cannot encode ErbB dimerization constraints.

---

## 6. The DCA/DCCM convergence claim (§2.7–2.8, Abstract, Figure 6)

**6.1 The published result reproduces.** Top-500 co-evolving distal pairs (|i−j| > 5) have mean |DCCM| = 0.287 against a sequence-separation null of 0.157 — **1.83×, p = 1.1 × 10⁻³³** (manuscript: 1.84×, 0.156, 4.56 × 10⁻³⁴).

**6.2 It does not survive a spatial control.** DCA couplings are enriched for three-dimensional contacts — this is the basis of contact prediction, and the manuscript cites Morcos et al. (2011) for exactly that. Contacting residues co-move for mechanical reasons. The top-500 co-evolving pairs have median Cα–Cα separation of **11.4 Å**; all other distal pairs, **24.3 Å**. Matching the null on Cα–Cα distance in 2 Å bins (coordinates from `egfr_clean.pdb`):

| Null model | Null mean \|DCCM\| | Enrichment | 95% CI | p |
|---|---|---|---|---|
| Sequence separation only (as published) | 0.157 | 1.83× | 1.67–2.01 | 1.1 × 10⁻³³ |
| Cα–Cα distance-matched | 0.281 | **1.02×** | 0.95–1.09 | **0.55** |

Partial Spearman correlation between coupling strength and |DCCM| controlling for Cα–Cα distance: **+0.028** (raw +0.060). The |i−j| > 5 filter removed backbone neighbours but not the spatial confound it was intended to address. The statement that the effect is "unlikely to arise from residual backbone proximity effects" is true of sequence proximity and silent about spatial proximity.

**6.3 The filter is inconsistently reported.** Every script uses |i−j| > 5 (`correlation_test.py`, `correlation_test_withneighbor.py`, `allostericcircuits.py`, `rank_couplings.py`). The manuscript states |i−j| > 10 in §2.7 and in the Figure 6D caption. No code uses 10.

**6.4 "Distance-matched" was not performed.** The abstract and Figure 6D describe a "distance-matched random baseline". A shared |i−j| threshold is not distribution matching. The term is now accurate only for the new analysis.

**6.5 `correlation_test.py` has no RNG seed**, so the published baseline and p-value are not reproducible as run.

**Replacement text for §2.8:**

> We compared the absolute dynamic correlation of the 500 most strongly co-evolving distal pairs (|i−j| > 5) against two null models. Against 500 random pairs satisfying the same sequence-separation filter, co-evolving pairs showed a mean |DCCM| of 0.287 versus 0.157, a 1.83-fold enrichment. However, co-evolving pairs are also far closer in space (median Cα–Cα 11.4 Å versus 24.3 Å), as expected given that direct couplings predict spatial contacts, and |DCCM| is itself strongly distance-dependent. When the null is matched on Cα–Cα distance in 2 Å bins, the enrichment falls to 1.02-fold (95% CI 0.95–1.09; p = 0.55), and the partial rank correlation between coupling strength and |DCCM| controlling for distance is +0.028. We therefore conclude that the apparent convergence between evolutionary coupling and dynamic correlation in this system is accounted for by three-dimensional proximity, and we do not interpret the molecular dynamics as independent support for the co-evolutionary predictions.

Remove from the abstract: "Molecular dynamics simulations provided **independent support** for these evolutionary predictions", the 1.84-fold enrichment, and the p-value. Remove from Figure 6D: "**strong evidence** that evolutionary constraints predict physical coordinated fluctuations".

---

## 7. L858R (§2.9, Figure 3A, community labels)

**7.1 Residue 853 is not the DFG motif.** Read from the alignment and confirmed against P00533: **I853**. The DFG motif is **D855-F856-G857**. The manuscript calls 853 "the DFG-switch" (§2.9), "DFG-motif (853)" (Figure 3A caption), and names Community 13 "Catalytic Switch (DFG)". The reported shortest path 858 → 857 → 856 → 855 → 854 → 853 *passes through* the actual DFG motif and terminates two residues beyond it at an isoleucine.

**7.2 The local-coupling-bias model is unsupported.** Across all 248 six-residue backbone windows in the domain (`backbone_window_cohesion.csv`), the 853–858 window sits at the **62nd percentile** by summed consecutive-bond coupling (14.71 vs median 13.66) and the **25th percentile** by summed coupling over all pairs in the window. The bond adjacent to L858 (857–858) is at the **34th percentile** — below median. The only elevated bond in the segment is D855–F856 (88th percentile), i.e. the conserved DFG motif itself. Residue 858 ranks **228 of 253** by node strength.

**7.3 There is a better result available.** L858 and the true DFG aspartate D855 are co-assigned to Community 12 in **93%** of 100 Louvain initialisations (858 and I853: 9%). This is a defensible module-level statement that the current analysis missed by targeting the wrong residue.

**7.4 The abstract/§2.9 tension.** The abstract calls the result "consistent with an allosteric mode of action" while §2.9 states the path "does not represent a long-range allosteric route". Both rest on different analyses, but the reversal reads as contradiction.

**Replacement text for §2.9, second paragraph:**

> The L858R activating mutation lies in Community 12, a nine-residue module spanning 854–862 that also contains the DFG aspartate D855; the two are co-assigned in 93% of 100 Louvain initialisations. Shortest-path analysis from residue 858 to D855 returns the sequence-adjacent backbone route, which reflects chain connectivity rather than a long-range allosteric channel. We tested whether this segment carries unusually strong co-evolutionary coupling by comparing it against all 248 six-residue backbone windows in the domain: the 853–858 window ranks at the 62nd percentile by summed consecutive-bond coupling and the 25th percentile by total within-window coupling, and residue 858 ranks 228th of 253 by node strength. The segment is therefore not distinguished by co-evolutionary cohesion, and we withdraw the local-coupling-bias model proposed in the preprint. The informative observation is the module-level one: L858 and the DFG motif belong to the same co-evolutionary community.

**Replacement text for the abstract's L858R clause:**

> …and placed the L858R activating mutation in the same co-evolutionary community as the DFG aspartate D855, while the minimal coupling route between them follows the backbone and does not constitute a long-range allosteric path.

---

## 8. Remaining corrections

**8.1 Threshold selection (Methods 3.3, §2.3, Figure 2).** Methods states the threshold maximizes Q subject to connectivity. Q is *maximal at the lowest threshold tested*: weighted Q = 0.8894 at 253 edges (21 communities, 50.20% connectivity) and falls monotonically. The selected threshold is the minimum density achieving full connectivity, nothing more. Remove "while maximizing Q" (Methods 3.3), "maximum modular resolution" (Figure 2), "mathematically optimal network scale" (Figure 5A) and "mathematically determined consensus maximum" (Figure 6).

**8.2 Figure 2 was computed with unweighted Louvain.** `generate_topology_data.py` calls `best_partition(G)` and `modularity(p, G)` without `weight='weight'`, while the production run (`1_8_egfrcommunity.py`) uses weighted Louvain. At 462 edges: weighted Q = 0.8036 (quoted in the text), unweighted Q = 0.7979 (plotted in Figure 2). Regenerate Figure 2 from `percolation_corrected.csv`.

**8.3 The percolation point is 458 edges, not 462.** At 457 edges the largest component holds 252 of 253 residues; at 458 it holds all 253. 462 = int(257 × 1.80) is four edges past the transition, because L was hard-coded as 257 while the network has 253 nodes. State the threshold as an edge count.

**8.4 The robustness claim is false (§2.3, lines 150–154).** Over 1.6L–2.0L the community count is **14, 15 or 16**, not stable. The top-5 betweenness residues are not stable either: at 1.6L and 1.7L the highest-betweenness residue is **849**, which does not appear in the published top five. Node strength *is* stable — 965, 717, 799, 851, 964 occupy the top five at every threshold tested. Report that asymmetry rather than a blanket stability claim.

**8.5 Residues 733 and 884 are not primary hubs.** §2.3 names them "the primary high-centrality residues" and §2.4 says 883 and 884 "rank among the primary infrastructure hubs". Actual node-strength ranks: **733 = 34/253, 884 = 98/253, 883 = 112/253**. Furthermore `fig8_perturbation_results.csv` — an analysis present in the data but absent from the manuscript — shows that removing 884 or 733 marginally *increases* global network efficiency (−0.01% and −0.02%), while 793 gives +2.26%, 858 +1.44% and 965 +1.00%. Either report this analysis or withdraw the hub claims.

**8.6 K745/E762 (§2.9).** The manuscript says K745 maps to Community 0 "consistent with its conserved salt-bridge interaction with Glu762". K745 is in Community 0 and E762 in Community 7; they are co-assigned in **0 of 100** Louvain runs. The partition separates the salt-bridge partners, so the observation does not support the claim. K745 also ranks 197/253 by node strength and 203/253 by betweenness.

**8.7 The RMSF gradient is 3.9-fold, not 3.4-fold.** §2.5 reports "3.4-fold … between the most dynamic and most rigid communities". Most dynamic C7 = 0.2135 nm; most rigid **C3 = 0.0545 nm** → **3.92×**. The 3.4× figure is C7/C14 (0.2135/0.0636 = 3.36), which is what Figure 4C's caption correctly states. Either say 3.9-fold or name Community 14 in the text.

**8.8 Residue 965 sits 5 residues from a truncated terminus.** The simulation spans 712–970 (`covar.log`: 259 Cα), so residue 965 is five residues from a free chain end and the C-terminal tail (971–1210) is absent entirely. The §2.4 caveat addresses crystallographic disorder, not the construct boundary. The §2.7 figure of 0.161 versus a 0.156 baseline is a 3% difference from a single unreplicated trajectory and should be reported as indistinguishable from noise or dropped. Reconcile the "αI-helix" label (Figure 3B) with the "disordered tail interface" description (§2.4).

**8.9 Connectivity at the lowest threshold.** §2.3 states "1.0L; 257 edges" and "50.2%". At 257 edges connectivity is 50.99%; the 50.2% figure corresponds to 253 edges. Pick one convention.

**8.10 1M17 is inhibitor-bound.** §2.3 and Methods 3.4 describe it as "the active EGFR crystal structure"; it is the erlotinib complex. State the bound ligand.

**8.11 References.** Two confirmed errors (`reference_audit.csv`):

- **Stamos et al. (2002)** is printed as "Crystal structure of the EGFR kinase domain in its active conformation. Cell 109, 711–719. doi:10.1016/S0092-8674(02)00792-5". The DOI **does not resolve**. Correct: Stamos J, Sliwkowski MX, Eigenbrot C (2002), *Structure of the epidermal growth factor receptor kinase domain alone and in complex with a 4-anilinoquinazoline inhibitor*, J Biol Chem 277:46265–46272, doi:10.1074/jbc.M207135200.
- **Zhang et al. (2006b)** is printed as "Structure of the epidermal growth factor receptor kinase domain in an inactive conformation. Nature 439, 379–386. doi:10.1038/nature04396". No such paper exists, and that DOI resolves to **McMahon & Gallop (2005), *Membrane curvature and mechanisms of dynamic cell membrane remodelling*, Nature 438:590–596** — an unrelated membrane-biology review. This reference is cited in the Introduction for the αC-helix/activation-loop/catalytic-loop switch. Replace with Zhang et al. (2006a) Cell 125:1137–1149, already in the list, or with Wood et al. (2004) Cancer Res 64:6652–6659 for the inactive conformation.
- Sutto et al. (2015) is missing its DOI (10.1073/pnas.1508584112).
- Yun et al. (2008) PNAS 105:2070–2075 (doi:10.1073/pnas.0709662105) must be added for the corrected T790M mechanism.

27 of 36 entries verified clean against CrossRef. Given that one reference is fabricated with a DOI belonging to a different paper, **the full bibliography should be re-verified by DOI before resubmission.**

**8.12 Title.** No conformational switch was simulated, observed or perturbed: one 100 ns trajectory from the active state, n = 1, monomeric, no inactive-state or mutant simulation. The paper's own closing sentence claims only a "correspondence". Suggested:

> *Co-evolutionary community architecture of the EGFR kinase domain and its relationship to conformational dynamics*

and strike "coordinate" from the hypothesis statement.

**8.13 Figure captions.** Bring every caption to the hedging level of the corresponding Results text. Specifically remove "mathematically optimal", "mathematically determined", "validated allosteric blueprint", and "strong evidence that evolutionary constraints predict physical coordinated fluctuations".

---

## 9. Items requiring compute beyond this workspace

Two items required resources unavailable in the audit workspace (8 GB RAM, no GPU, no compute provider). **9.1 has since been executed on HiPerGator and is addressed below; 9.2 remains open.**

**9.1 Replicate molecular dynamics. — ADDRESSED.** The original finding: the DCCM, the RMSF hierarchy and the enrichment analysis all rested on a single 100 ns trajectory (`covar.log`: 50,001 frames, 259 Cα, one replicate); §2.5 conceded that multi-replicate validation was needed, but that caveat was never repeated in §2.7 or §2.8, which depend on the same trajectory.

Four further 100 ns trajectories (replicates 2–5) were run on HiPerGator, branched from the same equilibrated NPT configuration — identical `npt.gro`, `npt.cpt` and `topol.top` by SHA256 — with fresh velocities from seeds 20250002–20250005 and production parameters identical to the published run. With replicate 1 this gives five replicates and 500 ns aggregate, at 400–418 ns/day per replicate. Results:

| Quantity | Result across five replicates |
|---|---|
| DCCM agreement (Pearson r, strict upper triangle, 253×253) | 0.746–0.899, median 0.840; replicate 1 vs others 0.798–0.888 |
| Enrichment, sequence-separation null | 1.45–1.88 (mean 1.68, SD 0.19); significant in every replicate |
| Enrichment, Cα–Cα distance-matched null | 0.95–1.02 (mean 0.98, SD 0.03); every bootstrap CI includes 1.0 |
| Replicate-mean DCCM | 1.73-fold and 0.97-fold under the two nulls |
| RMSF agreement (Pearson / Spearman) | 0.821–0.954 / 0.873–0.957 |
| Mechanical hierarchy | Community 7 most mobile and Community 3 most rigid in 5/5; C7/C3 gradient 2.8–4.0-fold |

The headline result is unchanged: the distance-matched enrichment is ~1.0 in every replicate, as expected, because the confound is geometric (median Cα–Cα 11.4 Å versus 24.3 Å) rather than stochastic. Two qualifications that the replicates introduce rather than remove. First, the *magnitude* of the RMSF gradient is replicate-dependent (2.8–4.0-fold against the published 3.9), and the rank order of individual RMSF peaks is not stable — residues 751–753 top the profile only in replicate 1 — so §2.5 and the Figure 4 caption now interpret the ordering of communities and the identity of the flexible regions, not individual peak residues. Second, replicate agreement establishes precision, not accuracy.

Reported in manuscript §2.9, Supplementary Figure S1 and Supplementary Table S1. Per-replicate matrices are in `data/processed/egfr_dccm_matrix_rep{2,3,4,5}.csv` with the replicate mean in `egfr_dccm_matrix_mean5.csv`; the replicate correlations and per-replicate enrichments are covered by `make verify`. Item 8.8's "single unreplicated trajectory" framing for residue 965 is superseded accordingly.

**9.2 RTK-only and ErbB-only alignments.** The decisive response to the scope problem. **Requirement:** subset PF07714 to metazoan sequences whose parent UniProt entry carries a transmembrane annotation, re-run plmc with identical parameters (8 cores, 64 GB, ~49 min for the current alignment), and repeat percolation, Louvain and the mutation mapping. Report Neff/L for each subset. The ErbB-restricted alignment will almost certainly be underpowered — the full family contains 17 ErbB sequences — and reporting its Neff/L is itself the justification for a superfamily-wide design. Exclude ERBB3 or report both ways, since it is a pseudokinase under relaxed selection at the catalytic residues.

They were listed in priority order: 9.1 was cheap and repaired an acknowledged weakness, and has now been done; 9.2 is the substantive new experiment and remains outstanding.

---

## 10. Supporting files

| File | Contents |
|---|---|
| `provenance_audit.md` | True pipeline chain, abandoned branches, five production-script bugs, confirmed-correct items |
| `corrected_statistics.csv` | Manuscript value vs verified value vs source file, 21 quantities |
| `alignment_composition.csv` | Taxonomic composition, raw vs retained, with enrichment factors |
| `fig_alignment_composition.png` | Composition, clustering bias, and the vertebrate/human/ErbB collapse |
| `percolation_corrected.csv` | Weighted and unweighted Q, connectivity and community count, 1.0L–4.0L |
| `fig_percolation_corrected.png` | Corrected percolation curve, the weighted/unweighted discrepancy, community instability |
| `intermodule_flux_corrected.csv` | All 105 pathways under four metrics with ranks |
| `fig_flux_matrix_corrected.png` | Score distribution, size artefact, and metric disagreement |
| `community_identity_corrected.csv` | All 15 communities with spans and observed secondary structure |
| `mutation_mapping_corrected.csv` | Per-residue community, node strength, betweenness and ranks |
| `fig_community_landmarks.png` | Linear community map, betweenness skyline, co-membership stability |
| `dccm_enrichment_corrected.csv` | Enrichment under both null models with CIs |
| `fig_enrichment_distance_matched.png` | Distance confound and the corrected enrichment |
| `backbone_window_cohesion.csv` | All 248 six-residue windows with percentiles |
| `fig_backbone_windows.png` | The 853–858 window against the domain-wide distribution |
| `reference_audit.csv` | Curated reference verdicts with corrected citations |
| `reference_audit_crossref_raw.csv` | Raw CrossRef match output for all 36 entries |

