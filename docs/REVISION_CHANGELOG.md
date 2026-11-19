# Revision changelog

Corrections from `CORRECTIONS.md` applied to `frontiers.tex`, `test.bib` and the figure set. Items §1–§8 applied in full; §9 (replicate MD and RTK-only DCA) deliberately not applied — it requires new compute — but is now stated as a limitation rather than left implicit.

**Validation after editing:** 37 bibliography entries, 37 distinct keys cited, no citation missing from the bibliography and no orphan entries; no `\ref` without a matching `\label`; braces balanced; `$` delimiters even; all six figure files resolve. No TeX installation was available in this environment, so the document was **not compiled** — please run `latexmk -pdf frontiers.tex` before submitting.

---

## 1. Title and framing

| Was | Now |
|---|---|
| *Evolutionary Hubs Coordinate the Conformational Switch of the Epidermal Growth Factor Receptor* | *Co-evolutionary Community Architecture of the EGFR Kinase Domain and its Relationship to Conformational Dynamics* |
| Running head: *Evolutionary Hubs in the EGFR Kinase Domain* | *Co-evolutionary Communities in the EGFR Kinase Domain* |

The hypothesis statement no longer asserts that hubs "coordinate the long-range allosteric transitions"; it now asks whether high-centrality residues correspond to structural elements and whether any dynamic correspondence survives a spatial control.

## 2. Alignment statistics (§2.1, Methods 3.1, Figure 1C)

Every reported alignment number replaced with the value verified from the pipeline that actually produced the couplings:

| Quantity | Was | Now |
|---|---|---|
| Raw family size | 1,051,876 | **400,321** |
| Alignment width | 3,667 columns | **3,839 columns** |
| Redundancy reduction | CD-HIT 80% then 62% | **single pass, `-c 0.62 -n 4`** (CD-HIT v4.7) |
| Final depth | 93,773 sequences | **32,912 sequences** (32,672 valid) |
| Focus positions | 257 | **253 match-state**; residues 712, 713, 967, 968 unscored; domain analysed is 714–966 |
| Depth justification | raw N/L ≈ 370 | **N_eff = 27,613; N_eff/L = 109** |

Added: the `-t 0.2` reweighting was effectively inert after 62% clustering (mean 80%-identity neighbourhood = 1.000), and PF07714 is identified under its current designation, *Protein tyrosine and serine/threonine kinase*.

## 3. Family scope (Abstract, Introduction, §2.1, Figure 1C)

Removed all claims of receptor-tyrosine-kinase or ErbB specificity. Added the measured composition of the retained alignment: vertebrate sequences fall from 35.8% of the raw family to 3.5% after clustering, leaving two human sequences and one ErbB-family sequence — the manually recovered EGFR reference itself. The Introduction now states explicitly that the couplings describe superfamily-conserved constraints, and points to the new Limitations section.

## 4. T790M (Abstract, Introduction, §2.9)

- Removed "placed the T790M resistance mutation at the network's primary information bottleneck".
- T790 is no longer described as acting "at a distance from the ATP-binding cleft"; it is the gatekeeper lining the cleft.
- §2.9 now reports T790's own centrality — betweenness 0.117 (rank 20/253), node strength 11.68 (rank 61/253) — and states that L792 and M793 fall in Communities 0 and 13, not Community 2 with T790, co-assigning in 8% of 100 Louvain initializations.
- The conventional mechanism is corrected to restoration of ATP affinity, citing **Yun et al. (2008)**, newly added.

## 5. Inter-module coupling (§2.6, Methods 3.4, Figure 5)

§2.6 rewritten and retitled *Inter-Module Coupling*:

- The coupling score is disclosed as average-product-corrected and therefore signed (mean 0.0004; 55.5% negative).
- Σ J_ij over all |A|×|B| pairs is no longer used to rank pathways: 75 of 105 module pairs return a negative value and Φ correlates negatively with pair count (r = −0.33).
- Coupling is now reported over the 462 edges defining the partition. Only 63 edges cross community boundaries, across 27 of 105 module pairs. Strongest: C3–C6 (8 edges, ΣJ = 18.01) and C3–C14 (7 edges, ΣJ = 12.95).
- "Top 0.1% of all inter-community pathway strengths" removed.
- The P-loop/αC-helix pair is disclosed as **one edge** of 462, ranking 22nd of 105, and is no longer presented as positive methodological validation.
- The dimer-interface/tail-anchor interpretation is withdrawn. Community 5 is identified as residues 909–922 (the αG–αH loop, not the αG helix); Community 10 flanks it in sequence; the three joining edges are 922–923 (|i−j| = 1), 908–909 (|i−j| = 1) and 918–926. The canonical activator-side dimer residues V948, M952, I953 are noted to lie in Communities 3, 11 and 11.
- The "external dimerization machinery … at the cell surface" sentence is deleted.
- Methods 3.4 now defines all three metrics actually computed and states that module pairs are keyed on the unordered pair (105, not 126).

## 6. DCA/DCCM convergence (§2.7, §2.8, Methods 3.6, Abstract, Figure 6)

§2.8 rewritten:

- The published result is retained and attributed: 1.83-fold against a sequence-separation null, p = 1.1 × 10⁻³³.
- Added the Cα–Cα distance-matched null: **1.02-fold (95% CI 0.95–1.09; p = 0.55)**, Mann–Whitney p = 0.62, with the median separations (11.4 Å vs 24.3 Å) that explain it, and the partial rank correlation controlling for distance (+0.028 vs +0.060 raw).
- Removed "independent support", "strong evidence", and the claim that the effect is "unlikely to arise from residual backbone proximity effects".
- Sequence-separation filter harmonized to |i−j| > 5 throughout; "distance-matched" now describes a null that is actually distance-matched.
- Methods 3.6 documents both nulls, the bootstrap, and a fixed RNG seed.

**A result changed as a consequence.** The consensus-maximum pair is filter-dependent. Under |i−j| > 10 the maximum is 799–841 (2.33) but leads the runner-up 728–741 (2.31) by **0.7%**; under the |i−j| > 5 filter used elsewhere the maximum is **717–725** (3.66), 49% clear. The previously reported consensus score of 0.6281 is not reproducible from the deposited data under any normalization tried. §2.7 now reports both filters and designates no single convergence locus.

## 7. L858R (§2.9, Abstract)

- Residue 853 is **isoleucine**; the DFG motif is D855–F856–G857. All references to 853 as a DFG position corrected.
- The local-coupling-bias model is withdrawn, with the quantification that retired it: the 853–858 window ranks at the 62nd percentile of 248 six-residue backbone windows by consecutive-bond coupling and the 25th by total within-window coupling; the 857–858 bond is at the 34th percentile; only D855–F856 is elevated (88th). Residue 858 ranks 228/253 by node strength.
- Replaced with the defensible result: L858 and D855 co-assign in **93%** of 100 Louvain initializations.
- The abstract no longer says "consistent with an allosteric mode of action".

## 8. Remaining corrections

- **8.1/8.3** Threshold selection no longer claims to maximize Q. Q is maximal at the lowest density tested (0.8894, 253 edges, 21 communities) and falls monotonically. Full connectivity is first reached at **458** edges; 462 is four edges past it.
- **8.2** Methods 3.3 now specifies weighted Louvain for the percolation sweep, matching the production partition.
- **8.4** The robustness claim is replaced with the measured asymmetry: node-strength hubs (965, 717, 799, 851, 964) are stable across 404–506 edges; community count varies 14–16; the top betweenness residue is 849 at 404 and 430 edges.
- **8.5** The hub claim for 883/884 is withdrawn (ranks 112 and 98 of 253), and the previously unreported knockout analysis is now included, including the negative efficiency changes for 884 and 733.
- **8.6** K745/E762 corrected: they are in different communities and co-assign in 0 of 100 runs, so the partition separates the salt-bridge pair.
- **8.7** RMSF gradient corrected to **3.9-fold** (C7 vs C3), with 3.4-fold identified as the C7:C14 ratio.
- **8.8** The simulated construct ends at residue 970; residue 965 lies five residues from a free terminus and the tail (971–1210) is absent. The 0.161-vs-0.156 result is now reported as indistinguishable from noise.
- **8.9/8.10** Connectivity at the lowest threshold stated consistently (253 edges, 50.2%); 1M17 described as erlotinib-bound.
- **8.13** All six captions brought to the hedging level of the Results text. Removed: "mathematically optimal", "mathematically determined consensus maximum", "validated allosteric blueprint", "maximum modular resolution", "strong evidence that evolutionary constraints predict physical coordinated fluctuations".

## Bibliography (§8.11)

| Entry | Action |
|---|---|
| `Stamos2002` | **Corrected.** Printed DOI `10.1016/S0092-8674(02)00792-5` does not resolve, and the title, journal, volume and pages were all wrong. Replaced with Stamos, Sliwkowski & Eigenbrot (2002), *Structure of the epidermal growth factor receptor kinase domain alone and in complex with a 4-anilinoquinazoline inhibitor*, J Biol Chem 277:46265–46272, doi:10.1074/jbc.M207135200. |
| `Zhang2006a` | **Removed as hallucinated.** No paper matches "Structure of the epidermal growth factor receptor kinase domain in an inactive conformation, Nature 439:379–386"; the printed DOI `10.1038/nature04396` resolves to McMahon & Gallop (2005), *Membrane curvature and mechanisms of dynamic cell membrane remodelling*, Nature 438:590–596. Its single citation — the αC-helix/activation-loop/catalytic-loop switch in the Introduction — now points to the genuine Zhang et al. (2006) *Cell* paper, which reports both conformations. |
| `Zhang2006b` → `Zhang2006` | Renamed; the a/b suffix is no longer needed. |
| `Yun2008` | **Added.** Yun et al. (2008) PNAS 105:2070–2075, doi:10.1073/pnas.0709662105, for the corrected T790M mechanism. |
| `JamesVerkhivker2014` | Full title restored (the subtitle had been dropped). |
| 19 entries | **DOIs added** after resolution against CrossRef: Astl & Verkhivker 2019, Behn 2025, Bowman 2024, Dixit & Verkhivker 2011, Ekambaram 2025, Ekeberg 2013, Hagberg 2008, Halabi 2009, Hilser 2012, James & Verkhivker 2014, Liu & Nussinov 2016, Jura 2009, Morcos 2011, Nussinov 2022, Nussinov & Tsai 2013, Sutto 2015, Xu 2025, Zhang 2006, Zhang 2016. |

**Every one of the 37 entries now carries a DOI that resolves to the cited paper.**

---

## Figures

All six figures were regenerated. The PyMOL renders (1A, 3C, 4B, 5A) carry no numerical
annotations and were reused unchanged; every panel with text or numbers baked into it was
redrawn from the deposited data.

| Figure | Was | Now |
|---|---|---|
| 1, panel B | DFG drawn at 853–855; αG at 905–920 | DFG at 855–857; αG at 893–908 (UniProt P00533 observed helix) |
| 1, panel C box I | "93,773 Tyrosine Kinases / Filtered at 62% Identity" | 32,912 sequences, PF07714, *N*<sub>eff</sub>/*L* = 109 |
| 1, panel C box III | "Optimal 1.80L Resolution" | 462-edge network, *Q* = 0.8036, full connectivity from 458 edges |
| 1, panel C box V | "Evolutionary-Dynamic Convergence" | enrichment absent under a distance-matched null (1.02×, *p* ≈ 0.5); no conformational transition simulated |
| 2 | unweighted *Q* plotted against a weighted annotation off the curve | weighted Louvain throughout; 458-edge transition marked; "mathematically optimal" dropped; monotonic decline of *Q* from its maximum at 253 edges annotated |
| 3, panel A | "Res 853 (DFG Motif)" | "Res 853 (Ile) — DFG is 855–857" |
| 3, side box | "Primary Bridge: Residue 792 … global bottleneck" | betweenness ranking stated as threshold-dependent (849 leads at 404–430 edges); node-strength top five stable across 404–506 |
| 4, panel C | "1.8x Rigidity Gradient" | 3.9× (C7 vs C3) and 3.4× (C7 vs C14), with the single-trajectory caveat on the panel |
| 5, panel B | "ALLOSTERIC INTEGRATION HUB" | hub label dropped; corrected module names; 63 of 462 edges inter-community stated on the panel |
| 5, panel C | signed flux matrix, values to −125 | coupling over the 462 network edges |
| 6, panel D | two bars, both means and both error bars written in as literals | three bars computed from the data, with both nulls and bootstrap intervals |

### Corrected scripts

The original scripts are in `EGFR_Final_Results/`. Corrected replacements are in
`figure_scripts/`, with a mapping table in `figure_scripts/README.md`. Each prints the
quantities it draws, so figure and manuscript can be checked against each other.

Defects fixed, beyond the figure text above:

- **`makefig2.py` (Figure 4)** selected the scaffold community with
  `summary_df['name'].str.contains('Scaffold')`, which matches "Scaffold-Switch Interface"
  (C8, 0.1209 nm) before "Main Engine Block (Scaffold)" (C14, 0.0636 nm) and returned the
  first. 0.2135/0.1209 = 1.766, printed as "1.8x". The replacement keys on community ID.
- **`makefig3.py` (Figure 5)** grouped on `('Name_i','Name_j')`, splitting each unordered
  community pair across two rows — 126 rows for 105 pairs — so the plotted values were
  partial sums. The replacement keys on the sorted pair.
- **`fig7spoof.py` (Figure 6)** hardcoded `means = [0.156, 0.287]` and
  `errors = [0.004, 0.007]`. The error bars were not derived from the data; the true
  standard errors are 0.005 and 0.009, and `makefig7.py` used a third, different pair
  (`[0.005, 0.008]`, commented "Standard error for a high-confidence plot"). The
  replacement computes both. Its heatmap comment,
  "Increased VMAX to 0.9 makes the 799/841 peak bright red", documents a colour scale
  chosen to emphasise a specific cell; the replacement leaves the scale untuned.
- **`backboneremove.py`** builds a graph with backbone edges excluded, then prints
  `path_remote = [910, 814, 813, 849, 817, 851, 820, 773, 823, 965]` as the
  "DIMER-TO-TAIL CONDUIT". All nine steps are genuine edges of the full 462-edge network,
  so the list is not fabricated — but step 814→813 has |i−j| = 1, so it is a backbone edge
  and the path cannot come from the graph the script constructs. It is also non-minimal in
  both graphs: the shortest 910→965 route is 6 steps in the full network and 8 steps with
  backbone edges removed. `compute_pathways.py` computes the path instead of printing one.
  The script's L858R result is correct and was honestly reported: residue 858's only edges
  in the 462-edge network are to 857 and 859, so removing backbone edges removes the node
  and no non-backbone path exists.
- **Edge count.** `int(257 * 1.80) = 462` uses *L* = 257, the number of focus columns
  requested, rather than *L* = 253, the number of match-state positions PLMC actually
  scored. The 462-edge network is retained for continuity with the deposited partition,
  and the percolation threshold is now reported as the directly located value (458 edges)
  rather than as a multiple of *L*.

All numbers appearing in the regenerated figures were recomputed from the primary data and
checked against the manuscript; the comparison is in `figure_number_audit.csv` (39 of 39
agree). The enrichment analysis is a resampling procedure, so the replacement script
averages over 200 draws: the ratios (1.83×, 1.02×) are stable across seeds, while *p* moves
between 0.50 and 0.56, which is why it is quoted to one decimal and the confidence interval
carries the inference.

---

## Not applied — §9

Deferred at your instruction; now disclosed in the new **Limitations** section (`\label{sec:limitations}`) rather than left unstated.

1. **Replicate molecular dynamics.** 3–5 × 100 ns from independent velocity seeds; `md.mdp`, `npt.gro`, `npt.cpt` and `topol.top` are all present in the project directory. Roughly 3–5 GPU-days. The distance-matched null should be applied per replicate.
2. **RTK-only and ErbB-only alignments.** Subset PF07714 to metazoan sequences with a transmembrane annotation, re-run PLMC with identical parameters (~49 min on 8 cores for the current alignment), and repeat percolation, Louvain and mutation mapping. Report N_eff/L for each subset; the ErbB-restricted set will likely be underpowered, which is itself the justification for a superfamily-wide design. Consider excluding ERBB3, a pseudokinase under relaxed selection at the catalytic residues.

---

## Replicate molecular dynamics (CORRECTIONS §9.1, added after the first revision round)

§9 is no longer deferred in full: four additional 100 ns trajectories were run on HiPerGator and
replicate MD is now reported in the manuscript. §9.2 (RTK-only and ErbB-only alignments) remains
the one outstanding item.

| Location | Was | Now |
|---|---|---|
| Abstract | "A 100 ns molecular dynamics simulation" | five independent 100 ns simulations, 500 ns aggregate; per-replicate ranges added for both nulls |
| §2.5 | "require multi-replicate validation before they can be considered statistically generalizable" | ordering reproduced in 5/5 replicates; magnitude reported as the replicate-dependent quantity (2.8–4.0-fold) |
| §2.7 | "a single unreplicated trajectory" | compared against the measured between-trajectory variation (0.280–0.422) |
| §2.8 | within-trajectory bootstrap only | per-replicate range 0.95–1.02 added, every CI spanning unity |
| §2.9 | — | **new subsection**: DCCM agreement, per-replicate enrichment, between-trajectory variability, RMSF reproducibility, and the two bounds (precision not accuracy; representativeness not existence) |
| §2.13 Limitations | "no replicate was run" | five replicates; 500 ns still short; agreement is precision, not accuracy |
| Methods 3.5 | single production run | five replicates, seeds, re-equilibration protocol, and the site-module vs conda-forge build difference stated plainly |
| Figure 4 caption | "require multi-replicate validation" | ordering reproduced in all five; individual RMSF peak ranks noted as trajectory-specific |
| Figure 6 caption | replicate-1 values only | replicate ranges added; **replicate 1 retained as the primary panel** |
| Supplementary | — | Figure S1 (per-replicate enrichment forest plot; 5×5 DCCM correlation matrix) and Table S1 (seed, build, ns/day, r vs rep 1, both enrichments) |

**Figure 6 decision.** Replicate 1 was kept as the primary panel rather than regenerating the
figure from the five-replicate mean. Its correlations to the other four (r = 0.798–0.888, median
0.840) lie inside the spread among replicates 2–5 (0.746–0.899, median 0.840), so the published
matrix is not an outlier, and the replicate-mean matrix reproduces both headline values
(1.73-fold and 0.97-fold). Keeping replicate 1 also keeps Figure 6 anchored to the same
trajectory as every other dynamic analysis in the paper.

**Distance-matched CI.** The upper bound is 1.12, not 1.11. `analysis/make_figure6.py` prints
`0.94-1.12` and Figure 6D is annotated with it; the bootstrap value is 1.1160. The 1.11 in
`results/tables/dccm_enrichment_corrected.csv` predated the repository pipeline and has been
corrected to agree with the figure and the text.

**Two citations added:** Bussi et al. (2007) for the velocity-rescaling thermostat and Bernetti &
Bussi (2020) for the C-rescale barostat used in the replicate re-equilibration.

**Validation after editing:** no `\ref` without a matching `\label`; every cited key present in
`test.bib`; all nine figure files resolve; `figure`/`table`/`tabular` environments balanced;
`make verify` passes 41/41 checks, 18 of them new and covering the replicate correlations,
per-replicate enrichments under both nulls, RMSF agreement and the community hierarchy. No TeX
installation was available in this environment, so the document was **not compiled** — run
`latexmk -pdf frontiers.tex` before submitting.
