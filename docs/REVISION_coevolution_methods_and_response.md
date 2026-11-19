# Beyond-DCA coevolution analyses for the EGFR allostery revision

Prepared for the preprint revision. All numbers below were computed in this session from
`egfr.plmc`, `kinase_final_aligned.fasta`, `egfr_dccm_matrix.csv` and 328 EGFR structures
retrieved from the PDB. Supporting tables are listed at the end.

---

## 1. Summary of what was decided and why

A literature assessment (37 CrossRef-verified references, see `method_assessment.md`) ranked six
candidate methods. Three were run; one was run as a supplement; two were skipped on
statistical grounds.

| Method | Verdict | Outcome here |
|---|---|---|
| Contact-masked DCA vs. an ErbB conformational ensemble | **Run** | Decisive — see §3.1, §3.2 |
| pySCA sector decomposition | **Run** | Positive result — see §3.3 |
| Site-anchored coupling profiles | **Run** | Null — see §3.4 |
| Hopfield-Potts spectral decomposition | Supplement | Uninformative — see §3.5 |
| ErbB-only subfamily refit | Skipped | Four paralogs cannot reach the effective-sequence floor |
| Restricted Boltzmann machines | Skipped | No pairwise read-out, no kinase precedent, no agreed null |

The governing prior from the literature is Anishchenko et al. (2017, PNAS,
10.1073/pnas.1702664114), who found little evidence for direct coupling between spatially
separated allosteric and functional sites across ~4,000 families, and reported that 91% of
coevolving pairs at 5–15 Å are in contact in at least one homologous structure. The analyses
below reproduce that conclusion in EGFR specifically.

---

## 2. Methods text (manuscript-ready)

### 2.1 Reconstruction of the coupling model

The plmc parameter file was parsed directly to recover the full coupling tensor
*J*<sub>ij</sub>(*a*,*b*) for the L = 253 modelled positions (EGFR 714–966; plmc drops the four
lowercase insertion columns corresponding to 712, 713, 967 and 968). Frobenius norms in the
zero-sum gauge followed by average-product correction reproduce the deposited coupling scores at
Pearson r = 0.9999997 (Spearman ρ = 0.9999995; maximum absolute deviation 0.018 on scores
spanning −1.78 to 6.04), confirming that the published score is an APC-corrected Frobenius norm
taken over the full 21-letter alphabet. Because inclusion of the gap state is a known confound,
all analyses were repeated with the norm restricted to the 20 amino-acid states.

### 2.2 Conformational contact mask

All EGFR entries in the PDB were retrieved by UniProt accession (P00533) and by sequence
similarity to the kinase domain (≥90% identity), giving 334 entries. Each polymer chain was
mapped onto EGFR 714–966 by global pairwise alignment with end-gap-free scoring (BLOSUM62,
gap open −11, extend −1) and retained when it aligned ≥150 positions at ≥85% identity. This
mapping step is essential: 39 chains, including 1M14, 1M17, 2GS2 and 2GS7, are numbered
according to the mature protein and are offset by +24 relative to the canonical sequence.
537 chains from 328 entries passed. For every residue pair we recorded the minimum heavy-atom
distance across all chains (capped at 15 Å) and the minimum Cα–Cα distance (uncapped), so that a
pair in contact in *any* solved conformational state can be excluded. Inter-chain distances were
computed separately from 50 biological assemblies containing two or more kinase-domain chains,
including the 2GS6 asymmetric activating dimer.

### 2.3 Matched-null enrichment testing

Coupling–dynamics concordance was assessed by comparing the mean absolute dynamic
cross-correlation of a coupling set against a null matched on both sequence separation
(bins 6–10, 11–20, 21–40, 41–80, 81–160, >160) and minimum Cα–Cα distance (2 Å bins). For each
test pair, up to 20 control pairs were drawn from the same joint bin and their mean taken,
giving a paired design assessed by Wilcoxon signed-rank test. Because residue pairs drawn from
253 positions share endpoints and are therefore not independent, every headline difference was
additionally assessed by a residue-level cluster bootstrap: positions were resampled with
replacement and each pair weighted by the product of its endpoints' multiplicities
(4,000 resamples). Pair-level p-values are reported for comparability with the original analysis,
but the residue-level intervals are the ones we interpret.

### 2.4 Statistical coupling analysis

SCA was run with pySCA 6.1. The 253 focus columns were extracted from the 32,912-sequence
alignment and sequences with more than 20% gaps over those columns were removed, leaving 17,150
sequences (effective sequence count 16,985 at an 80% identity threshold). The SCA positional
correlation matrix was computed with Frobenius-norm dimension reduction and regularisation
λ = 0.03. Significance of eigenmodes was assessed against 10 randomised alignments preserving
per-position amino-acid frequencies; six eigenmodes exceeded the mean second eigenvalue of the
randomised trials (5.32). Sectors were defined by independent component analysis of the six
significant eigenvectors with a 95% cutoff on the fitted t-distribution of each component.
Because single-sector proteins can have their SCA signal dominated by conservation, two controls
were computed: the rank correlation between positional relative entropy and maximum independent
component loading, and a conservation-matched permutation null in which random position sets were
drawn to reproduce the sector's distribution over relative-entropy deciles (2,000 draws).

### 2.5 Site-anchored coupling profiles

Binding sites were defined from bound ligands rather than by hand. The orthosteric site is the
union of residues within 4.5 Å of erlotinib in 1M17 and of AMP-PNP in 5D41 and 6DUK (29 residues).
The allosteric site is the union of residues within 4.5 Å of EAI045 in 5D41 and
JBJ-04-125-02 in 6DUK, with residues shared with the orthosteric site removed (16 residues).
For each position, coupling to the orthosteric site was summarised as the mean of its five
largest couplings to site residues more than 8 Å away, so that the score cannot be produced by a
direct contact. Allosteric-pocket residues were compared against non-site residues matched on
minimum Cα–Cα distance to the orthosteric site.

### 2.6 Spectral decomposition of the coupling tensor (supplementary)

The coupling tensor was assembled into a 5,060 × 5,060 symmetric matrix in the zero-sum gauge
with the gap state removed and self-couplings set to zero, and diagonalised. For each eigenmode
we recorded the eigenvalue, the participation ratio over positions, and the enrichment of its
pairwise contribution on contacts relative to non-contacts. This is a descriptive decomposition
of an already-fitted plmDCA model and is **not** a Hopfield-Potts fit in the sense of
Cocco et al. (2013) or Shimagaki & Weigt (2019); the λ<sub>J</sub> = 1.0 regularisation used in the
original fit preferentially shrinks small couplings and therefore distorts the low-eigenvalue tail
those papers identify as contact-like.

---

## 3. Results

### 3.1 The published "long-range" set is dominated by contacts

Of the top 500 couplings at |i−j| > 5, **43.0% are within 5 Å in at least one of the 328 solved
EGFR structures** and 59.0% are within 8 Å; the median minimum heavy-atom distance is 6.17 Å.
Only 28.6% are further than 12 Å in every solved state. Restating the census at the CASP
long-range convention (|i−j| ≥ 24) improves this but does not fix it: 32.0% are still contacts.

Of the four pairs highlighted in the current allosteric census, three — 717–725 (|i−j| = 8),
846–852 (6) and 728–741 (13) — are short-range by the CASP convention, and only 799–841 (42)
qualifies as long-range.

### 3.2 The evolution–dynamics concordance does not survive a matched null

The published comparison reproduces exactly: mean |DCCM| of 0.287 for the top-500 set versus
0.154–0.159 for a uniformly drawn null, a 1.80–1.87× enrichment at p ≈ 10⁻²⁵. However, the top-500
set has median sequence separation 32 and median minimum Cα–Cα distance 9.9 Å, against 21.2 Å for
all pairs at |i−j| > 5, so the null is not comparable.

| Comparison | Mean \|DCCM\| | Enrichment | Pair-level p | Residue-level 95% CI on the paired difference |
|---|---|---|---|---|
| Top-500 vs. uniform null (as published) | 0.287 vs. 0.158 | 1.80× | 1.7 × 10⁻²⁵ | — (unmatched) |
| Top-500 vs. separation + distance matched null | 0.287 vs. 0.242 | **1.17–1.19×** | 6.2 × 10⁻⁸ | **[+0.012, +0.076]** — excludes zero |
| Top-500 distal only (>8 Å, all states) vs. matched null | 0.146 vs. 0.153 | **0.95×** | 2.1 × 10⁻³ | **[−0.025, +0.014]** — includes zero |
| Top-500 at >12 Å vs. matched null | 0.133 vs. 0.140 | 0.95× | 6.8 × 10⁻⁴ | includes zero |

A modest residual enrichment of about 1.2× survives both the matched null and residue-level
resampling, and is carried by pairs in physical contact. For couplings that are distal in every
solved EGFR structure, **there is no detectable concordance with the molecular dynamics in either
direction**: the pair-level test appears to put the distal set significantly *below* the matched
null (p = 2.1 × 10⁻³), but that significance is an artefact of treating 500 pairs drawn from 203
distinct residues as independent observations — the residue-level interval comfortably includes
zero. The correct statement is a null result, not a reversal. Ranking by SCA instead of DCA does
not change this (0.90× on distal pairs), nor does excluding the gap state (0.95×).

### 3.3 SCA sectors: the one positive result

Six significant eigenmodes were recovered (λ₁ = 34.9, λ₂ = 13.8, versus 5.32 for the randomised
second eigenvalue), decomposing into six sectors covering 101 of 253 positions.

**Sector 2 (22 positions)** is the catalytic/regulatory spine: 718, 719, 721, 723, 724, 726, 727,
729, 742, 743, **745** (catalytic lysine), **762** (αC glutamate), **835–837** (HRD motif), 842,
844, **855–857** (DFG motif), **858** (the oncogenic L858R position) and 881. It spans a median
sequence separation of 82 residues, is fully contiguous in the 8 Å contact graph, and has a mean
pairwise Cα–Cα distance of 10.9 Å against 18.0 Å for conservation-matched random position sets
(**p < 0.0005**, 0 of 2,000 draws). Sectors 5 and 6 are also significantly compact
(p < 0.0005 and p = 0.004); sectors 1 and 4 are not.

The compactness is not an artefact of the sector containing contiguous stretches of sequence.
Sector 2 comprises 13 runs (lengths 4, 3, 2, 2, 2, 2 and seven singletons), and against a null
that reproduces exactly that run-length structure at random positions the mean pairwise distance
is 20.4 ± 2.0 Å versus the observed 10.9 Å (**p < 0.0005**). Sector 2 is therefore spatially
clustered beyond what either its conservation profile or its sequence-block structure predicts.

**The necessary caveat.** Sector membership is strongly predicted by conservation alone
(Spearman ρ = 0.832 between positional relative entropy and maximum independent component
loading; sector positions have mean relative entropy 1.507 versus 0.492 for unassigned
positions). The defensible claim is therefore about *spatial organisation*, not about
membership: the sector's three-dimensional compactness is not reproduced by conservation-matched
controls, whereas which positions are in it largely is.

### 3.4 Site-anchored coupling: null

The allosteric pocket is not preferentially coupled to the ATP site relative to
distance-matched non-site residues, by either method — 1.05× (p = 0.81) for DCA and 1.02×
(p = 0.88) for SCA. Note that the EAI045/JBJ pocket lies 3.2–8.5 Å from the ATP site, so the
distance matching is stringent and "long-range" is not really the right description of this
pocket in any case.

### 3.5 Spectral decomposition: uninformative

Contact enrichment ranges from 0.83 to 1.18 across eigenmodes, i.e. essentially none, and every
mode is delocalised (participation ratio 186–218 of 253 positions; median 216). There is no
separation into localised contact-like and extended collective modes. Sampling only the extremes
of |λ| gives an apparently strong correlation between eigenvalue and contact enrichment
(ρ = −0.63), but sampling across the full eigenvalue range gives ρ = −0.10 (p = 0.27). We
attribute the absence of structure to the λ<sub>J</sub> = 1.0 regularisation. Report as a
supplementary negative result or omit.

### 3.6 Gap-state check: immaterial

Restricting the Frobenius norm to the 20 amino-acid states changes little:
Spearman ρ = 0.986 with the published scores overall and 0.987 at |i−j| ≥ 24, with 462 of 500
top long-range pairs shared. No conclusion depends on this choice.

---

## 4. What this means for the manuscript

**Figure 7 cannot be presented as a "smoking gun" and should be rebuilt.** The honest version of
the panel is panel b of the new figure: the enrichment under the published null, under a matched
null, and for distal pairs only. The reframed claim — that evolutionary couplings and dynamic
correlations agree where residues touch, and that neither identifies long-range functional
coupling in this system — is defensible, interesting, and consistent with the published
literature on the limits of DCA.

**Replace the pairwise long-range claim with the sector result.** Sector 2 is a genuine
long-range, spatially contiguous functional network linking the αC glutamate, the catalytic
lysine, the HRD and DFG motifs and L858, and its compactness is not explained by conservation.
That is the strongest coevolutionary statement the data support, and it is a statement about
collective organisation rather than about individual residue pairs.

**Do not claim DCA detects the allosteric pocket.** Both the pairwise and the site-anchored tests
are null.

---

## 5. Specific corrections to the existing pipeline

1. **Invented error bars.** `makefig7.py` and `fig7spoof.py` hardcode `errors = [0.005, 0.008]`
   with the comment "Standard error for a high-confidence plot". These are not computed from the
   data. The bootstrap intervals are in `enrichment_tests.csv`. (The *means*, 0.287 and 0.156,
   are real and reproduce.)
2. **Confounded null.** `correlation_test.py` draws control pairs uniformly, which does not match
   the sequence-separation or distance distribution of the test set. See §3.2.
3. **"Long-range" defined as |i−j| > 5.** The field convention is |i−j| ≥ 24. Three of the four
   highlighted census pairs do not meet it.
4. **No 3D masking anywhere in the pipeline.** 43% of the top-500 set are contacts.
5. **Signed coupling scores summed over module pairs.** `allosteric_highway_quantification.csv`
   sums APC-corrected scores, which can be negative, producing negative "highway" sums and
   negative densities (e.g. community 0↔3: sum −140.5, density −0.209). A sum over signed
   APC values is not interpretable as a connection strength; use positive-part sums, rank-based
   aggregation, or the matched-null framework instead.
6. **Inconsistent community names** between `quantify_allostery.py` and
   `allosteric_highway_quantification.csv` (e.g. community 6 is "Central Allosteric Hub
   (αD/Linker)" in one and "Central Linker Hub" in the other; community 3 likewise). Fix before
   submission.
7. **Alignment coverage.** `hhfilter -cov 50` admits sequences covering only half the domain;
   47.9% of the 32,912 sequences carry more than 20% gaps across the 253 focus columns. This
   inflates the nominal effective sequence count. It does not change the conclusions above
   (§3.6), but a reviewer may raise it.

---

## 6. Anticipated reviewer objections and responses

**"DCA cannot detect allosteric couplings; your long-range claims are artefacts."**
Agreed, and now demonstrated quantitatively in this system rather than asserted. We report the
contact composition of the coupling set, the collapse of the enrichment under a matched null, and
the absence of any distal signal, and we have moved the coevolutionary claim to sector
organisation.

**"Your sectors are just conserved positions."**
Partly — we report ρ = 0.832 between conservation and component loading explicitly. The claim we
make is about three-dimensional compactness, tested against conservation-matched controls, where
sector 2 is far outside the null (0 of 2,000 draws).

**"Why should the dynamic cross-correlation matrix be the arbiter of functional coupling?"**
It should not, and we no longer treat it as one. A single ~100 ns trajectory gives a noisy
estimate dominated by local and collective rigid-body motion. This is a further reason to
present §3.2 as a negative control rather than as validation.

**"Did you use a conformational ensemble or a single structure?"**
328 EGFR entries spanning active, αC-out inactive, ligand-bound and asymmetric-dimer states; a
pair is excluded if it contacts in any of them.

**"Would a protein language model do better?"**
Methods that extract couplings from language models are benchmarked on contact recovery, and a
fine-tuned ESM2 binding-site predictor collapses on kinase allosteric sites specifically
(AUPR 0.06 versus 0.64–0.76 for orthosteric sites). We did not pursue this.

---

## 7. Software and data

pySCA 6.1 (`ranganathanlab/pySCA`, master); Python 3.11.15; NumPy 2.4.6; SciPy 1.17.1;
Biopython 1.88; pandas 2.3.3; NetworkX 3.6.1. Structures retrieved from RCSB PDB
(search and data APIs, October 2026).

### Supporting files

- `Figure_R1_coevolution_long_range.png` — the four-panel revision figure
- `results_summary.csv` — all headline numbers in one table
- `pair_table.csv.gz` — all 31,878 residue pairs with couplings, cross-state distances, SCA, DCCM, sector and community assignments
- `sca_sectors.csv`, `sector_spatial_tests.csv` — sector membership and the conservation-matched spatial tests
- `enrichment_tests.csv`, `head_to_head_distal.csv` — matched-null enrichment results
- `site_anchored_test.csv`, `site_definitions.json` — site-anchored analysis and the ligand-derived site definitions
- `hopfield_patterns.csv` — eigenmode table for the supplementary spectral analysis
- `egfr_sca_sectors.pdb` — EGFR structure with the SCA sector index in the B-factor column
- `contact_qc.json` — every chain used in the contact mask, with its alignment identity and the 18 rejected chains
- `method_assessment.md`, `method_assessment.csv` — the literature assessment behind the method choice
