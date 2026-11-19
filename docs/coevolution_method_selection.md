# Beyond-DCA methods for the EGFR allostery revision

The literature does not support the claim that pairwise DCA couplings mark allosteric
communication between spatially separated sites, and the strongest evidence against it is
exactly what a hostile reviewer will cite. Across ~4,000 families, 25% of directly coevolving
pairs are separated by more than 5 Å and 3% by more than 15 Å, but 91% of the 5–15 Å
exceptions are in contact in at least one homologous structure, 35% of the >15 Å exceptions sit
at homo-oligomeric interfaces, and the authors state plainly that they find little evidence for
direct coupling between spatially separated allosteric and functional sites or for increased
coupling along putative allosteric pathways
([Anishchenko et al. 2017](https://doi.org/10.1073/pnas.1702664114)). In evolved in-silico
allosteric materials, where the ground truth is known by construction, DCA predicts point-mutation
costs well and recovers short-range epistasis but fails on long-range epistasis, a failure the
authors argue is generic whenever function requires subparts to act in concert
([Bravi et al. 2020](https://doi.org/10.1371/journal.pcbi.1007630);
[Yan et al. 2017](https://doi.org/10.1073/pnas.1615536114)). A decade earlier the same
conclusion was reached from the opposite direction: because allosteric mechanisms diverge rapidly
even among close homologs, there is little reason to expect MSA-derived correlations to recover
them at all ([Livesay et al. 2011](https://doi.org/10.1007/978-1-61779-334-9_21)).

The implication for this manuscript is that no beyond-DCA method will rescue the original
pairwise long-range claim, and chasing one is the wrong revision strategy. Two moves are
defensible. The first is subtractive: mask against a conformational ensemble and ask honestly
whether any distal coupling survives a matched null — a negative answer is publishable and
pre-empts the obvious attack. The second is a change of statistic: shift the claim from individual
long-range pairs to collective sector organisation and to site-aggregated coupling, which is where
the functional validation and the kinase precedent actually are.

## Two pipeline-specific problems to fix first

Two features of the current analysis are independently attackable regardless of which method is
added. First, the sequence-separation filter: |i−j| > 5 is the field's standard *minimum*
exclusion, not a definition of long range. CASP assessment defines long-range contacts as pairs
separated by at least 24 residues along the sequence
([Monastyrskyy et al. 2013](https://doi.org/10.1002/prot.24340)). The manuscript's top "allosteric
census" pairs — 717–725 (separation 8), 846–852 (6), 728–741 (13) — are short-range by that
convention, and only 799–841 (42) qualifies. Re-stating the census at |i−j| ≥ 24 costs nothing and
removes an easy objection. Second, the published score is an APC-corrected Frobenius norm over the
full 21-letter alphabet with the gap symbol included. Long stretches of gaps are a data feature
atypical of independent samples drawn from a Potts model, and explicitly modelling them improves
contact prediction ([Feinauer et al. 2014](https://doi.org/10.1371/journal.pcbi.1003847)); with
3,839 alignment columns mapped onto 253 model positions, the gap channel carries indel and
phylogenetic structure rather than interaction. Recomputing the norm over the 20 amino acid states
and reporting the rank correlation against the published scores is a cheap robustness check that
closes this off. Neither issue requires a new method.

## Methods worth running

### Contact-masked DCA against a conformational ensemble

This is the highest-payoff addition and it is a subtraction rather than a new inference, so its
evidence base is the structural literature itself. The AlphaFold-era restatement of the
Anishchenko analysis searched 2,012 families for couplings unsatisfied in *any* AF model of family
members and found that ~12% of coevolving pairs remain spatially distant, with those pairs more
likely to be disordered than their contacting counterparts
([Karamanos 2023](https://doi.org/10.1002/bip.23530)). That ~12% is the realistic ceiling on what a
distal set can contain, and it sets the expectation for this manuscript: if a far larger fraction of
the top-500 survives masking, the mask is too narrow. The design constraint follows directly from
Anishchenko: a single-structure mask is indefensible, because most apparently distal pairs are
contacts in some homolog or oligomeric state — which for ErbB means the mask must span active and
inactive conformers, the asymmetric activator/receiver dimer interface, and AF models, not 1M17
alone. Two nulls are then mandatory. Coupling magnitude decays with both sequence separation and
3D distance, so the null must be matched on both; and phylogeny-induced spurious couplings are
comparable in magnitude to the first false-positive contact predictions, which is precisely the
intermediate-magnitude band where a top-500 set with median separation 32 lives
([Rodriguez Horta & Weigt 2021](https://doi.org/10.1371/journal.pcbi.1008957)). The
phylogeny-preserving resampling scheme in that paper — randomising sequence data while retaining
conservation and phylogenetic relations but destroying structural and functional coupling — is the
right second null, and it is cheap. Positive precedent for combining structure with coevolution
exists and includes a kinase: a three-parameter model integrating local pocket descriptors,
coevolutionary coupling and dynamic allostery ranked every known allosteric pocket in the top three
across five proteins including p38-α, and led to two experimentally confirmed novel sites
([La Sala et al. 2023](https://doi.org/10.1039/d2sc06272k)).

### SCA sector decomposition

SCA has the best experimental validation of any method on this list, and all of it is against
functional rather than structural coupling. The original method was validated against double-mutant
thermodynamic cycles in the PDZ domain, where it predicted long-range energetic coupling to a
binding-site residue that mutagenesis confirmed
([Lockless & Ranganathan 1999](https://doi.org/10.1126/science.286.5438.295)); sectors were then
defined as quasi-independent, physically connected groups with distinct functional roles
([Halabi et al. 2009](https://doi.org/10.1016/j.cell.2009.07.038)). A complete single-mutation scan
of PSD95-PDZ3 showed sector positions are functionally sensitive while non-sector positions tolerate
substitution, and that adaptation to a new binding specificity initiates exclusively within sector
residues ([McLaughlin et al. 2012](https://doi.org/10.1038/nature11500)). The decisive test came
from a deep coupling scan measuring thousands of pairwise amino acid couplings across several
homologs, which found cooperative interactions loaded in a sparse, conserved, spatially contiguous
network that SCA captured quantitatively
([Salinas & Ranganathan 2018](https://doi.org/10.7554/eLife.34300)). Generative tests point the same
way, from artificial WW domains with class-specific peptide recognition
([Russ et al. 2005](https://doi.org/10.1038/nature03990)) to synthetic chorismate mutases with
natural-like catalytic function ([Russ et al. 2020](https://doi.org/10.1126/science.aba3304)), and
the method is packaged for direct use ([Rivoire et al. 2016](https://doi.org/10.1371/journal.pcbi.1004817)).

For this manuscript the governing fact is that the kinase domain has already been decomposed this
way. pySCA applied to a kinome-wide alignment identified 79 significantly coupled positions
resolving into five independent components and three quasi-independent sectors encoding catalytic
activity, substrate specificity and regulation, with the catalytic sector enriched for somatic
cancer mutations and regulatory-sector residues contacting allosteric kinase inhibitors
([Creixell et al. 2018](https://doi.org/10.1016/j.cels.2018.08.008)). That study used 389 human
kinases, which after redundancy correction gave 265 *effective* sequences. The present alignment has
N_eff = 27,613 — roughly a hundredfold deeper — so a tyrosine-kinase SCA here is better powered than
the published kinase precedent, and mapping the resulting sectors onto EGFR 714–966 is a direct
extension rather than a novel application.

The standing critique must be addressed head-on rather than ignored. For proteins in which SCA finds
a single sector, conservation is the dominating factor in the SCA matrix and conservation alone makes
statistically equivalent functional predictions; the recommendation from that analysis was to move
experimental attention to proteins with several sectors
([Tesileanu et al. 2015](https://doi.org/10.1371/journal.pcbi.1004091)). The kinase domain is exactly
such a protein, which is why the critique bites least here — but the conservation-only control should
be computed and reported alongside the sectors, not asserted away. A second caution is that
phylogenetic clades can masquerade as sectors; separating inter-clade from intra-clade coupling
yields sectors that are more spatially compact, phylogenetically distinct, and more strongly enriched
for activity-disrupting mutations ([Colavin et al. 2022](https://doi.org/10.1038/s41598-021-04260-1)),
and is worth applying if the first-pass sectors look clade-structured.

### Site-anchored coupling profiles

The cheapest useful addition changes the statistic from "which pairs rank highest" to "how strongly
is this pocket coupled to the orthosteric site", and the evidence says the relevant signal lives in
the weak couplings that top-N pair ranking discards. Allosteric sites were found to be strongly
coupled to the orthosteric site relative to non-functional sites, key allosteric residues were
inferred by comparing each pocket residue's coupling to the functional site, and the predictions
matched prior experimental work for BCR-ABL1, Tar and PDZ3 as well as key cancer mutations
([Xie et al. 2023](https://doi.org/10.7554/eLife.81850)). BCR-ABL1 is a tyrosine kinase with the
canonical myristoyl allosteric pocket, so this is close precedent, and the required input is only
the existing APC-corrected score matrix plus site definitions. The method's weakness is that it is
entirely determined by how anchor and control sets are drawn: conserved anchors inflate coupling to
everything, so decoy pockets matched for size and mean conservation are not optional. Site
definitions should be fixed from the erlotinib contacts in 1M17 and from independent structural and
pharmacological evidence for the αC-out pocket and the asymmetric-dimer interface, and stated in
Methods before the numbers are computed.

## Methods to downgrade or skip

### Hopfield-Potts spectral decomposition — supplementary only

The eigen-decomposition story is attractive and cheap, and the literature firmly establishes the
split between contact-like and collective structure: introducing the Hopfield-Potts model to
interpolate between PCA and DCA showed that low-eigenvalue correlation modes, which PCA discards,
are highly localised and concentrated on few sites that are in close 3D contact, while the top modes
are collective ([Cocco et al. 2013](https://doi.org/10.1371/journal.pcbi.1003176)). Independent
random-matrix theory shows phylogeny produces a power-law tail in the covariance spectrum that is
essentially independent of tree topology, and that truncating the largest eigenvectors improves
contact prediction ([Qin & Colwell 2018](https://doi.org/10.1073/pnas.1711913115)) — a principled
criterion for where to cut. What the literature does *not* establish is that the collective modes are
allosteric. Twenty to forty Hopfield patterns suffice for a close-to-generative model and form
interpretable motifs, but their distributed collective nature intrinsically limits contact-map
prediction, and the patterns cluster sequences into functional *subfamilies*
([Shimagaki & Weigt 2019](https://doi.org/10.1103/PhysRevE.100.032128)) — meaning a top mode may be a
clade label rather than a communication pathway. Validation throughout is against contacts and
generative capacity, never against measured functional coupling.

Three technical hazards apply specifically to decomposing a plmc J tensor. The eigenspectrum of J is
gauge-dependent, so an explicit zero-sum (Ising) gauge fixing must precede the decomposition.
λ_J = 1.0 shrinks small couplings preferentially and therefore distorts exactly the low-eigenvalue
tail identified as contact-like, and regularisation strength is known to reshape inferred Potts
parameters ([Haldane & Levy 2019](https://doi.org/10.1103/PhysRevE.99.032405)). And
eigen-decomposing an already-fitted plmDCA coupling matrix is a descriptive spectral analysis, not
the reduced-rank inference performed in either source paper; presenting it as a Hopfield-Potts fit
would be inaccurate. Cost is negligible — eigendecomposition of the 5,313 × 5,313 coupling matrix
took 11.2 s at 0.75 GB peak on this machine — so run it, label it descriptive, and report each top
mode's correlation with clade identity.

The closest kinase precedent for decomposing a Potts model is energetic rather than spectral. A Potts
model of kinase residue covariation, combined with an alignment of 2,896 PDB structures, predicted
DFG-out propensity and identified activation-loop, αC-helix and HRD interactions as the dominant
contributors to DFG-in stability ([Haldane et al. 2016](https://doi.org/10.1002/pro.2954)); the
pairwise terms were shown necessary and sufficient to capture higher-than-pairwise kinase mutation
patterns ([Haldane et al. 2018](https://doi.org/10.1016/j.bpj.2017.10.028)); and the approach was
extended to a tyrosine-kinase-specific result — that TKs carry a 4–6 kcal/mol smaller free-energy
penalty than STKs for the folded-activation-loop DFG-out state — validated by 108 protein-mutation
free-energy perturbation simulations
([Gizzio et al. 2022](https://doi.org/10.7554/eLife.83368);
[Gizzio et al. 2024](https://doi.org/10.1038/s41467-024-50812-0)). If this manuscript wants a
statistical-energy claim about tyrosine-kinase allostery, that is the template, and it is a larger
undertaking than a spectral supplement.

### Subfamily-differential coupling — cite, do not compute

The kinase evidence here is strong but answers a different question. Determinants of specificity in
the human kinome form sparse networks of non-conserved residues spanning distant regions, three of
which (αC1, αC3, APE-7) were experimentally validated, implying a role for inter-residue allostery
in specificity and an evolutionary decoupling of activity from specificity
([Creixell et al. 2015](https://doi.org/10.1016/j.cell.2015.08.057)). The blocker is statistical, not
computational. Four human ErbB paralogs, and even a metazoa-wide ErbB alignment, yield at most a few
thousand sequences; model quality degrades with alignment depth most severely for contact prediction
and for double-mutant-cycle epistasis predictions, and overfitting becomes negligible only above a
few thousand effective sequences ([Haldane & Levy 2019](https://doi.org/10.1103/PhysRevE.99.032405)).
Differences between two separately regularised fits at different N_eff are dominated by sampling
noise, phylogenetic stratification is confounded with subfamily by construction, and depth-matched
bootstrap controls would likely absorb the whole signal. A TK-versus-STK split would be adequately
powered, but that analysis has substantially been done (Gizzio et al. 2022, 2024) and is a different
paper. Cite this line of work in the Discussion as the natural next step.

### Restricted Boltzmann machines — skip

RBM features on protein families are interpretable and map to tertiary contacts, extended secondary
motifs, disordered regions, activity and ligand specificity — and also to phylogenetic identity —
across 20 families with detailed results for Kunitz, WW, Hsp70 and synthetic lattice proteins
([Tubiana et al. 2019](https://doi.org/10.7554/eLife.39397)). None of that is a benchmark showing
RBMs recover long-range functional coupling better than DCA, and extra model capacity is not
automatically an improvement: comparing generative capacity through higher-order sequence
statistics ranked the pairwise Potts Hamiltonian above a variational autoencoder, with the
site-independent model last ([McGee et al. 2021](https://doi.org/10.1038/s41467-021-26529-9)). No
kinase application surfaced in these searches. The decisive practical objection is structural: a
hidden unit is a weighted pattern over positions, not a residue pair, so an RBM cannot answer whether
pair i–j is allosterically coupled; there is no agreed significance null; and each persistent
contrastive divergence training run takes hours on CPU with hyperparameter selection on top. This is
the only candidate here with a material compute burden, and it buys the least.

## Verdict

**Run, in this order.** (1) Contact-masked DCA against an ErbB conformational ensemble with
separation- and distance-matched nulls plus the phylogeny-preserving resampling null — it answers the
reviewer directly, costs hours, and is robust to its own negative result. (2) pySCA sector
decomposition of the tyrosine-kinase alignment mapped to EGFR 714–966, with the conservation-only
control reported alongside — it has the only genuine functional validation on this list and a direct
kinase precedent at a hundredth of the present alignment depth. (3) Site-anchored coupling profiles
for the ATP site versus candidate allosteric pockets with matched decoys — minutes of work, kinase
precedent in BCR-ABL1, and it supplies the orthosteric-versus-allosteric contrast the figure needs.

**Run as a labelled supplement if budget allows.** (4) Hopfield-Potts spectral decomposition of the
gauge-fixed J tensor, presented as a descriptive separation of localised from collective structure
with clade correlations shown, not as evidence of allostery.

**Skip.** (5) An ErbB-only subfamily-differential refit — too shallow to defend, and the hostile
reading ("small-sample artifact") would probably be correct. (6) RBMs — highest cost, no pairwise
read-out, no null, no kinase precedent.

Alongside these, restate the long-range census at |i−j| ≥ 24 and report the gap-excluded
Frobenius-norm robustness check. Both are free and both remove objections that no added method
would address.

## References

- Anishchenko I, Ovchinnikov S, Kamisetty H, Baker D (2017) Origins of coevolution between residues distant in protein 3D structures. *PNAS*. [10.1073/pnas.1702664114](https://doi.org/10.1073/pnas.1702664114) · PMID 28784799
- Bravi B, Ravasio R, Brito C, Wyart M (2020) Direct coupling analysis of epistasis in allosteric materials. *PLoS Comput Biol*. [10.1371/journal.pcbi.1007630](https://doi.org/10.1371/journal.pcbi.1007630) · PMID 32119660
- Cocco S, Monasson R, Weigt M (2013) From principal component to direct coupling analysis of coevolution in proteins. *PLoS Comput Biol*. [10.1371/journal.pcbi.1003176](https://doi.org/10.1371/journal.pcbi.1003176) · PMID 23990764
- Colavin A, Atolia E, Bitbol A-F, Huang KC (2022) Extracting phylogenetic dimensions of coevolution reveals hidden functional signals. *Sci Rep*. [10.1038/s41598-021-04260-1](https://doi.org/10.1038/s41598-021-04260-1) · PMID 35039514
- Creixell P, Palmeri A, Miller CJ, et al. (2015) Unmasking determinants of specificity in the human kinome. *Cell*. [10.1016/j.cell.2015.08.057](https://doi.org/10.1016/j.cell.2015.08.057) · PMID 26388442
- Creixell P, Pandey JP, Palmeri A, et al. (2018) Hierarchical organization endows the kinase domain with regulatory plasticity. *Cell Syst*. [10.1016/j.cels.2018.08.008](https://doi.org/10.1016/j.cels.2018.08.008) · PMID 30243563
- Dunn SD, Wahl LM, Gloor GB (2008) Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. *Bioinformatics*. [10.1093/bioinformatics/btm604](https://doi.org/10.1093/bioinformatics/btm604) · PMID 18057019
- Ekeberg M, Lövkvist C, Lan Y, Weigt M, Aurell E (2013) Improved contact prediction in proteins: using pseudolikelihoods to infer Potts models. *Phys Rev E*. [10.1103/PhysRevE.87.012707](https://doi.org/10.1103/PhysRevE.87.012707)
- Feinauer C, Skwark MJ, Pagnani A, Aurell E (2014) Improving contact prediction along three dimensions. *PLoS Comput Biol*. [10.1371/journal.pcbi.1003847](https://doi.org/10.1371/journal.pcbi.1003847) · PMID 25299132
- Gizzio J, Thakur A, Haldane A, Levy RM (2022) Evolutionary divergence in the conformational landscapes of tyrosine vs serine/threonine kinases. *eLife*. [10.7554/eLife.83368](https://doi.org/10.7554/eLife.83368) · PMID 36562610
- Gizzio J, Thakur A, Haldane A, Post CB, Levy RM (2024) Evolutionary sequence and structural basis for the distinct conformational landscapes of Tyr and Ser/Thr kinases. *Nat Commun*. [10.1038/s41467-024-50812-0](https://doi.org/10.1038/s41467-024-50812-0) · PMID 39095350
- Halabi N, Rivoire O, Leibler S, Ranganathan R (2009) Protein sectors: evolutionary units of three-dimensional structure. *Cell*. [10.1016/j.cell.2009.07.038](https://doi.org/10.1016/j.cell.2009.07.038) · PMID 19703402
- Haldane A, Flynn WF, He P, Vijayan RSK, Levy RM (2016) Structural propensities of kinase family proteins from a Potts model of residue co-variation. *Protein Sci*. [10.1002/pro.2954](https://doi.org/10.1002/pro.2954) · PMID 27241634
- Haldane A, Flynn WF, He P, Levy RM (2018) Coevolutionary landscape of kinase family proteins. *Biophys J*. [10.1016/j.bpj.2017.10.028](https://doi.org/10.1016/j.bpj.2017.10.028) · PMID 29320688
- Haldane A, Levy RM (2019) Influence of multiple-sequence-alignment depth on Potts statistical models of protein covariation. *Phys Rev E*. [10.1103/PhysRevE.99.032405](https://doi.org/10.1103/PhysRevE.99.032405) · PMID 30999494
- Hopf TA, Ingraham JB, Poelwijk FJ, et al. (2017) Mutation effects predicted from sequence co-variation. *Nat Biotechnol*. [10.1038/nbt.3769](https://doi.org/10.1038/nbt.3769) · PMID 28092658
- Hopf TA, Green AG, Schubert B, et al. (2019) The EVcouplings Python framework for coevolutionary sequence analysis. *Bioinformatics*. [10.1093/bioinformatics/bty862](https://doi.org/10.1093/bioinformatics/bty862) · PMID 30304492
- Karamanos TK (2023) Chasing long-range evolutionary couplings in the AlphaFold era. *Biopolymers*. [10.1002/bip.23530](https://doi.org/10.1002/bip.23530) · PMID 36752285 · PMC10909459
- La Sala G, Pfleger C, Käck H, et al. (2023) Combining structural and coevolution information to unveil allosteric sites. *Chem Sci*. [10.1039/d2sc06272k](https://doi.org/10.1039/d2sc06272k) · PMID 37389247
- Livesay DR, Kreth KE, Fodor AA (2011) A critical evaluation of correlated mutation algorithms and coevolution within allosteric mechanisms. *Methods Mol Biol*. [10.1007/978-1-61779-334-9_21](https://doi.org/10.1007/978-1-61779-334-9_21) · PMID 22052502
- Lockless SW, Ranganathan R (1999) Evolutionarily conserved pathways of energetic connectivity in protein families. *Science*. [10.1126/science.286.5438.295](https://doi.org/10.1126/science.286.5438.295) · PMID 10514373
- Marks DS, Colwell LJ, Sheridan R, et al. (2011) Protein 3D structure computed from evolutionary sequence variation. *PLoS One*. [10.1371/journal.pone.0028766](https://doi.org/10.1371/journal.pone.0028766) · PMID 22163331
- McGee F, Hauri S, Novinger Q, et al. (2021) The generative capacity of probabilistic protein sequence models. *Nat Commun*. [10.1038/s41467-021-26529-9](https://doi.org/10.1038/s41467-021-26529-9) · PMID 34728624
- McLaughlin RN Jr, Poelwijk FJ, Raman A, Gosal WS, Ranganathan R (2012) The spatial architecture of protein function and adaptation. *Nature*. [10.1038/nature11500](https://doi.org/10.1038/nature11500) · PMID 23041932
- Monastyrskyy B, D'Andrea D, Fidelis K, Tramontano A, Kryshtafovych A (2013) Evaluation of residue-residue contact prediction in CASP10. *Proteins*. [10.1002/prot.24340](https://doi.org/10.1002/prot.24340) · PMID 23760879
- Morcos F, Pagnani A, Lunt B, et al. (2011) Direct-coupling analysis of residue coevolution captures native contacts across many protein families. *PNAS*. [10.1073/pnas.1111471108](https://doi.org/10.1073/pnas.1111471108) · PMID 22106262
- Qin C, Colwell LJ (2018) Power law tails in phylogenetic systems. *PNAS*. [10.1073/pnas.1711913115](https://doi.org/10.1073/pnas.1711913115) · PMID 29311320
- Rivoire O, Reynolds KA, Ranganathan R (2016) Evolution-based functional decomposition of proteins. *PLoS Comput Biol*. [10.1371/journal.pcbi.1004817](https://doi.org/10.1371/journal.pcbi.1004817) · PMID 27254668
- Rodriguez Horta E, Weigt M (2021) On the effect of phylogenetic correlations in coevolution-based contact prediction in proteins. *PLoS Comput Biol*. [10.1371/journal.pcbi.1008957](https://doi.org/10.1371/journal.pcbi.1008957) · PMID 34029316
- Russ WP, Lowery DM, Mishra P, Yaffe MB, Ranganathan R (2005) Natural-like function in artificial WW domains. *Nature*. [10.1038/nature03990](https://doi.org/10.1038/nature03990) · PMID 16177795
- Russ WP, Figliuzzi M, Stocker C, et al. (2020) An evolution-based model for designing chorismate mutase enzymes. *Science*. [10.1126/science.aba3304](https://doi.org/10.1126/science.aba3304) · PMID 32703877
- Salinas VH, Ranganathan R (2018) Coevolution-based inference of amino acid interactions underlying protein function. *eLife*. [10.7554/eLife.34300](https://doi.org/10.7554/eLife.34300) · PMID 30024376
- Shimagaki K, Weigt M (2019) Selection of sequence motifs and generative Hopfield-Potts models for protein families. *Phys Rev E*. [10.1103/PhysRevE.100.032128](https://doi.org/10.1103/PhysRevE.100.032128) · PMID 31639992
- Tesileanu T, Colwell LJ, Leibler S (2015) Protein sectors: statistical coupling analysis versus conservation. *PLoS Comput Biol*. [10.1371/journal.pcbi.1004091](https://doi.org/10.1371/journal.pcbi.1004091) · PMID 25723535
- Tubiana J, Cocco S, Monasson R (2019) Learning protein constitutive motifs from sequence data. *eLife*. [10.7554/eLife.39397](https://doi.org/10.7554/eLife.39397) · PMID 30857591
- Xie J, Zhang W, Zhu X, et al. (2023) Coevolution-based prediction of key allosteric residues for protein function regulation. *eLife*. [10.7554/eLife.81850](https://doi.org/10.7554/eLife.81850) · PMID 36799896
- Yan L, Ravasio R, Brito C, Wyart M (2017) Architecture and coevolution of allosteric materials. *PNAS*. [10.1073/pnas.1615536114](https://doi.org/10.1073/pnas.1615536114) · PMID 28223497
