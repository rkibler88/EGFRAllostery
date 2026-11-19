"""Replaces 1bschematic.py and 1c.py (submitted as FIGURE_1_overview.png panels B and C).

Changes:
  * 1bschematic.py drew the DFG motif at residues 853-855 and the alphaG helix
    at 905-920. The DFG motif is D855-F856-G857, and the observed alphaG helix
    in UniProt P00533 spans 893-908. Both corrected, and the DFG position is
    asserted against the residue identities in egfr_couplings.txt.
  * 1c.py hardcoded "93,773 Tyrosine Kinases" in box I and an "Evolutionary-
    Dynamic Convergence" label on the arrow. The alignment actually used was
    32,912 sequences, and the convergence does not survive a spatial control."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import matplotlib.pyplot as plt, matplotlib.patches as patches

N_SEQ, N_PAIRS, N_POS, EDGES, Q, PERC = 32912, 31878, 253, 462, 0.8036, 458
NEFF_OVER_L, ENRICH_MATCHED, CI_MATCHED = 109, 1.02, (0.94, 1.12)

ref = "".join(open(DATA / "egfr_reference_aligned.fasta").read().split("\n")[1:]).strip()
col2res, cur = {}, 712
for i, ch in enumerate(ref):
    if ch.isalpha(): col2res[i + 1] = cur; cur += 1
aa = {}
for line in open(DATA / "egfr_couplings.txt"):
    p = line.split()
    for c, x in ((int(p[0]), p[1]), (int(p[2]), p[3])):
        if c in col2res: aa[col2res[c]] = x
assert (aa[855], aa[856], aa[857]) == ("D", "F", "G"), "DFG is not at 855-857"
assert aa[858] == "L" and aa[790] == "T"
print("residue identity checks passed: DFG = D855-F856-G857, L858, T790")

C = {'ECD':'#D3D3D3','Kinase':'#708090','P-loop':'#FFFF00','alphaC':'#FF69B4','Hinge':'#FFA500',
     'A-loop':'#FF0000','DFG':'#FF8C00','alphaG':'#228B22','Tail':'#800080'}
fig, (a1, a2) = plt.subplots(2, 1, figsize=(14, 7), gridspec_kw={'height_ratios':[1, 2]})
plt.subplots_adjust(hspace=0.7)
a1.set_xlim(0, 1210); a1.set_ylim(0, 10); a1.axis('off')
a1.add_patch(patches.Rectangle((1, 4), 1210, 2, color=C['ECD']))
a1.add_patch(patches.Rectangle((712, 4), 256, 2, color=C['Kinase']))
a1.text(605, 7.5, 'Full EGFR (1-1210)', ha='center', fontweight='bold', fontsize=12)
a1.text(840, 2.5, 'Kinase Domain', ha='center', fontsize=10)
a2.set_xlim(712, 970); a2.set_ylim(0, 12)
a2.set_xlabel('Residue Number (UniProt P00533)', fontweight='bold', fontsize=11)
a2.get_yaxis().set_visible(False)
for s in ['top','right','left']: a2.spines[s].set_visible(False)
a2.add_patch(patches.Rectangle((712, 4), 258, 3, color='#F5F5F5', ec='black', lw=1, zorder=1))
for ms, me, mc, mn, yp in [(718,724,C['P-loop'],'P-loop',8.5),(752,768,C['alphaC'],'αC-helix',8.5),
                           (790,795,C['Hinge'],'Hinge',8.5),(853,874,C['A-loop'],'Activation Loop',10.5),
                           (855,857,C['DFG'],'DFG (855–857)',8.5),(893,908,C['alphaG'],'αG-helix',8.5),
                           (960,968,C['Tail'],'Tail Base',8.5)]:
    w = me - ms
    a2.add_patch(patches.Rectangle((ms, 4), w, 3, color=mc, zorder=2))
    a2.plot([ms + w/2, ms + w/2], [7, yp - 0.2], color='black', lw=0.5, zorder=0)
    a2.text(ms + w/2, yp, mn, ha='center', fontsize=9, fontweight='bold')
a2.plot(858, 3.5, marker='^', color='black', markersize=10)
a2.text(858, 1.5, 'L858R\n(Oncogenic)', ha='center', fontsize=9, fontweight='bold')
plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_1_overview_panelB.png", dpi=300, bbox_inches='tight'); plt.close()

fig, ax = plt.subplots(figsize=(16, 14)); ax.set_xlim(0, 100); ax.set_ylim(-15, 110); ax.axis('off')
bs = "round,pad=1.0"; fh = {'fontweight':'bold','fontsize':13}; ft = {'fontsize':11,'linespacing':1.4}
for x, y, w, h, t, b, fc, ec in [
 (4,75,26,25,"I. EVOLUTIONARY INPUT",
  f"{N_SEQ:,} sequences (Pfam PF07714,\nprotein Tyr and Ser/Thr kinase)\n62% identity clustering, $N_{{eff}}/L$ = {NEFF_OVER_L}",'#E3F2FD','#1565C0'),
 (37,75,26,25,"II. CO-EVOLUTION (DCA)",
  f"PLMC Model Inference\n{N_PAIRS:,} APC-corrected couplings\nover {N_POS} match-state positions",'#E3F2FD','#1565C0'),
 (70,75,26,25,"III. MODULAR NETWORK",
  f"15 Communities\n{EDGES}-edge network ($Q$ = {Q})\nFull connectivity from {PERC} edges",'#FFF3E0','#E65100'),
 (37,35,26,25,"IV. PHYSICS COMPARISON",
  "100 ns Explicit Solvent MD\nof the monomeric domain\nDynamic Correlation (DCCM)",'#E8F5E9','#2E7D32')]:
    ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle=bs, facecolor=fc, edgecolor=ec, lw=2))
    ax.text(x + w/2, y + 21, t, ha='center', **fh, color=ec)
    ax.text(x + w/2, y + 5, b, ha='center', **ft)
ax.add_patch(patches.FancyBboxPatch((25, 0), 50, 25, boxstyle=bs, facecolor='#FFEBEE', edgecolor='#C62828', lw=3))
ax.text(50, 21, "V. WHAT THE ANALYSIS SUPPORTS", ha='center', fontweight='bold', fontsize=16, color='#C62828')
ax.text(50, 5, "Community-level description of co-evolutionary organization\n"
               f"Enrichment in correlated motion absent under a Cα–Cα\n"
               f"distance-matched null ({ENRICH_MATCHED}×, 95% CI {CI_MATCHED[0]}–{CI_MATCHED[1]})\n"
               "No conformational transition was simulated", ha='center', fontsize=11.5, linespacing=1.6)
apr = dict(arrowstyle='simple', color='#424242', mutation_scale=30)
ax.annotate('', xy=(37, 87.5), xytext=(30, 87.5), arrowprops=apr)
ax.annotate('', xy=(70, 87.5), xytext=(63, 87.5), arrowprops=apr)
ax.annotate('', xy=(75, 15), xytext=(83, 75),
            arrowprops=dict(arrowstyle='simple', color='#E65100', connectionstyle="arc3,rad=-0.4", mutation_scale=30))
ax.annotate('', xy=(50, 25), xytext=(50, 35),
            arrowprops=dict(arrowstyle='simple', color='#2E7D32', mutation_scale=30))
ax.text(84, 45, "Compared against\nconformational dynamics\n(not independent validation)",
        color='#C62828', fontweight='bold', ha='center', fontsize=11,
        bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'))
plt.title("Computational Framework for Mapping EGFR Co-evolutionary Architecture",
          fontsize=19, fontweight='bold', y=1.02)
plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_1_overview_panelC.png", dpi=300, bbox_inches='tight')
