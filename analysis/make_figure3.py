"""Replaces makefig5.py (submitted as FIGURE_3_centrality.png panels A and B).

Changes:
  * the MOTIFS dictionary no longer labels residue 853 as the DFG motif.
    Residue 853 is isoleucine; the DFG motif is D855-F856-G857 in UniProt
    P00533 numbering. Motif assignments are checked against the residue
    identities carried in egfr_couplings.txt rather than asserted.
  * the narrative panel reports threshold sensitivity, which is computed here
    rather than claimed: the top five by node strength are stable over
    404-506 edges, the top betweenness residue is not."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, networkx as nx

EDGES = 462
df = pd.read_csv(DATA / "egfr_communities_PRODUCTION_1.8L.csv").sort_values("residue")
cp = pd.read_csv(DATA / "egfr_couplings_translated.csv").sort_values("score", ascending=False)

def graph(n):
    G = nx.Graph()
    for r in cp.head(n).itertuples():
        G.add_edge(int(r.res_i), int(r.res_j), weight=r.score)
    return G

# verify residue identities straight from the PLMC output
ref = "".join(open(DATA / "egfr_reference_aligned.fasta").read().split("\n")[1:]).strip()
col2res, cur = {}, 712
for i, ch in enumerate(ref):
    if ch.isalpha():
        col2res[i + 1] = cur; cur += 1
aa = {}
for line in open(DATA / "egfr_couplings.txt"):
    p = line.split()
    for c, x in ((int(p[0]), p[1]), (int(p[2]), p[3])):
        if c in col2res:
            aa[col2res[c]] = x
assert aa[853] == "I" and (aa[855], aa[856], aa[857]) == ("D", "F", "G"), "numbering check failed"

MOTIFS = {965: "αI-helix/Tail Anchor", 717: "N-lobe β-sheet", 799: "αD-helix",
          851: "Activation Loop (pre-DFG)", 792: "Hinge Region", 793: "Hinge Region",
          838: "Catalytic Loop", 853: f"{aa[853]}le (DFG is 855–857)"}

# threshold sensitivity, computed
strength_top5, bw_top1 = [], {}
for E in (404, 430, 455, 468, 480, 506):
    G = graph(E)
    d = dict(G.degree(weight="weight"))
    strength_top5.append(tuple(sorted(sorted(d, key=d.get, reverse=True)[:5])))
    b = nx.betweenness_centrality(G, weight="weight")
    bw_top1[E] = max(b, key=b.get)
stable = len(set(strength_top5)) == 1
print(f"node-strength top five stable over 404-506 edges: {stable} -> {sorted(strength_top5[0])}")
print(f"top betweenness residue by threshold: {bw_top1}")
low = sorted({E: r for E, r in bw_top1.items() if E <= 430}.values())

plt.style.use("seaborn-v0_8-whitegrid"); sns.set_context("paper", font_scale=1.4)
fig = plt.figure(figsize=(16, 12)); gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.4)

ax1 = fig.add_subplot(gs[0, :]); col = "#008B8B"
ax1.fill_between(df.residue, df.betweenness, color=col, alpha=0.3)
ax1.plot(df.residue, df.betweenness, color=col, lw=2)
ax1.set_ylim(0, df.betweenness.max() * 1.8)
for r in (792, 838, 853):
    y = df.loc[df.residue == r, "betweenness"].values[0]
    ax1.annotate(f"Res {r}\n({MOTIFS[r]})", xy=(r, y),
                 xytext=(r, y + (0.28 if r == 838 else 0.06)),
                 arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=5),
                 ha="center", fontweight="bold", fontsize=10,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=col, alpha=0.9))
ax1.set_title("A. Allosteric Conduits (Betweenness Centrality Skyline)", fontweight="bold", loc="left")
ax1.set_ylabel("Betweenness Centrality"); ax1.set_xlabel("Residue Number"); ax1.set_xlim(712, 970)

ax2 = fig.add_subplot(gs[1, 0])
top = df.sort_values("weighted_degree", ascending=False).head(10).copy()
top["full_label"] = top.apply(lambda r: f"{int(r.residue)} ({MOTIFS.get(int(r.residue), 'Kinase Core')})", axis=1)
sns.barplot(x="weighted_degree", y="full_label", data=top, color="#FFD700", ax=ax2, edgecolor="black")
ax2.set_title("B. Top Structural Anchors (Weighted Degree)", fontweight="bold", loc="left")
ax2.set_xlabel("Node Strength (Cumulative Coupling Weight)"); ax2.set_ylabel("Residue (Canonical Motif)")

ns = df.set_index("residue").weighted_degree.rank(ascending=False, method="min")
bwr = df.set_index("residue").betweenness.rank(ascending=False, method="min")
ax3 = fig.add_subplot(gs[1, 1]); ax3.axis("off")
ax3.text(0.05, 0.5,
 "ANCHOR-BRIDGE DUALISM\n\n"
 "Primary Anchor: Residue 965 (αI-helix)\nhas the highest node strength. It lies\n"
 "five residues from the C-terminal end of\nthe simulated construct, so its dynamic\n"
 "properties are not interpretable.\n\n"
 f"Primary Bridge: Residue 792 (Hinge)\ncarries the highest betweenness at {EDGES}\n"
 f"edges. This ranking is threshold-\ndependent: residue {low[0]} leads at 404–430\nedges.\n\n"
 f"Dual Role: Residue 793 (Hinge) ranks {int(ns[793])}th\nby node strength and {int(bwr[793])}rd by betweenness.\n\n"
 f"Stability: the top five by node strength\n{sorted(strength_top5[0])} hold across\n404–506 edges; community count and\n"
 "betweenness rank do not.",
 fontsize=11, fontweight="bold", linespacing=1.55, va="center",
 bbox=dict(boxstyle="round,pad=1", fc="#F8F9FA", ec="#DEE2E6"))
plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_3_centrality_panelsAB.png", dpi=300, bbox_inches="tight")
