"""Replaces makefig4.py (submitted as FIGURE_2_percolation.png).

Changes:
  * modularity is computed with WEIGHTED Louvain, matching the production
    partition. makefig4.py plotted the unweighted curve from
    fig4_percolation_data.csv while the manuscript quoted the weighted value,
    so the annotated point did not lie on the plotted line.
  * the percolation point is located by scanning edge counts rather than
    assumed at 1.80L, and L is the true node count (253), not 257.
  * the title no longer claims a "mathematically optimal" scale: Q is maximal
    at the lowest density tested and decreases monotonically."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns
import networkx as nx, community as community_louvain

L, SEEDS = 253, 20
cp = pd.read_csv(DATA / "egfr_couplings_translated.csv").sort_values("score", ascending=False)

def graph(n):
    G = nx.Graph()
    for r in cp.head(n).itertuples():
        G.add_edge(int(r.res_i), int(r.res_j), weight=r.score)
    return G

rows = []
for m in np.round(np.arange(1.00, 4.001, 0.05), 2):
    E = int(L * m); G = graph(E)
    best = max(community_louvain.modularity(
                 community_louvain.best_partition(G, weight="weight", random_state=s),
                 G, weight="weight") for s in range(SEEDS))
    cc = len(max(nx.connected_components(G), key=len))
    rows.append(dict(edges=E, connectivity=100 * cc / L, modularity=best))
df = pd.DataFrame(rows)
df.to_csv(RESULTS / "tables" / "percolation_weighted.csv", index=False)

percolation = next(E for E in range(int(L), int(4 * L))
                   if graph(E).number_of_nodes() == L and nx.is_connected(graph(E)))
qmax_i = df.modularity.idxmax()
print(f"full connectivity first reached at {percolation} edges")
print(f"Q maximal ({df.modularity[qmax_i]:.4f}) at {int(df.edges[qmax_i])} edges")

plt.style.use("seaborn-v0_8-whitegrid"); sns.set_context("paper", font_scale=1.5)
fig, ax1 = plt.subplots(figsize=(10, 7))
ax1.set_xlabel("Edges Retained in Network", fontweight="bold")
ax1.set_ylabel("Giant Component Connectivity (%)", color="#1f77b4", fontweight="bold")
l1 = ax1.plot(df.edges, df.connectivity, color="#1f77b4", lw=3, label="Connectivity", zorder=3)
ax1.tick_params(axis="y", labelcolor="#1f77b4"); ax1.set_ylim(40, 105)

ax2 = ax1.twinx(); ax2.grid(False)
l2 = ax2.plot(df.edges, df.modularity, color="#d62728", lw=3, ls="--",
              label="Modularity (Q), weighted", zorder=3)
ax2.set_ylabel("Modularity Score (Q)", color="#d62728", fontweight="bold")
ax2.tick_params(axis="y", labelcolor="#d62728"); ax2.set_ylim(0.6, 0.95)

ax1.axvline(percolation, color="black", ls=":", lw=2, zorder=4)
ax1.annotate(f"Full connectivity first\nreached at {percolation} edges",
             xy=(percolation, 100), xytext=(percolation + 100, 88),
             arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=8),
             fontsize=12, fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85))
ax2.annotate(f"Q is maximal here ({df.modularity[qmax_i]:.4f})\nand falls monotonically",
             xy=(df.edges[qmax_i], df.modularity[qmax_i]), xytext=(330, 0.925),
             arrowprops=dict(facecolor="#8E2A20", shrink=0.05, width=1, headwidth=7),
             fontsize=10, color="#8E2A20", fontweight="bold")
ax1.text(300, 45, "FRAGMENTED\n(Disconnected)", fontsize=10, fontstyle="italic", ha="center", color="gray")
ax1.text(900, 45, "SATURATED\n(Weak edges)",   fontsize=10, fontstyle="italic", ha="center", color="gray")
ax1.legend(l1 + l2, [l.get_label() for l in l1 + l2], loc="center right", frameon=True)
plt.title("Determination of the Network Threshold by Percolation Analysis",
          fontweight="bold", size=16, pad=20)
plt.tight_layout(); plt.savefig(RESULTS / "figures" / "FIGURE_2_percolation.png", dpi=300)
