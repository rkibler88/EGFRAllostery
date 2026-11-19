"""Replaces makefig3.py (submitted as FIGURE_5_flux.png panels B and C).

Changes:
  * community pairs are keyed on the UNORDERED pair. makefig3.py grouped on
    ('Name_i','Name_j'), which splits a pair across two rows depending on
    residue order; its input table had 126 rows for 105 community pairs, so
    the plotted values were partial sums.
  * the plotted quantity is coupling over the 462 edges that define the
    partition, not the signed sum over all |A|x|B| residue pairs. The scores
    are average-product-corrected and mean-zero, so the all-pairs sum is
    negative for 75 of 105 pairs and is not interpretable as pathway strength.
  * the "ALLOSTERIC INTEGRATION HUB" callout is removed; no module is
    designated a hub, since the ranking depends on which metric is used."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import pandas as pd, numpy as np, networkx as nx, matplotlib.pyplot as plt, seaborn as sns

EDGES = 462
CN = {0:"β1–β3 N-lobe",1:"P-loop / N-lobe core",2:"β4–β5 / hinge",3:"C-lobe core (αE/αH/αI)",
      4:"β4–β5 / αC distal",5:"αG–αH loop",6:"αD / catalytic-loop linker",7:"Regulatory αC-helix",
      8:"Activation loop / αEF",9:"αE-helix",10:"αH region",11:"αI / tail base",
      12:"Activation loop (DFG-proximal)",13:"Hinge–αE–activation loop",14:"Catalytic loop / αF"}

comm = pd.read_csv(DATA / "egfr_communities_PRODUCTION_1.8L.csv")
cmap = dict(zip(comm.residue, comm.community)); sizes = comm.community.value_counts().to_dict()
cp = pd.read_csv(DATA / "egfr_couplings_translated.csv").sort_values("score", ascending=False)

net = cp.head(EDGES).copy()
net["ci"] = net.res_i.map(cmap); net["cj"] = net.res_j.map(cmap)
inter = net[net.ci != net.cj].copy()
inter["pair"] = [tuple(sorted((a, b))) for a, b in zip(inter.ci, inter.cj)]
agg = inter.groupby("pair")["score"].agg(coupling="sum", edges="count").reset_index()
agg[["A", "B"]] = pd.DataFrame(agg.pair.tolist(), index=agg.index)
agg = agg.drop(columns="pair").sort_values("coupling", ascending=False)
agg.to_csv(RESULTS / "tables" / "intermodule_coupling_networkedges.csv", index=False)
print(f"{int(agg.edges.sum())} of {EDGES} edges are inter-community, "
      f"across {len(agg)} of {15*14//2} module pairs")
print(agg.head(4).to_string(index=False))
strongest = [(int(r.A), int(r.B)) for r in agg.head(2).itertuples()]

sns.set_context("paper", font_scale=1.4)
M = pd.DataFrame(np.nan, index=[CN[i] for i in range(15)], columns=[CN[i] for i in range(15)])
for r in agg.itertuples():
    M.iloc[int(r.A), int(r.B)] = M.iloc[int(r.B), int(r.A)] = r.coupling
plt.figure(figsize=(14, 11))
sns.heatmap(M, cmap="YlGnBu", annot=True, fmt=".2f", annot_kws={"size": 7.5},
            cbar_kws={"label": "$\\Sigma J_{ij}$ over inter-community network edges"},
            linewidths=.5, linecolor="white")
plt.title("C. Inter-Module Coupling Matrix (restricted to the 462 network edges)",
          fontweight="bold", loc="left", pad=20)
plt.tight_layout(); plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_5_flux_panelC.png", dpi=300); plt.close()

G = nx.Graph()
for c, n in sizes.items(): G.add_node(CN[c], size=n)
for r in agg.itertuples(): G.add_edge(CN[int(r.A)], CN[int(r.B)], weight=float(r.coupling))
plt.figure(figsize=(14, 14))
pos = nx.spring_layout(G, k=4.0, iterations=200, seed=42, weight="weight")
nx.draw_networkx_edges(G, pos, width=[G[u][v]["weight"]/4 for u, v in G.edges],
                       alpha=0.35, edge_color="gray")
hl = [(CN[a], CN[b]) for a, b in strongest]
nx.draw_networkx_edges(G, pos, edgelist=hl, width=[G[u][v]["weight"]/4 for u, v in hl], edge_color="#E53935")
nx.draw_networkx_nodes(G, pos, node_size=[G.nodes[n]["size"]*200 for n in G.nodes],
                       node_color="#1E88E5", edgecolors="white", linewidths=2)
for nd, (x, y) in pos.items():
    plt.text(x, y, s=nd, fontsize=10, fontweight="bold", ha="center", va="center",
             bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", boxstyle="round,pad=0.2"))
plt.title("B. Inter-Module Coupling Network", fontweight="bold", loc="left", size=18, pad=30)
plt.figtext(0.5, 0.045,
    f"Edge width $\\propto$ $\\Sigma J_{{ij}}$ over the {EDGES} network edges; red = the two strongest pairs "
    f"(C{strongest[0][0]}–C{strongest[0][1]}, C{strongest[1][0]}–C{strongest[1][1]}).\n"
    f"Only {int(agg.edges.sum())} of {EDGES} edges are inter-community, across {len(agg)} of 105 possible module pairs.",
    ha="center", fontsize=11, color="0.25")
plt.axis("off"); plt.margins(0.2)
plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_5_flux_panelB.png", dpi=300, bbox_inches="tight")
