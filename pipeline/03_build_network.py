"""
Stage 3 - residue interaction network and Louvain partition.

Input : egfr_couplings_translated.csv  (couplings mapped to UniProt P00533 numbering)
Output: egfr_communities_PRODUCTION_1.8L.csv

L is ambiguous in this project and the distinction matters:
  257 = EGFR residues spanning 712-968, the focus region requested from plmc
  253 = match-state (uppercase) columns actually scored; 712, 713, 967, 968 are
        insert-state and were dropped, so the analysed domain is 714-966
The deposited network uses int(257 * 1.80) = 462 edges. The percolation threshold
located directly is 458 edges (at 457 the giant component holds 252 of 253 nodes).
462 is retained here for continuity with the deposited partition.
"""
import pandas as pd, networkx as nx, community as community_louvain

L_FOCUS, L_SCORED, MULTIPLIER, SEEDS = 257, 253, 1.80, 100
EDGES = int(L_FOCUS * MULTIPLIER)

df = pd.read_csv("egfr_couplings_translated.csv").sort_values("score", ascending=False)
G = nx.Graph()
for r in df.head(EDGES).itertuples():
    G.add_edge(int(r.res_i), int(r.res_j), weight=r.score)
assert G.number_of_nodes() == L_SCORED, G.number_of_nodes()

best_q, best = -1, None
for seed in range(SEEDS):
    p = community_louvain.best_partition(G, weight="weight", random_state=seed)
    q = community_louvain.modularity(p, G, weight="weight")
    if q > best_q:
        best_q, best = q, p
print(f"{EDGES} edges, {G.number_of_nodes()} nodes, "
      f"{len(set(best.values()))} communities, weighted Q = {best_q:.4f}")

pd.DataFrame({
    "community": pd.Series(best),
    "weighted_degree": pd.Series(dict(G.degree(weight="weight"))),
    "betweenness": pd.Series(nx.betweenness_centrality(G, weight="weight")),
}).rename_axis("residue").reset_index().to_csv("egfr_communities_PRODUCTION_1.8L.csv", index=False)
