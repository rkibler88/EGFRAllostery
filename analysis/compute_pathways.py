"""Replaces backboneremove.py.

Two changes:
  * the dimer-to-tail path is COMPUTED, not written in as a literal list.
    The literal in backboneremove.py, [910,814,813,849,817,851,820,773,823,965],
    is a valid walk through the 462-edge network, but it is not the shortest
    path between its endpoints (the shortest is 6 steps, not 9).
  * the L858R result is reported as found rather than narrated. Residue 858's
    only edges in the 462-edge network are backbone edges, so deleting |i-j|=1
    removes the node entirely and no non-backbone path exists."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import networkx as nx, pandas as pd

EDGES = 462
cp = pd.read_csv(DATA / "egfr_couplings_translated.csv").sort_values("score", ascending=False)
top = cp.head(EDGES)

G, G_nb = nx.Graph(), nx.Graph()
for r in top.itertuples():
    i, j = int(r.res_i), int(r.res_j)
    G.add_edge(i, j, weight=1.0 / r.score)
    if abs(i - j) > 1:
        G_nb.add_edge(i, j, weight=1.0 / r.score)

def path(graph, a, b, label):
    if a not in graph or b not in graph:
        missing = [x for x in (a, b) if x not in graph]
        print(f"{label}: no path -- residue(s) {missing} have no qualifying edges")
        return None
    try:
        p = nx.shortest_path(graph, a, b, weight="weight")
        print(f"{label}: {p}  ({len(p)-1} steps)")
        return p
    except nx.NetworkXNoPath:
        print(f"{label}: no path")
        return None

print(f"network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
      f"connected={nx.is_connected(G)}\n")

# L858R to the DFG aspartate. 853 is isoleucine; the DFG motif is D855-F856-G857.
path(G,    858, 855, "858 -> D855, full network          ")
path(G_nb, 858, 855, "858 -> D855, backbone edges removed")

# alphaG-alphaH loop to the tail base, computed rather than asserted.
path(G, 910, 965, "910 -> 965, full network          ")
path(G_nb, 910, 965, "910 -> 965, backbone edges removed")
