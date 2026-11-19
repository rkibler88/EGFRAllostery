"""Network knockout: change in global efficiency when a residue is removed.

Quoted in Results ("removing residue 793 reduces global network efficiency by
2.26%..."). Rewritten from run_perturbation.py in the working directory, which
hardcoded the edge count and printed to screen; here the 462-edge network is
built the same way as everywhere else and the table is written out.

Output: results/tables/network_knockout.csv
"""
import pathlib
import pandas as pd, networkx as nx

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, TABLES = ROOT / "data" / "processed", ROOT / "results" / "tables"
EDGES = 462
TARGETS = [793, 858, 965, 884, 733, 717, 799, 851, 790, 745]

df = pd.read_csv(DATA / "egfr_couplings_translated.csv").sort_values("score", ascending=False)
G = nx.Graph()
for r in df.head(EDGES).itertuples():
    G.add_edge(int(r.res_i), int(r.res_j), weight=r.score)

eff_wt = nx.global_efficiency(G)
rows = []
for node in TARGETS:
    if node not in G:
        rows.append({"residue": node, "in_network": False, "efficiency_drop_pct": float("nan")})
        continue
    H = G.copy(); H.remove_node(node)
    rows.append({"residue": node, "in_network": True,
                 "efficiency_drop_pct": 100 * (eff_wt - nx.global_efficiency(H)) / eff_wt})

out = pd.DataFrame(rows).sort_values("efficiency_drop_pct", ascending=False)
TABLES.mkdir(parents=True, exist_ok=True)
out.to_csv(TABLES / "network_knockout.csv", index=False)
print(f"wild-type global efficiency = {eff_wt:.6f} ({G.number_of_nodes()} nodes, {EDGES} edges)")
print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
