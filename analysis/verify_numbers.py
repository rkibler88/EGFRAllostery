"""
Regression test: every number quoted in the manuscript and drawn in the figures,
recomputed from data/processed. Exits non-zero if any check fails.

    python analysis/verify_numbers.py        # or: make verify
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
import numpy as np, pandas as pd, networkx as nx, community as community_louvain

EDGES, SEEDS = 462, 100
cp = pd.read_csv(DATA / "egfr_couplings_translated.csv").sort_values("score", ascending=False)

def graph(n):
    G = nx.Graph()
    for r in cp.head(n).itertuples():
        G.add_edge(int(r.res_i), int(r.res_j), weight=r.score)
    return G

G = graph(EDGES)
bw = nx.betweenness_centrality(G, weight="weight")
ns = dict(G.degree(weight="weight"))
bwr = pd.Series(bw).rank(ascending=False, method="min")
nsr = pd.Series(ns).rank(ascending=False, method="min")
qbest, sbest = max((community_louvain.modularity(
    community_louvain.best_partition(G, weight="weight", random_state=s), G, weight="weight"), s)
    for s in range(SEEDS))
part = community_louvain.best_partition(G, weight="weight", random_state=sbest)

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

cds = pd.read_csv(DATA / "community_dynamics_summary.csv")
rmsf = lambda c: float(cds.loc[cds.community == c, "mean_rmsf"].iloc[0])

comm = pd.read_csv(DATA / "egfr_communities_PRODUCTION_1.8L.csv")
cmap = dict(zip(comm.residue, comm.community))
net = cp.head(EDGES).copy()
net["ci"], net["cj"] = net.res_i.map(cmap), net.res_j.map(cmap)
inter = net[net.ci != net.cj].copy()
inter["pair"] = [tuple(sorted((a, b))) for a, b in zip(inter.ci, inter.cj)]
agg = inter.groupby("pair")["score"].agg(coupling="sum", edges="count").sort_values("coupling", ascending=False)

# --- replicate molecular dynamics (CORRECTIONS 9.1, manuscript 2.9) ---------
# Five 100 ns trajectories. Replicate 1 is the published run. The enrichment
# method below is the one in make_figure6.py, re-run per replicate with the
# same fixed seed so differences reflect the trajectory and not the draw.
from itertools import combinations
from scipy.stats import pearsonr, spearmanr

MASK, N_SAMPLE, MDSEED, DRAWS, NBOOT = 5, 500, 20261013, 200, 2000

def _load_dccm(p):
    m = pd.read_csv(p, index_col=0)
    m.index = m.index.astype(int); m.columns = m.columns.astype(int)
    return m

DCCM = {1: _load_dccm(DATA / "egfr_dccm_matrix.csv")}
for _r in (2, 3, 4, 5):
    DCCM[_r] = _load_dccm(DATA / f"egfr_dccm_matrix_rep{_r}.csv")
reslist = list(DCCM[1].index)

_iu = np.triu_indices(len(reslist), k=1)
_vec = {r: m.values[_iu] for r, m in DCCM.items()}
rep_r = {(a, b): pearsonr(_vec[a], _vec[b])[0] for a, b in combinations(sorted(DCCM), 2)}
rep_r_all = np.array(sorted(rep_r.values()))
rep_r_vs1 = np.array([v for (a, b), v in rep_r.items() if a == 1])

def _ca(path):
    xyz, ri = [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            ri.append(int(line[22:26]))
    ri = np.array(ri)
    return np.array(xyz), (ri + 711 if ri.min() == 1 else ri)

_xyz, _ri = _ca(DATA / "egfr_clean.pdb")
_P = np.array([_xyz[list(_ri).index(r)] for r in reslist])
_DIST = np.linalg.norm(_P[:, None] - _P[None, :], axis=-1)
_dmap = {(reslist[i], reslist[j]): _DIST[i, j]
         for i in range(len(reslist)) for j in range(len(reslist))}
# file order, not score order: make_figure6.py reads the table unsorted and the
# null draws are index-order dependent, so sorting here would shift the ratio
# in the third decimal and disagree with the figure.
_cp_raw = pd.read_csv(DATA / "egfr_couplings_translated.csv")
_base = _cp_raw[(abs(_cp_raw.res_i - _cp_raw.res_j) > MASK) &
                _cp_raw.res_i.isin(reslist) & _cp_raw.res_j.isin(reslist)].copy()
_base["dist"] = [_dmap[(a, b)] for a, b in zip(_base.res_i, _base.res_j)]
_bins = np.arange(0, 65, 2.0)

def _enrich(dccm):
    """(sequence-separation, distance-matched) enrichment and bootstrap CIs."""
    rng = np.random.default_rng(MDSEED)
    pairs = _base.copy()
    pairs["absdccm"] = [abs(dccm.loc[a, b]) for a, b in zip(pairs.res_i, pairs.res_j)]
    top = pairs.nlargest(N_SAMPLE, "score"); rest = pairs.drop(top.index)
    test = top.absdccm.values
    tb, rb = np.digitize(top.dist, _bins), np.digitize(rest.dist, _bins)
    bc = pd.Series(tb).value_counts()
    pools = {b: rest.index[rb == b].to_numpy() for b in bc.index}
    seq, dst = [], []
    for _ in range(DRAWS):
        seq.append(rest.absdccm.sample(N_SAMPLE, random_state=int(rng.integers(1 << 31))).values)
        picks = [rng.choice(pools[b], size=k, replace=len(pools[b]) < k)
                 for b, k in bc.items() if len(pools[b])]
        dst.append(rest.loc[np.concatenate(picks)].absdccm.values)
    out = []
    for draws in (seq, dst):
        pooled = np.concatenate(draws)
        boot = [rng.choice(test, test.size, replace=True).mean() /
                rng.choice(pooled, N_SAMPLE, replace=True).mean() for _ in range(NBOOT)]
        out.append((test.mean() / pooled.mean(), *np.percentile(boot, [2.5, 97.5])))
    return out

_enr = {r: _enrich(DCCM[r]) for r in sorted(DCCM)}
enr_seq = np.array([_enr[r][0][0] for r in sorted(_enr)])
enr_dst = np.array([_enr[r][1][0] for r in sorted(_enr)])
dst_ci_spans_one = all(_enr[r][1][1] < 1.0 < _enr[r][1][2] for r in _enr)

rmsf_prof = {}
for _r, _p in [(1, DATA / "rmsf.xvg")] + [(r, DATA / f"rmsf_rep{r}.xvg") for r in (2, 3, 4, 5)]:
    _a = np.array([l.split() for l in open(_p) if l[0] not in "#@"], float)
    rmsf_prof[_r] = pd.Series(_a[:, 1], index=_a[:, 0].astype(int)).reindex(reslist)
rmsf_r = np.array([pearsonr(rmsf_prof[a], rmsf_prof[b])[0]
                   for a, b in combinations(sorted(rmsf_prof), 2)])
_cmap_res = dict(zip(comm.residue, comm.community))
_hier = []
for _r, _s in rmsf_prof.items():
    _d = pd.DataFrame({"rmsf": _s}); _d["c"] = _d.index.map(_cmap_res)
    _m = _d.dropna().groupby("c").rmsf.mean()
    _hier.append((int(_m.idxmax()), int(_m.idxmin())))

checks = [
    ("scored coupling pairs",            31878, len(cp)),
    ("match-state positions",            253,   len(set(cp.res_i) | set(cp.res_j))),
    ("communities at 462 edges",         15,    len(set(part.values()))),
    ("weighted Q at 462 edges",          0.8036, round(qbest, 4)),
    ("percolation threshold (edges)",    458,   next(E for E in range(400, 520)
                                                  if graph(E).number_of_nodes() == 253
                                                  and nx.is_connected(graph(E)))),
    ("betweenness rank of 792",          1,     int(bwr[792])),
    ("betweenness rank of 793",          3,     int(bwr[793])),
    ("betweenness of 792",               0.3118, round(bw[792], 4)),
    ("node-strength rank of 793",        7,     int(nsr[793])),
    ("top-5 node strength",              [965, 717, 799, 851, 964],
                                         sorted(nsr[nsr <= 5].index.tolist(), key=lambda r: nsr[r])),
    ("most rigid community",             3,     int(cds.loc[cds.mean_rmsf.idxmin(), "community"])),
    ("RMSF gradient C7/C3",              3.92,  round(rmsf(7) / rmsf(3), 2)),
    ("RMSF gradient C7/C14",             3.36,  round(rmsf(7) / rmsf(14), 2)),
    ("inter-community edges",            63,    int(agg.edges.sum())),
    ("module pairs with an edge",        27,    len(agg)),
    ("strongest pair C3-C6 coupling",    18.00, round(float(agg.coupling.iloc[0]), 2)),
    ("strongest pair C3-C6 edges",       8,     int(agg.edges.iloc[0])),
    ("second pair C3-C14 coupling",      12.95, round(float(agg.coupling.iloc[1]), 2)),
    ("residue 853 identity",             "I",   aa[853]),
    ("DFG motif at 855-857",             "DFG", aa[855] + aa[856] + aa[857]),
    ("residue 790 identity (gatekeeper)","T",   aa[790]),
    ("residue 858 identity",             "L",   aa[858]),
    ("858 has only backbone edges",      True,  all(abs(858 - j) == 1 for j in G.neighbors(858))),
    # --- replicate MD (CORRECTIONS 9.1, manuscript 2.9) ---
    ("replicate DCCM r, min",            0.746, round(float(rep_r_all.min()), 3)),
    ("replicate DCCM r, median",         0.840, round(float(np.median(rep_r_all)), 3)),
    ("replicate DCCM r, max",            0.899, round(float(rep_r_all.max()), 3)),
    ("replicate DCCM r vs rep1, min",    0.798, round(float(rep_r_vs1.min()), 3)),
    ("replicate DCCM r vs rep1, max",    0.888, round(float(rep_r_vs1.max()), 3)),
    ("replicate pairs compared",         10,    len(rep_r_all)),
    ("rep1 enrichment, sequence null",   1.83,  round(float(enr_seq[0]), 2)),
    ("rep1 enrichment, distance null",   1.02,  round(float(enr_dst[0]), 2)),
    ("enrichment seq null, min",         1.45,  round(float(enr_seq.min()), 2)),
    ("enrichment seq null, max",         1.88,  round(float(enr_seq.max()), 2)),
    ("enrichment distance null, min",    0.95,  round(float(enr_dst.min()), 2)),
    ("enrichment distance null, max",    1.02,  round(float(enr_dst.max()), 2)),
    ("enrichment distance null, SD",     0.03,  round(float(enr_dst.std(ddof=1)), 2)),
    ("all distance-null CIs span 1.0",   True,  dst_ci_spans_one),
    ("replicate RMSF r, min",            0.821, round(float(rmsf_r.min()), 3)),
    ("replicate RMSF r, max",            0.954, round(float(rmsf_r.max()), 3)),
    ("C7 most mobile in all replicates", True,  all(h[0] == 7 for h in _hier)),
    ("C3 most rigid in all replicates",  True,  all(h[1] == 3 for h in _hier)),
]

fails = []
for name, expect, got in checks:
    ok = (expect == got) if isinstance(expect, (str, bool, list)) else abs(expect - got) <= 0.02 * max(1, abs(expect))
    print(f"{'PASS' if ok else 'FAIL'}  {name:<34} expected {expect}  got {got}")
    if not ok:
        fails.append(name)

print(f"\n{len(checks) - len(fails)}/{len(checks)} checks passed")
if fails:
    print("FAILED:", ", ".join(fails)); sys.exit(1)
