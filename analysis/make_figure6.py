"""Replaces fig7spoof.py / makefig7.py (submitted as FIGURE_6_dccm.png).

Changes:
  * panel D is COMPUTED. fig7spoof.py hardcoded means = [0.156, 0.287] and
    errors = [0.004, 0.007]; the error bars were not derived from anything
    (makefig7.py used a different invented pair, [0.005, 0.008]). Here the
    means and standard errors come from the samples.
  * a Ca-Ca distance-matched null is added alongside the sequence-separation
    null. DCA preferentially recovers spatial contacts, and contacting
    residues co-move for mechanical reasons, so a sequence filter alone does
    not isolate allostery.
  * the mask used for normalisation is the same |i-j| > 5 filter used
    everywhere else, not the |i-j| > 10 used only here.
  * the annotation states that the 799-841 maximum is filter-dependent."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from matplotlib.patches import Rectangle
from scipy.stats import ttest_ind

MASK, N_SAMPLE, SEED, REPS = 5, 500, 20261013, 200
rng = np.random.default_rng(SEED)

dca = np.load(DATA / "fig7_dca_aligned.npy")
dccm = pd.read_csv(DATA / "egfr_dccm_matrix.csv", index_col=0)
dccm.index = dccm.index.astype(int); dccm.columns = dccm.columns.astype(int)
res = list(dccm.index); n = len(res)
cp = pd.read_csv(DATA / "egfr_couplings_translated.csv")

def ca_coords(path):
    xyz, ri = [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            ri.append(int(line[22:26]))
    ri = np.array(ri)
    return np.array(xyz), (ri + 711 if ri.min() == 1 else ri)

xyz, ri = ca_coords(DATA / "egfr_clean.pdb")
P = np.array([xyz[list(ri).index(r)] for r in res])
DIST = np.linalg.norm(P[:, None] - P[None, :], axis=-1)
dmap = {(res[i], res[j]): DIST[i, j] for i in range(n) for j in range(n)}

pairs = cp[(abs(cp.res_i - cp.res_j) > MASK) &
           cp.res_i.isin(res) & cp.res_j.isin(res)].copy()
pairs["absdccm"] = [abs(dccm.loc[a, b]) for a, b in zip(pairs.res_i, pairs.res_j)]
pairs["dist"]    = [dmap[(a, b)] for a, b in zip(pairs.res_i, pairs.res_j)]
top  = pairs.nlargest(N_SAMPLE, "score")
rest = pairs.drop(top.index)
test = top.absdccm.values

# Both nulls are averaged over REPS draws. A single draw is seed-dependent:
# the distance-matched enrichment moves between about 0.98 and 1.08 across seeds,
# so one draw cannot be quoted to two decimals.
bins = np.arange(0, 65, 2.0)
tb, rb = np.digitize(top.dist, bins), np.digitize(rest.dist, bins)
bin_counts = pd.Series(tb).value_counts()
pools = {b: rest.index[rb == b].to_numpy() for b in bin_counts.index}

seq_draws, dist_draws = [], []
for _ in range(REPS):
    seq_draws.append(rest.absdccm.sample(N_SAMPLE, random_state=int(rng.integers(1 << 31))).values)
    picks = [rng.choice(pools[b], size=k, replace=len(pools[b]) < k)
             for b, k in bin_counts.items() if len(pools[b])]
    dist_draws.append(rest.loc[np.concatenate(picks)].absdccm.values)

def summarise(draws, label):
    """Enrichment with a bootstrap CI over pairs. Resampling both the test set and
    the pooled null propagates sampling noise in both, which the spread across null
    draws alone does not. p is reported to one decimal: the draw-to-draw spread makes
    further digits unstable, so the interval is the inferential statement."""
    pooled = np.concatenate(draws)
    ratio = test.mean() / pooled.mean()
    boot = [rng.choice(test, test.size, replace=True).mean() /
            rng.choice(pooled, N_SAMPLE, replace=True).mean() for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    pv = ttest_ind(test, pooled, equal_var=False).pvalue
    print(f"{label:<34} null {pooled.mean():.4f}  enrichment {ratio:.2f}x  "
          f"bootstrap 95% CI {lo:.2f}-{hi:.2f}  p {pv:.2g}")
    return pooled.mean(), pooled, ratio, pv, (lo, hi)

seq_mean,  seq_pool,  e_seq, p_seq, ci_seq = summarise(seq_draws,  "sequence-separation null")
dist_mean, dist_pool, e_dst, p_dst, ci_dst = summarise(dist_draws, "Ca-Ca distance-matched null")
print(f"co-evolving distal pairs           mean {test.mean():.4f}  n {len(test)}")

sem = lambda v: v.std(ddof=1) / np.sqrt(N_SAMPLE)
groups = [("Random distal pairs\n(sequence-separation null)", seq_mean,  sem(seq_pool),  "#D3D3D3"),
          ("Distance-matched null\n(Ca-Ca, 2 A bins)",        dist_mean, sem(dist_pool), "#90A4AE"),
          ("Co-evolving\ndistal pairs",                        test.mean(), sem(test),   "#C62828")]

def normalise(M, mask=MASK):
    masked = M.copy()
    for i in range(n):
        masked[i, max(0, i - mask):min(n, i + mask + 1)] = 0
    return np.clip(M / np.percentile(masked, 99.5), 0, 1)

mirror = np.where(np.arange(n)[:, None] < np.arange(n)[None, :],
                  normalise(dca), normalise(np.abs(dccm.values))).astype(float)
np.fill_diagonal(mirror, np.nan)

plt.style.use("seaborn-v0_8-whitegrid")
fig = plt.figure(figsize=(24, 12))
gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 0.7], wspace=0.3)
ax1 = fig.add_subplot(gs[0]); sns.set_context("paper", font_scale=1.7)
sns.heatmap(mirror, cmap="YlOrRd", vmin=0.1, vmax=0.9, ax=ax1,
            cbar_kws={"label": "Normalized Interaction Strength"})
i799, i841 = res.index(799), res.index(841)
for x, y in [(i841, i799), (i799, i841)]:
    ax1.add_patch(Rectangle((x - 5, y - 5), 10, 10, fill=False, edgecolor="black", lw=4, zorder=10))
ax1.annotate("799-841 pair\nmaximum only under |i-j|>10,\nand by a 0.7% margin",
             xy=(i841, i799), xytext=(i841 + 34, i799 - 40),
             arrowprops=dict(arrowstyle="->", color="black", lw=3, connectionstyle="arc3,rad=-0.2"),
             fontweight="bold", fontsize=14,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.9))
ax1.set_title("A-C. Evolutionary Coupling (DCA) vs. Dynamic Correlation (DCCM)",
              fontweight="bold", loc="left", pad=30, size=20)
ax1.set_xlabel("Residue Number (Top: Evolutionary Couplings)", fontweight="bold", size=16)
ax1.set_ylabel("Residue Number (Bottom: Dynamic Correlations)", fontweight="bold", size=16)
tk = np.arange(0, n, 25)
ax1.set_xticks(tk); ax1.set_xticklabels([res[t] for t in tk])
ax1.set_yticks(tk); ax1.set_yticklabels([res[t] for t in tk])

ax2 = fig.add_subplot(gs[1])
ax2.bar([g[0] for g in groups], [g[1] for g in groups],
        yerr=[g[2] for g in groups], capsize=12,
        color=[g[3] for g in groups], edgecolor="black", lw=3)
ax2.text(0.5, 0.368, f"{e_seq:.2f}x vs sequence-separation null",
         ha="center", fontsize=15, fontweight="bold", color="#C62828")
ax2.text(1.5, 0.330, f"{e_dst:.2f}x vs distance-matched null\n95% CI {ci_dst[0]:.2f}–{ci_dst[1]:.2f}",
         ha="center", fontsize=14, fontweight="bold", color="#37474F")
ax2.set_ylim(0, 0.42)
ax2.set_ylabel("Mean Absolute Coordination (|C$_{ij}$|)", fontweight="bold", size=16)
ax2.set_title("D. Enrichment Collapses Under a Spatial Control",
              fontweight="bold", loc="left", pad=30, size=20)
ax2.tick_params(axis="x", labelsize=11)
plt.savefig(RESULTS / "figures" / "FIGURE_6_dccm.png", dpi=300, bbox_inches="tight")
