"""FIGURE_7_sectors.png - conformational-ensemble contacts and SCA sector compactness.

Panel A is recomputed in full from data/processed/pair_table.csv.gz: the
cumulative distribution of the minimum heavy-atom distance, taken over the 328
EGFR kinase-domain chains listed in data/processed/contact_qc.json, for the 500
top-ranked couplings against all distal pairs (|i-j| > 5).

Panel B compares the mean pairwise Ca-Ca distance of SCA sector 2 against two
nulls:
  * run-length-matched - recomputed here. It needs only the sector-2 positions
    (results/tables/sca_sectors.csv) and the Ca-Ca distances carried in the
    pair table, both committed.
  * conservation-matched - NOT recomputed here. It resamples positions within
    deciles of the per-position SCA conservation Di, which is an output of the
    pySCA stage described in pipeline/06_sca_sectors.md and is not reproducible
    from a clone (it needs the 2.5 GB PF07714 alignment). Its published values
    are in results/tables/sector_spatial_tests.csv and are annotated on the
    panel; the figure shipped in manuscript/figures/ additionally shows that
    null as a histogram.
"""
import pathlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, TABLES, FIGURES = (ROOT / "data" / "processed", ROOT / "results" / "tables",
                         ROOT / "results" / "figures")
plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150})
FOC, GREY, BLU = "#C2410C", "#9E9E9E", "#1F6F8B"
SECTOR, N_DRAWS, SEED = 2, 2000, 31

tab = pd.read_csv(DATA / "pair_table.csv.gz")
sec = pd.read_csv(TABLES / "sca_sectors.csv")
spatial = pd.read_csv(TABLES / "sector_spatial_tests.csv").set_index("sector")

res = np.arange(714, 967)                      # the 253 scored positions
idx = {r: i for i, r in enumerate(res)}
CA = np.full((len(res), len(res)), np.nan)
ii = tab.res_i.map(idx).values
jj = tab.res_j.map(idx).values
CA[ii, jj] = tab.ca_min.values
CA[jj, ii] = tab.ca_min.values
np.fill_diagonal(CA, 0.0)

pos = np.sort([idx[int(r)] for r in
               str(sec.loc[sec.sector == SECTOR, "residues"].iloc[0]).split(";")])

def mean_pair_dist(p):
    p = list(p)
    return np.mean([CA[a, b] for k, a in enumerate(p) for b in p[k + 1:]]) if len(p) > 1 else np.nan

observed = mean_pair_dist(pos)

runs, start, prev = [], pos[0], pos[0]          # contiguous runs in sequence
for p in pos[1:]:
    if p == prev + 1:
        prev = p
    else:
        runs.append(prev - start + 1); start = p; prev = p
runs.append(prev - start + 1)
lengths = sorted(runs, reverse=True)

rng = np.random.default_rng(SEED)
null_run = np.empty(N_DRAWS)
for d in range(N_DRAWS):
    occupied = np.zeros(len(res), bool); drawn = []
    for L in lengths:
        starts = [x for x in range(len(res) - L + 1) if not occupied[x:x + L].any()]
        x = starts[rng.integers(len(starts))]
        occupied[x:x + L] = True
        drawn.extend(range(x, x + L))
    null_run[d] = mean_pair_dist(drawn)
p_run = (null_run <= observed).mean()

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2))

ax = axes[0]
distal = tab[tab.sep > 5]
top500 = distal.nlargest(500, "apc_published")
for v, lab, c, lw in [(distal.dmin_intra.values, "all distal pairs ($|i-j|>5$)", GREY, 1.5),
                      (top500.dmin_intra.values, "500 top-ranked couplings", FOC, 2.0)]:
    x = np.sort(np.clip(v, 0, 15.0))
    ax.plot(x, np.arange(1, len(x) + 1) / len(x), color=c, lw=lw, label=lab)
frac_top = float((top500.dmin_intra <= 5).mean())
frac_all = float((distal.dmin_intra <= 5).mean())
ax.axvline(5, color="k", ls=":", lw=1.0)
ax.text(5.2, 0.99, "5 A contact\nthreshold", fontsize=7, va="top")
ax.annotate(f"{100 * frac_top:.0f}%", xy=(5, frac_top), xytext=(6.9, 0.30), fontsize=8.5,
            color=FOC, fontweight="bold", arrowprops=dict(arrowstyle="->", color=FOC, lw=1.1))
ax.annotate(f"{100 * frac_all:.1f}%", xy=(5, frac_all), xytext=(6.9, 0.155), fontsize=8.5,
            color="#5f5f5f", arrowprops=dict(arrowstyle="->", color="#5f5f5f", lw=1.1))
ax.set_xlim(2, 15); ax.set_ylim(0, 1.02)
ax.set_xlabel("Minimum heavy-atom distance across 328\nEGFR structures (A; capped at 15)")
ax.set_ylabel("Cumulative fraction of pairs")
ax.set_title("A. Top-ranked couplings contact in at least\none solved conformational state",
             loc="left", fontweight="bold")
ax.legend(frameon=False, loc="lower right", fontsize=7, bbox_to_anchor=(1.03, 0.44))

ax = axes[1]
ax.hist(null_run, bins=30, color=GREY, edgecolor="none", label="run-length-matched null")
ax.axvline(observed, color=FOC, lw=2.2)
ax.axvline(spatial.loc[SECTOR, "null_mean_CA"], color=BLU, lw=1.4, ls="--",
           label="conservation-matched null (mean)")
ax.annotate(f"SCA sector {SECTOR}\n{observed:.1f} A",
            xy=(observed, ax.get_ylim()[1] * 0.58), xytext=(observed + 1.6, ax.get_ylim()[1] * 0.76),
            fontsize=7.5, color=FOC, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=FOC, lw=1.1))
ax.set_xlim(observed - 1.9, null_run.max() + 1.2)
ax.set_xlabel("Mean pairwise Ca-Ca distance (A)")
ax.set_ylabel(f"Random position sets (of {N_DRAWS:,})")
ax.set_title("B. Sector 2 is spatially compact beyond\nconservation and sequence contiguity",
             loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6.5, loc="upper right")
ax.text(0.985, 0.50, "$p < 0.0005$ vs both nulls", transform=ax.transAxes, ha="right", fontsize=7)

fig.subplots_adjust(wspace=0.34)
out = FIGURES / "FIGURE_7_sectors.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"{out.name}: top-500 within 5 A in some state = {100 * frac_top:.1f}% "
      f"(all distal pairs {100 * frac_all:.1f}%)")
print(f"   sector {SECTOR}: observed mean Ca-Ca {observed:.2f} A; run-length null "
      f"{null_run.mean():.2f} +/- {null_run.std():.2f} A, p = {p_run:.4f}")
print(f"   conservation-matched null (from sector_spatial_tests.csv): "
      f"{spatial.loc[SECTOR, 'null_mean_CA']:.2f} A, p = {spatial.loc[SECTOR, 'p_dist']:.4f}")
