"""FIGURE_S1_replicates.png - between-replicate agreement for the five 100 ns trajectories.

Panel A: enrichment of the top-500 co-evolving pairs in correlated motion, per
         replicate, under the sequence-separation null and the distance-matched
         null (results/tables/replicate_enrichment.csv).
Panel B: Pearson r between replicate DCCMs, strict upper triangle
         (results/tables/replicate_dccm_corr_matrix.csv).

Both inputs are produced by pipeline/05_md_replicates/aggregate_replicates.py.
Replicate 1 is the trajectory used throughout the main text.
"""
import pathlib
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
TABLES, FIGURES = ROOT / "results" / "tables", ROOT / "results" / "figures"
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150})

E = pd.read_csv(TABLES / "replicate_enrichment.csv")
E = E[E.replicate.astype(str) != "mean5"].copy()
E["replicate"] = E.replicate.astype(int)
M = pd.read_csv(TABLES / "replicate_dccm_corr_matrix.csv", index_col=0).values

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3),
                               gridspec_kw=dict(width_ratios=[1.15, 1.0], wspace=0.42))

cols = {"sequence-separation": "#90A4AE", "distance-matched": "#C62828"}
off = {"sequence-separation": +0.16, "distance-matched": -0.16}
for null, g in E.groupby("null"):
    axA.errorbar(g.enrichment, g.replicate.values + off[null],
                 xerr=[g.enrichment - g.ci_low, g.ci_high - g.enrichment],
                 fmt="o", ms=4.4, lw=1.3, capsize=2.6, color=cols[null], zorder=3)
axA.axvline(1.0, color="0.25", lw=1.0, ls="--", zorder=1)
axA.text(1.0, 5.6, " no enrichment", rotation=90, va="bottom", ha="left", fontsize=6, color="0.25")
axA.set_yticks([1, 2, 3, 4, 5]); axA.set_yticklabels([f"rep {i}" for i in range(1, 6)])
axA.set_ylim(5.7, 0.45); axA.set_xlim(0.80, 2.15)
axA.set_xlabel("Enrichment in correlated motion (fold)")
axA.set_title("Enrichment survives the sequence null\nand vanishes under distance matching", loc="left")
axA.text(1.05, 2.4, "distance-matched null", color=cols["distance-matched"],
         fontsize=6, ha="left", fontweight="bold")
axA.text(1.95, 0.95, "sequence-separation null", color=cols["sequence-separation"],
         fontsize=6, ha="center", fontweight="bold")
axA.text(2.12, 5.55, "rep 1 = published", fontsize=6, color="0.35", ha="right")

im = axB.imshow(M, cmap="YlGnBu", vmin=0.70, vmax=1.0)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axB.text(j, i, "1" if i == j else f"{M[i, j]:.2f}", ha="center", va="center",
                 fontsize=6, color="white" if M[i, j] > 0.90 else "0.15")
axB.set_xticks(range(5)); axB.set_xticklabels([str(i) for i in range(1, 6)])
axB.set_yticks(range(5)); axB.set_yticklabels([str(i) for i in range(1, 6)])
axB.set_xlabel("Replicate"); axB.set_ylabel("Replicate")
axB.set_title("DCCMs agree across replicates\n(Pearson r, strict upper triangle)", loc="left")
cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04)
cb.set_label("Pearson r", fontsize=6); cb.ax.tick_params(labelsize=6)
for ax, letter in ((axA, "A"), (axB, "B")):
    ax.text(-0.14, 1.12, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

out = FIGURES / "FIGURE_S1_replicates.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
iu = np.triu_indices(M.shape[0], 1)
print(f"{out.name}: inter-replicate r = {M[iu].min():.2f}-{M[iu].max():.2f} "
      f"(mean {M[iu].mean():.2f})")
