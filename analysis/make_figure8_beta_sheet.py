"""FIGURE_8_beta_sheet.png - correlated motion tracks beta-sheet pairing.

Everything here is recomputed from data/processed/pair_table_ss.csv.gz, which
carries, for each of the 31,878 residue pairs: the published APC coupling, the
dynamic cross-correlation from replicate 1, the minimum Ca-Ca distance in the
AlphaFold model, the per-residue secondary-structure class (DSSP over 501
trajectory frames, >=50% of frames) and the backbone H-bond occupancy.

Panel A: mean |DCCM| against Ca-Ca distance, split by secondary-structure pair
         class - the beta-beta advantage is confined to contact range.
Panel B: the 15 cross-strand beta registers against distance-matched non-beta
         pairs (Mann-Whitney U).
Panel C: composition of the most strongly correlated pairs.

Distal pairs only (|i-j| > 5). Statistics are written to
results/tables/ss_test_summary.csv by the same analysis; the values printed
here should match it.
"""
import pathlib
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, FIGURES = ROOT / "data" / "processed", ROOT / "results" / "figures"
plt.rcParams.update({"font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5,
                     "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150})
FOC, GREY, BLU, ORA = "#C2410C", "#9E9E9E", "#1F6F8B", "#E58B5A"

tab = pd.read_csv(DATA / "pair_table_ss.csv.gz")
t = tab[tab.sep > 5].copy()
cs = t[t.cross_strand]
lo, hi = cs.ca_min.min(), cs.ca_min.max()
ctrl = t[(t.pair_class != "EE") & (t.ca_min >= lo) & (t.ca_min <= hi)]
mwu = mannwhitneyu(cs.absdccm, ctrl.absdccm)

fig = plt.figure(figsize=(7.4, 3.1))
gs = fig.add_gridspec(1, 3, wspace=0.46, width_ratios=[1.1, 0.8, 1.05])

ax = fig.add_subplot(gs[0, 0])
t2 = t.assign(db=pd.cut(t.ca_min, np.arange(4, 40, 4), labels=False))
for pc, lab, c in [("EE", "beta-beta", FOC), ("EH", "beta-alpha", ORA),
                   ("HH", "alpha-alpha", BLU), ("LL", "loop-loop", GREY)]:
    s = (t2[t2.pair_class == pc].groupby("db")
         .agg(m=("absdccm", "mean"), x=("ca_min", "mean"), n=("absdccm", "size")))
    s = s[s.n >= 20]
    ax.plot(s.x, s.m, "o-", color=c, ms=3.5, lw=1.6, label=lab)
ax.axvspan(4, 12, color=FOC, alpha=0.07, zorder=0)
ax.set_xlabel("Minimum Ca-Ca distance (A)"); ax.set_ylabel("Mean |dynamic correlation|")
ax.set_title("A. The beta-beta advantage is confined\nto contact range", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6, loc="upper right")

ax = fig.add_subplot(gs[0, 1])
rng = np.random.default_rng(1)
for k, (v, c) in enumerate([(cs.absdccm.values, FOC), (ctrl.absdccm.values, GREY)]):
    ax.scatter(k + rng.uniform(-0.14, 0.14, len(v)), v, s=9, color=c, alpha=0.7, linewidths=0)
    ax.plot([k - 0.27, k + 0.27], [np.median(v)] * 2, color="k", lw=1.7)
ax.set_xticks([0, 1])
ax.set_xticklabels(["cross-strand\nregisters", "distance-\nmatched"], fontsize=5.6)
ax.set_ylabel("|dynamic correlation|")
ax.set_title(f"B. {cs.absdccm.mean() / ctrl.absdccm.mean():.1f}x higher at the\nsame distance",
             loc="left", fontweight="bold")
_mant, _exp = f"{mwu.pvalue:.1e}".split("e")
ax.text(0.5, 0.03, rf"$p = {_mant}\times10^{{{int(_exp)}}}$",
        transform=ax.transAxes, ha="center", fontsize=6.5)
ax.set_xlim(-0.75, 1.75); ax.margins(y=0.12)

ax = fig.add_subplot(gs[0, 2])
ns = [50, 100, 200, 500, 1000, 2000, 5000]
ee = [100 * (t.nlargest(n, "absdccm").pair_class == "EE").mean() for n in ns]
hh = [100 * (t.nlargest(n, "absdccm").pair_class == "HH").mean() for n in ns]
bg_e, bg_h = 100 * (t.pair_class == "EE").mean(), 100 * (t.pair_class == "HH").mean()
ax.semilogx(ns, ee, "o-", color=FOC, ms=3.5, lw=1.6, label="beta-beta")
ax.semilogx(ns, hh, "s-", color=BLU, ms=3.5, lw=1.6, label="alpha-alpha")
ax.axhline(bg_e, color=FOC, ls=":", lw=1.0); ax.axhline(bg_h, color=BLU, ls=":", lw=1.0)
ax.annotate(f"beta-beta background {bg_e:.1f}%", xy=(60, bg_e), xytext=(62, 24), fontsize=5.8,
            color=FOC, arrowprops=dict(arrowstyle="->", color=FOC, lw=0.8))
ax.annotate(f"alpha-alpha background {bg_h:.1f}%", xy=(1500, bg_h), xytext=(140, 45), fontsize=5.8,
            color=BLU, arrowprops=dict(arrowstyle="->", color=BLU, lw=0.8))
ax.set_xlabel("Top $N$ pairs by |dynamic correlation|"); ax.set_ylabel("Percent of the top $N$")
ax.set_title("C. The most correlated pairs\nare beta-sheet pairs", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6, loc="center right")
ax.set_xticks([50, 500, 5000]); ax.set_xticklabels(["50", "500", "5k"]); ax.set_ylim(-3, 80)

out = FIGURES / "FIGURE_8_beta_sheet.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"{out.name}: cross-strand registers n={len(cs)} mean |DCCM|={cs.absdccm.mean():.3f}; "
      f"distance-matched non-beta n={len(ctrl)} mean={ctrl.absdccm.mean():.3f}; "
      f"ratio {cs.absdccm.mean() / ctrl.absdccm.mean():.2f}x, MWU p={mwu.pvalue:.2e}")
print(f"   top-50 beta-beta {ee[0]:.0f}% vs background {bg_e:.1f}%")
