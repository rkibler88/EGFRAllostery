"""Replaces makefig2.py (submitted as FIGURE_4_rmsf.png panels A and C).

Change: the rigidity gradient is computed from community IDs, not from a
substring match on community names. makefig2.py used
summary_df['name'].str.contains('Scaffold'), which matches both
"Scaffold-Switch Interface" (C8, 0.1209 nm) and "Main Engine Block (Scaffold)"
(C14, 0.0636 nm) and selected the first, giving 0.2135/0.1209 = 1.8x. The
correct most-mobile / most-rigid ratio is C7/C3 = 3.9x; against C14 it is 3.4x."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / "data" / "processed", ROOT / "results"
sys.path.insert(0, str(ROOT / "analysis"))

import pandas as pd, matplotlib.pyplot as plt, seaborn as sns

CLUTCH, SCAFFOLD = 7, 14
cc, sc = "#FF1493", "#2E8B57"

def parse_xvg(fn):
    rows = [[float(x) for x in l.split()] for l in open(fn) if not l.startswith(("@", "#"))]
    return pd.DataFrame(rows, columns=["residue", "rmsf"])

rmsf = parse_xvg(DATA / "rmsf.xvg")
summary = pd.read_csv(DATA / "community_dynamics_summary.csv").sort_values("mean_rmsf", ascending=False)

val = lambda c: float(summary.loc[summary.community == c, "mean_rmsf"].iloc[0])
most_mobile = int(summary.community.iloc[0]); most_rigid = int(summary.community.iloc[-1])
r_extreme = val(most_mobile) / val(most_rigid)
r_scaffold = val(CLUTCH) / val(SCAFFOLD)
print(f"most mobile C{most_mobile} = {val(most_mobile):.4f} nm; most rigid C{most_rigid} = {val(most_rigid):.4f} nm")
print(f"gradient C{most_mobile}/C{most_rigid} = {r_extreme:.2f}x ; C{CLUTCH}/C{SCAFFOLD} = {r_scaffold:.2f}x")

plt.style.use("seaborn-v0_8-whitegrid"); sns.set_context("paper", font_scale=1.5)
fig, ax1 = plt.subplots(figsize=(16, 6))
ax1.plot(rmsf.residue, rmsf.rmsf, color="black", lw=1.5, zorder=3)
ax1.fill_between(rmsf.residue, rmsf.rmsf, color="gray", alpha=0.2, zorder=2)
ax1.axvspan(748, 765, color=cc, alpha=0.3, label=f"Regulatory αC-helix (C{CLUTCH})")
ax1.axvspan(875, 896, color=sc, alpha=0.3, label=f"Catalytic loop / αF (C{SCAFFOLD})")
ax1.axvline(970, color="0.45", lw=1.2, ls="--")
ax1.annotate("construct ends at 970;\nresidue 965 lies 5 residues\nfrom a free terminus",
             xy=(969, 0.22), xytext=(900, 0.40), fontsize=10, color="0.25",
             arrowprops=dict(arrowstyle="->", lw=1, color="0.45"))
ax1.set_title("A. Residue-Level Fluctuation (RMSF)", fontweight="bold", loc="left")
ax1.set_ylabel("RMSF (nm)"); ax1.set_xlabel("Residue Number"); ax1.set_xlim(712, 975)
ax1.legend(loc="upper center", frameon=True, ncol=2)
plt.tight_layout(); plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_4_rmsf_panelA.png", dpi=300, bbox_inches="tight"); plt.close()

CN = {0:"β1–β3 N-lobe",1:"P-loop / N-lobe core",2:"β4–β5 / hinge",3:"C-lobe core (αE/αH/αI)",
      4:"β4–β5 / αC distal",5:"αG–αH loop",6:"αD / catalytic-loop linker",7:"Regulatory αC-helix",
      8:"Activation loop / αEF",9:"αE-helix",10:"αH region",11:"αI / tail base",
      12:"Activation loop (DFG-proximal)",13:"Hinge–αE–activation loop",14:"Catalytic loop / αF"}
summary = summary.reset_index(drop=True)
summary["label"] = [f"{CN[int(c)]} (C{int(c)})" for c in summary.community]
colors = ["#D3D3D3"] * len(summary)
colors[int(summary.index[summary.community == CLUTCH][0])] = cc
colors[int(summary.index[summary.community == SCAFFOLD][0])] = sc

fig, ax2 = plt.subplots(figsize=(16, 6))
ax2.bar(summary.label, summary.mean_rmsf, color=colors, edgecolor="black", lw=1)
hi, lo = val(most_mobile), val(most_rigid)
ax2.annotate("", xy=(0, hi + 0.010), xytext=(len(summary) - 1, lo + 0.010),
             arrowprops=dict(arrowstyle="<->", lw=1.6, color="black"))
ax2.text((len(summary) - 1) / 2, hi + 0.030,
         f"{r_extreme:.1f}× gradient (C{most_mobile} vs C{most_rigid});  {r_scaffold:.1f}× against C{SCAFFOLD}",
         ha="center", fontsize=14, fontweight="bold")
ax2.text((len(summary) - 1) / 2, hi + 0.016,
         "single 100 ns trajectory of the monomeric domain — requires multi-replicate validation",
         ha="center", fontsize=10, color="0.3", style="italic")
ax2.set_title("C. Community-Level Mechanical Hierarchy", fontweight="bold", loc="left")
ax2.set_ylabel("Mean RMSF (nm)"); ax2.set_ylim(0, hi * 1.26)
plt.xticks(rotation=45, ha="right")
plt.tight_layout(); plt.savefig(RESULTS / "figures" / "panels" / "FIGURE_4_rmsf_panelC.png", dpi=300, bbox_inches="tight")
