#!/usr/bin/env python3
"""Revised Figure 4: GC bias and equal-GC contrasts (replaces old Fig 4A-C).

Purpose  Show (A) IP vs INPUT coverage depends on GC differently, (B) log2U rises with GC, (C) effects that remain at equal GC.
Inputs   qc_library_bias/output/coverage_by_gc.csv, fig4B_gc_correlation/output/gc_vs_log2U_100bp.csv.gz (merged track),
         genic_vs_intergenic/output, subtelomere_position/output (all produced by run_all.py).
Output   output/figure_gcbias.svg (+ png), COMPARISON.md
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
A = HERE.parent
OUT = HERE / "output"; OUT.mkdir(exist_ok=True)
# Style follows the original figures: Arial (Liberation Sans is metric-identical and is renamed to Arial in the SVG),
# plain uppercase panel letters (not bold, #1b1918) at the top-left of each panel, black 0.8 pt axes, dark-grey text, no titles/text boxes.
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 8, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "axes.labelsize": 9, "axes.labelcolor": "#333333",
                     "xtick.color": "#333333", "ytick.color": "#333333", "xtick.labelcolor": "#333333", "ytick.labelcolor": "#333333",
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "text.color": "#333333"})


def letter(fig, ax, txt, x_in):
    """Panel letter: plain Arial 22 pt, top of the letter level with the top of the plot; x_in = distance from the figure's left edge (inches)."""
    pos = ax.get_position()
    fig.text(x_in / fig.get_figwidth(), pos.y1, txt, ha="left", va="top", fontsize=22, fontweight="normal", color="#1b1918")


COL = {"3d7_rep1": "#1b6ca8", "3d7_rep2": "#e07b00", "DMSO": "#6a994e", "RAPA": "#b5179e", "3d7_merged": "#444444"}
LAB = {"3d7_rep1": "3D7 rep 1", "3d7_rep2": "3D7 rep 2", "DMSO": "DMSO", "RAPA": "RAPA", "3d7_merged": "3D7 merged"}
ORDER = ["3d7_rep1", "3d7_rep2", "DMSO", "RAPA"]

fig = plt.figure(figsize=(7.2, 5.0), dpi=150)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.6, wspace=1.0, left=0.17, right=0.98, top=0.95, bottom=0.1)
gsc = fig.add_gridspec(1, 1, left=0.29, right=0.98, top=0.95 - 0.95 * 0.0, bottom=0.1)
axa, axb, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1:]), None

# A: IP / INPUT mean coverage by GC class
cg = pd.read_csv(A / "qc_library_bias/output/coverage_by_gc.csv")
for s in ORDER:
    d = cg[cg["sample"] == s]
    axa.plot(range(len(d)), (d.ip_mean / d.input_mean).values, "-o", ms=3, color=COL[s], label=LAB[s])
axa.axhline(1, color="grey", lw=0.5, ls=":")
axa.set_xticks(range(5)); axa.set_xticklabels(["<15", "15-25", "25-35", "35-50", ">50"], rotation=30)
axa.set_xlabel("100 bp bin GC (%)"); axa.set_ylabel("mean IP / mean INPUT\ncoverage"); axa.legend(frameon=False, loc="upper left")


# B: merged log2U vs GC
g = pd.read_csv(A / "fig4B_gc_correlation/output/gc_vs_log2U_100bp.csv.gz")
g["gcb"] = (g.gc * 100).round().clip(0, 60)
m = g.groupby("gcb").log2U.agg(["mean", "sem", "size"]); m = m[m["size"] >= 20]
axb.errorbar(m.index, m["mean"], yerr=1.96 * m["sem"], fmt="o", ms=3, color=COL["3d7_merged"], lw=0.6)
axb.axhline(0, color="grey", lw=0.5, ls=":")
axb.set_xlabel("100 bp bin GC (%)"); axb.set_ylabel("mean log2(U-enrichment)\n(3D7 merged)")


# C: forest of equal-GC contrasts
gv = pd.read_csv(A / "genic_vs_intergenic/output/genic_contrast.csv").set_index("dataset")
st = pd.read_csv(A / "subtelomere_position/output/subtelomere_contrast.csv")
rows = [("Genic vs intergenic", lambda s: (gv.loc[s, "equalGC_genic_minus_inter"], gv.loc[s, "B_lo"], gv.loc[s, "B_hi"])),
        ("VSA vs other genic", lambda s: (gv.loc[s, "equalGC_VSA_minus_otherGenic"], gv.loc[s, "C_lo"], gv.loc[s, "C_hi"])),
        ("VSA, subtelomeric", lambda s: tuple(st[(st.dataset == s) & (st.group == "VSA/subtelomeric")][["equalGC_vs_other_internal", "lo", "hi"]].iloc[0])),
        ("VSA, internal", lambda s: tuple(st[(st.dataset == s) & (st.group == "VSA/internal")][["equalGC_vs_other_internal", "lo", "hi"]].iloc[0])),
        ("Non-VSA subtelomeric", lambda s: tuple(st[(st.dataset == s) & (st.group == "other/subtelomeric")][["equalGC_vs_other_internal", "lo", "hi"]].iloc[0]))]
axc = fig.add_subplot(gs[1, :]); axc.set_position([0.29, gs[1, :].get_position(fig).y0, 0.98 - 0.29, gs[1, :].get_position(fig).height])
y, ticks, labs, data = 0, [], [], []
for name, f in rows:
    for i, s in enumerate(ORDER):
        v, lo, hi = f(s)
        axc.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", ms=3.5, color=COL[s], lw=1, capsize=0, label=LAB[s] if name == rows[0][0] else None)
        data.append(dict(contrast=name, dataset=s, estimate=v, lo=lo, hi=hi)); y += 1
    ticks.append(y - 2.5); labs.append(name); y += 1
axc.axvline(0, color="grey", lw=0.6); axc.set_yticks(ticks); axc.set_yticklabels(labs); axc.invert_yaxis()
axc.set_xlabel("difference in log2(U-enrichment) at equal GC (95 % CI, 100 kb block bootstrap)")
axc.legend(frameon=False, ncol=4, loc="lower right", bbox_to_anchor=(1.0, 1.0))

letter(fig, axa, "A", 0.05); letter(fig, axb, "B", axb.get_position().x0 * fig.get_figwidth() - 1.0); letter(fig, axc, "C", 0.05)
fig.savefig(OUT / "figure_gcbias.svg", bbox_inches="tight"); fig.savefig(OUT / "figure_gcbias.png", bbox_inches="tight", dpi=600)
svg = (OUT / "figure_gcbias.svg"); svg.write_text(svg.read_text().replace("Liberation Sans", "Arial"))
pd.DataFrame(data).round(4).to_csv(OUT / "figure_data.csv", index=False)
(HERE / "COMPARISON.md").write_text("# Revised Figure 4 (GC bias)\n\nPanel data in output/figure_data.csv; figure_gcbias.svg is the vector file.\n\n" + pd.DataFrame(data).round(3).to_markdown(index=False) + "\n")
