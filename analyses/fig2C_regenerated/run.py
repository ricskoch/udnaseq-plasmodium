#!/usr/bin/env python3
"""Fig 2C regenerated from the pipeline tables.

A  per-gene mean log2(U-enrichment), var / rif / stevor vs other genes, 4 datasets; BH-adjusted Welch p above each bracket (raw, no GC control).
B  the same VSA-vs-other-genes question at equal GC (bin-level contrast, 100 kb block-bootstrap 95 % CI), from fig_gcbias/output/figure_data.csv.
Outputs  output/figure_2C_regenerated.svg/.png, output/figure_data.csv, COMPARISON.md
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent; A = HERE.parent
OUT = HERE / "output"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 8, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "axes.labelsize": 9, "axes.labelcolor": "#333333",
                     "xtick.color": "#333333", "ytick.color": "#333333", "xtick.labelcolor": "#333333", "ytick.labelcolor": "#333333",
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "text.color": "#333333"})
CATS = ["Others", "RIF", "STEVOR", "VAR"]
COL = {"Others": "#7f7f7f", "RIF": "#1f77b4", "STEVOR": "#17becf", "VAR": "#2ca02c"}
DS = [("3d7_rep1", "3D7 rep1"), ("3d7_rep2", "3D7 rep2"), ("DMSO", "DMSO"), ("RAPA", "RAPA")]
t = pd.read_csv(A / "fig2C_3A_S1D_S2_gene_enrichment/output/table_S1_gene_log2U.csv")
st = pd.read_csv(A / "fig2C_3A_S1D_S2_gene_enrichment/output/stats_fig2C_VSA_vs_others_exact_annotation.csv")
eq = pd.read_csv(A / "fig_gcbias/output/figure_data.csv"); eq = eq[eq.contrast == "VSA vs other genic"].set_index("dataset")


def pfmt(p):
    return "p<1e-99" if p < 1e-99 else f"p={p:.0e}".replace("e-0", "e-")


fig = plt.figure(figsize=(8.4, 3.0), dpi=150)
gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1, 0.1, 0.9], wspace=0.12, left=0.09, right=0.985, top=0.86, bottom=0.2)
axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
for i, (ax, (k, lab)) in enumerate(zip(axes, DS)):
    c = f"{k}_exact"
    data = [t.loc[t.category == g, c].dropna().values for g in CATS]
    parts = ax.violinplot(data, showextrema=False, showmedians=True, widths=0.8)
    for b, g in zip(parts["bodies"], CATS):
        b.set_facecolor(COL[g]); b.set_edgecolor("k"); b.set_linewidth(0.5); b.set_alpha(0.7)
    parts["cmedians"].set_color("k"); parts["cmedians"].set_linewidth(0.8)
    for j, fam in enumerate(CATS[1:]):
        p = st[(st.dataset == c) & (st.family == fam)].p_BH.iloc[0]
        y = 0.62 + 0.17 * j
        ax.plot([1, j + 2], [y, y], c="k", lw=0.6)
        ax.text((j + 3) / 2, y + 0.01, pfmt(p), ha="center", va="bottom", fontsize=6.5)
    ax.axhline(0, c="grey", lw=0.5, ls=":")
    ax.set_xticks(range(1, 5)); ax.set_xticklabels([g.lower() if g != "Others" else "other" for g in CATS], fontsize=7.5)
    ax.set_xlabel(lab, labelpad=2); ax.set_ylim(-0.9, 1.25)
    if i: ax.tick_params(labelleft=False); ax.spines["left"].set_visible(False); ax.tick_params(left=False)
axes[0].set_ylabel("Mean log2(U-enrichment) per gene", fontsize=8.5)
axb = fig.add_subplot(gs[0, 5])
y = np.arange(4)[::-1]
est = [eq.loc[k, "estimate"] for k, _ in DS]; lo = [eq.loc[k, "lo"] for k, _ in DS]; hi = [eq.loc[k, "hi"] for k, _ in DS]
axb.errorbar(est, y, xerr=[np.array(est) - lo, np.array(hi) - est], fmt="o", color="#222222", ms=4, lw=1, capsize=2)
axb.axvline(0, color="grey", lw=0.6, ls=":")
axb.set_yticks(y); axb.set_yticklabels([l for _, l in DS]); axb.set_ylim(-0.6, 3.6)
axb.set_xlabel("VSA - other genes,\nequal GC (log2U)", labelpad=2)
axa0 = axes[0].get_position(); axb0 = axb.get_position()
fig.text(0.0, axa0.y1 + 0.11, "C", ha="left", va="top", fontsize=22, color="#1b1918")
fig.text(axb0.x0 - 0.1, axa0.y1 + 0.11, "D", ha="left", va="top", fontsize=22, color="#1b1918")
for ext in ("svg", "png"):
    fig.savefig(OUT / f"figure_2C_regenerated.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
s = OUT / "figure_2C_regenerated.svg"; s.write_text(s.read_text().replace("Liberation Sans", "Arial"))
rows = []
for k, _ in DS:
    for fam in CATS[1:]:
        r = st[(st.dataset == f"{k}_exact") & (st.family == fam)].iloc[0]
        rows.append(dict(dataset=k, family=fam, n_family=r.n_family, n_others=r.n_others, mean_diff_raw=r.mean_diff, p_BH=r.p_BH,
                         equalGC_VSA_vs_other=eq.loc[k, "estimate"], equalGC_lo=eq.loc[k, "lo"], equalGC_hi=eq.loc[k, "hi"]))
d = pd.DataFrame(rows); d.to_csv(OUT / "figure_data.csv", index=False)
(HERE / "COMPARISON.md").write_text("# Fig 2C regenerated\n\nOriginal: violins per dataset with significance stars (Welch, BH). Regenerated: same data from table_S1_gene_log2U.csv (exact annotation), BH p values printed, plus equal-GC VSA contrast (panel D) from fig_gcbias.\n\n" + d.to_markdown(index=False, floatfmt=".3g") + "\n")
