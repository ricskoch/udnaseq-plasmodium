#!/usr/bin/env python3
"""Supplementary items for the revised manuscript.

Fig S4  RAPA vs DMSO: paired gene differences (RAPA - DMSO) for var / rif / stevor, raw and GC-corrected (from gene_ranking_go).
Table S2  G4 proximity of the top-4,034 hotspots: observed vs random and GC-matched null (from fig4C_gc_matched_null).
Outputs  output/figure_s4_rapa_dmso.svg, output/table_s2_g4.csv, COMPARISON.md
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent; A = HERE.parent
OUT = HERE / "output"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
d = pd.read_csv(A / "gene_ranking_go/output/differential_summary.csv")
d = d[(d.a == "RAPA") & (d.b == "DMSO")].set_index("score")
SC = [("raw", "raw", "#888888"), ("gc_bin", "GC-normalised (per track)", "#1b6ca8"), ("gc_gene", "gene-GC residual", "#e07b00")]
fig, ax = plt.subplots(figsize=(4.6, 2.8), dpi=150)
w = 0.26
for i, (k, lab, c) in enumerate(SC):
    v = [d.loc[k, f"{f}_mean_diff"] for f in ("VAR", "RIF", "STEVOR")]
    ax.bar(np.arange(3) + (i - 1) * w, v, w, color=c, label=lab)
ax.axhline(0, color="k", lw=0.6); ax.set_xticks(range(3)); ax.set_xticklabels(["var", "rif", "stevor"])
ax.set_ylabel("RAPA - DMSO\nmean log2(U-enrichment)"); ax.legend(frameon=False, fontsize=6)
ax.set_title("Fig S4  var/rif difference between RAPA and DMSO is GC-linked", loc="left", fontsize=7.5, fontweight="bold")
fig.savefig(OUT / "figure_s4_rapa_dmso.svg", bbox_inches="tight"); fig.savefig(OUT / "figure_s4_rapa_dmso.png", bbox_inches="tight", dpi=600)
g = pd.read_csv(A / "fig4C_gc_matched_null/output/matched_null_summary.csv")
t = g[["dataset", "matching", "n_hotspots", "obs_frac_within_100bp", "null_frac_within_100bp", "fold_frac_within_100bp", "p_frac_within_100bp"]].copy()
t["matching"] = t.matching.map({"gc10": "random bins, 10 bp GC matched", "gc10+500": "10 bp + 500 bp window GC matched"})
t.columns = ["dataset", "null model", "hotspots (n)", "observed within 100 bp of G4", "null expectation", "fold", "p (permutation)"]
t.round(4).to_csv(OUT / "table_s2_g4.csv", index=False)
(HERE / "COMPARISON.md").write_text("# Supplementary items\n\nFig S4 data (RAPA - DMSO):\n\n" + d[["VAR_mean_diff", "RIF_mean_diff", "STEVOR_mean_diff"]].round(4).to_markdown() + "\n\nTable S2:\n\n" + t.round(4).to_markdown(index=False) + "\n")
