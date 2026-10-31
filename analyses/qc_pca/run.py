#!/usr/bin/env python3
"""QC - PCA of the eight RPGC coverage tracks (not a manuscript figure).

Purpose
    Exploratory check of whether samples separate by library type (IP vs INPUT)
    or by experiment/batch (3D7 rep1, rep2, PfUNG-iKO DMSO/RAPA).

Method (identical to legacy/pca_analysis/pca_analysis.ipynb)
    features = all 100 bp bins of the RPGC coverage tracks (samples = rows),
    z-scored per bin across the 8 samples (StandardScaler), PCA keeping the
    components that explain >= 80 % of the variance.
Outputs
    output/pca_results.csv, pca_variance.csv, pca_plot.png/svg, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report, savefig, status  # noqa: E402

OUT = HERE / "output"
MARKERS = {"3d7_rep1": "o", "3d7_rep2": "s", "DMSO": "^", "RAPA": "D"}
COLORS = {"INPUT": "C0", "UDNAIP": "C1"}


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("QC - PCA of RPGC coverage tracks",
                 "Principal component analysis of the 8 IP/INPUT RPGC coverage tracks (100 bp bins).")
    T = io.load_tracks_npz("tracks_100bp")
    data = {}
    for s in config.SAMPLE_ORDER:
        data[f"{s}_INPUT"] = T[f"input_{s}"]
        data[f"{s}_UDNAIP"] = T[f"ip_{s}"]
    keys = sorted(data)
    X = StandardScaler().fit_transform(np.vstack([data[k] for k in keys]))
    full = PCA().fit(X)
    n = int(np.argmax(np.cumsum(full.explained_variance_ratio_) >= 0.80) + 1)
    pca = PCA(n_components=n)
    S = pca.fit_transform(X)
    res = pd.DataFrame(S, columns=[f"PC{i + 1}" for i in range(n)]).assign(Sample=keys)[["Sample"] + [f"PC{i + 1}" for i in range(n)]]
    res.round(6).to_csv(OUT / "pca_results.csv", index=False)
    pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(full.explained_variance_ratio_))],
                  "explained_variance_ratio": full.explained_variance_ratio_}).to_csv(OUT / "pca_variance.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 4.2))
    for _, r in res.iterrows():
        exp, cond = r.Sample.rsplit("_", 1)
        ax.scatter(r.PC1, r.PC2, marker=MARKERS[exp], c=COLORS[cond], s=90, edgecolors="k", label=f"{exp} {cond}")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, OUT, "pca_plot")

    leg = pd.read_csv(config.LEGACY / "pca_analysis" / "pca_results.csv").set_index("Sample")
    new = res.set_index("Sample").reindex(leg.index)
    for pc in leg.columns:
        sign = np.sign(np.dot(leg[pc], new[pc]))  # PCA sign is arbitrary
        d = float(np.max(np.abs(leg[pc] - sign * new[pc])))
        rep.check(f"{pc} scores", f"{leg[pc].abs().max():.1f} (max |score|)", f"max |diff| = {d:.2e}",
                  status(d / leg[pc].abs().max()), "sign-aligned")
    rep.check("PCs kept (>= 80 % variance)", 3, n, "MATCH" if n == 3 else "DIFFERS")
    rep.note("Explained variance: " + ", ".join(f"PC{i + 1} {v:.1%}" for i, v in enumerate(full.explained_variance_ratio_[:4])))
    rep.note("PC1 separates the 3D7 libraries (both INPUTs and rep1 IP, positive) from the PfUNG-iKO IP libraries "
             "(negative); PC2 isolates the PfUNG-iKO INPUT libraries; PC3 is driven by 3d7_rep2 IP alone. The structure "
             "mixes library type (IP/INPUT) with experiment batch (3D7 vs iKO), so it is a QC view, not evidence for "
             "uracil biology. The interpretation text in the legacy notebook (e.g. 'genomic rearrangement in "
             "knockout parasites') is not supported by this analysis and was not carried over.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
