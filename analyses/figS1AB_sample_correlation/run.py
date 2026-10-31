#!/usr/bin/env python3
"""Fig. S1A/B - Correlation of log2(U-enrichment) between U-DNA-Seq datasets.

Purpose
    Quantify replicate agreement (3D7 rep1 vs rep2) and DMSO vs RAPA similarity
    ("only a moderate correlation between biological replicates"; "DMSO/RAPA
    correlate more strongly than the 3D7 replicates").

Method (replaces deepTools multiBigwigSummary + plotCorrelation)
    Input: GEO per-sample 1 bp log2 tracks (the legacy matrices were built from these).
    Three region sets, each with Spearman and Pearson correlation:
      genome_100bp  100 bp bins, blacklisted bins removed
      genome_10kb   10 kb bins, whole genome
      var_genes     the 62 var genes (data/BED_files/var_genes/var_genes.bed)
    Two ways of averaging the signal per region:
      deeptools  pyBigWig zoom-level summaries, as deepTools does -> reproduces the legacy figures
      exact      true per-base mean (recommended; zoom summaries are approximations)
Outputs
    output/correlation_<set>_<method>.csv/.png/.svg, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402
from udnaseq.report import Report, savefig, status  # noqa: E402

OUT = HERE / "output"
LEG = config.LEGACY / "heatmap_correlation"
LEGACY_NPZ = {"genome_100bp": "heatmap_data_all_bin100bp_blacklist.npz",
              "genome_10kb": "heatmap_data_all.npz",
              "var_genes": "heatmap_data_var.npz"}


def corr(X: np.ndarray, method: str) -> np.ndarray:
    X = X[~np.isnan(X).any(axis=1)]
    return spearmanr(X).statistic if method == "spearman" else np.corrcoef(X.T)


def plot_corr(C: pd.DataFrame, title: str):
    n = len(C)
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    im = ax.imshow(C.values, cmap="RdYlBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(n), C.columns, rotation=45, ha="right")
    ax.set_yticks(range(n), C.index)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{C.values[i, j]:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046)
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. S1A/B - sample correlation",
                 "Spearman/Pearson correlation of log2(U-enrichment) between 3d7_rep1, 3d7_rep2, DMSO and RAPA.")
    labels = config.SAMPLE_ORDER

    # --- region sets ---------------------------------------------------------
    b100 = genome.add_blacklist_flag(io.bins_frame(100))
    regions = {"genome_100bp": b100[~b100.blacklisted].reset_index(drop=True),
               "genome_10kb": io.bins_frame(10_000),
               "var_genes": io.read_bed(config.VAR_BED)}
    paths = {s: io.geo_path(config.SAMPLES[s]["geo_1bp"]) for s in labels}

    for name, reg in regions.items():
        leg = np.load(LEG / LEGACY_NPZ[name], allow_pickle=True)
        L = leg["matrix"]
        assert list(leg["labels"]) == labels, leg["labels"]
        X = {mode: np.column_stack([io.bigwig_stats(paths[s], reg, exact=(mode == "exact")) for s in labels])
             for mode in ("deeptools", "exact")}
        Xd = X["deeptools"]
        # legacy rows are in deepTools' parallel-chunk order -> compare sorted columns (same multiset of values)
        d = max(float(np.max(np.abs(np.sort(Xd[~np.isnan(Xd[:, j]), j]) - np.sort(L[~np.isnan(L[:, j]), j]))))
                if (~np.isnan(Xd[:, j])).sum() == (~np.isnan(L[:, j])).sum() else np.inf for j in range(4))
        rep.check(f"{name}: region values (deepTools mode)", f"{L.shape[0]} rows", f"{Xd.shape[0]} rows",
                  status(d), f"max |diff| of sorted columns = {d:.1e}")
        de = float(np.nanmax(np.abs(X["exact"] - Xd)))
        rep.check(f"{name}: exact vs zoom-level means", "", f"max |diff| = {de:.3f}", "INFO",
                  "error introduced by deepTools' approximate bigWig summaries")
        if name != "genome_100bp":  # per-bin table for 227k bins would be ~18 MB
            pd.DataFrame(np.column_stack([reg[["chrom", "start", "end"]].values, Xd, X["exact"]]),
                         columns=["chrom", "start", "end"] + [f"{s}_deeptools" for s in labels] +
                         [f"{s}_exact" for s in labels]).to_csv(OUT / f"values_{name}.csv.gz", index=False)
        for method in ("spearman", "pearson"):
            Cl = corr(L, method)
            iu = np.triu_indices(4, 1)
            for mode in ("deeptools", "exact"):
                C = pd.DataFrame(corr(X[mode], method), index=labels, columns=labels)
                C.round(4).to_csv(OUT / f"correlation_{name}_{method}_{mode}.csv")
                savefig(plot_corr(C, f"{method.capitalize()} - {name.replace('_', ' ')} ({mode})"), OUT,
                        f"correlation_{name}_{method}_{mode}")
                d = float(np.max(np.abs(Cl - C.values)))
                rep.check(f"{name}: {method} ({mode})", " ".join(f"{v:.2f}" for v in Cl[iu]),
                          " ".join(f"{v:.2f}" for v in C.values[iu]),
                          status(d, 1e-3, 0.03) if mode == "deeptools" else "INFO", f"max |diff| = {d:.3f}")

    rep.note("Upper-triangle order: rep1-rep2, rep1-DMSO, rep1-RAPA, rep2-DMSO, rep2-RAPA, DMSO-RAPA.")
    rep.note("The legacy .npz matrices (and therefore the published S1 heatmaps) were computed from the GEO 1 bp "
             "tracks, not from the 100 bp MAPQ-filtered tracks referenced in deepTools_correlation_plotting.md.")
    rep.note("Conclusions are unchanged with exact means: replicate correlation is moderate, DMSO-RAPA is the "
             "highest pair.")
    rep.note("DMSO/RAPA labels are in the correct order in the legacy files (the documented multiBigwigSummary "
             "command lists RAPA before DMSO but the matrices agree with the labels).")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
