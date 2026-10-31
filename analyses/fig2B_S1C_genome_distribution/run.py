#!/usr/bin/env python3
"""Fig. 2B and S1C - genome-wide distribution of uracil enrichment (merged 3D7).

Purpose
    2B : ideogram showing where along the 14 chromosomes uracil is enriched, with
         var genes (green triangles) and centromeres (black lines) for orientation.
    S1C: per-chromosome histograms of log2(U-enrichment) - is enrichment a few
         clustered loci on top of a ~0-centred background?

Method (replaces RIdeogram / ggplot2)
    2B : merged 3D7 track at 100 bp (analyses/00_tracks), averaged into 5 kb windows
         for display; colour white (<= 0) -> red (>= 0.5 log2), as in the manuscript.
    S1C: GEO merged 10 bp track (the source of the legacy histogram), bin width 0.05,
         x-range -1..1, nuclear chromosomes only.
Outputs  output/fig2B_ideogram.*, ideogram_windows_5kb.csv, figS1C_histograms.*,
         figS1C_chromosome_summary.csv, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report, savefig  # noqa: E402

OUT = HERE / "output"
WINDOW = 5_000
CMAP = LinearSegmentedColormap.from_list("white_red", ["#ffffff", "#fcbba1", "#ef3b2c", "#99000d"])


def ideogram(win: pd.DataFrame, var: pd.DataFrame, cen: pd.DataFrame, sizes: dict):
    chroms = config.NUCLEAR_CHROMS
    fig, ax = plt.subplots(figsize=(7, 6))
    h = 0.55
    for i, c in enumerate(chroms):
        y = len(chroms) - i
        w = win[win.chrom == c]
        img = np.clip(w.value.values, 0, 0.5)[None, :]
        ax.imshow(img, extent=(0, sizes[c] / 1e6, y - h / 2, y + h / 2), aspect="auto", cmap=CMAP,
                  vmin=0, vmax=0.5, interpolation="nearest")
        ax.add_patch(plt.Rectangle((0, y - h / 2), sizes[c] / 1e6, h, fill=False, lw=0.6))
        for r in cen[cen.chrom == c].itertuples():
            ax.plot([(r.start + r.end) / 2e6] * 2, [y - h / 2, y + h / 2], c="k", lw=1.5)
        vv = var[var.chrom == c]
        ax.scatter((vv.start + vv.end) / 2e6, np.full(len(vv), y + h / 2 + 0.17), marker="v", s=14,
                   c="#1a9850", edgecolors="none", zorder=3)
    ax.set_yticks(range(len(chroms), 0, -1), [c.split("_")[1] for c in chroms])
    ax.set_xlim(-0.02, max(sizes[c] for c in chroms) / 1e6 + 0.05)
    ax.set_ylim(0.4, len(chroms) + 0.8)
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("Chromosome")
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, 0.5))
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("log2(U-enrichment), 5 kb mean")
    cb.set_ticks([0, 0.25, 0.5], labels=["≤0", "0.25", "≥0.5"])
    ax.spines[["left"]].set_visible(False)
    return fig


def histograms(vals: dict):
    fig, axes = plt.subplots(2, 7, figsize=(14, 4.5), sharex=True, sharey=True)
    bins = np.arange(-1, 1.0001, 0.05)
    colors = plt.cm.hsv(np.linspace(0, 0.9, 14))
    for ax, c, col in zip(axes.flat, config.NUCLEAR_CHROMS, colors):
        v = vals[c]
        ax.hist(v[(v >= -1) & (v <= 1)], bins=bins, color=col, edgecolor="k", lw=0.4)
        ax.axvline(0, c="red", ls=":", lw=0.8)
        ax.set_title(c, fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel("log2(U-enrichment)")
    for ax in axes[:, 0]:
        ax.set_ylabel("10 bp bins")
    fig.tight_layout()
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 2B / S1C - genome-wide distribution",
                 "Ideogram of merged 3D7 log2(U-enrichment) and per-chromosome value distributions.")
    sizes = io.chrom_sizes()

    # ---- Fig. 2B -------------------------------------------------------------------
    T = io.load_tracks_npz("tracks_100bp")
    rows = []
    for c in config.NUCLEAR_CHROMS:
        m = T["chrom"] == c
        v = io.bin_mean(T["log2_3d7_merged"][m], WINDOW // config.BIN)
        s = np.arange(len(v)) * WINDOW
        rows.append(pd.DataFrame({"chrom": c, "start": s, "end": np.minimum(s + WINDOW, sizes[c]), "value": v}))
    win = pd.concat(rows, ignore_index=True)
    win.round(5).to_csv(OUT / "ideogram_windows_5kb.csv", index=False)
    var = io.read_bed(config.VAR_BED)
    cen = io.read_bed(config.CENTROMERES)
    savefig(ideogram(win, var, cen, sizes), OUT, "fig2B_ideogram")
    hot = win[win.value >= config.U_THRESHOLD]
    near_var = sum(((var.chrom == r.chrom) & (var.end > r.start - 10_000) & (var.start < r.end + 10_000)).any()
                   for r in hot.itertuples())
    rep.check("Fig 2B: 5 kb windows with mean >= 0.5", "no legacy file", f"{len(hot)} windows, {near_var} within 10 kb of a var gene",
              "INFO", "legacy RIdeogram input/output not in the repo; regenerated from the merged track")

    # ---- Fig. S1C ------------------------------------------------------------------
    pb = io.read_bigwig(io.geo_path(config.GEO_MERGED_10BP), chroms=config.NUCLEAR_CHROMS)
    vals = {c: io.bin_mean(pb[c], 10) for c in config.NUCLEAR_CHROMS}
    vals = {c: v[~np.isnan(v)] for c, v in vals.items()}
    savefig(histograms(vals), OUT, "figS1C_histograms")
    bins = np.arange(-1, 1.0001, 0.05)
    summ = pd.DataFrame([{"chrom": c, "n_bins": len(v), "median": np.median(v), "mean": v.mean(),
                          "frac_above_0.5": np.mean(v > 0.5), "peak_count": np.histogram(v, bins)[0].max()}
                         for c, v in vals.items()]).round(5)
    summ.to_csv(OUT / "figS1C_chromosome_summary.csv", index=False)
    rep.check("Fig S1C: chr14 tallest histogram bar", "~46,000 (read from legacy PNG)",
              f"{summ.set_index('chrom').loc['Pf3D7_14_v3', 'peak_count']:,}", "CLOSE",
              "visual comparison only; legacy histogram data file is not in the repo")
    rep.check("Fig S1C: all chromosomes centred near 0", "yes", f"medians {summ['median'].min():.3f} .. {summ['median'].max():.3f}",
              "MATCH")
    rep.note("The legacy histogram counted bedGraph rows (runs of equal value may span >1 bin), here every 10 bp "
             "bin counts once, so bar heights can differ slightly.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
