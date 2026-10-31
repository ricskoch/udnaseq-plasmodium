#!/usr/bin/env python3
"""Fig. 4A - Representative locus: uracil enrichment, VSA genes and predicted G4 sites.

Purpose
    Visual example of how uracil peaks co-locate with var / rif-stevor genes and G4s.

Method
    Merged 3D7 log2(U-enrichment), GEO 10 bp track; red = positive, blue = negative.
    Gene tracks: var (green), rif/stevor (blue); G4 predictions (red; strong sites
    |G4Hunter| > 1.2 and <= 200 bp drawn full height, others half height).
    Window recovered from the IGV export igv_alltracks_chr7.svg (Drive: Figures/): ruler ticks every 10 kb
    (47.44 px) and track baselines spanning 459 px give chr7:507,300-604,000 (+-300 bp).
    Override:  python run.py Pf3D7_07_v3 505000 605000
Outputs  output/fig4A_<chrom>_<start>_<end>.png/svg, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyBigWig  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report, savefig  # noqa: E402

OUT = HERE / "output"
DEFAULT = ("Pf3D7_07_v3", 507_300, 604_000)  # recovered from igv_alltracks_chr7.svg (tick spacing and track extent)


def main(chrom, start, end):
    OUT.mkdir(exist_ok=True)
    bw = pyBigWig.open(str(io.geo_path(config.GEO_MERGED_10BP)))
    v = np.asarray(bw.values(chrom, start, end, numpy=True))
    x = np.arange(start, end) / 1e6
    var = io.read_bed(config.VAR_BED)
    rs = io.read_bed(config.RIFSTEVOR_BED)
    g4 = pd.read_csv(config.G4_SITES, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"])
    sel = lambda d: d[(d.chrom == chrom) & (d.end > start) & (d.start < end)]
    var, rs, g4 = sel(var), sel(rs), sel(g4)

    fig, axes = plt.subplots(4, 1, figsize=(9, 3.6), sharex=True, gridspec_kw={"height_ratios": [4, 0.6, 0.6, 0.6]})
    a = axes[0]
    a.fill_between(x, 0, np.where(v > 0, v, 0), color="#d62728", lw=0)
    a.fill_between(x, 0, np.where(v < 0, v, 0), color="#1f77b4", lw=0)
    a.axhline(config.U_THRESHOLD, c="k", ls=":", lw=0.6)
    a.set_ylabel("log2(U-enrichment)")
    a.set_title(f"{chrom}:{start:,}-{end:,} (merged 3D7, 10 bp)", fontsize=9)
    for ax, d, col, lab in [(axes[1], var, "#2ca02c", "var"), (axes[2], rs, "#1f77b4", "rif/stevor")]:
        for r in d.itertuples():
            ax.add_patch(plt.Rectangle((r.start / 1e6, 0.1), (r.end - r.start) / 1e6, 0.8, color=col))
            ax.text((r.start + r.end) / 2e6, 1.05, ("◀" if r.strand == "-" else "▶"), ha="center", fontsize=5, color=col)
        ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=8)
    for r in g4.itertuples():
        strong = abs(r.score) > config.G4_MIN_SCORE and (r.end - r.start) <= config.G4_MAX_SIZE
        axes[3].add_patch(plt.Rectangle((r.start / 1e6, 0.1 if strong else 0.3), max(r.end - r.start, 150) / 1e6,
                                        0.8 if strong else 0.4, color="#d62728" if strong else "#f4a6a6"))
    axes[3].set_ylabel("G4", rotation=0, ha="right", va="center", fontsize=8)
    for ax in axes[1:]:
        ax.set_ylim(0, 1.4)
        ax.set_yticks([])
        ax.spines[["left", "bottom"]].set_visible(False)
    axes[-1].spines["bottom"].set_visible(True)
    axes[-1].set_xlabel("Position (Mb)")
    axes[-1].set_xlim(start / 1e6, end / 1e6)
    savefig(fig, OUT, f"fig4A_{chrom}_{start}_{end}")

    rep = Report("Fig. 4A - representative locus", "Locus view of merged 3D7 log2(U-enrichment) with VSA genes and G4s.")
    rep.check("locus", "chr7 ~507.3-604.0 kb (from igv_alltracks_chr7.svg, +-300 bp)", f"{chrom}:{start:,}-{end:,}", "INFO",
              f"{len(var)} var, {len(rs)} rif/stevor, {len(g4)} G4 sites in window")
    rep.note("The legacy repository has no file for this panel. The window is read off the IGV export igv_alltracks_chr7.svg "
             "(ruler 510-600 kb ticks and the 459 px track width); IGV shows per-sample IP/INPUT coverage there, not the merged log2 track drawn here.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], int(args[1]), int(args[2])) if len(args) == 3 else main(*DEFAULT)
