#!/usr/bin/env python3
"""Fig. 4B per dataset - GC% vs log2(U-enrichment) for 3d7_rep1, 3d7_rep2, DMSO, RAPA, merged.

Purpose
    Check whether the GC-uracil correlation (r ~ 0.42-0.44 in the merged track, which is ~89 % rep1)
    is reproduced by the independent replicate and by the DMSO/RAPA samples.
Method
    Same as analyses/fig4B_gc_correlation, but per sample: 100 bp bins, blacklist removed, bins with
    log2 == 0 exactly dropped. GC% from the genome FASTA. Pearson and Spearman with 95 % bootstrap CI
    (1000 resamples of bins, fixed seed; bins are autocorrelated, so CIs are optimistic) and the slope per 10 GC% points.
Outputs  output/per_dataset_stats.csv, fig4B_per_dataset.png/svg, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE.parents[1] / "analyses" / "fig4B_gc_correlation"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report, savefig  # noqa: E402
from run import read_fasta, gc_fraction  # noqa: E402  (fig4B_gc_correlation/run.py)

OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
N_BOOT, SEED = 1000, 20240209


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 4B per dataset - GC% vs log2(U-enrichment)",
                 "Does the GC-uracil correlation hold in each dataset separately? 100 bp bins, blacklist removed.")
    T = io.load_tracks_npz("tracks_100bp")
    seqs = read_fasta(config.FASTA)
    gc = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gc[m] = gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    rng = np.random.default_rng(SEED)
    rows, fig, axs = [], *plt.subplots(1, len(SETS), figsize=(3 * len(SETS), 3), sharex=True, sharey=True)
    for ax, s in zip(axs, SETS):
        y = T[f"log2_{s}"]
        keep = (~T["blacklisted"]) & (y != 0)
        x, v = gc[keep], y[keep]
        rp, rs = stats.pearsonr(x, v).statistic, stats.spearmanr(x, v).statistic
        fit = stats.linregress(x * 100, v)
        idx = rng.integers(0, len(x), (N_BOOT, len(x)))
        bp = np.array([np.corrcoef(x[i], v[i])[0, 1] for i in idx])
        rows.append(dict(dataset=s, n_bins=len(x), pearson=rp, pearson_ci_lo=np.quantile(bp, .025), pearson_ci_hi=np.quantile(bp, .975),
                         spearman=rs, slope_log2U_per_10GCpct=fit.slope * 10, sd_log2U=v.std(), frac_log2U_gt_0p5=np.mean(v > config.U_THRESHOLD)))
        ax.scatter(x * 100, v, s=1, lw=0, c="grey", alpha=.3, rasterized=True)
        xx = np.array([0, 80]); ax.plot(xx, fit.intercept + fit.slope * xx, c="r", lw=1)
        ax.set_title(f"{s}\nr={rp:.2f} ρ={rs:.2f}", fontsize=8); ax.set_xlim(0, 80); ax.set_ylim(-1.5, 2)
        ax.set_xlabel("GC%")
    axs[0].set_ylabel("log2(U-enrichment)")
    savefig(fig, OUT, "fig4B_per_dataset")
    df = pd.DataFrame(rows)
    df.round(4).to_csv(OUT / "per_dataset_stats.csv", index=False)
    for r in rows:
        rep.check(f"{r['dataset']}: Pearson r (95% CI)", "merged: 0.436", f"{r['pearson']:.3f} ({r['pearson_ci_lo']:.3f}, {r['pearson_ci_hi']:.3f})",
                  "INFO", f"Spearman {r['spearman']:.3f}; slope {r['slope_log2U_per_10GCpct']:+.3f} log2U per 10 GC%; n = {r['n_bins']:,}; "
                  f"{r['frac_log2U_gt_0p5']:.2%} bins log2U > 0.5")
    rep.note("Bootstrap over bins ignores spatial autocorrelation; read the CI as a lower bound on uncertainty.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
