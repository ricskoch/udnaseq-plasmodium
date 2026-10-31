#!/usr/bin/env python3
"""00b - How good is "average the GEO 1 bp log2 track" as a stand-in for a binned log2 ratio?

Purpose
    Per-sample 10 bp tracks do not exist (only 1 bp tracks on GEO). The per-dataset analyses (3B, 4C, 3rd-party
    reuse) therefore average the 1 bp log2 tracks into 10 bp bins. Before trusting that, measure the error where
    the exact answer exists: at 100 bp, the exact value is log2((IP+1)/(INPUT+1)) of the RPGC coverage
    (analyses/00_tracks); the approximation is the mean of the GEO 1 bp log2 values inside the same bin.
    mean(log2 x) != log2(mean x), so some bias is expected; the question is how large it is relative to the signal.
Method
    For each of 3d7_rep1, 3d7_rep2, DMSO, RAPA: approx100 = mean of 1 bp GEO values over covered bases per 100 bp bin.
    For 3d7_merged: GEO merged 10 bp track averaged to 100 bp vs the exact (read-share weighted) 100 bp merged track.
    Compare on non-blacklisted bins: Pearson/Spearman, bias (approx - exact), RMSE, 99th percentile |error|, slope,
    and whether the log2U > 0.5 hotspot call (the threshold used in Fig. 4C) agrees.
Outputs  output/approximation_error.csv, approximation_by_signal_level.csv, approximation_error.png/svg, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report, savefig  # noqa: E402

OUT = HERE / "output"


def approx_100bp(path, T):
    """Mean of the GEO track over covered bases in the 100 bp bins of T (NaN if no covered base)."""
    pb = io.read_bigwig(path)
    out = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        out[m] = io.bin_mean(pb[str(c)], config.BIN)[T["start"][m] // config.BIN]
    return out


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("00b - 1 bp to 100 bp averaging error",
                 "Exact 100 bp log2U (from RPGC coverage) vs the mean of the GEO log2 track inside the same bin. "
                 "Non-blacklisted bins.")
    T = io.load_tracks_npz("tracks_100bp")
    rows, strat = [], []
    fig, axs = plt.subplots(1, 5, figsize=(15, 3.2), sharex=True, sharey=True)
    for ax, s in zip(axs, config.SAMPLE_ORDER + ["3d7_merged"]):
        f = config.GEO_MERGED_10BP if s == "3d7_merged" else config.SAMPLES[s]["geo_1bp"]
        a = approx_100bp(io.geo_path(f), T)
        e = T[f"log2_{s}"]
        ok = (~T["blacklisted"]) & ~np.isnan(a)
        n_nan = int(((~T["blacklisted"]) & np.isnan(a)).sum())
        x, y = e[ok], a[ok]
        err = y - x
        fit = stats.linregress(x, y)
        hx, hy = x > config.U_THRESHOLD, y > config.U_THRESHOLD
        rows.append(dict(dataset=s, source="GEO 10 bp merged" if s == "3d7_merged" else "GEO 1 bp", n_bins=len(x),
                         bins_without_geo_value=n_nan,
                         pearson=np.corrcoef(x, y)[0, 1], spearman=stats.spearmanr(x, y).statistic,
                         bias_mean=err.mean(), bias_median=np.median(err), rmse=np.sqrt(np.mean(err ** 2)),
                         p99_abs_err=np.quantile(np.abs(err), .99), max_abs_err=np.abs(err).max(),
                         slope=fit.slope, intercept=fit.intercept, sd_exact=x.std(),
                         rmse_over_sd=np.sqrt(np.mean(err ** 2)) / x.std(),
                         hot_exact=int(hx.sum()), hot_approx=int(hy.sum()), hot_both=int((hx & hy).sum()),
                         hot_jaccard=(hx & hy).sum() / max((hx | hy).sum(), 1)))
        q = pd.qcut(x, 5, duplicates="drop")
        for k, g in pd.DataFrame({"x": x, "y": y, "q": q}).groupby("q", observed=True):
            strat.append(dict(dataset=s, exact_range=str(k), n=len(g), exact_mean=g.x.mean(), approx_mean=g.y.mean(),
                              bias=(g.y - g.x).mean()))
        ax.scatter(x, y, s=1, lw=0, c="grey", alpha=.25, rasterized=True)
        ax.plot([-1.5, 2], [-1.5, 2], c="r", lw=.8)
        ax.set_title(f"{s}\nr={rows[-1]['pearson']:.3f}, bias={rows[-1]['bias_mean']:+.3f}", fontsize=8)
        ax.set_xlim(-1.5, 2); ax.set_ylim(-1.5, 2); ax.set_xlabel("exact 100 bp log2U")
    axs[0].set_ylabel("mean of GEO track in bin")
    savefig(fig, OUT, "approximation_error")
    df = pd.DataFrame(rows)
    df.round(4).to_csv(OUT / "approximation_error.csv", index=False)
    pd.DataFrame(strat).round(4).to_csv(OUT / "approximation_by_signal_level.csv", index=False)
    for r in rows:
        rep.check(f"{r['dataset']} ({r['source']}): approx vs exact", "exact RPGC log2 ratio",
                  f"r = {r['pearson']:.4f}, rho = {r['spearman']:.4f}", "INFO",
                  f"bias {r['bias_mean']:+.4f} (median {r['bias_median']:+.4f}); RMSE {r['rmse']:.4f} = "
                  f"{r['rmse_over_sd']:.2f} SD of signal; 99th pct |err| {r['p99_abs_err']:.3f}; slope {r['slope']:.3f}; "
                  f"hotspots (>0.5) exact {r['hot_exact']}, approx {r['hot_approx']}, shared {r['hot_both']} "
                  f"(Jaccard {r['hot_jaccard']:.2f}); {r['bins_without_geo_value']} bins without GEO value")
    rep.note("Error is a per-bin systematic shrinkage/offset plus noise; see approximation_by_signal_level.csv for how it "
             "depends on signal level.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
