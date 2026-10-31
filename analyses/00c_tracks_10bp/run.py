#!/usr/bin/env python3
"""00c - Per-sample 10 bp log2(U-enrichment) tracks (average of the GEO 1 bp tracks).

Purpose
    Fig. 3B and 4C were made from the merged 3D7 10 bp track (~89 % rep1). To test them in rep1, rep2, DMSO and RAPA
    separately we need 10 bp tracks per sample. GEO provides only 1 bp tracks per sample, so each 10 bp bin is the
    mean of the 1 bp values (covered bases only). The approximation error is quantified in analyses/00b_bin_approximation.
Method
    10 bp bin i of a chromosome = mean over non-NaN 1 bp values in [10 i, 10 i + 10).
    Sanity check: the same procedure on the merged GEO 10 bp track is the identity; rep1 (the dominant replicate of the
    merged track) should correlate with it far better than rep2.
Outputs  output/tracks_10bp.npz (git-ignored), COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report  # noqa: E402

OUT = HERE / "output"


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("00c - per-sample 10 bp tracks",
                 "10 bp log2U per sample = mean of the GEO 1 bp track in each bin; merged = GEO merged 10 bp track as is.")
    tr = {}
    for s in config.SAMPLE_ORDER:
        tr[s] = io.read_binned_bigwig(io.geo_path(config.SAMPLES[s]["geo_1bp"]), 10, config.NUCLEAR_CHROMS)
    tr["3d7_merged"] = io.read_binned_bigwig(io.geo_path(config.GEO_MERGED_10BP), 10, config.NUCLEAR_CHROMS)
    flat = {f"{s}|{c}": tr[s][c] for s in tr for c in config.NUCLEAR_CHROMS}
    np.savez_compressed(OUT / "tracks_10bp.npz", **flat)
    cat = lambda s: np.concatenate([tr[s][c] for c in config.NUCLEAR_CHROMS])
    m = cat("3d7_merged")
    for s in config.SAMPLE_ORDER:
        v = cat(s)
        ok = ~np.isnan(v) & ~np.isnan(m)
        rep.check(f"{s} vs GEO merged (10 bp, nuclear)", "-", f"r = {np.corrcoef(v[ok], m[ok])[0, 1]:.3f}", "INFO",
                  f"Spearman {stats.spearmanr(v[ok], m[ok]).statistic:.3f}; n = {ok.sum():,} bins; "
                  f"bins > 0.5: {np.sum(v[ok] > .5):,} (merged: {np.sum(m[ok] > .5):,})")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
