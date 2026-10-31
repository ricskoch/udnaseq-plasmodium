#!/usr/bin/env python3
"""Fig. 3B - Uracil enrichment profiles along var and rif/stevor genes.

Purpose
    Show where within the VSA genes uracil accumulates (var: first exon,
    ~first 5-6 kb; rif/stevor: second exon).

Method (replaces deepTools computeMatrix scale-regions + plotHeatmap)
    Input: GEO merged 3D7 log2 track at 10 bp (the legacy matrices were built from it).
    Every gene is scaled to a fixed body length and split into 10 bp bins; flanks are
    taken unscaled; minus-strand genes are reversed so that all run TSS -> TES.
      var        : 2 kb upstream | body scaled to 10 kb | 2 kb downstream  (62 genes)
      rif/stevor : 1 kb upstream | body scaled to 1.5 kb | 1 kb downstream (212 genes)
                   + a 2 kb-body version, because the manuscript legend states 2,000 bp
    Bin value = mean of the bigWig over the bin (pyBigWig stats, as deepTools).
    Heatmap rows are sorted by the row sum (descending); NaN (no mapped reads) is white.
Outputs  output/<set>_matrix.csv.gz, <set>_profile.csv, fig3B_<set>.png/svg, COMPARISON.md
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
from udnaseq.report import Report, savefig, status  # noqa: E402

OUT = HERE / "output"
LEG = config.LEGACY / "var_heatmaps"
SETS = {
    "var": dict(bed=config.VAR_BED, up=2000, body=10000, down=2000, legacy="var_scaled_matrix.tab"),
    "rifstevor": dict(bed=config.RIFSTEVOR_BED, up=1000, body=1500, down=1000, legacy="rifinstevor_scaled_matrix.tab"),
    "rifstevor_body2kb": dict(bed=config.RIFSTEVOR_BED, up=1000, body=2000, down=1000, legacy=None),
}
BS = 10


def scale_regions(bw, regions: pd.DataFrame, up, body, down, bs=BS) -> np.ndarray:
    f = lambda x: np.array([np.nan if y is None else y for y in x], float)
    rows = []
    for r in regions.itertuples():
        s, e, L = int(r.start), int(r.end), bw.chroms()[r.chrom]
        if r.strand == "-":
            left, right = down, up
        else:
            left, right = up, down
        lv = f(bw.stats(r.chrom, max(s - left, 0), s, nBins=left // bs)) if s - left >= 0 else np.full(left // bs, np.nan)
        bv = f(bw.stats(r.chrom, s, e, nBins=body // bs))
        rv = f(bw.stats(r.chrom, e, min(e + right, L), nBins=right // bs)) if e + right <= L else np.full(right // bs, np.nan)
        v = np.concatenate([lv, bv, rv])
        rows.append(v[::-1] if r.strand == "-" else v)
    return np.vstack(rows)


def heatmap(M, names, up, body, down, title, highlight=None):
    order = np.argsort(-np.nansum(M, axis=1), kind="stable")
    M, names = M[order], np.asarray(names)[order]
    x = np.arange(M.shape[1])
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(4.2, 2.2 + 0.075 * len(M) if len(M) < 80 else 9),
                                 gridspec_kw={"height_ratios": [1, 4]}, sharex=True)
    mu, sd = np.nanmean(M, 0), np.nanstd(M, 0)
    a1.plot(x, mu, c="#08306b", lw=1.2)
    a1.fill_between(x, mu - sd, mu + sd, color="#9ecae1", alpha=0.6, lw=0)
    a1.set_ylabel("log2 U")
    a1.set_title(title, fontsize=9)
    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")
    im = a2.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=-0.5, vmax=0.8, interpolation="nearest")
    ticks = [0, up // BS, (up + body) // BS, M.shape[1] - 1]
    a2.set_xticks(ticks, [f"-{up / 1000:g} kb", "TSS", "TES", f"+{down / 1000:g} kb"])
    for t in ticks[1:3]:
        a2.axvline(t, c="k", lw=0.5, ls="--")
        a1.axvline(t, c="k", lw=0.5, ls="--")
    if len(M) <= 80:
        a2.set_yticks(range(len(M)), names, fontsize=5)
        for lab, n in zip(a2.get_yticklabels(), names):
            if highlight and n in highlight:
                lab.set_color(highlight[n]); lab.set_fontweight("bold")
    else:
        a2.set_yticks([])
        a2.set_ylabel(f"{len(M)} genes")
    fig.colorbar(im, ax=a2, fraction=0.04, pad=0.02, label="log2(U-enrichment)")
    fig.tight_layout()
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 3B - VSA metaprofiles",
                 "Scaled-region matrices (computeMatrix scale-regions equivalent) of the merged 3D7 10 bp track "
                 "over var and rif/stevor genes.")
    bw = pyBigWig.open(str(io.geo_path(config.GEO_MERGED_10BP)))
    hl = {g: "#1f3fbf" for g in config.VAR3} | {g: "#c00000" for g in config.VAR2CSA}
    for name, p in SETS.items():
        reg = io.read_bed(p["bed"])
        M = scale_regions(bw, reg, p["up"], p["body"], p["down"])
        pd.DataFrame(M.round(6), index=reg.name).to_csv(OUT / f"{name}_matrix.csv.gz", compression={"method": "gzip", "mtime": 0})
        x = (np.arange(M.shape[1]) * BS) - p["up"]
        pd.DataFrame({"position_bp_scaled": x, "mean": np.nanmean(M, 0), "sd": np.nanstd(M, 0),
                      "n": np.sum(~np.isnan(M), 0)}).round(6).to_csv(OUT / f"{name}_profile.csv", index=False)
        savefig(heatmap(M, reg.name.values, p["up"], p["body"], p["down"],
                        f"{name}: body scaled to {p['body'] / 1000:g} kb, {BS} bp bins", hl), OUT, f"fig3B_{name}")
        if p["legacy"]:
            params, lreg, L = io.read_deeptools_matrix(LEG / p["legacy"])
            assert list(lreg.name) == list(reg.name)
            same_nan = np.array_equal(np.isnan(L), np.isnan(M))
            err = np.nanmax(np.abs(L - M), axis=1)
            rep.check(f"{name}: matrix {L.shape}", f"body {params['body'][0]} bp, bin {params['bin size'][0]}",
                      f"max |diff| = {np.nanmax(err):.1e}", status(np.nanmax(err), 1e-4, 0.01) if same_nan else "DIFFERS",
                      f"{(err > 1e-4).sum()} of {len(err)} rows differ by > 1e-4 (bin-boundary rounding); "
                      f"NaN pattern {'identical' if same_nan else 'differs'}")
            dprof = float(np.nanmax(np.abs(np.nanmean(L, 0) - np.nanmean(M, 0))))
            rep.check(f"{name}: mean profile", "", f"max |diff| = {dprof:.1e}", status(dprof, 1e-4, 0.01))
        mu = np.nanmean(M, 0)
        b0, b1 = p["up"] // BS, (p["up"] + p["body"]) // BS
        peak = (np.nanargmax(mu[b0:b1]) * BS) * 1.0
        rep.note(f"{name}: mean log2U upstream {np.nanmean(mu[:b0]):.3f}, body {np.nanmean(mu[b0:b1]):.3f}, "
                 f"downstream {np.nanmean(mu[b1:]):.3f}; body maximum at {peak:.0f} of {p['body']} scaled bp.")
    v = pd.read_csv(OUT / "var_profile.csv")
    body = v[(v.position_bp_scaled >= 0) & (v.position_bp_scaled < 10000)]
    rep.check("var: enrichment concentrated in first ~5-6 kb (exon 1)", "yes",
              f"mean 0-6 kb {body[body.position_bp_scaled < 6000]['mean'].mean():.3f} vs 6-10 kb "
              f"{body[body.position_bp_scaled >= 6000]['mean'].mean():.3f}", "INFO")
    rep.note("The legacy rif/stevor matrix used a 1,500 bp scaled body (and var_heatmaps/deeptools.sh says 2,000). "
             "The manuscript legend says 2,000 bp. Both versions are produced; please align the legend with the "
             "version used in the figure.")
    rep.note("Both legacy heatmap SVGs carry the title 'Rif/stevor gene U-enrichment with 10bp binsize', i.e. the "
             "var heatmap (UDNA_heatmap_over_var_scaled.svg) is mislabelled - check the panel title in the figure.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
