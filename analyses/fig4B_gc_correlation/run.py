#!/usr/bin/env python3
"""Fig. 4B - Correlation between GC content and uracil enrichment.

Purpose
    Test whether more uracilated 100 bp bins are more GC-rich (expected if uracil
    arises from cytosine deamination rather than dUTP misincorporation).

Method (replaces bedtools subtract + bedtools nuc + R cor())
    Bins: merged 3D7 track, 100 bp, blacklisted bins removed (analyses/00_tracks).
          Bins with log2 == 0 exactly (empty in the deepTools output) are dropped, as in legacy.
    GC%:  fraction of G/C among the bin's bases (bedtools nuc 'pct_gc' definition: G+C / bin length).
    Statistics: Pearson and Spearman correlation, least-squares line.
Outputs  output/gc_vs_log2U_100bp.csv.gz, fig4B_gc_correlation.png/svg, stats.csv, COMPARISON.md
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
from udnaseq.report import Report, savefig, status  # noqa: E402

OUT = HERE / "output"
LEGACY_GC = config.LEGACY / "GCcontent" / "GCcorrelation" / "averageGCcontent_100bin_merged_blacklisted.bed"


def read_fasta(path) -> dict[str, np.ndarray]:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = np.frombuffer("".join(buf).upper().encode(), dtype="S1")
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    seqs[name] = np.frombuffer("".join(buf).upper().encode(), dtype="S1")
    return seqs


def gc_fraction(seq: np.ndarray, starts, ends) -> np.ndarray:
    gc = np.concatenate([[0], np.cumsum((seq == b"G") | (seq == b"C"))])
    return (gc[ends] - gc[starts]) / (ends - starts)


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 4B - GC content vs log2(U-enrichment)",
                 "Correlation of GC% and merged 3D7 log2(U-enrichment) across 100 bp bins (blacklist removed).")
    T = io.load_tracks_npz("tracks_100bp")
    keep = (~T["blacklisted"]) & (T["log2_3d7_merged"] != 0)
    df = pd.DataFrame({"chrom": T["chrom"][keep], "start": T["start"][keep], "end": T["end"][keep],
                       "log2U": T["log2_3d7_merged"][keep]})
    seqs = read_fasta(config.FASTA)
    df["gc"] = np.nan
    for c, idx in df.groupby("chrom").groups.items():
        d = df.loc[idx]
        df.loc[idx, "gc"] = gc_fraction(seqs[c], d.start.values, d.end.values)
    df.round(5).to_csv(OUT / "gc_vs_log2U_100bp.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})

    r_p = stats.pearsonr(df.gc, df.log2U)
    r_s = stats.spearmanr(df.gc, df.log2U)
    fit = stats.linregress(df.gc * 100, df.log2U)
    nuc = df[df.chrom.isin(config.NUCLEAR_CHROMS)]
    st = pd.DataFrame([{"set": "all bins", "n": len(df), "pearson_r": r_p.statistic, "pearson_p": r_p.pvalue,
                        "spearman_rho": r_s.statistic, "slope_per_GCpct": fit.slope, "intercept": fit.intercept},
                       {"set": "nuclear chromosomes", "n": len(nuc), "pearson_r": stats.pearsonr(nuc.gc, nuc.log2U).statistic,
                        "pearson_p": np.nan, "spearman_rho": stats.spearmanr(nuc.gc, nuc.log2U).statistic,
                        "slope_per_GCpct": np.nan, "intercept": np.nan}])
    st.to_csv(OUT / "stats.csv", index=False)

    # ---- figure (legacy style: jet colour by log2U, y from -0.1) --------------------------
    fig, ax = plt.subplots(figsize=(5.5, 4.4))
    o = np.argsort(df.log2U.values)
    sc = ax.scatter(df.gc.values[o] * 100, df.log2U.values[o], c=df.log2U.values[o], cmap="jet", vmin=-0.8, vmax=0.5,
                    s=2, lw=0, rasterized=True)
    xx = np.array([0, 80])
    ax.plot(xx, fit.intercept + fit.slope * xx, c="k", lw=1)
    ax.set_xlim(0, 80)
    ax.set_ylim(-1, 1.6)
    ax.axhline(0, c="grey", lw=0.5)
    ax.set_xlabel("GC%")
    ax.set_ylabel("log2(U-enrichment)")
    ax.text(0.03, 0.95, f"Pearson r = {r_p.statistic:.2f}\nSpearman ρ = {r_s.statistic:.2f}\nn = {len(df):,} bins",
            transform=ax.transAxes, va="top", fontsize=8)
    fig.colorbar(sc, ax=ax, fraction=0.04, label="log2(U-enrichment)")
    savefig(fig, OUT, "fig4B_gc_correlation")

    # ---- compare with legacy -------------------------------------------------------------
    leg = pd.read_csv(LEGACY_GC, sep="\t")
    leg.columns = ["chrom", "start", "end", "log2U", "gc"]
    lr_p, lr_s = leg.gc.corr(leg.log2U), leg.gc.corr(leg.log2U, method="spearman")
    rep.check("Pearson r", f"{lr_p:.3f} (legacy table); manuscript 'r = 0.42'", f"{r_p.statistic:.3f}",
              status(abs(lr_p - r_p.statistic), 1e-3, 0.02))
    rep.check("Spearman rho", f"{lr_s:.3f} (legacy table)", f"{r_s.statistic:.3f}", status(abs(lr_s - r_s.statistic), 1e-3, 0.02))
    rep.check("bins", f"{len(leg):,}", f"{len(df):,}", "DIFFERS" if abs(len(leg) - len(df)) > 100 else "MATCH",
              "legacy table stops at Pf3D7_13_v3:2,249,900 (truncated output; 2.25-2.89 Mb of chr13 missing). Correlations are unaffected")
    l100 = leg[(leg.end - leg.start) == 100]
    m = l100.merge(df, on=["chrom", "start", "end"], suffixes=("_leg", ""))
    rep.check("per-bin GC% (shared 100 bp bins)", f"{len(l100):,} bins", f"{len(m):,} matched",
              status(float(np.abs(m.gc_leg - m.gc).max())), f"max |diff| = {np.abs(m.gc_leg - m.gc).max():.1e}")
    rep.check("per-bin log2U (shared 100 bp bins)", "", "", status(float(np.abs(m.log2U_leg - m.log2U).max())),
              f"max |diff| = {np.abs(m.log2U_leg - m.log2U).max():.1e}")
    rep.note("The manuscript's r = 0.42 is closest to the Spearman coefficient; the Pearson coefficient is "
             f"{r_p.statistic:.2f}. Suggest reporting both, or stating which one is used.")
    rep.note(f"Nuclear chromosomes only (excluding apicoplast/mitochondrion): Pearson {st.pearson_r[1]:.3f}, "
             f"Spearman {st.spearman_rho[1]:.3f}.")
    rep.note("With >200,000 bins any r is 'significant' (p ~ 0); the effect size is moderate: "
             f"+{fit.slope * 10:.3f} log2U per 10 GC% points.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
