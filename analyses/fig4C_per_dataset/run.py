#!/usr/bin/env python3
"""Fig. 4C per dataset - are uracil hotspots near predicted G4s in rep1, rep2, DMSO, RAPA (and merged)?

Purpose
    Fig. 4C used the merged 3D7 10 bp track (~89 % rep1). Check that the hotspot-G4 proximity is present in each dataset.
Method
    Per-dataset 10 bp tracks (analyses/00c_tracks_10bp: mean of GEO 1 bp tracks; merged = GEO merged track).
    Hotspot definitions per dataset (nuclear chromosomes, blacklisted bins removed):
      fixed    : 10 bp bins with log2U > 0.5 (the published cut-off). Counts differ 3-25-fold between datasets because the
                 1 bp-average tracks have different scale/noise (see 00b), so this is NOT like-for-like.
      matched  : the N top-ranked bins, N = number of merged-track hotspots (like-for-like across datasets).
    G4s: |G4Hunter| > 1.2 and size <= 200 bp filtered first (the 'prefiltered', recommended version of Fig. 4C).
    Statistics per dataset and definition: signed distance median / IQR, fraction overlapping and within 100 bp, and a
    random-position null (same number of 10 bp bins, uniformly on nuclear chromosomes, blacklist removed, 200
    permutations). Enrichment = observed / null mean; one-sided permutation p = (1 + #null >= obs) / (1 + 200).
    The null is not GC-matched (see analyses/fig4C_gc_matched_null for that).
Outputs  output/per_dataset_summary.csv, fig4C_per_dataset.png/svg, COMPARISON.md
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from udnaseq import config, genome  # noqa: E402
from udnaseq.report import SAMPLE_COLORS, Report, savefig  # noqa: E402

spec = importlib.util.spec_from_file_location("fig4c", HERE.parents[0] / "fig4C_g4_distance" / "run.py")
f4c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f4c)

spec2 = importlib.util.spec_from_file_location("fig4b", HERE.parents[0] / "fig4B_gc_correlation" / "run.py")
f4b = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(f4b)

OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
N_PERM, SEED = 200, 20240209


def dataset_bins(z, s):
    """All 10 bp bins of the nuclear genome with a value, blacklist removed -> DataFrame(chrom, start, end, value)."""
    bl = genome.blacklist()
    parts = []
    for c in config.NUCLEAR_CHROMS:
        v = z[f"{s}|{c}"]
        st = np.arange(len(v)) * 10
        ok = ~np.isnan(v) & ~genome.blacklist_mask(c, st, st + 10, bl)
        parts.append(pd.DataFrame({"chrom": c, "start": st[ok], "end": st[ok] + 10, "value": v[ok]}))
    return pd.concat(parts, ignore_index=True)


def exact_100bp_top(strong, n_top=1000, seed=SEED):
    """Robustness check at 100 bp with the EXACT RPGC-based log2U (no 1 bp averaging): what do the top-1000 bins look like?"""
    from udnaseq import io
    T = io.load_tracks_npz("tracks_100bp")
    seqs = f4b.read_fasta(config.FASTA)
    gc = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gc[m] = f4b.gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    genes = genome.protein_coding_genes()
    vsa = genes[genes.category != "Others"]
    rng = np.random.default_rng(seed)
    ok = np.isin(T["chrom"], config.NUCLEAR_CHROMS) & ~T["blacklisted"]
    out = []
    by_chrom = lambda ix: {c: T["start"][ix][T["chrom"][ix] == c] for c in np.unique(T["chrom"][ix])}
    for s in SETS:
        y = T[f"log2_{s}"]
        idx = np.where(ok & (y != 0))[0]
        top = idx[np.argsort(-y[idx], kind="stable")[:n_top]]
        rnd = idx[rng.integers(0, len(idx), 20000)]
        d = f4c.abs_distance_fast(by_chrom(top), strong, width=100)
        d0 = f4c.abs_distance_fast(by_chrom(rnd), strong, width=100)
        inv = 0
        for c, g in vsa.groupby("chrom"):
            st = T["start"][top][T["chrom"][top] == c]
            inv += sum(np.sum((st < r.end) & (st + 100 > r.start)) for r in g.itertuples())
        out.append(dict(dataset=s, n_top=n_top, min_log2U=y[top].min(), mean_GC=np.nanmean(gc[top]), genome_mean_GC=np.nanmean(gc[idx]),
                        frac_within_500bp_strong_G4=np.mean(d <= 500), random_frac_within_500bp=np.mean(d0 <= 500),
                        frac_in_VSA_genes=inv / n_top))
    return pd.DataFrame(out)


def null_stats(n, strong, rng, n_perm):
    old = f4c.N_PERM
    f4c.N_PERM = n_perm
    try:
        return f4c.permutation_null(n, strong, rng)
    finally:
        f4c.N_PERM = old


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 4C per dataset - hotspot to G4 distance",
                 "Signed distance from uracil hotspots (10 bp bins) to the closest strong predicted G4, per dataset; "
                 "fixed cut-off (log2U > 0.5) and count-matched (top-N) definitions; random-position null.")
    z = np.load(config.ANALYSES / "00c_tracks_10bp" / "output" / "tracks_10bp.npz")
    g4 = f4c.g4_sites()
    strong = g4[(g4.score.abs() > config.G4_MIN_SCORE) & (g4["size"] <= config.G4_MAX_SIZE)]
    bins = {s: dataset_bins(z, s) for s in SETS}
    n_match = int((bins["3d7_merged"].value > config.U_THRESHOLD).sum())
    rng = np.random.default_rng(SEED)
    rows, dists = [], {}
    for defn in ("fixed", "matched"):
        for s in SETS:
            b = bins[s]
            H = b[b.value > config.U_THRESHOLD] if defn == "fixed" else b.nlargest(n_match, "value", keep="all").iloc[:n_match]
            H = H.sort_values(["chrom", "start"]).reset_index(drop=True)
            cl = genome.closest_signed(H, strong)
            sm = f4c.summarise(cl.distance)
            obs = f4c.abs_distance_fast({c: g.start.values for c, g in H.groupby("chrom")}, strong)
            null = null_stats(len(H), strong, rng, N_PERM)
            row = dict(definition=defn, dataset=s, n_hotspots=len(H), min_log2U=H.value.min(), median_dist=sm["median"],
                       q25=sm["q25"], q75=sm["q75"])
            for k, o in (("frac_overlap", np.mean(obs == 0)), ("frac_within_100bp", np.mean(obs <= 100))):
                row[f"obs_{k}"] = o
                row[f"null_{k}"] = null[k].mean()
                row[f"fold_{k}"] = o / null[k].mean()
                row[f"p_{k}"] = (1 + (null[k] >= o).sum()) / (1 + N_PERM)
            rows.append(row)
            dists[(defn, s)] = cl.distance.values
            print(defn, s, len(H), flush=True)
    D = pd.DataFrame(rows)
    D.round(5).to_csv(OUT / "per_dataset_summary.csv", index=False)

    X = exact_100bp_top(strong)
    X.round(4).to_csv(OUT / "exact_100bp_top1000_profile.csv", index=False)
    for r in X.itertuples():
        rep.check(f"{r.dataset}: top {r.n_top} exact 100 bp bins", "", f"GC {r.mean_GC:.1%} (genome {r.genome_mean_GC:.1%}); "
                  f"{r.frac_within_500bp_strong_G4:.1%} within 500 bp of strong G4 (random {r.random_frac_within_500bp:.1%}); "
                  f"{r.frac_in_VSA_genes:.0%} in var/rif/stevor genes", "INFO", f"lowest log2U {r.min_log2U:.2f}")

    fig, axs = plt.subplots(1, 2, figsize=(9, 3))
    for ax, k in zip(axs, ("frac_within_100bp", "frac_overlap")):
        for j, defn in enumerate(("fixed", "matched")):
            d = D[D.definition == defn]
            x = np.arange(len(SETS)) + (j - .5) * .38
            ax.bar(x, d[f"obs_{k}"], .36, color=["#8B0000", "#00008B", "#008B8B", "#008B00", "k"], alpha=1 if j else .45,
                   label=defn)
            ax.scatter(x, d[f"null_{k}"], c="orange", s=14, zorder=3, label="random null" if (j == 0) else None)
        ax.set_xticks(range(len(SETS)), SETS, rotation=30, fontsize=7)
        ax.set_title(k.replace("frac_", "fraction ").replace("_", " "), fontsize=9)
    axs[0].legend(fontsize=7, frameon=False)
    savefig(fig, OUT, "fig4C_per_dataset")

    for r in D.itertuples():
        holds = r.fold_frac_within_100bp > 1 and r.p_frac_within_100bp < 0.05
        rep.check(f"{r.dataset} [{r.definition}, n = {r.n_hotspots:,}]: within 100 bp of strong G4",
                  "merged published: 25.5 % vs 1.1 % null", f"{r.obs_frac_within_100bp:.1%} vs {r.null_frac_within_100bp:.1%} null "
                  f"({r.fold_frac_within_100bp:.1f}x)", "HOLDS" if holds else "NOT SUPPORTED",
                  f"overlap {r.obs_frac_overlap:.1%} vs {r.null_frac_overlap:.1%} ({r.fold_frac_overlap:.1f}x); median distance "
                  f"{r.median_dist:g}, IQR ({r.q25:g}, {r.q75:g}); permutation p = {r.p_frac_within_100bp:.3g}; "
                  f"lowest hotspot log2U {r.min_log2U:.2f}")
    rep.note(f"Count-matched N = {n_match:,} (merged-track hotspots, nuclear, blacklist removed). The fixed log2U > 0.5 cut-off "
             "is shown because it is the published definition, but its hotspot count is dataset-dependent.")
    rep.note("Hotspot bins are not independent (adjacent 10 bp bins of one peak); permutations place bins independently, so "
             "p-values are optimistic. Compare effect sizes (fold enrichment), not p-values, across datasets.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
