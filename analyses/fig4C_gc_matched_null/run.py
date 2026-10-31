#!/usr/bin/env python3
"""Fig. 4C GC-matched null and a GC + G4 model - does G4 proximity add anything beyond GC content?

Purpose
    Fig. 4C's random null (analyses/fig4C_g4_distance) is not GC-matched. G4 sites are GC-rich and uracil enrichment
    correlates with GC (Fig. 4B), so part of the hotspot-G4 excess could be a GC effect. Test that per dataset.
Part A - GC-matched null (10 bp hotspots, count-matched to the merged track: top N = 4,034 bins per dataset)
    Random 10 bp bins (nuclear, blacklist removed, bins with data) drawn so the hotspot GC distribution is reproduced exactly:
      match 'gc10'      : GC count of the 10 bp bin itself (0-10)
      match 'gc10+500'  : GC count of the bin x GC% of the 500 bp window centred on it (classes <15, 15-25, 25-35, 35-45, >=45 %)
    Same statistics and 200 permutations as fig4C_per_dataset; strata without candidates are dropped and reported.
Part B - linear model at 100 bp (exact RPGC log2U, per dataset, nuclear, blacklist removed, log2U != 0)
      M0: log2U ~ GC decile (100 bp) + GC decile (1 kb window)
      M1: M0 + G4  (bin within 100 bp of a strong G4)
      M2: M1 + var/rif/stevor gene (bin overlaps a VSA gene)
    G4 coefficient in log2U units with cluster-robust SE (clusters = 1 Mb blocks), and the R2 gain of G4 over GC alone.
Outputs  output/matched_null_summary.csv, model_coefficients.csv, COMPARISON.md
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402
from udnaseq.report import Report  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


f4c = load("f4c", "fig4C_g4_distance/run.py")
f4b = load("f4b", "fig4B_gc_correlation/run.py")
pd_ = load("f4cpd", "fig4C_per_dataset/run.py")

OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
N_PERM, SEED = 200, 20240209
W500 = 250


def gc_arrays(seqs):
    """Per chromosome: cumulative GC count (for window GC)."""
    return {c: np.concatenate([[0], np.cumsum((s == b"G") | (s == b"C"))]) for c, s in seqs.items()}


def window_gc(cum, starts, ends):
    L = len(cum) - 1
    a, b = np.clip(starts, 0, L), np.clip(ends, 0, L)
    return (cum[b] - cum[a]) / np.maximum(b - a, 1)


def add_gc(df, cum):
    gc10 = np.zeros(len(df), int)
    gc500 = np.zeros(len(df))
    for c, ix in df.groupby("chrom").indices.items():
        s, e = df.start.values[ix], df.end.values[ix]
        gc10[ix] = cum[c][np.minimum(e, len(cum[c]) - 1)] - cum[c][s]
        gc500[ix] = window_gc(cum[c], s + 5 - W500, s + 5 + W500)
    cls500 = np.digitize(gc500, [.15, .25, .35, .45])
    return df.assign(gc10=gc10, gc500=gc500, cls500=cls500)


def draw(pool, keys_obs, key_col, rng):
    """Draw one random sample with the same stratum counts as the observed hotspots."""
    idx = []
    for k, n in keys_obs.items():
        cand = pool.groups.get(k)
        if cand is None:
            continue
        idx.append(rng.choice(cand, size=n, replace=len(cand) < n))
    return np.concatenate(idx)


def matched_null(H, bins, strong, mode, rng):
    cols = ["gc10"] if mode == "gc10" else ["gc10", "cls500"]
    keyf = lambda d: list(zip(*[d[c].values for c in cols]))
    hk = pd.Series(keyf(H)).value_counts()
    bk = pd.Series(keyf(bins))
    groups = {k: np.asarray(v) for k, v in bk.groupby(bk).groups.items()}
    miss = {k: n for k, n in hk.items() if k not in groups}
    hk = hk[[k in groups for k in hk.index]]
    pool = type("P", (), {"groups": groups})
    res = []
    chrom, start = bins.chrom.values, bins.start.values
    for _ in range(N_PERM):
        ix = draw(pool, hk.to_dict(), cols, rng)
        d = f4c.abs_distance_fast({c: start[ix][chrom[ix] == c] for c in np.unique(chrom[ix])}, strong)
        res.append(dict(frac_overlap=np.mean(d == 0), frac_within_100bp=np.mean(d <= 100)))
    return pd.DataFrame(res), sum(miss.values())


def per_bin_abs_dist(T, strong, width=100):
    out = np.full(len(T["chrom"]), np.nan)
    for c in config.NUCLEAR_CHROMS:
        m = np.where(T["chrom"] == c)[0]
        out[m] = f4c.abs_distance_fast({c: T["start"][m]}, strong, width=width)
    return out


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 4C GC-matched null and GC + G4 model",
                 "Does hotspot-G4 proximity survive matching for GC content, per dataset?")
    seqs = f4b.read_fasta(config.FASTA)
    cum = gc_arrays(seqs)
    g4 = f4c.g4_sites()
    strong = g4[(g4.score.abs() > config.G4_MIN_SCORE) & (g4["size"] <= config.G4_MAX_SIZE)]
    z = np.load(config.ANALYSES / "00c_tracks_10bp" / "output" / "tracks_10bp.npz")
    rng = np.random.default_rng(SEED)

    # ------------------------------------------------------------ Part A
    rows = []
    allb = {s: pd_.dataset_bins(z, s) for s in SETS}
    n_match = int((allb["3d7_merged"].value > config.U_THRESHOLD).sum())
    for s in SETS:
        b = add_gc(allb[s], cum)
        H = b.nlargest(n_match, "value").sort_values(["chrom", "start"]).reset_index(drop=True)
        obs = f4c.abs_distance_fast({c: g.start.values for c, g in H.groupby("chrom")}, strong)
        row0 = dict(dataset=s, n_hotspots=len(H), hotspot_mean_gc10=H.gc10.mean() / 10, bins_mean_gc10=b.gc10.mean() / 10,
                    hotspot_mean_gc500=H.gc500.mean(), bins_mean_gc500=b.gc500.mean())
        for mode in ("gc10", "gc10+500"):
            null, miss = matched_null(H, b, strong, mode, rng)
            row = dict(row0, matching=mode, hotspots_without_stratum=miss)
            for k, o in (("frac_overlap", np.mean(obs == 0)), ("frac_within_100bp", np.mean(obs <= 100))):
                row[f"obs_{k}"], row[f"null_{k}"] = o, null[k].mean()
                row[f"null_{k}_q975"] = null[k].quantile(.975)
                row[f"fold_{k}"] = o / null[k].mean()
                row[f"p_{k}"] = (1 + (null[k] >= o).sum()) / (1 + N_PERM)
            rows.append(row)
            print(s, mode, flush=True)
    A = pd.DataFrame(rows)
    A.round(5).to_csv(OUT / "matched_null_summary.csv", index=False)
    for r in A.itertuples():
        keeps = r.fold_frac_within_100bp > 1 and r.p_frac_within_100bp < 0.05
        rep.check(f"{r.dataset} [{r.matching}]: within 100 bp of strong G4", "unmatched null: see fig4C_per_dataset",
                  f"{r.obs_frac_within_100bp:.1%} vs {r.null_frac_within_100bp:.1%} GC-matched null ({r.fold_frac_within_100bp:.1f}x)",
                  "SURVIVES" if keeps else "LOST",
                  f"overlap {r.obs_frac_overlap:.1%} vs {r.null_frac_overlap:.1%} ({r.fold_frac_overlap:.1f}x); p = {r.p_frac_within_100bp:.3g}; "
                  f"hotspot GC10 {r.hotspot_mean_gc10:.0%} vs all bins {r.bins_mean_gc10:.0%}; 500 bp-window GC {r.hotspot_mean_gc500:.0%} vs {r.bins_mean_gc500:.0%}; "
                  f"{r.hotspots_without_stratum} hotspots without candidates")

    # ------------------------------------------------------------ Part B
    T = io.load_tracks_npz("tracks_100bp")
    gc100 = np.full(len(T["chrom"]), np.nan)
    gc1k = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gc100[m] = window_gc(cum[str(c)], T["start"][m], T["end"][m])
        gc1k[m] = window_gc(cum[str(c)], T["start"][m] - 450, T["start"][m] + 550)
    near = per_bin_abs_dist(T, strong) <= 100  # NaN (non-nuclear contigs) -> False
    genes = genome.protein_coding_genes()
    vsa = genes[genes.category != "Others"]
    in_vsa = np.zeros(len(T["chrom"]), bool)
    for c, g in vsa.groupby("chrom"):
        m = np.where(T["chrom"] == c)[0]
        for r in g.itertuples():
            in_vsa[m] |= (T["start"][m] < r.end) & (T["end"][m] > r.start)
    ok0 = np.isin(T["chrom"], config.NUCLEAR_CHROMS) & ~T["blacklisted"]
    cluster = pd.factorize(pd.Series(T["chrom"]) + "_" + pd.Series(T["start"] // 1_000_000).astype(str))[0]
    mrows = []
    for s in SETS:
        y = T[f"log2_{s}"]
        ok = ok0 & (y != 0) & ~np.isnan(gc100)
        d100 = pd.qcut(gc100[ok], 10, labels=False, duplicates="drop")
        d1k = pd.qcut(gc1k[ok], 10, labels=False, duplicates="drop")
        gcm = pd.concat([pd.get_dummies(d100, prefix="g100", drop_first=True), pd.get_dummies(d1k, prefix="g1k", drop_first=True)],
                        axis=1).astype(float)
        base = sm.add_constant(gcm.reset_index(drop=True))
        g4v = pd.Series(near[ok].astype(float), name="G4_within_100bp")
        vv = pd.Series(in_vsa[ok].astype(float), name="in_VSA_gene")
        yy = y[ok]
        cl = cluster[ok]
        fits = {}
        for name, X in (("gc_only", base), ("gc+G4", pd.concat([base, g4v], axis=1)), ("gc+G4+VSA", pd.concat([base, g4v, vv], axis=1)),
                        ("G4_only", sm.add_constant(g4v))):
            fits[name] = sm.OLS(yy, X).fit(cov_type="cluster", cov_kwds={"groups": cl})
        for name, f in fits.items():
            r = dict(dataset=s, model=name, n_bins=int(ok.sum()), r2=f.rsquared, frac_G4_bins=g4v.mean())
            for term in ("G4_within_100bp", "in_VSA_gene"):
                if term in f.params:
                    r[f"{term}_coef"], r[f"{term}_se"] = f.params[term], f.bse[term]
            mrows.append(r)
    M = pd.DataFrame(mrows)
    M.round(5).to_csv(OUT / "model_coefficients.csv", index=False)
    for s in SETS:
        m = M[M.dataset == s].set_index("model")
        g, g1, g2, go = m.loc["gc_only"], m.loc["gc+G4"], m.loc["gc+G4+VSA"], m.loc["G4_only"]
        rep.check(f"{s}: G4 within 100 bp, log2U effect", "-",
                  f"alone {go.G4_within_100bp_coef:+.3f}; after GC {g1.G4_within_100bp_coef:+.3f} (SE {g1.G4_within_100bp_se:.3f}); "
                  f"after GC + VSA {g2.G4_within_100bp_coef:+.3f} (SE {g2.G4_within_100bp_se:.3f})", "INFO",
                  f"R2 GC only {g.r2:.3f} -> +G4 {g1.r2:.3f} -> +VSA {g2.r2:.3f}; {g.frac_G4_bins:.1%} of bins near a strong G4")
    rep.note("Part B uses 100 bp bins with the exact RPGC log2U. Bins are spatially autocorrelated; SEs are clustered by 1 Mb block "
             "but remain approximate.")
    rep.note("'G4_only' is the unadjusted difference in mean log2U between bins near/not near a strong G4; 'after GC' is the "
             "same contrast within GC deciles of the 100 bp bin and its 1 kb window.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
