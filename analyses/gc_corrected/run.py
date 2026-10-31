#!/usr/bin/env python3
"""GC-corrected U-enrichment: do the VSA and G4 findings survive removing GC-dependent coverage bias?

Purpose
    The FASTQ rebuild (alignment/README.md) shows 2-3x stronger GC bias in DMSO/RAPA IP and INPUT than in rep1/rep2, and
    the 100 bp log2U correlates with GC in every dataset (Fig. 4B). Remove the GC-dependent coverage expectation from IP and
    INPUT separately and repeat the two central tests at 100 bp (exact RPGC tracks).
Method
    For each sample and for IP and INPUT separately: coverage divided by the mean coverage of nuclear, non-blacklisted
    bins in the same GC class (40 quantile classes of 100 bp GC fraction), rescaled to the overall mean
    (explicit GC normalisation, the same idea as deepTools correctGCBias, but on binned tracks).
    log2U_gc = log2((IP_gc + 1) / (INPUT_gc + 1)).
    Test 1 (Fig 2C): mean log2U_gc of var / rif / stevor genes vs all other genes (Welch), per dataset.
        Gene value = mean over 100 bp bins whose centre lies in the gene (nuclear genes).
    Test 1b (gene-level GC control, no bin-level correction): raw gene mean log2U of "other" genes regressed on gene GC
        (quadratic); residual of each VSA family relative to that expectation (Welch vs the residuals of other genes).
    Test 2 (Fig 4C at 100 bp): the 1,000 highest log2U_gc bins (nuclear, non-blacklisted) per dataset; fraction within 500 bp
        of a strong G4 vs random bins and vs random bins matched to the same GC decile distribution; GC and VSA overlap.
Outputs  output/gene_enrichment_gc_corrected.csv, hotspot_g4_gc_corrected.csv, COMPARISON.md
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402
from udnaseq.report import Report  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


f4c = load("f4c", "fig4C_g4_distance/run.py")
f4b = load("f4b", "fig4B_gc_correlation/run.py")
OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
N_TOP, SEED, NCLASS = 1000, 20240209, 40


def gc_normalise(cov, gcf, ref):
    edges = np.quantile(gcf[ref], np.linspace(0, 1, NCLASS + 1)[1:-1])
    cls = np.digitize(gcf, edges)
    mean = np.array([cov[ref & (cls == k)].mean() for k in range(NCLASS)])
    return cov / mean[cls] * cov[ref].mean()


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("GC-corrected U-enrichment", "Fig 2C (VSA genes) and Fig 4C-like (G4 proximity) tests after explicit GC normalisation of IP and INPUT.")
    T = io.load_tracks_npz("tracks_100bp")
    seqs = f4b.read_fasta(config.FASTA)
    gcf = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gcf[m] = f4b.gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    nuc = np.isin(T["chrom"], config.NUCLEAR_CHROMS)
    ref = nuc & ~T["blacklisted"] & ~np.isnan(gcf)
    g4 = f4c.g4_sites()
    strong = g4[(g4.score.abs() > config.G4_MIN_SCORE) & (g4["size"] <= config.G4_MAX_SIZE)]
    genes = genome.protein_coding_genes()
    genes = genes[genes.chrom.isin(config.NUCLEAR_CHROMS)]
    ctr = T["start"] + config.BIN // 2
    rng = np.random.default_rng(SEED)
    grows, hrows, rrows = [], [], []
    ggc = np.array([seqs[r.chrom][r.start:r.end] for r in genes.itertuples()], dtype=object)
    ggc = np.array([np.mean((x == b"G") | (x == b"C")) for x in ggc])
    deciles = pd.qcut(gcf[ref], 10, labels=False, duplicates="drop")
    dec_full = np.full(len(gcf), -1)
    dec_full[np.where(ref)[0]] = deciles
    vsa = genes[genes.category != "Others"]
    for s in SETS:
        if s == "3d7_merged":
            ip, inp = T["ip_3d7_merged"], T["input_3d7_merged"]
        else:
            ip, inp = T[f"ip_{s}"], T[f"input_{s}"]
        ipc, inc = gc_normalise(ip, np.nan_to_num(gcf), ref), gc_normalise(inp, np.nan_to_num(gcf), ref)
        y = np.log2((ipc + 1) / (inc + 1))
        y0 = T[f"log2_{s}"]
        r_gc = stats.pearsonr(gcf[ref], y[ref]).statistic
        r_gc0 = stats.pearsonr(gcf[ref], y0[ref]).statistic
        # Test 1: gene means
        vals = {}
        for c, g in genes.groupby("chrom"):
            idx = np.where(T["chrom"] == c)[0]
            cc = ctr[idx]
            for r in g.itertuples():
                m = (cc >= r.start) & (cc < r.end)
                vals[r.id] = (y[idx][m].mean() if m.any() else np.nan, y0[idx][m].mean() if m.any() else np.nan)
        gv = genes.assign(u_gc=[vals[i][0] for i in genes.id], u_raw=[vals[i][1] for i in genes.id])
        others = gv[gv.category == "Others"]
        for fam in ("VAR", "RIF", "STEVOR"):
            f = gv[gv.category == fam]
            a, b = f.u_gc.dropna(), others.u_gc.dropna()
            w = stats.ttest_ind(a, b, equal_var=False)
            grows.append(dict(dataset=s, family=fam, n=len(a), mean_family_gc=a.mean(), mean_others_gc=b.mean(), diff_gc=a.mean() - b.mean(),
                              welch_p_gc=w.pvalue, diff_raw=f.u_raw.mean() - others.u_raw.mean(), pearson_gc_vs_log2U_raw=r_gc0,
                              pearson_gc_vs_log2U_gc=r_gc))
        # Test 1b: residual vs gene-GC expectation (raw tracks)
        oth = (gv.category == "Others").values & ~np.isnan(gv.u_raw.values)
        co = np.polyfit(ggc[oth], gv.u_raw.values[oth], 2)
        resid = gv.u_raw.values - np.polyval(co, ggc)
        for fam in ("VAR", "RIF", "STEVOR"):
            mf = (gv.category == fam).values & ~np.isnan(resid)
            w = stats.ttest_ind(resid[mf], resid[oth], equal_var=False)
            rrows.append(dict(dataset=s, family=fam, n=int(mf.sum()), mean_gene_GC=ggc[mf].mean(), others_mean_gene_GC=ggc[oth].mean(),
                              raw_diff=gv.u_raw.values[mf].mean() - gv.u_raw.values[oth].mean(), resid_diff=resid[mf].mean() - resid[oth].mean(),
                              welch_p=w.pvalue))
        # Test 2: top bins
        idx = np.where(ref & (y0 != 0))[0]
        top = idx[np.argsort(-y[idx], kind="stable")[:N_TOP]]
        by = lambda ix: {c: T["start"][ix][T["chrom"][ix] == c] for c in np.unique(T["chrom"][ix])}
        d = f4c.abs_distance_fast(by(top), strong, width=100)
        rnd = idx[rng.integers(0, len(idx), 20000)]
        d0 = f4c.abs_distance_fast(by(rnd), strong, width=100)
        # GC-decile matched random bins
        want = np.bincount(dec_full[top], minlength=10)
        pool = [idx[dec_full[idx] == k] for k in range(10)]
        mm = []
        for _ in range(200):
            ix = np.concatenate([rng.choice(pool[k], size=want[k]) for k in range(10) if want[k]])
            mm.append(np.mean(f4c.abs_distance_fast(by(ix), strong, width=100) <= 500))
        inv = 0
        for c, g in vsa.groupby("chrom"):
            st = T["start"][top][T["chrom"][top] == c]
            inv += sum(np.sum((st < r.end) & (st + 100 > r.start)) for r in g.itertuples())
        obs = np.mean(d <= 500)
        hrows.append(dict(dataset=s, n_top=N_TOP, min_log2U_gc=y[top].min(), mean_GC=np.nanmean(gcf[top]), frac_within_500bp_G4=obs,
                          random_frac=np.mean(d0 <= 500), gc_matched_random_frac=np.mean(mm), gc_matched_p=(1 + np.sum(np.array(mm) >= obs)) / 201,
                          frac_in_VSA_genes=inv / N_TOP))
        print(s, flush=True)
    G, H, R = pd.DataFrame(grows), pd.DataFrame(hrows), pd.DataFrame(rrows)
    R.round(5).to_csv(OUT / "gene_enrichment_vs_gene_gc_expectation.csv", index=False)
    for r in R.itertuples():
        rep.check(f"{r.dataset}: {r.family} vs gene-GC expectation", f"raw diff {r.raw_diff:+.3f}",
                  f"{r.resid_diff:+.3f} (Welch p = {r.welch_p:.1e})", "HOLDS" if r.resid_diff > 0 and r.welch_p < 1e-3 else "NOT SUPPORTED",
                  f"n = {r.n}; gene GC {r.mean_gene_GC:.1%} vs {r.others_mean_gene_GC:.1%} for other genes")
    G.round(5).to_csv(OUT / "gene_enrichment_gc_corrected.csv", index=False)
    H.round(5).to_csv(OUT / "hotspot_g4_gc_corrected.csv", index=False)
    for r in G.itertuples():
        rep.check(f"{r.dataset}: {r.family} vs other genes after GC correction", f"raw diff {r.diff_raw:+.3f}",
                  f"{r.diff_gc:+.3f} (Welch p = {r.welch_p_gc:.1e})", "HOLDS" if r.diff_gc > 0 and r.welch_p_gc < 1e-3 else "NOT SUPPORTED",
                  f"n = {r.n}; Pearson(GC, log2U) raw {r.pearson_gc_vs_log2U_raw:.2f} -> corrected {r.pearson_gc_vs_log2U_gc:.2f}")
    for r in H.itertuples():
        rep.check(f"{r.dataset}: top {r.n_top} GC-corrected 100 bp bins near strong G4", "-",
                  f"{r.frac_within_500bp_G4:.1%} within 500 bp vs {r.random_frac:.1%} random, {r.gc_matched_random_frac:.1%} GC-matched (p = {r.gc_matched_p:.3g})",
                  "INFO", f"GC {r.mean_GC:.1%}; {r.frac_in_VSA_genes:.0%} in var/rif/stevor genes; lowest log2U_gc {r.min_log2U_gc:.2f}")
    rep.note("GC correction is applied to binned 100 bp tracks (40 GC quantile classes), not to reads. Because uracil genuinely "
             "correlates with GC, correction removes part of the real signal as well as the technical bias; compare it with the raw result, "
             "do not replace it.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
