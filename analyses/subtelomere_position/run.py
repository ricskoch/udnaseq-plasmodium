#!/usr/bin/env python3
"""Does chromosome position (subtelomeric vs internal) explain the VSA signal at equal GC?

Purpose  VSA genes cluster in subtelomeres. Separate 'VSA identity' from 'subtelomeric location' without extra data.
Method   100 bp genic bins (nuclear, non-blacklisted). Subtelomeric = within 100 kb of either chromosome end (config below).
         Four groups: VSA/subtelomeric, VSA/internal, other/subtelomeric, other/internal. Equal-GC contrast = mean over 20 GC
         classes of group-mean log2U minus the reference group (other/internal), weights min(n_a, n_b) (classes with >=5 bins each).
         95 % CI by 100 kb block bootstrap (200 resamples). Raw tracks.
Outputs  output/subtelomere_contrast.csv, COMPARISON.md
"""
import importlib.util, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from udnaseq import config, genome, io  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


f4b = load("f4b", "fig4B_gc_correlation/run.py")
OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
NC, NB, SEED, SUBTEL = 20, 200, 20240209, 100_000


def main():
    OUT.mkdir(exist_ok=True)
    T = io.load_tracks_npz("tracks_100bp")
    seqs = f4b.read_fasta(config.FASTA)
    clen = {c: len(s) for c, s in seqs.items()}
    gcf = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gcf[m] = f4b.gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    ref = np.isin(T["chrom"], config.NUCLEAR_CHROMS) & ~T["blacklisted"] & ~np.isnan(gcf)
    genes = genome.protein_coding_genes()
    genes = genes[genes.chrom.isin(config.NUCLEAR_CHROMS)]
    ctr = T["start"] + config.BIN // 2
    genic = np.zeros(len(ctr), bool); vsa = np.zeros(len(ctr), bool)
    for c, g in genes.groupby("chrom"):
        idx = np.where(T["chrom"] == c)[0]
        for r in g.itertuples():
            m = (ctr[idx] >= r.start) & (ctr[idx] < r.end)
            genic[idx[m]] = True
            if r.category != "Others":
                vsa[idx[m]] = True
    L = np.array([clen.get(str(c), 0) for c in T["chrom"]])
    sub = (ctr < SUBTEL) | (ctr > L - SUBTEL)
    G = {"VSA/subtelomeric": genic & vsa & sub, "VSA/internal": genic & vsa & ~sub,
         "other/subtelomeric": genic & ~vsa & sub, "other/internal": genic & ~vsa & ~sub}
    base = G["other/internal"]
    edges = np.quantile(gcf[ref], np.linspace(0, 1, NC + 1)[1:-1]); cls = np.digitize(gcf, edges)
    block = (ctr // 100_000).astype(int) + 1000 * pd.factorize(T["chrom"])[0]
    ub = np.unique(block[ref]); rng = np.random.default_rng(SEED)

    def contrast(y, a, b, wgt):
        d, ww = [], []
        for k in range(NC):
            ia, ib = ref & a & (cls == k) & (wgt > 0), ref & b & (cls == k) & (wgt > 0)
            sa, sb = wgt[ia].sum(), wgt[ib].sum()
            if ia.sum() >= 5 and ib.sum() >= 5:
                d.append(np.average(y[ia], weights=wgt[ia]) - np.average(y[ib], weights=wgt[ib])); ww.append(min(sa, sb))
        return np.average(d, weights=ww) if d else np.nan

    rows = []
    ones = np.ones(len(ctr))
    boots = []
    for _ in range(NB):
        cnt = pd.Series(rng.choice(ub, len(ub))).value_counts()
        boots.append(pd.Series(block).map(cnt).fillna(0).values.astype(float))
    for s in SETS:
        y = T[f"log2_{s}"]
        for name, a in G.items():
            if name == "other/internal":
                continue
            c0 = contrast(y, a, base, ones)
            bs = [contrast(y, a, base, w) for w in boots]
            rows.append(dict(dataset=s, group=name, n_bins=int((ref & a).sum()), mean_gc=round(float(gcf[ref & a].mean()), 3),
                             equalGC_vs_other_internal=c0, lo=np.nanpercentile(bs, 2.5), hi=np.nanpercentile(bs, 97.5)))
        print(s, flush=True)
    R = pd.DataFrame(rows).round(4); R.to_csv(OUT / "subtelomere_contrast.csv", index=False)
    piv = R.pivot(index="dataset", columns="group", values="equalGC_vs_other_internal")
    (HERE / "COMPARISON.md").write_text("\n".join(["# Subtelomeric position vs VSA identity at equal GC", "",
        "Contrast = log2U difference to non-VSA genic bins in internal regions, equal-GC, raw tracks. Subtelomeric = within 100 kb of a chromosome end.", "",
        R.to_markdown(index=False), "", "Point estimates:", "", piv.round(3).to_markdown(), ""]))


if __name__ == "__main__":
    main()
