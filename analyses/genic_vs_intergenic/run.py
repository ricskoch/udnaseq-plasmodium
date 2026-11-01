#!/usr/bin/env python3
"""Genic vs intergenic log2U: how much is GC, how much is left at equal GC?

Purpose
    Genes are GC-richer than intergenic DNA, so a GC correction can remove a genuine genic/intergenic contrast. Quantify it.
Method (100 bp bins, nuclear, non-blacklisted; genic = bin centre inside a protein-coding gene span, else intergenic)
    A. Genic minus intergenic mean log2U: raw, and after per-track GC normalisation (gc_bin).
    B. Equal-GC contrast: 20 quantile classes of bin GC; genic-minus-intergenic mean log2U within each class, averaged with weights
       min(n_genic, n_intergenic). Uses raw tracks, so the GC confound is removed without touching real differences at equal GC.
    C. Same equal-GC contrast for VSA-gene bins vs other-gene bins (both genic).
    Bootstrap by 100 kb blocks (200 resamples) for the 95 % CI of B and C.
Outputs  output/genic_contrast.csv, COMPARISON.md
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


gcc = load("gcc", "gc_corrected/run.py")
f4b = load("f4b", "fig4B_gc_correlation/run.py")
OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
NC, NB, SEED = 20, 200, 20240209


def contrast(y, a, b, cls, w=None):
    """weighted mean over GC classes of mean(y[a]) - mean(y[b])."""
    d, ww = [], []
    for k in range(NC):
        ia, ib = a & (cls == k), b & (cls == k)
        na, nb = ia.sum(), ib.sum()
        if na >= 5 and nb >= 5:
            d.append(y[ia].mean() - y[ib].mean())
            ww.append(min(na, nb))
    return np.average(d, weights=ww) if d else np.nan


def main():
    OUT.mkdir(exist_ok=True)
    T = io.load_tracks_npz("tracks_100bp")
    seqs = f4b.read_fasta(config.FASTA)
    gcf = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gcf[m] = f4b.gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    nuc = np.isin(T["chrom"], config.NUCLEAR_CHROMS)
    ref = nuc & ~T["blacklisted"] & ~np.isnan(gcf)
    genes = genome.protein_coding_genes()
    genes = genes[genes.chrom.isin(config.NUCLEAR_CHROMS)]
    ctr = T["start"] + config.BIN // 2
    genic = np.zeros(len(ctr), bool)
    vsa = np.zeros(len(ctr), bool)
    for c, g in genes.groupby("chrom"):
        idx = np.where(T["chrom"] == c)[0]
        for r in g.itertuples():
            m = (ctr[idx] >= r.start) & (ctr[idx] < r.end)
            genic[idx[m]] = True
            if r.category != "Others":
                vsa[idx[m]] = True
    inter = ref & ~genic
    gen = ref & genic
    oth = gen & ~vsa
    edges = np.quantile(gcf[ref], np.linspace(0, 1, NC + 1)[1:-1])
    cls = np.digitize(gcf, edges)
    block = (T["start"] // 100_000).astype(int) + 1000 * pd.factorize(T["chrom"])[0]
    ub = np.unique(block[ref])
    rng = np.random.default_rng(SEED)
    rows = []
    print(f"bins: genic {gen.sum()}, intergenic {inter.sum()}, VSA-gene {(gen & vsa).sum()}; GC genic {gcf[gen].mean():.3f}, intergenic {gcf[inter].mean():.3f}")
    for s in SETS:
        ip, inp = (T["ip_3d7_merged"], T["input_3d7_merged"]) if s == "3d7_merged" else (T[f"ip_{s}"], T[f"input_{s}"])
        ipc, inc = gcc.gc_normalise(ip, np.nan_to_num(gcf), ref), gcc.gc_normalise(inp, np.nan_to_num(gcf), ref)
        yr, yc = T[f"log2_{s}"], np.log2((ipc + 1) / (inc + 1))
        A_raw = yr[gen].mean() - yr[inter].mean()
        A_gc = yc[gen].mean() - yc[inter].mean()
        B = contrast(yr, gen, inter, cls)
        C = contrast(yr, gen & vsa, oth, cls)
        bs_B, bs_C = [], []
        for _ in range(NB):
            pick = rng.choice(ub, len(ub))
            cnt = pd.Series(pick).value_counts()
            wgt = pd.Series(block).map(cnt).fillna(0).values  # resample via weights
            wv = np.repeat(np.arange(len(block)), 0)  # placeholder to keep memory low
            sel = np.isin(block, cnt.index)
            # weighted class means
            def wc(a, b):
                d, ww = [], []
                for k in range(NC):
                    ia, ib = a & (cls == k) & sel, b & (cls == k) & sel
                    sa, sb = wgt[ia].sum(), wgt[ib].sum()
                    if sa >= 5 and sb >= 5:
                        d.append(np.average(yr[ia], weights=wgt[ia]) - np.average(yr[ib], weights=wgt[ib]))
                        ww.append(min(sa, sb))
                return np.average(d, weights=ww) if d else np.nan
            bs_B.append(wc(gen, inter))
            bs_C.append(wc(gen & vsa, oth))
        rows.append(dict(dataset=s, genic_minus_inter_raw=A_raw, genic_minus_inter_gcnorm=A_gc, retained_pct=100 * A_gc / A_raw if A_raw else np.nan,
                         equalGC_genic_minus_inter=B, B_lo=np.nanpercentile(bs_B, 2.5), B_hi=np.nanpercentile(bs_B, 97.5),
                         equalGC_VSA_minus_otherGenic=C, C_lo=np.nanpercentile(bs_C, 2.5), C_hi=np.nanpercentile(bs_C, 97.5)))
        print(rows[-1], flush=True)
    R = pd.DataFrame(rows).round(4)
    R.to_csv(OUT / "genic_contrast.csv", index=False)
    md = ["# Genic vs intergenic log2U: raw, GC-normalised and at equal GC", "",
          __doc__.split("Method")[1].split("Outputs")[0].strip(), "",
          f"Bins: genic {gen.sum()}, intergenic {inter.sum()}; mean GC genic {100*gcf[gen].mean():.1f} %, intergenic {100*gcf[inter].mean():.1f} %.", "",
          R.to_markdown(index=False), ""]
    (HERE / "COMPARISON.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
