#!/usr/bin/env python3
"""Is the VSA excess in U-DNA-Seq 'heterochromatin' (H3K9me3) at equal GC?

Purpose  The IFA shows U-DNA / PfUNG / H3K9me3 colocalisation. Genomic counterpart: uracil vs H3K9me3 state, and VSA identity vs state.
Data     Public 3D7 wild-type H3K9me3 ChIP-seq and INPUT bigWigs, trophozoite 30-40 hpi (GEO GSE157709, WT rep2: GSM5353855 / GSM5353856).
         Rep1 (GSM4773592) was not downloadable from GEO (HTTP 500 at time of run). Stage is earlier than the Percoll-synchronised
         late schizonts of the U-DNA-Seq samples (last 4 h of the cycle); var/rif/stevor heterochromatin is stable over the cycle,
         but this is an assumption, not a measurement.
Method   H3K9me3 score per 100 bp bin = log2((IP/mean IP + 0.1)/(INPUT/mean INPUT + 0.1)). Heterochromatin = score >= HET (1.0;
         sensitivity 0.5 and 1.5). Equal-GC contrasts as in genic_vs_intergenic (20 GC classes, weights min(n_a, n_b), 100 kb block bootstrap, raw tracks):
         A. het vs non-het, all nuclear bins; B. same, genic only; C. VSA vs other genic bins within het bins and within non-het bins.
Outputs  output/h3k9me3_contrast.csv, COMPARISON.md
"""
import importlib.util, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyBigWig  # noqa: E402
from udnaseq import config, genome, io  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


f4b = load("f4b", "fig4B_gc_correlation/run.py")
OUT = HERE / "output"
D = config.DATA / "h3k9me3"
IP, INP = D / "GSM5353855_WT_H3K9me3_T_rep2_HM_ChIP-seq.bw", D / "GSM5353856_WT_input_T_rep2_HM_ChIP-seq.bw"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
NC, NB, SEED = 20, 200, 20240209
HETS = [0.5, 1.0, 1.5]


def bin_values(path, T):
    b = pyBigWig.open(str(path)); chroms = b.chroms()
    out = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        c = str(c)
        if c not in chroms:
            continue
        idx = np.where(T["chrom"] == c)[0]
        for i0 in range(0, len(idx), 2000):
            ii = idx[i0:i0 + 2000]
            s, e = T["start"][ii], np.minimum(T["start"][ii] + config.BIN, chroms[c])
            out[ii] = [b.stats(c, int(x), int(y), exact=False)[0] if y > x else np.nan for x, y in zip(s, e)]
    return np.array([np.nan if v is None else v for v in out], float)


def main():
    OUT.mkdir(exist_ok=True)
    T = io.load_tracks_npz("tracks_100bp")
    seqs = f4b.read_fasta(config.FASTA)
    gcf = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gcf[m] = f4b.gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    ref = np.isin(T["chrom"], config.NUCLEAR_CHROMS) & ~T["blacklisted"] & ~np.isnan(gcf)
    ip, inp = bin_values(IP, T), bin_values(INP, T)
    ok = ref & ~np.isnan(ip) & ~np.isnan(inp)
    h3 = np.log2((ip / np.nanmean(ip[ok]) + 0.1) / (inp / np.nanmean(inp[ok]) + 0.1))
    ref = ok
    genes = genome.protein_coding_genes(); genes = genes[genes.chrom.isin(config.NUCLEAR_CHROMS)]
    ctr = T["start"] + config.BIN // 2
    genic = np.zeros(len(ctr), bool); vsa = np.zeros(len(ctr), bool)
    for c, g in genes.groupby("chrom"):
        idx = np.where(T["chrom"] == c)[0]
        for r in g.itertuples():
            m = (ctr[idx] >= r.start) & (ctr[idx] < r.end)
            genic[idx[m]] = True
            if r.category != "Others":
                vsa[idx[m]] = True
    edges = np.quantile(gcf[ref], np.linspace(0, 1, NC + 1)[1:-1]); cls = np.digitize(gcf, edges)
    block = (ctr // 100_000).astype(int) + 1000 * pd.factorize(T["chrom"])[0]
    ub = np.unique(block[ref]); rng = np.random.default_rng(SEED)
    boots = []
    for _ in range(NB):
        cnt = pd.Series(rng.choice(ub, len(ub))).value_counts()
        boots.append(pd.Series(block).map(cnt).fillna(0).values.astype(float))
    ones = np.ones(len(ctr))

    def contrast(y, a, b, w):
        d, ww = [], []
        for k in range(NC):
            ia, ib = ref & a & (cls == k) & (w > 0), ref & b & (cls == k) & (w > 0)
            if ia.sum() >= 5 and ib.sum() >= 5:
                d.append(np.average(y[ia], weights=w[ia]) - np.average(y[ib], weights=w[ib])); ww.append(min(w[ia].sum(), w[ib].sum()))
        return np.average(d, weights=ww) if d else np.nan

    print(f"bins with H3K9me3 data: {ref.sum()}; fraction het(>=1): {np.mean(h3[ref] >= 1):.3f}; VSA-bin het fraction {np.mean(h3[ref & vsa] >= 1):.2f}; other genic {np.mean(h3[ref & genic & ~vsa] >= 1):.3f}", flush=True)
    rows = []
    for het_cut in HETS:
        het = h3 >= het_cut
        tests = {"A all bins: het - nonhet": (het, ~het),
                 "B genic: het - nonhet": (genic & het, genic & ~het),
                 "C1 within het: VSA - other genic": (genic & het & vsa, genic & het & ~vsa),
                 "C2 within non-het: VSA - other genic": (genic & ~het & vsa, genic & ~het & ~vsa)}
        for s in SETS:
            y = T[f"log2_{s}"]
            for name, (a, b) in tests.items():
                c0 = contrast(y, a, b, ones)
                bs = [contrast(y, a, b, w) for w in boots]
                rows.append(dict(het_cutoff=het_cut, dataset=s, test=name, n_a=int((ref & a).sum()), n_b=int((ref & b).sum()),
                                 equalGC_diff=c0, lo=np.nanpercentile(bs, 2.5), hi=np.nanpercentile(bs, 97.5)))
        print("done", het_cut, flush=True)
    R = pd.DataFrame(rows).round(4); R.to_csv(OUT / "h3k9me3_contrast.csv", index=False)
    piv = R[R.het_cutoff == 1.0].pivot(index="dataset", columns="test", values="equalGC_diff")
    (HERE / "COMPARISON.md").write_text("\n".join(["# U-DNA-Seq vs H3K9me3 state at equal GC", "", __doc__.split("Data")[1].split("Outputs")[0].strip(), "",
        f"Het fraction of nuclear bins (score >= 1): {np.mean(h3[ref] >= 1):.3f}; of VSA-gene bins: {np.mean(h3[ref & vsa] >= 1):.2f}; of other genic bins: {np.mean(h3[ref & genic & ~vsa] >= 1):.3f}.", "",
        "Point estimates at HET = 1.0:", "", piv.round(3).to_markdown(), "", "All cut-offs with 95 % CI:", "", R.to_markdown(index=False), ""]))


if __name__ == "__main__":
    main()
