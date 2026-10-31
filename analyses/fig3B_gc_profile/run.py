#!/usr/bin/env python3
"""Fig. 3B vs GC - does the position-along-gene pattern of uracil simply follow the GC profile of the genes?

Purpose
    var genes: uracil 0-6 kb > 6-10 kb. rif/stevor: body > flanks. If GC content shows the same position pattern and
    log2U tracks GC along the gene, the profile is not independent evidence of a uracil-specific mechanism.
Method
    Same gene layout and scaling as fig3B_per_dataset. A 10 bp GC-fraction track is built from the genome FASTA.
    For each dataset: mean-profile Pearson r between log2U and GC (scaled bins), and per-gene log2U ~ GC slope within the
    var body (pooled over genes: bins' GC and log2U, gene-wise demeaned). Also GC difference (0-6 kb minus 6-10 kb;
    body minus flanks) for comparison with the log2U difference.
Outputs  output/gc_profile_summary.csv, COMPARISON.md
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import Report  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


f3 = load("f3", "fig3B_per_dataset/run.py")
f4b = load("f4b", "fig4B_gc_correlation/run.py")
BS = 10


def main():
    (HERE / "output").mkdir(exist_ok=True)
    rep = Report("Fig. 3B vs GC profile", "Does log2U along var and rif/stevor genes follow the GC profile of the genes?")
    seqs = f4b.read_fasta(config.FASTA)
    gc10 = {}
    for c in config.NUCLEAR_CHROMS:
        s = seqs[c]
        n = len(s) // BS + (len(s) % BS > 0)
        is_gc = ((s == b"G") | (s == b"C")).astype(float)
        pad = np.zeros(n * BS)
        pad[:len(s)] = is_gc
        cnt = np.full(n * BS, 0.0)
        cnt[:len(s)] = 1
        gc10[c] = pad.reshape(n, BS).sum(1) / np.maximum(cnt.reshape(n, BS).sum(1), 1)
    z = np.load(config.ANALYSES / "00c_tracks_10bp" / "output" / "tracks_10bp.npz")
    rows = []
    for gs, p in f3.LAYOUT.items():
        reg = io.read_bed(p["bed"])
        reg = reg[reg.chrom.isin(config.NUCLEAR_CHROMS)].reset_index(drop=True)
        u, b = p["up"] // BS, p["body"] // BS
        G = f3.matrix(gc10, reg, p["up"], p["body"], p["down"])
        if gs == "var":
            split = lambda M: (np.nanmean(M[:, u:u + 600], 1), np.nanmean(M[:, u + 600:u + b], 1))
        else:
            split = lambda M: (np.nanmean(M[:, u + int(.3 * b):u + int(.9 * b)], 1),
                               np.nanmean(np.hstack([M[:, :u], M[:, u + b:]]), 1))
        ga, gb = split(G)
        for s in config.SAMPLE_ORDER + ["3d7_merged"]:
            M = f3.matrix({c: z[f"{s}|{c}"] for c in config.NUCLEAR_CHROMS}, reg, p["up"], p["body"], p["down"])
            mu, gm = np.nanmean(M, 0), np.nanmean(G, 0)
            ok = ~np.isnan(mu) & ~np.isnan(gm)
            ua, ub = split(M)
            # within-gene association: demean both by gene, pool bins
            Md, Gd = M - np.nanmean(M, 1, keepdims=True), G - np.nanmean(G, 1, keepdims=True)
            o2 = ~np.isnan(Md) & ~np.isnan(Gd)
            r_within = np.corrcoef(Md[o2], Gd[o2])[0, 1]
            rows.append(dict(geneset=gs, dataset=s, r_mean_profile_log2U_vs_GC=np.corrcoef(mu[ok], gm[ok])[0, 1],
                             r_within_gene_bins=r_within, gc_diff=np.nanmean(ga - gb), log2U_diff=np.nanmean(ua - ub),
                             gc_A=np.nanmean(ga), gc_B=np.nanmean(gb)))
    D = pd.DataFrame(rows)
    D.round(4).to_csv(HERE / "output" / "gc_profile_summary.csv", index=False)
    for r in D.itertuples():
        rep.check(f"{r.dataset}: {r.geneset}", "-", f"mean-profile r(log2U, GC) = {r.r_mean_profile_log2U_vs_GC:.2f}; within-gene bins r = {r.r_within_gene_bins:.2f}",
                  "INFO", f"GC region A {r.gc_A:.1%} vs B {r.gc_B:.1%} (diff {r.gc_diff:+.1%}); log2U diff {r.log2U_diff:+.2f}")
    rep.note("Regions: var A = scaled 0-6 kb, B = 6-10 kb; rif/stevor A = 30-90 % of body, B = flanks.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
