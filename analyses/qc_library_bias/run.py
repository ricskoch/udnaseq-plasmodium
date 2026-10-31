#!/usr/bin/env python3
"""QC - library depth, evenness and GC bias of IP and INPUT per sample (why do the datasets disagree on G4/GC hotspots?).

Purpose
    DMSO shows no hotspot-G4 enrichment while rep1 does (fig4C_per_dataset). Check whether depth, coverage evenness or a
    GC-dependent coverage bias in IP/INPUT could explain it.
Method
    From the 100 bp RPGC coverage tracks (nuclear, blacklist removed): mean, CV (sd/mean) and fraction of zero bins of IP and
    INPUT; correlation IP vs INPUT; slope of coverage vs GC% (per 10 GC points) for IP and INPUT; mean coverage by GC class.
    SRA read-pair counts (PRJNA1076816, NovaSeq 6000, 2 x 150, taken from NCBI SRA run info) give the sequencing depth
    before duplicate removal; depth = pairs x 301 bp / 23.3 Mb.
Outputs  output/library_qc.csv, output/coverage_by_gc.csv, COMPARISON.md
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

spec = importlib.util.spec_from_file_location("f4b", HERE.parents[0] / "fig4B_gc_correlation" / "run.py")
f4b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f4b)

# NCBI SRA run info, BioProject PRJNA1076816 (read pairs before duplicate removal)
SRA = {"3d7_rep1": dict(ip=("SRR27975802", 31566218), input=("SRR27975806", 23935477)),
       "3d7_rep2": dict(ip=("SRR27975801", 9140590), input=("SRR27975805", 21068631)),
       "DMSO": dict(ip=("SRR27975800", 19094719), input=("SRR27975804", 10733351)),
       "RAPA": dict(ip=("SRR27975799", 20317860), input=("SRR27975803", 14003413))}
GENOME_MB, FRAG = 23.3, 301
EDGES = [0, .15, .25, .35, .5, 1.01]


def main():
    rep = Report("QC - depth, evenness and GC bias per library", "Does depth or GC-dependent coverage bias differ between datasets?")
    T = io.load_tracks_npz("tracks_100bp")
    seqs = f4b.read_fasta(config.FASTA)
    gc = np.full(len(T["chrom"]), np.nan)
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        gc[m] = f4b.gc_fraction(seqs[str(c)], T["start"][m], T["end"][m])
    nuc = np.isin(T["chrom"], config.NUCLEAR_CHROMS) & ~T["blacklisted"]
    g = gc[nuc]
    rows, byg = [], []
    slope = lambda v: np.polyfit(g * 100, v, 1)[0] * 10
    for s in config.SAMPLE_ORDER:
        ip, inp = T[f"ip_{s}"][nuc], T[f"input_{s}"][nuc]
        r = dict(sample=s, ip_pairs_M=SRA[s]["ip"][1] / 1e6, input_pairs_M=SRA[s]["input"][1] / 1e6,
                 ip_depth_x=SRA[s]["ip"][1] * FRAG / GENOME_MB / 1e6, input_depth_x=SRA[s]["input"][1] * FRAG / GENOME_MB / 1e6,
                 ip_cv=ip.std() / ip.mean(), input_cv=inp.std() / inp.mean(), ip_zero_pct=np.mean(ip == 0) * 100,
                 input_zero_pct=np.mean(inp == 0) * 100, r_ip_input=np.corrcoef(ip, inp)[0, 1],
                 input_gc_slope_per10=slope(inp), ip_gc_slope_per10=slope(ip), log2U_sd=T[f"log2_{s}"][nuc].std())
        rows.append(r)
        for a, b in zip(EDGES[:-1], EDGES[1:]):
            m = (g >= a) & (g < b)
            byg.append(dict(sample=s, gc_class=f"{a:.0%}-{min(b, 1):.0%}", n_bins=int(m.sum()), input_mean=inp[m].mean(), ip_mean=ip[m].mean()))
        rep.check(f"{s}: depth / evenness", "-",
                  f"IP {r['ip_pairs_M']:.1f} M pairs ({r['ip_depth_x']:.0f}x), INPUT {r['input_pairs_M']:.1f} M ({r['input_depth_x']:.0f}x)", "INFO",
                  f"CV IP {r['ip_cv']:.2f} / INPUT {r['input_cv']:.2f}; zero bins IP {r['ip_zero_pct']:.1f} % / INPUT {r['input_zero_pct']:.1f} %; "
                  f"r(IP, INPUT) {r['r_ip_input']:.2f}; coverage slope per 10 GC points: INPUT {r['input_gc_slope_per10']:+.2f}, IP {r['ip_gc_slope_per10']:+.2f}")
    pd.DataFrame(rows).round(4).to_csv(HERE / "output" / "library_qc.csv", index=False) if (HERE / "output").exists() else None
    (HERE / "output").mkdir(exist_ok=True)
    pd.DataFrame(rows).round(4).to_csv(HERE / "output" / "library_qc.csv", index=False)
    pd.DataFrame(byg).round(4).to_csv(HERE / "output" / "coverage_by_gc.csv", index=False)
    rep.note("RPGC coverage is scaled to 1x genome content, so means are ~1.5 for every library; depth comes from the SRA counts.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
