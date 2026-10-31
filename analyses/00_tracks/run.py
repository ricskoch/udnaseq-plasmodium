#!/usr/bin/env python3
"""00 - Build the log2(U-enrichment) tracks used by all figure analyses.

Purpose
    Regenerate every derived U-DNA-Seq track from the raw-most data in the repo,
    the RPGC-normalised IP and INPUT coverage bigWigs (100 bp bins).

Method (replaces deepTools bigwigCompare + bedtools subtract)
    log2U = log2((IP + 1) / (INPUT + 1))                    per sample, per 100 bp bin
    merged 3D7 = log2((w*IP1 + (1-w)*IP2 + 1) / (u*IN1 + (1-u)*IN2 + 1))
        The legacy merged track came from merged BAMs, i.e. replicates weighted by
        their read share.  BAM read counts are not in the repo, so w and u are
        estimated by least squares against the legacy merged bedGraph.
    Blacklisted bins = bins overlapping data/BED_files/blacklist_telomer_pf11.bed.

Outputs (output/)
    tracks_100bp.npz                 all bins x {IP, INPUT, log2} per sample + merged  (git-ignored)
    <sample>.log2.bin100bp.bw        genome browser tracks                             (git-ignored)
    <sample>.log2.bin100bp.blacklisted.bedGraph                                        (git-ignored)
    merge_weights.json               estimated replicate weights of the merged track
    COMPARISON.md                    agreement with legacy/data/log2_bw/*
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402
from udnaseq.report import Report, status  # noqa: E402

OUT = HERE / "output"
LEG = config.LEGACY / "data" / "log2_bw"
LEGACY_SAMPLE = {
    "3d7_rep1": ("3d7_rep1.MAPQfiltered.IPvsINPUT.bin100bp.log2.bw", "3d7_rep1.IP.vs.INPUT.blacklisted.bin100bp.log2.bedgraph"),
    "3d7_rep2": ("3d7_rep2.MAPQfiltered.IPvsINPUT.bin100bp.log2.bw", "3d7_rep2.IP.vs.INPUT.blacklisted.bin100bp.log2.bedgraph"),
    "DMSO": ("DMSO.MAPQfiltered.IPvsINPUT.bin100bp.log2.bw", "DMSO.IP.vs.INPUT.blacklisted.bin100bp.log2.bedgraph"),
    "RAPA": ("RAPA.MAPQfiltered.IPvsINPUT.bin100bp.log2.bw", "RAPA.IP.vs.INPUT.blacklisted.bin100bp.log2.bedgraph"),
}


def log2ratio(ip, inp):
    return np.log2((ip + config.PSEUDOCOUNT) / (inp + config.PSEUDOCOUNT))


def expand_to_bins(bg: pd.DataFrame, binsize=config.BIN) -> pd.Series:
    """Split merged bedGraph runs back into 100 bp bins -> Series indexed by (chrom, start)."""
    n = ((bg.end - bg.start + binsize - 1) // binsize).astype(int)
    rep = bg.loc[bg.index.repeat(n)].copy()
    rep["start"] = rep.start + rep.groupby(level=0).cumcount() * binsize
    return pd.Series(rep.value.values, index=pd.MultiIndex.from_arrays([rep.chrom, rep.start]))


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("00 - log2(U-enrichment) tracks",
                 "Per-sample and merged 3D7 log2(IP/INPUT) tracks rebuilt from the RPGC coverage bigWigs "
                 "(data/*.RPGC.bw) and compared bin-by-bin with the legacy tracks in legacy/data/log2_bw/.")

    # ---- load coverage and build the bin table --------------------------------
    cov = {}
    for s in config.SAMPLE_ORDER:
        cov[s] = {k: io.read_binned_bigwig(config.SAMPLES[s][k]) for k in ("ip", "input")}
    frame = io.tracks_to_frame(cov[config.SAMPLE_ORDER[0]]["ip"], config.BIN)[["chrom", "start", "end"]]
    frame = genome.add_blacklist_flag(frame)
    flat = lambda d: np.concatenate([d[c] for c in sorted(d)])
    for s in config.SAMPLE_ORDER:
        frame[f"ip_{s}"] = flat(cov[s]["ip"])
        frame[f"input_{s}"] = flat(cov[s]["input"])
        assert not frame[[f"ip_{s}", f"input_{s}"]].isna().any().any(), "RPGC tracks should cover all bins"
        frame[f"log2_{s}"] = log2ratio(frame[f"ip_{s}"], frame[f"input_{s}"])

    # ---- merged 3D7 track: fit replicate weights ------------------------------
    leg_m = io.read_bedgraph(LEG / "3d7_merged.IP.vs.INPUT.bin100bp.log2.bedGraph")
    leg_m = leg_m[(leg_m.end - leg_m.start) == config.BIN]
    key = pd.MultiIndex.from_arrays([frame.chrom, frame.start])
    j = key.get_indexer(pd.MultiIndex.from_arrays([leg_m.chrom, leg_m.start]))
    ok = j >= 0
    y = leg_m.value.values[ok]
    ip1, ip2 = frame["ip_3d7_rep1"].values[j[ok]], frame["ip_3d7_rep2"].values[j[ok]]
    in1, in2 = frame["input_3d7_rep1"].values[j[ok]], frame["input_3d7_rep2"].values[j[ok]]
    model = lambda p: log2ratio(p[0] * ip1 + (1 - p[0]) * ip2, p[1] * in1 + (1 - p[1]) * in2)
    fit = minimize(lambda p: np.mean((model(p) - y) ** 2), x0=[0.5, 0.5], bounds=[(0, 1), (0, 1)])
    w, u = fit.x
    resid = np.abs(model(fit.x) - y)
    weights = {"ip_weight_rep1": round(float(w), 4), "input_weight_rep1": round(float(u), 4),
               "n_bins_fitted": int(ok.sum()), "pearson_r": float(np.corrcoef(model(fit.x), y)[0, 1]),
               "median_abs_residual": float(np.median(resid)), "p99_abs_residual": float(np.quantile(resid, .99))}
    (OUT / "merge_weights.json").write_text(json.dumps(weights, indent=2))
    frame["ip_3d7_merged"] = w * frame["ip_3d7_rep1"] + (1 - w) * frame["ip_3d7_rep2"]
    frame["input_3d7_merged"] = u * frame["input_3d7_rep1"] + (1 - u) * frame["input_3d7_rep2"]
    frame["log2_3d7_merged"] = log2ratio(frame["ip_3d7_merged"], frame["input_3d7_merged"])

    simple = log2ratio((ip1 + ip2) / 2, (in1 + in2) / 2)
    rep.check("merged 3D7 (100 bp): equal-weight average of replicates", "legacy bedGraph",
              f"r = {np.corrcoef(simple, y)[0, 1]:.3f}", "DIFFERS",
              "the legacy merged track is not a simple average")
    rep.check("merged 3D7 (100 bp): read-share weighted merge", "legacy bedGraph",
              f"r = {weights['pearson_r']:.6f}, 99th pct |diff| = {weights['p99_abs_residual']:.1e}",
              status(weights["p99_abs_residual"]), f"IP weight rep1 = {w:.3f}, INPUT weight rep1 = {u:.3f}")
    rep.note(f"The merged 3D7 track is dominated by replicate 1 on the IP side ({w:.0%} of IP reads vs "
             f"{u:.0%} of INPUT reads). Results labelled 'merged' therefore mostly reflect replicate 1.")

    # ---- compare per-sample tracks with legacy --------------------------------
    for s, (bwf, bgf) in LEGACY_SAMPLE.items():
        leg = io.read_binned_bigwig(LEG / bwf)
        lv = np.concatenate([leg[c] for c in sorted(leg)])
        new = frame[f"log2_{s}"].values
        m = ~np.isnan(lv)
        d = float(np.max(np.abs(lv[m] - new[m])))
        # deepTools drops bins whose log2 ratio is exactly 0 (IP == INPUT, mostly both empty)
        dropped = ~m
        frac_equal = np.mean(frame[f"ip_{s}"].values[dropped] == frame[f"input_{s}"].values[dropped]) if dropped.any() else 1
        rep.check(f"{s}: log2 bigWig (all bins)", bwf, f"max |diff| = {d:.1e}", status(d),
                  f"{dropped.sum()} bins absent in legacy bigWig; {frac_equal:.0%} of them have IP == INPUT (log2 = 0)")
        lb = expand_to_bins(io.read_bedgraph(LEG / bgf))
        nb = pd.Series(new, index=key)[~frame.blacklisted.values]
        common = lb.index.intersection(nb.index)
        d2 = float(np.max(np.abs(lb[common] - nb[common])))
        rep.check(f"{s}: blacklisted bedGraph", f"{len(lb)} bins", f"{len(nb)} bins, {len(common)} shared",
                  status(d2), f"max |diff| = {d2:.1e}; legacy omits log2 = 0 bins")

    # ---- blacklist sanity --------------------------------------------------------
    bl = genome.blacklist()
    sizes = io.chrom_sizes()
    missing = []
    for c in config.NUCLEAR_CHROMS:
        b = bl[bl.chrom == c]
        if not (b.start == 0).any():
            missing.append(f"{c} 5' end")
        if not (b.end >= sizes[c] - 1).any():
            missing.append(f"{c} 3' end")
    rep.note("Blacklist (telomeres + Pf11-1) does NOT cover: " + ", ".join(missing) +
             ". These subtelomeric regions are therefore included in all 'blacklisted' analyses, as in the legacy work. "
             "Kept unchanged here for reproducibility; consider extending the blacklist.")

    # ---- write outputs ---------------------------------------------------------
    cols = ["chrom", "start", "end", "blacklisted"] + [c for c in frame.columns if c.split("_")[0] in ("ip", "input", "log2")]
    np.savez_compressed(OUT / "tracks_100bp.npz", **{c: frame[c].to_numpy() if c != "chrom" else frame[c].to_numpy().astype("U16")
                                                     for c in cols})
    for s in config.SAMPLE_ORDER + ["3d7_merged"]:
        df = frame[["chrom", "start", "end"]].assign(value=frame[f"log2_{s}"])
        io.write_bigwig(df, OUT / f"{s}.log2.bin100bp.bw")
        io.write_bedgraph(df[~frame.blacklisted.values], OUT / f"{s}.log2.bin100bp.blacklisted.bedGraph")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
