#!/usr/bin/env python3
"""Fig. 4C - Distance between highly uracilated regions and predicted G-quadruplexes.

Purpose
    Test whether uracil hotspots (log2U > 0.5) sit on/near predicted G4-forming
    sequences, supporting cytosine deamination opposite G4s as the uracil source.

Method (replaces bedtools closest -D ref + R/tidyverse)
    U hotspots : 10 bp bins of the GEO merged 3D7 track with log2U > 0.5.
    G4 sites   : data/BED_files/G4_coordinates.bed (G4Hunter, Gazanion et al. 2020).
    Distance   : signed distance to the closest G4 relative to the reference strand
                 (0 = overlap, ties all reported) - bedtools closest -D ref semantics.
    published  : closest G4 among ALL sites, then keep pairs with |G4Hunter| > 1.2 and
                 G4 size <= 200 bp  (reproduces the manuscript numbers)
    prefiltered: filter G4 sites first, then find the closest one (recommended: every
                 hotspot keeps a partner; no hotspot is silently dropped)
    null       : the same number of 10 bp bins placed at random on the nuclear chromosomes
                 (blacklist excluded), 200 permutations; fraction overlapping / within 100 bp.
    Also reproduces the legacy 100 bp version (legacy/g4coordinates/.../0p5_to_G4_coordinates.bed).
Outputs  output/*.csv, fig4C_*.png/svg, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402
from udnaseq.report import Report, savefig, status  # noqa: E402

OUT = HERE / "output"
LEG100 = config.LEGACY / "g4coordinates" / "g4coordinates" / "0p5_to_G4_coordinates.bed"
N_PERM = 200
SEED = 20240209


def g4_sites():
    g = pd.read_csv(config.G4_SITES, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"])
    return g.assign(size=g.end - g.start)[["chrom", "start", "end", "score", "size"]]


def hotspots_10bp():
    pb = io.read_bigwig(io.geo_path(config.GEO_MERGED_10BP))
    rows = []
    for c, v in pb.items():
        b = io.bin_mean(v, 10)
        s = np.arange(len(b)) * 10
        m = b > config.U_THRESHOLD
        rows.append(pd.DataFrame({"chrom": c, "start": s[m], "end": np.minimum(s[m] + 10, len(v)), "value": b[m]}))
    return pd.concat(rows, ignore_index=True)


def summarise(d) -> dict:
    d = np.asarray(d, float)
    q25, q75 = np.quantile(d, [0.25, 0.75])  # numpy 'linear' == R quantile type 7
    return {"n": len(d), "median": np.median(d), "q25": q25, "q75": q75, "IQD": q75 - q25, "mean": d.mean(),
            "frac_overlap": np.mean(d == 0), "frac_within_100bp": np.mean(np.abs(d) <= 100)}


def abs_distance_fast(chrom_starts: dict, g4: pd.DataFrame, width=10) -> np.ndarray:
    """Unsigned bedtools-style distance of [s, s+width) to the nearest G4 (vectorised, for permutations)."""
    out = []
    for c, s in chrom_starts.items():
        g = g4[g4.chrom == c].sort_values("start")
        if g.empty:
            continue
        gs, ge = g.start.values, np.maximum.accumulate(g.end.values)
        e = s + width
        i = np.searchsorted(gs, e, side="left")  # G4s starting before bin end: candidates for overlap/upstream
        left_end = np.where(i > 0, ge[np.maximum(i - 1, 0)], -10**12)
        d_left = np.where(left_end > s, 0, s - left_end + 1)
        right_start = np.where(i < len(gs), gs[np.minimum(i, len(gs) - 1)], 10**12)
        d_right = right_start - e + 1
        out.append(np.minimum(d_left, d_right))
    return np.concatenate(out)


def permutation_null(n_bins, g4, rng):
    sizes = io.chrom_sizes()
    bl = genome.blacklist()
    chroms = config.NUCLEAR_CHROMS
    L = np.array([sizes[c] for c in chroms], float)
    res = []
    for _ in range(N_PERM):
        counts = rng.multinomial(n_bins, L / L.sum())
        starts = {}
        for c, k in zip(chroms, counts):
            s = rng.integers(0, sizes[c] // 10, size=k * 2) * 10
            s = s[~genome.blacklist_mask(c, s, s + 10, bl)][:k]
            starts[c] = s
        d = abs_distance_fast(starts, g4)
        res.append({"frac_overlap": np.mean(d == 0), "frac_within_100bp": np.mean(d <= 100), "median_abs": np.median(d)})
    return pd.DataFrame(res)


def hist(d, title, xlim=1000):
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.hist(d, bins=np.arange(-xlim - 15, xlim + 16, 30), color="grey", edgecolor="white", lw=0.2)
    ax.set_xlim(-xlim, xlim)
    ax.set_xlabel("Distance to closest predicted G4 (bp)")
    ax.set_ylabel("Count")
    s = summarise(d)
    ax.set_title(f"{title}\nn = {s['n']}, median = {s['median']:g}, IQR = ({s['q25']:g}, {s['q75']:g})", fontsize=8)
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 4C - uracil hotspots vs predicted G4 sites",
                 "Signed distance from 10 bp bins with log2U > 0.5 (merged 3D7) to the closest predicted G4.")
    g4 = g4_sites()
    U = hotspots_10bp()
    U.round(5).to_csv(OUT / "hotspots_10bp_log2U_gt0.5.csv", index=False)

    # ---- published procedure -----------------------------------------------------------
    allc = genome.closest_signed(U, g4)
    pub = allc[(allc.score_b.abs() > config.G4_MIN_SCORE) & (allc.size_b <= config.G4_MAX_SIZE)]
    pub.to_csv(OUT / "closest_published_filter_after.csv.gz", index=False)
    sp = summarise(pub.distance)
    rep.check("hotspot-G4 pairs kept", "-", f"{sp['n']} of {len(allc)} pairs ({len(U)} hotspot bins)", "INFO")
    rep.check("median distance", "0", f"{sp['median']:g}", "MATCH" if sp["median"] == 0 else "DIFFERS")
    rep.check("quartiles (IQD)", "(-206.5, 389.5), IQD 596", f"({sp['q25']:g}, {sp['q75']:g}), IQD {sp['IQD']:g}",
              "MATCH" if (sp["q25"], sp["q75"]) == (-206.5, 389.5) else "DIFFERS")
    rep.check("filtered G4 size mean / median", "33.13 / 29", f"{pub.size_b.mean():.2f} / {pub.size_b.median():g}",
              "MATCH" if round(pub.size_b.mean(), 2) == 33.13 and pub.size_b.median() == 29 else "DIFFERS")
    savefig(hist(pub.distance, "Published procedure (closest G4, then filter pairs)"), OUT, "fig4C_published")

    # ---- recommended: filter G4 first ------------------------------------------------------
    strong = g4[(g4.score.abs() > config.G4_MIN_SCORE) & (g4["size"] <= config.G4_MAX_SIZE)]
    pre = genome.closest_signed(U, strong)
    pre.to_csv(OUT / "closest_prefiltered_G4.csv.gz", index=False)
    spr = summarise(pre.distance)
    rep.check("prefiltered G4 (recommended): median, quartiles", "-",
              f"{spr['median']:g}, ({spr['q25']:g}, {spr['q75']:g}), n = {spr['n']}", "INFO",
              f"published filter drops {len(U) - pub[['chrom_a', 'start_a']].drop_duplicates().shape[0]} of {len(U)} hotspot bins")
    savefig(hist(pre.distance, "Recommended (strong G4s first, then closest)"), OUT, "fig4C_prefiltered")

    # ---- permutation null ------------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    Un = U[U.chrom.isin(config.NUCLEAR_CHROMS)]
    Un = Un[~genome.add_blacklist_flag(Un).blacklisted.values]
    obs = abs_distance_fast({c: g.start.values for c, g in Un.groupby("chrom")}, strong)
    null = permutation_null(len(Un), strong, rng)
    null.to_csv(OUT / "permutation_null.csv", index=False)
    for k, o in [("frac_overlap", np.mean(obs == 0)), ("frac_within_100bp", np.mean(obs <= 100))]:
        p = (1 + (null[k] >= o).sum()) / (1 + N_PERM)
        rep.check(f"null model: {k} (strong G4s, nuclear, blacklist removed)",
                  f"random: {null[k].mean():.3f} (95% {null[k].quantile(.025):.3f}-{null[k].quantile(.975):.3f})",
                  f"observed: {o:.3f} ({o / null[k].mean():.1f}x)", "INFO", f"permutation p = {p:.3g} ({N_PERM} perms)")
    pd.DataFrame([dict(set="published", **sp), dict(set="prefiltered", **spr)]).round(4).to_csv(OUT / "distance_summary.csv", index=False)

    # ---- legacy 100 bp version -------------------------------------------------------------
    T = io.load_tracks_npz("tracks_100bp")
    m = (~T["blacklisted"]) & (T["log2_3d7_merged"] > config.U_THRESHOLD)
    U100 = pd.DataFrame({"chrom": T["chrom"][m], "start": T["start"][m], "end": T["end"][m], "value": T["log2_3d7_merged"][m]})
    c100 = genome.closest_signed(U100, g4)
    leg = pd.read_csv(LEG100, sep="\t", header=None)
    same = (len(leg) == len(c100)) and np.array_equal(leg[10].values, c100.distance.values) \
        and np.array_equal(leg[1].values, c100.start_a.values)
    rep.check("legacy 100 bp closest table", f"{len(leg)} rows", f"{len(c100)} rows",
              "MATCH" if same else "DIFFERS", "same hotspot bins, partners and signed distances" if same else "")
    rep.note("Order matters: the published procedure finds the closest G4 among ALL predictions and then discards "
             "hotspots whose closest G4 is weak (|score| <= 1.2) or large (> 200 bp), instead of re-matching them to "
             "the closest strong G4. Filtering G4s first keeps every hotspot; the median stays 0 but the IQR widens. "
             "Suggest describing/using the prefiltered version.")
    rep.note("The manuscript reports no background expectation. The permutation null above (random 10 bp bins, "
             "same chromosomes, blacklist removed) gives that reference; random bins are not GC-matched, and G4s "
             "are GC-rich, so part of the excess is expected from the GC correlation in Fig. 4B.")
    rep.note("'Interquartile distance' in the manuscript is the IQR of the signed distances.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
