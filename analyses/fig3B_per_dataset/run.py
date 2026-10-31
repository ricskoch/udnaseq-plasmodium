#!/usr/bin/env python3
"""Fig. 3B per dataset - uracil along var and rif/stevor genes in rep1, rep2, DMSO, RAPA (and merged).

Purpose
    The merged 3D7 track behind Fig. 3B is ~89 % rep1. Test whether the two claims hold in each dataset separately:
      (a) var genes: uracil is concentrated in the first ~5-6 kb (exon 1) rather than the 3' end of the gene
      (b) rif/stevor: uracil is enriched in the gene body (exon 2) relative to the flanks
Method
    Per-dataset 10 bp tracks from analyses/00c_tracks_10bp (mean of the GEO 1 bp tracks; merged = GEO merged track).
    Gene layout as in analyses/fig3B_vsa_metaprofiles: var 2 kb | 10 kb scaled body | 2 kb; rif/stevor 1 kb | 1.5 kb | 1 kb;
    minus-strand genes reversed. Bin value = mean of non-NaN base-pair values per bin edge (equal-width split of the gene).
    Per gene: mean log2U in region A and region B; effect = mean over genes of (A - B); 95 % CI by bootstrapping genes
    (2000 resamples, fixed seed); Wilcoxon signed-rank on the paired gene values; fraction of genes with A > B.
      var        A = scaled 0-6 kb   B = scaled 6-10 kb
      rif/stevor A = scaled 30-90 % of body (exon 2 side)   B = both 1 kb flanks (mean)
                 A2 = scaled 30-90 % vs first 30 % of body
    Also: peak position of the mean profile, and Pearson correlation of each dataset's mean profile with merged / with the others.
Outputs  output/per_dataset_effects.csv, per_dataset_profiles.csv, fig3B_per_dataset_profiles.png/svg, COMPARISON.md
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from udnaseq import config, io  # noqa: E402
from udnaseq.report import SAMPLE_COLORS, Report, savefig  # noqa: E402

OUT = HERE / "output"
warnings.filterwarnings("ignore", message="Mean of empty slice")
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
COL = SAMPLE_COLORS | {"3d7_merged": "#000000"}
BS, N_BOOT, SEED = 10, 2000, 20240209
LAYOUT = {"var": dict(bed=config.VAR_BED, up=2000, body=10000, down=2000),
          "rifstevor": dict(bed=config.RIFSTEVOR_BED, up=1000, body=1500, down=1000)}


def region_bins(tr10, c, a, b, n):
    """n equal-width bins over [a, b) from a 10 bp track; each 10 bp bin is assigned to the output bin of its centre."""
    v = tr10[c]
    if a < 0 or b > len(v) * BS:
        return np.full(n, np.nan)
    ctr = (np.arange(a // BS, -(-b // BS)) * BS + BS / 2)
    ctr = ctr[(ctr >= a) & (ctr < b)]
    idx = np.minimum(((ctr - a) / (b - a) * n).astype(int), n - 1)
    vals = v[(ctr // BS).astype(int)]
    out = np.full(n, np.nan)
    for k in np.unique(idx):
        x = vals[idx == k]
        x = x[~np.isnan(x)]
        if len(x):
            out[k] = x.mean()
    return out


def matrix(tr10, reg, up, body, down):
    rows = []
    for r in reg.itertuples():
        s, e = int(r.start), int(r.end)
        l, rt = (down, up) if r.strand == "-" else (up, down)
        v = np.concatenate([region_bins(tr10, r.chrom, s - l, s, l // BS), region_bins(tr10, r.chrom, s, e, body // BS),
                            region_bins(tr10, r.chrom, e, e + rt, rt // BS)])
        rows.append(v[::-1] if r.strand == "-" else v)
    return np.vstack(rows)


def boot_ci(d, rng):
    d = d[~np.isnan(d)]
    idx = rng.integers(0, len(d), (N_BOOT, len(d)))
    m = d[idx].mean(1)
    return np.quantile(m, .025), np.quantile(m, .975)


def effect(name, A, B, rng, extra):
    d = A - B
    ok = ~np.isnan(d)
    lo, hi = boot_ci(d, rng)
    w = stats.wilcoxon(d[ok])
    return dict(name=name, n_genes=int(ok.sum()), A_mean=np.nanmean(A), B_mean=np.nanmean(B), diff=np.nanmean(d),
                diff_ci_lo=lo, diff_ci_hi=hi, frac_genes_A_gt_B=float(np.mean(d[ok] > 0)), wilcoxon_p=w.pvalue,
                cohen_dz=np.nanmean(d) / np.nanstd(d, ddof=1), **extra)


def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 3B per dataset - var and rif/stevor metaprofiles",
                 "Do the var (first ~6 kb) and rif/stevor (exon 2) enrichments hold in each dataset separately? "
                 "10 bp tracks; genes scaled as in the merged figure.")
    z = np.load(io.config.ANALYSES / "00c_tracks_10bp" / "output" / "tracks_10bp.npz")
    rng = np.random.default_rng(SEED)
    eff, prof = [], []
    fig, axs = plt.subplots(1, 2, figsize=(9, 3))
    for ax, (gs, p) in zip(axs, LAYOUT.items()):
        reg = io.read_bed(p["bed"])
        reg = reg[reg.chrom.isin(config.NUCLEAR_CHROMS)].reset_index(drop=True)
        u, b, d = p["up"] // BS, p["body"] // BS, p["down"] // BS
        mus = {}
        for s in SETS:
            tr = {c: z[f"{s}|{c}"] for c in config.NUCLEAR_CHROMS}
            M = matrix(tr, reg, p["up"], p["body"], p["down"])
            mu, n = np.nanmean(M, 0), np.sum(~np.isnan(M), 0)
            sem = np.nanstd(M, 0, ddof=1) / np.sqrt(n)
            mus[s] = mu
            x = np.arange(M.shape[1]) * BS - p["up"]
            prof.append(pd.DataFrame({"geneset": gs, "dataset": s, "position_bp_scaled": x, "mean": mu, "sem": sem, "n": n}))
            ax.plot(x, mu, c=COL[s], lw=1.1 if s != "3d7_merged" else 1.4, ls="-" if s != "3d7_merged" else "--", label=s)
            body = M[:, u:u + b]
            sm = np.convolve(np.nan_to_num(np.nanmean(body, 0)), np.ones(5) / 5, "same")
            ex = dict(geneset=gs, dataset=s, peak_pos_frac_of_body=float(np.argmax(sm[2:-2]) + 2) / b)
            if gs == "var":
                A, B = np.nanmean(body[:, :600], 1), np.nanmean(body[:, 600:], 1)
                eff.append(effect("var: 0-6 kb minus 6-10 kb", A, B, rng, ex))
                eff.append(effect("var: 0-6 kb body minus flanks", A, np.nanmean(np.hstack([M[:, :u], M[:, u + b:]]), 1), rng, ex))
            else:
                A = np.nanmean(body[:, int(.3 * b):int(.9 * b)], 1)
                fl = np.nanmean(np.hstack([M[:, :u], M[:, u + b:]]), 1)
                eff.append(effect("rif/stevor: body 30-90 % minus flanks", A, fl, rng, ex))
                eff.append(effect("rif/stevor: body 30-90 % minus first 30 %", A, np.nanmean(body[:, :int(.3 * b)], 1), rng, ex))
        for t in (0, p["body"]):
            ax.axvline(t, c="k", lw=.5, ls=":")
        ax.set_title(f"{gs} (n = {len(reg)})", fontsize=9)
        ax.set_xlabel("bp (body scaled)"); ax.axhline(0, c="grey", lw=.4)
        mm = pd.DataFrame(mus)
        for s in config.SAMPLE_ORDER:
            rep.note(f"{gs}: mean-profile correlation {s} vs merged r = {mm[s].corr(mm['3d7_merged']):.2f}; "
                     f"vs rep1 {mm[s].corr(mm['3d7_rep1']):.2f}; vs rep2 {mm[s].corr(mm['3d7_rep2']):.2f}")
    axs[0].set_ylabel("mean log2(U-enrichment)"); axs[1].legend(fontsize=7, frameon=False)
    savefig(fig, OUT, "fig3B_per_dataset_profiles")
    E = pd.DataFrame(eff)
    E.round(5).to_csv(OUT / "per_dataset_effects.csv", index=False)
    pd.concat(prof).round(5).to_csv(OUT / "per_dataset_profiles.csv", index=False)
    for r in E.itertuples():
        holds = r.diff_ci_lo > 0 and r.wilcoxon_p < 0.05
        rep.check(f"{r.dataset}: {r.name}", "merged claim: positive" if "flanks" in r.name or "first" in r.name or "kb minus" in r.name else "",
                  f"{r.diff:+.3f} ({r.diff_ci_lo:+.3f}, {r.diff_ci_hi:+.3f})", "HOLDS" if holds else "NOT SUPPORTED",
                  f"n genes {r.n_genes}; A {r.A_mean:+.3f} vs B {r.B_mean:+.3f}; {r.frac_genes_A_gt_B:.0%} of genes A > B; "
                  f"Wilcoxon p = {r.wilcoxon_p:.1e}; dz = {r.cohen_dz:.2f}; mean-profile peak at {r.peak_pos_frac_of_body:.0%} of body")
    rep.note("Genes are the unit of replication (bootstrap over genes); rep1/rep2/DMSO/RAPA are single libraries, so a dataset "
             "'holds' only means the pattern is present in that library.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
