#!/usr/bin/env python3
"""Fig. 2C, 3A, S1D, S2 / Table S1 - gene-level log2(U-enrichment).

These panels all derive from one table: the mean log2(U-enrichment) of each of
the 5318 protein-coding genes, so they live in one folder.

Purpose
    Table S1 / Fig. S2 : rank genes by uracilation; top 200 -> ShinyGO GO enrichment
    Fig. 2C            : VSA families (VAR, RIF, STEVOR) vs all other genes, per dataset;
                         var genes RAPA vs DMSO
    Fig. S1D           : the same gene means grouped by family, compared across datasets
    Fig. 3A            : per-var-gene mean +/- SD of the two 3D7 replicates, by var group

Method
    Per-sample gene means: GEO 1 bp log2 tracks averaged over each gene.
      legacy mode (reproduces published numbers): pyBigWig zoom-level mean over
          [GFF start, end) - i.e. the 1-based GFF start used as a 0-based start
      exact mode (recommended): exact per-base mean over [start-1, end)
    Merged 3D7 gene means (ranking for GO): 100 bp merged track from 00_tracks,
      sum over gene / gene length (bigWigAverageOverBed 'mean0', as in the legacy xlsx,
      which used the same 1 bp start offset; *_exact uses correct coordinates).
    Tests: Welch two-sample t-test (scipy, equal_var=False), Benjamini-Hochberg
      across all comparisons shown in a panel. RAPA vs DMSO for var genes also as
      a paired t-test, because both values come from the same genes.
Outputs  output/: table_S1_gene_log2U.csv, top200_for_GO.txt, stats_*.csv, figures, COMPARISON.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402
from udnaseq.report import CATEGORY_COLORS, SAMPLE_COLORS, Report, savefig, status  # noqa: E402

OUT = HERE / "output"
LEG = config.LEGACY / "gene_enrichment"
CATS = ["Others", "RIF", "STEVOR", "VAR"]
VSA = ["RIF", "STEVOR", "VAR"]
# p-values printed on legacy/gene_enrichment/UDNAenrichment_violinplot_pvalue.eps
# (unadjusted Welch t-tests; 2.2e-16 = R's display floor).  Panel -> dataset mapping
# was not recorded; it is inferred below by matching the STEVOR p-values.
LEGACY_EPS_P = {1: {"RIF": 2.2e-16, "STEVOR": 3.5e-14, "VAR": 2.2e-16},
                2: {"RIF": 2.2e-16, "STEVOR": 1.3e-15, "VAR": 2.2e-16},
                3: {"RIF": 2.2e-16, "STEVOR": 6.5e-13, "VAR": 2.2e-16}}


# ----------------------------------------------------------------- tables ---
def gene_table(genes: pd.DataFrame) -> pd.DataFrame:
    t = genes[["id", "symbol", "description", "chrom", "start", "end", "strand", "category", "category_legacy"]].copy()
    legacy_coords = genes.assign(start=genes.start + 1)  # replicate the 1 bp offset of the legacy table
    for s in config.SAMPLE_ORDER:
        p = io.geo_path(config.SAMPLES[s]["geo_1bp"])
        t[s] = io.bigwig_stats(p, legacy_coords, exact=False)
        t[f"{s}_exact"] = io.bigwig_stats(p, genes, exact=True)
    T = io.load_tracks_npz("tracks_100bp")
    merged = {}
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        merged[c] = np.repeat(T["log2_3d7_merged"][m], 100)[: int(T["end"][m][-1])]
    t["3d7_merged_mean0"] = genome.interval_means(merged, legacy_coords, empty_as_zero=True)
    t["3d7_merged_mean0_exact"] = genome.interval_means(merged, genes, empty_as_zero=True)
    return t


# ------------------------------------------------------------------ stats ---
def welch(a, b):
    a, b = np.asarray(a)[~np.isnan(a)], np.asarray(b)[~np.isnan(b)]
    r = stats.ttest_ind(a, b, equal_var=False)
    return r.statistic, r.pvalue, len(a), len(b), a.mean() - b.mean()


def r_pooled_pairwise(values, groups, g1, g2):
    """R pairwise.t.test(pool.sd = TRUE) as used in legacy gene_enrichment.R."""
    df = pd.DataFrame({"v": values, "g": groups}).dropna()
    s2 = sum(((x - x.mean()) ** 2).sum() for _, x in df.groupby("g").v) / (len(df) - df.g.nunique())
    a, b = df[df.g == g1].v, df[df.g == g2].v
    t = (a.mean() - b.mean()) / np.sqrt(s2 * (1 / len(a) + 1 / len(b)))
    return 2 * stats.t.sf(abs(t), len(df) - df.g.nunique())


def vsa_vs_others(t: pd.DataFrame, cols, catcol):
    rows = []
    for c in cols:
        for fam in VSA:
            tt, p, n1, n2, d = welch(t.loc[t[catcol] == fam, c], t.loc[t[catcol] == "Others", c])
            rows.append({"dataset": c, "family": fam, "n_family": n1, "n_others": n2, "mean_diff": d,
                         "welch_t": tt, "p": p,
                         "p_pooled_sd_unadjusted": r_pooled_pairwise(t[c], t[catcol], fam, "Others")})
    r = pd.DataFrame(rows)
    r["p_BH"] = multipletests(r.p, method="fdr_bh")[1]
    return r


def var_across_datasets(t, cols, catcol):
    v = t[t[catcol] == "VAR"]
    rows = []
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            _, p, *_ = welch(v[a], v[b])
            ok = v[[a, b]].dropna()
            rows.append({"a": a, "b": b, "n_pairs": len(ok), "mean_a": v[a].mean(), "mean_b": v[b].mean(),
                         "welch_p": p, "paired_p": stats.ttest_rel(ok[a], ok[b]).pvalue,
                         "wilcoxon_signed_rank_p": stats.wilcoxon(ok[a], ok[b]).pvalue})
    r = pd.DataFrame(rows)
    for k in ("welch_p", "paired_p", "wilcoxon_signed_rank_p"):
        r[k + "_BH"] = multipletests(r[k], method="fdr_bh")[1]
    return r


# ------------------------------------------------------------------ plots ---
def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def violin_by_dataset(t, cols, catcol, st, title):
    fig, axes = plt.subplots(1, len(cols), figsize=(2.6 * len(cols), 3.6), sharey=True)
    for ax, c in zip(axes, cols):
        data = [t.loc[t[catcol] == k, c].dropna().values for k in CATS]
        parts = ax.violinplot(data, showextrema=False, showmedians=True, widths=0.8)
        for body, k in zip(parts["bodies"], CATS):
            body.set_facecolor(CATEGORY_COLORS[k]); body.set_edgecolor("k"); body.set_alpha(0.7)
        parts["cmedians"].set_color("k")
        ax.set_xticks(range(1, 5), [f"{k}\n(n={len(d)})" for k, d in zip(CATS, data)], fontsize=7)
        ymax = 0.75
        for i, fam in enumerate(VSA):
            p = st[(st.dataset == c) & (st.family == fam)].p_BH.iloc[0]
            y = ymax + 0.1 * i
            ax.plot([1, i + 2], [y, y], c="k", lw=0.6)
            ax.text((i + 3) / 2, y, stars(p), ha="center", va="bottom", fontsize=8)
        ax.axhline(0, c="grey", lw=0.5, ls=":")
        ax.set_title(c.replace("_exact", ""), fontsize=9)
        ax.set_ylim(-1.1, 1.1)
    axes[0].set_ylabel("Mean log2(U-enrichment) per gene")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    return fig


def violin_by_family(t, cols, catcol):
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.4), sharey=True)
    for ax, k in zip(axes, CATS):
        data = [t.loc[t[catcol] == k, c].dropna().values for c in cols]
        parts = ax.violinplot(data, showextrema=False, showmedians=True)
        for body, c in zip(parts["bodies"], cols):
            body.set_facecolor(SAMPLE_COLORS[c]); body.set_alpha(0.6); body.set_edgecolor("k")
        parts["cmedians"].set_color("k")
        ax.set_xticks(range(1, len(cols) + 1), cols, rotation=45, ha="right", fontsize=8)
        ax.set_title(k, fontsize=9)
        ax.axhline(0, c="grey", lw=0.5, ls=":")
    axes[0].set_ylabel("Mean log2(U-enrichment) per gene")
    fig.tight_layout()
    return fig


def var_groups_plot(v: pd.DataFrame):
    order = ["Group A", "Group A/B", "Group B", "Group B/C", "Group C", "var2csa", "unassigned"]
    pal = dict(zip(order, ["#d62728", "#ff7f0e", "#2ca02c", "#17becf", "#1f77b4", "#e377c2", "#7f7f7f"]))
    v = v.assign(o=v.group.map({g: i for i, g in enumerate(order)})).sort_values(["o", "mean"], ascending=[True, False])
    fig, ax = plt.subplots(figsize=(12, 4))
    x = np.arange(len(v))
    ax.bar(x, v["mean"], yerr=v["sd"], color=v.group.map(pal), edgecolor="k", lw=0.3, capsize=1.5, error_kw={"lw": 0.6})
    ax.set_xticks(x, v.id, rotation=90, fontsize=6)
    for lab, gid in zip(ax.get_xticklabels(), v.id):
        if gid in config.VAR3:
            lab.set_color("#1f3fbf"); lab.set_fontweight("bold")
        if gid in config.VAR2CSA:
            lab.set_color("#c00000"); lab.set_fontweight("bold")
    ax.axhline(0, c="k", lw=0.5)
    ax.set_ylabel("log2(U-enrichment)\nmean ± SD (3D7 rep1, rep2)")
    for g in order:
        if (v.group == g).any():
            ax.bar([0], [0], color=pal[g], label=g)
    ax.legend(ncol=7, fontsize=7, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_xlim(-0.7, len(v) - 0.3)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------- main ---
def main():
    OUT.mkdir(exist_ok=True)
    rep = Report("Fig. 2C, 3A, S1D, S2 / Table S1 - gene-level uracil enrichment",
                 "Mean log2(U-enrichment) per protein-coding gene, VSA family comparisons and var gene groups.")
    genes = genome.protein_coding_genes()
    t = gene_table(genes)
    num = t.select_dtypes("number").columns
    t[num.drop(["start", "end"])] = t[num.drop(["start", "end"])].round(6)
    t.to_csv(OUT / "table_S1_gene_log2U.csv", index=False)

    # ---- compare with legacy tables ------------------------------------------------
    leg = pd.read_csv(LEG / "list_of_average_uenrichment.csv", index_col=0).set_index("id")
    ti = t.set_index("id").reindex(leg.index)
    for s in config.SAMPLE_ORDER:
        d = float(np.nanmax(np.abs(ti[s] - leg[s])))
        same_nan = (ti[s].isna() == leg[s].isna()).all()
        rep.check(f"gene means {s} (legacy mode)", f"{leg[s].notna().sum()} genes", f"{ti[s].notna().sum()} genes",
                  status(d) if same_nan else "DIFFERS", f"max |diff| = {d:.1e}")
        de = float(np.nanmax(np.abs(ti[s] - ti[f"{s}_exact"])))
        rep.check(f"gene means {s}: exact vs legacy mode", "", f"max |diff| = {de:.3f}", "INFO",
                  f"median |diff| = {np.nanmedian(np.abs(ti[s] - ti[f'{s}_exact'])):.4f}")
    xl = pd.read_excel(LEG / "average_uenrichment" / "3d7_merged.log2uenrichment.xlsx").set_index("ID")
    tm = t.set_index("id").reindex(xl.index)
    d = float(np.nanmax(np.abs(tm["3d7_merged_mean0"] - xl["Mean w non-covered"])))
    rep.check("merged 3D7 gene means (xlsx)", "5318 genes", f"{tm['3d7_merged_mean0'].notna().sum()} genes",
              status(d, 1e-4, 0.02), f"max |diff| = {d:.1e}; {(np.abs(tm['3d7_merged_mean0'] - xl['Mean w non-covered']) > 1e-4).sum()} "
              "genes > 1e-4, caused by ~176 near-empty bins that the legacy merged bedGraph leaves blank")

    # ---- Table S1 ranking & top 200 for GO -------------------------------------------
    ranked = t.sort_values("3d7_merged_mean0", ascending=False)
    top = ranked.head(200)
    top.id.to_csv(OUT / "top200_for_GO.txt", index=False, header=False)
    leg_top = set(xl.sort_values("Mean w non-covered", ascending=False).head(200).index)
    ov = len(leg_top & set(top.id))
    rep.check("top-200 genes for ShinyGO", "legacy xlsx top 200", f"{ov}/200 shared", "MATCH" if ov == 200 else "CLOSE")
    comp = top.category.value_counts()
    rep.note("Top-200 composition (annotation-based): " + ", ".join(f"{k} {v}" for k, v in comp.items()) +
             ". GO enrichment itself (Fig. S2) is run on the ShinyGO web server and is not scripted here; "
             "submit output/top200_for_GO.txt with all 5318 genes as background.")

    # ---- Fig. 2C ---------------------------------------------------------------------
    cols = config.SAMPLE_ORDER
    st = vsa_vs_others(t, cols, "category_legacy")
    st.to_csv(OUT / "stats_fig2C_VSA_vs_others_legacy.csv", index=False)
    st_ex = vsa_vs_others(t, [f"{c}_exact" for c in cols], "category")
    st_ex.to_csv(OUT / "stats_fig2C_VSA_vs_others_exact_annotation.csv", index=False)
    savefig(violin_by_dataset(t, cols, "category_legacy", st, "Legacy definition (gene symbol), Welch t-test, BH"),
            OUT, "fig2C_violin")
    savefig(violin_by_dataset(t, [f"{c}_exact" for c in cols], "category", st_ex,
                              "Exact means, annotation-based families (var2csa counted as VAR)"), OUT, "fig2C_violin_exact")
    rep.check("Fig 2C: every VSA family > Others in every dataset (BH p < 0.001)", "yes (**)",
              "yes" if (st.p_BH < 0.001).all() and (st.mean_diff > 0).all() else "NO", "MATCH" if (st.p_BH < 0.001).all() else "DIFFERS",
              f"largest BH p = {st.p_BH.max():.1e} ({st.loc[st.p_BH.idxmax(), 'family']}, {st.loc[st.p_BH.idxmax(), 'dataset']})")
    stev = st[st.family == "STEVOR"].set_index("dataset").p
    rep.check("legacy EPS (3 panels): STEVOR vs Others p, unadjusted",
              ", ".join(f"{v['STEVOR']:.1e}" for v in LEGACY_EPS_P.values()),
              ", ".join(f"{k} {v:.1e}" for k, v in stev.items()), "UNRESOLVED",
              "no dataset/coordinate/resolution variant reproduces these; the EPS was likely made from an earlier "
              "data version. RIF and VAR p < 2.2e-16 in both. Conclusion unchanged.")

    va = var_across_datasets(t, cols, "category_legacy")
    va.to_csv(OUT / "stats_fig2C_var_across_datasets.csv", index=False)
    rd = va[(va.a == "DMSO") & (va.b == "RAPA")].iloc[0]
    rep.check("Fig 2C: var genes RAPA vs DMSO (Welch, BH over 6 pairs)", "p < 0.001",
              f"p_BH = {rd.welch_p_BH:.2g} (mean {rd.mean_b:.3f} vs {rd.mean_a:.3f})",
              "MATCH" if rd.welch_p_BH < 0.001 else "DIFFERS")
    rep.check("var genes RAPA vs DMSO, paired t-test (recommended)", "-",
              f"p_BH = {rd.paired_p_BH:.2g}; Wilcoxon p_BH = {rd.wilcoxon_signed_rank_p_BH:.2g}", "INFO",
              "same 61 genes in both samples -> paired test is the appropriate design")

    # ---- Fig. S1D --------------------------------------------------------------------
    savefig(violin_by_family(t, cols, "category_legacy"), OUT, "figS1D_violin_by_family")

    # ---- Fig. 3A ---------------------------------------------------------------------
    var = io.read_bed(config.VAR_BED)
    grp = io.read_bed(config.VAR_GROUPS_BED, names=["chrom", "start", "end", "name", "group"])
    gmap = dict(zip(grp.name, grp.group.fillna("").str.replace("Group", "Group ")))
    v = t.set_index("id").loc[var.name, ["3d7_rep1", "3d7_rep2"]].copy()
    v["mean"] = v.mean(axis=1)
    v["sd"] = v[["3d7_rep1", "3d7_rep2"]].std(axis=1, ddof=1)
    v["group"] = [("var2csa" if i in config.VAR2CSA else gmap.get(i) or "unassigned") for i in v.index]
    v["subclass"] = ["var3" if i in config.VAR3 else "var2csa" if i in config.VAR2CSA else "" for i in v.index]
    v = v.reset_index().rename(columns={"name": "id"})
    v.round(6).to_csv(OUT / "fig3A_var_gene_means.csv", index=False)
    savefig(var_groups_plot(v), OUT, "fig3A_var_groups")
    gsum = v.groupby("group")["mean"].agg(["count", "mean", "median"]).round(3)
    gsum.to_csv(OUT / "fig3A_group_summary.csv")
    rep.note("Fig. 3A group means: " + "; ".join(f"{g}: n={r['count']}, mean={r['mean']:.3f}" for g, r in gsum.iterrows()))
    rep.note("var3: " + ", ".join(f"{r.id} {r['mean']:.3f}±{r.sd:.3f}" for _, r in v[v.subclass == "var3"].iterrows()) +
             "; var2csa: " + ", ".join(f"{r['mean']:.3f}±{r.sd:.3f}" for _, r in v[v.subclass == "var2csa"].iterrows()) +
             f"; all other var genes: median {v[v.subclass == ''].loc[:, 'mean'].median():.3f}.")
    rep.note("Var group annotation issues (data/BED_files/var_genes/): vargene_groups.txt maps MAL7P1.1 to "
             "PF3D7_0425800 (duplicate of PFD1235w); PF3D7_0700100 is therefore ungrouped. PF3D7_0600600, "
             "PF3D7_0833500 and PF3D7_1240900 have no group. Shown as 'unassigned' - please check.")
    rep.note("Legacy Fig. 2C family definition used the gene symbol only, so var2csa (symbol VAR2CSA) and three "
             "PfEMP1 exon-2 fragments were counted as 'Others'. The *_exact_annotation outputs count them as VAR; "
             "conclusions do not change.")
    rep.note("Legacy R code ran pairwise.t.test with the default pool.sd = TRUE (pooled-SD t-tests), while the "
             "manuscript states Welch's t-test. Both are in stats_fig2C_VSA_vs_others_*.csv (column p_pooled_sd_unadjusted); all are p < 1e-9.")
    rep.write(HERE / "COMPARISON.md")


if __name__ == "__main__":
    main()
