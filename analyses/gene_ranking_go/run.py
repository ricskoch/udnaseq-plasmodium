#!/usr/bin/env python3
"""Gene ranking, GO enrichment, feature patterns and differential patterns - raw vs GC-corrected.

Purpose
    Repeat the original logic on the validated tracks: rank genes by mean log2U, take the top n, run GO enrichment, look at
    where within genes uracil sits, and compare conditions (RAPA vs DMSO, rep1 vs rep2). Everything is done on the raw
    100 bp log2U and on two GC-controlled versions, to see what survives.
Scores (per dataset, nuclear protein-coding genes; gene value = mean over 100 bp bins whose centre lies in the gene, or the bin
containing the gene midpoint if none)
    raw        : log2((IP+1)/(INPUT+1)) from the RPGC coverage (analyses/00_tracks)
    gc_bin     : same after dividing IP and INPUT by the mean coverage of their 100 bp GC class (analyses/gc_corrected)
    gc_gene    : raw gene mean minus the value expected from the gene's own GC (quadratic fit over all nuclear genes)
Analyses
    1 top-n composition (n = 50, 100, 200, 400; 200 is the manuscript value): % var/rif/stevor, mean gene GC, overlap (Jaccard)
    2 GO enrichment of the top n (Fisher exact, BH; GO terms propagated through is_a; background = annotated nuclear genes;
      terms with 5-1000 background genes and >= 3 hits). Annotations: UniProt UP000001450 (scripts/fetch_go.py).
    3 feature patterns: first exon, last exon, introns, 1 kb up- and downstream, for var, rif/stevor, other multi-exon genes
    4 differential: RAPA vs DMSO and rep1 vs rep2 - paired gene differences, Spearman between datasets, top-200 differential
      genes (z-scored scores), their composition, GC and GO terms
Outputs  output/gene_scores.csv.gz, topn_composition.csv, go_enrichment.csv, feature_patterns.csv, differential_*.csv, COMPARISON.md
"""
import importlib.util
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from udnaseq import config, genome, io  # noqa: E402


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[0] / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


gcc = load("gcc", "gc_corrected/run.py")
f4b = load("f4b", "fig4B_gc_correlation/run.py")
gm = load("gm", "fig4C_gc_matched_null/run.py")
OUT = HERE / "output"
SETS = config.SAMPLE_ORDER + ["3d7_merged"]
SCORES = ["raw", "gc_bin", "gc_gene"]
NS = [50, 100, 200, 400]
GO_DIR = config.DATA / "go"


# ---------------------------------------------------------------- helpers ---
def bin_prefix(T, y):
    """per-chromosome (starts, prefix sums) for fast interval means over 100 bp bins."""
    out = {}
    for c in np.unique(T["chrom"]):
        m = T["chrom"] == c
        out[str(c)] = (T["start"][m] + config.BIN // 2, np.r_[0, np.cumsum(y[m])])
    return out


def interval_mean(pref, chrom, s, e):
    if chrom not in pref or e <= s:
        return np.nan
    ctr, cs = pref[chrom]
    i, j = np.searchsorted(ctr, s), np.searchsorted(ctr, e)
    if j > i:
        return (cs[j] - cs[i]) / (j - i)
    k = min(max(np.searchsorted(ctr, (s + e) / 2) - 1, 0), len(ctr) - 1)
    return cs[k + 1] - cs[k]


def read_go():
    parents, names, ns, cur = defaultdict(set), {}, {}, None
    for line in open(GO_DIR / "go-basic.obo"):
        line = line.strip()
        if line == "[Term]":
            cur = {}
        elif line.startswith("[") or not line:
            if line.startswith("["):
                cur = None
        elif cur is not None:
            if line.startswith("id: "):
                cur["id"] = line[4:]
            elif line.startswith("name: "):
                names[cur["id"]], ns[cur["id"]] = line[6:], None
            elif line.startswith("namespace: "):
                ns[cur["id"]] = line[11:]
            elif line.startswith("is_a: "):
                parents[cur["id"]].add(line[6:].split(" ")[0])
    anc = {}

    def up(t):
        if t not in anc:
            anc[t] = {t}.union(*[up(p) for p in parents.get(t, ())]) if parents.get(t) else {t}
        return anc[t]

    ann = pd.read_csv(GO_DIR / "pf3d7_uniprot_go.tsv", sep="\t").dropna()
    gene_go = defaultdict(set)
    for orf, go in zip(ann["Gene Names (ORF)"], ann["Gene Ontology IDs"]):
        terms = {t.strip() for t in go.split(";") if t.strip() in names}
        full = set().union(*[up(t) for t in terms]) if terms else set()
        for g in re.split(r"[\s;]+", orf):
            if g.startswith("PF3D7_"):
                gene_go[g] |= full
    return gene_go, names, ns


def bh(p):
    p = np.asarray(p)
    o = np.argsort(p)
    q = np.empty(len(p))
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def go_enrich(top_ids, bg_ids, gene_go, names, ns, min_bg=5, max_bg=1000, min_hit=3):
    bg = [g for g in bg_ids if g in gene_go]
    top = [g for g in top_ids if g in gene_go]
    if len(top) < 5:
        return pd.DataFrame()
    cnt_bg, cnt_top = defaultdict(int), defaultdict(int)
    for g in bg:
        for t in gene_go[g]:
            cnt_bg[t] += 1
    for g in top:
        for t in gene_go[g]:
            cnt_top[t] += 1
    rows = [(t, cnt_top[t], cnt_bg[t]) for t in cnt_top if cnt_top[t] >= min_hit and min_bg <= cnt_bg[t] <= max_bg
            and ns.get(t) in ("biological_process", "molecular_function", "cellular_component")]
    if not rows:
        return pd.DataFrame()
    d = pd.DataFrame(rows, columns=["go", "hits", "bg"])
    d["p"] = stats.hypergeom.sf(d.hits - 1, len(bg), d.bg, len(top))
    d["fdr"] = bh(d.p.values)
    d["fold"] = (d.hits / len(top)) / (d.bg / len(bg))
    d["term"], d["ns"] = d.go.map(names), d.go.map(ns)
    d["n_top_annotated"], d["n_bg"] = len(top), len(bg)
    return d.sort_values("p")


def gff_structure(genes):
    g = pd.read_csv(config.GFF, sep="\t", comment="#", header=None)
    ex = g[g[2] == "exon"].copy()
    ex["tid"] = ex[8].str.extract(r"Parent=([^;]+)")[0]
    ex["gid"] = ex[8].str.extract(r"gene_id=([^;]+)")[0]
    # keep the first transcript (.1) of each gene
    ex = ex[ex.tid.str.endswith(".1") | ~ex.tid.str.contains(r"\.\d+$")]
    out = {}
    for gid, d in ex.groupby("gid"):
        d = d.sort_values(3)
        out[gid] = list(zip(d[3] - 1, d[4]))
    return out


# ------------------------------------------------------------------- main ---
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
    genes = genes[genes.chrom.isin(config.NUCLEAR_CHROMS)].reset_index(drop=True)
    cum = gm.gc_arrays(seqs)
    genes["gc"] = [gm.window_gc(cum[r.chrom], np.array([r.start]), np.array([r.end]))[0] for r in genes.itertuples()]
    genes["length"] = genes.end - genes.start
    genes["is_vsa"] = genes.category != "Others"

    # ---- scores
    tracks, prefs = {}, {}
    for s in SETS:
        ip, inp = (T["ip_3d7_merged"], T["input_3d7_merged"]) if s == "3d7_merged" else (T[f"ip_{s}"], T[f"input_{s}"])
        ipc, inc = gcc.gc_normalise(ip, np.nan_to_num(gcf), ref), gcc.gc_normalise(inp, np.nan_to_num(gcf), ref)
        tracks[s] = {"raw": T[f"log2_{s}"], "gc_bin": np.log2((ipc + 1) / (inc + 1))}
        prefs[s] = {k: bin_prefix(T, v) for k, v in tracks[s].items()}
    S = genes[["id", "symbol", "category", "chrom", "start", "end", "gc", "length", "is_vsa"]].copy()
    for s in SETS:
        for k in ("raw", "gc_bin"):
            S[f"{s}|{k}"] = [interval_mean(prefs[s][k], r.chrom, r.start, r.end) for r in genes.itertuples()]
        co = np.polyfit(S.gc, S[f"{s}|raw"].fillna(S[f"{s}|raw"].mean()), 2)
        S[f"{s}|gc_gene"] = S[f"{s}|raw"] - np.polyval(co, S.gc)
    S.round(5).to_csv(OUT / "gene_scores.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})

    md = []
    # ---- 1 top-n composition
    comp = []
    for s in SETS:
        for k in SCORES:
            v = S[f"{s}|{k}"]
            for n in NS:
                top = S.loc[v.nlargest(n).index]
                comp.append(dict(dataset=s, score=k, n=n, pct_vsa=100 * top.is_vsa.mean(), pct_var=100 * (top.category == "VAR").mean(),
                                 mean_gc=top.gc.mean(), median_len=top.length.median()))
    C = pd.DataFrame(comp)
    C.round(3).to_csv(OUT / "topn_composition.csv", index=False)
    jac = []
    for k in SCORES:
        for n in (200,):
            tops = {s: set(S.id[S[f"{s}|{k}"].nlargest(n).index]) for s in SETS}
            for a in config.SAMPLE_ORDER:
                for b in config.SAMPLE_ORDER:
                    if a < b:
                        jac.append(dict(score=k, n=n, a=a, b=b, jaccard=len(tops[a] & tops[b]) / len(tops[a] | tops[b])))
    J = pd.DataFrame(jac)
    J.round(3).to_csv(OUT / "topn_overlap.csv", index=False)
    md += ["## 1. Top-n genes: composition", "", "Percent var/rif/stevor among the top n (all genes: %.1f %% VSA, mean gene GC %.1f %%):" % (100 * S.is_vsa.mean(), 100 * S.gc.mean()), "",
           C[C.n.isin([100, 200, 400])].pivot_table(index=["dataset", "score"], columns="n", values="pct_vsa").round(1).to_markdown(), "",
           "Mean gene GC (%) of the top 200:", "", C[C.n == 200].assign(gc=lambda d: (100 * d.mean_gc).round(1)).pivot(index="dataset", columns="score", values="gc").to_markdown(), "",
           "Overlap of the top-200 gene lists between datasets (Jaccard):", "", J.pivot_table(index=["a", "b"], columns="score", values="jaccard").round(2).to_markdown(), ""]

    # ---- 2 GO
    gene_go, names, ns = read_go()
    bg = S.id[S.id.isin(gene_go)].tolist()
    rows, summ = [], []
    for s in SETS:
        for k in SCORES:
            for n in (100, 200, 400):
                for excl in (False, True):
                    pool = S[~S.is_vsa] if excl else S
                    top = pool.loc[pool[f"{s}|{k}"].nlargest(n).index].id.tolist()
                    d = go_enrich(top, bg if not excl else [g for g in bg if g in set(pool.id)], gene_go, names, ns)
                    if len(d):
                        d = d.assign(dataset=s, score=k, n=n, vsa_excluded=excl)
                        rows.append(d[d.fdr < 0.05])
                    sig = d[d.fdr < 0.05] if len(d) else d
                    summ.append(dict(dataset=s, score=k, n=n, vsa_excluded=excl, n_annotated=int(d.n_top_annotated.iloc[0]) if len(d) else 0,
                                     n_sig_terms=len(sig), top_terms="; ".join(f"{t} ({f:.1e})" for t, f in zip(sig.term.head(4), sig.fdr.head(4)))))
    G = pd.concat(rows) if rows else pd.DataFrame()
    G.round(6).to_csv(OUT / "go_enrichment.csv", index=False)
    GS = pd.DataFrame(summ)
    GS.to_csv(OUT / "go_summary.csv", index=False)
    md += ["## 2. GO enrichment of the top n (BH FDR < 0.05; terms with 5-1000 background genes)", "",
           "Number of significant terms (VSA genes included / VSA genes removed before ranking):", ""]
    pv = GS.assign(v=lambda d: d.n_sig_terms.astype(str)).pivot_table(index=["dataset", "score", "n"], columns="vsa_excluded", values="n_sig_terms", aggfunc="first")
    pv.columns = ["incl. VSA", "VSA removed"]
    md += [pv.to_markdown(), "", "Top terms, n = 200, VSA genes included:", ""]
    t200 = GS[(GS.n == 200) & (~GS.vsa_excluded)][["dataset", "score", "n_annotated", "top_terms"]]
    md += [t200.to_markdown(index=False), "", "Top terms, n = 200, VSA genes removed before ranking:", ""]
    t200x = GS[(GS.n == 200) & (GS.vsa_excluded)][["dataset", "score", "n_annotated", "top_terms"]]
    md += [t200x.to_markdown(index=False), ""]

    # ---- 3 features
    struct = gff_structure(genes)
    frows = []
    for s in SETS:
        for k in ("raw", "gc_bin"):
            pref = prefs[s][k]
            for r in genes.itertuples():
                ex = struct.get(r.id, [])
                if len(ex) < 2:
                    continue
                cls = "VAR" if r.category == "VAR" else ("RIF/STEVOR" if r.category in ("RIF", "STEVOR") else "other multi-exon")
                feats = {"exon_first": ex[:1], "exon_last": ex[-1:], "introns": [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]}
                if r.strand == "-":
                    feats = {"exon_first": ex[-1:], "exon_last": ex[:1], "introns": feats["introns"]}
                    up, dn = (r.end, r.end + 1000), (r.start - 1000, r.start)
                else:
                    up, dn = (r.start - 1000, r.start), (r.end, r.end + 1000)
                feats["upstream_1kb"], feats["downstream_1kb"] = [up], [dn]
                for fn, iv in feats.items():
                    iv = [(a, b) for a, b in iv if b > a and a >= 0]
                    if not iv:
                        continue
                    v = np.nanmean([interval_mean(pref, r.chrom, a, b) for a, b in iv])
                    frows.append((s, k, cls, r.id, fn, v))
    F = pd.DataFrame(frows, columns=["dataset", "score", "class", "id", "feature", "value"])
    FP = F.groupby(["dataset", "score", "class", "feature"]).value.agg(["mean", "count"]).reset_index()
    FP.round(4).to_csv(OUT / "feature_patterns.csv", index=False)
    # feature GC
    ggc = []
    for r in genes.itertuples():
        ex = struct.get(r.id, [])
        if len(ex) < 2:
            continue
        cls = "VAR" if r.category == "VAR" else ("RIF/STEVOR" if r.category in ("RIF", "STEVOR") else "other multi-exon")
        feats = {"exon_first": ex[:1] if r.strand != "-" else ex[-1:], "exon_last": ex[-1:] if r.strand != "-" else ex[:1],
                 "introns": [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]}
        for fn, iv in feats.items():
            v = [gm.window_gc(cum[r.chrom], np.array([a]), np.array([b]))[0] for a, b in iv if b > a]
            if v:
                ggc.append((cls, fn, np.mean(v)))
    GCF = pd.DataFrame(ggc, columns=["class", "feature", "gc"]).groupby(["class", "feature"]).gc.mean().reset_index()
    GCF.round(4).to_csv(OUT / "feature_gc.csv", index=False)
    md += ["## 3. Where within genes: mean log2U by gene feature (raw / GC-corrected), multi-exon genes", "",
           "Feature GC (%): " + "; ".join(f"{r['class']} {r.feature} {100 * r.gc:.0f}" for _, r in GCF.iterrows()), ""]
    for s in SETS:
        sub = FP[FP.dataset == s].pivot_table(index=["class", "feature"], columns="score", values="mean").round(3)
        n_ = FP[(FP.dataset == s) & (FP.score == "raw")].set_index(["class", "feature"])["count"]
        sub["n_genes"] = n_
        md += [f"### {s}", "", sub.to_markdown(), ""]

    # ---- 4 differential
    drows, dgo = [], []
    for a, b in (("RAPA", "DMSO"), ("3d7_rep1", "3d7_rep2")):
        for k in SCORES:
            va, vb = S[f"{a}|{k}"], S[f"{b}|{k}"]
            ok = va.notna() & vb.notna()
            rho = stats.spearmanr(va[ok], vb[ok]).statistic
            za, zb = (va - va.mean()) / va.std(), (vb - vb.mean()) / vb.std()
            delta = za - zb
            top_up, top_dn = S.loc[delta.nlargest(200).index], S.loc[delta.nsmallest(200).index]
            row = dict(a=a, b=b, score=k, spearman_between=rho, delta_vs_gc_rho=stats.spearmanr(delta[ok], S.gc[ok]).statistic,
                       up_pct_vsa=100 * top_up.is_vsa.mean(), down_pct_vsa=100 * top_dn.is_vsa.mean(), up_gc=top_up.gc.mean(), down_gc=top_dn.gc.mean())
            for fam in ("VAR", "RIF", "STEVOR"):
                m = ((S.category == fam) & ok).values
                d = (va - vb)[m]
                row[f"{fam}_mean_diff"] = d.mean()
                row[f"{fam}_paired_p"] = stats.ttest_1samp(d, 0).pvalue
            drows.append(row)
            for lab, top in (("up_in_" + a, top_up), ("up_in_" + b, top_dn)):
                d = go_enrich(top.id.tolist(), bg, gene_go, names, ns)
                sig = d[d.fdr < 0.05] if len(d) else d
                dgo.append(dict(a=a, b=b, score=k, direction=lab, n_sig=len(sig), top_terms="; ".join(f"{t} ({f:.1e})" for t, f in zip(sig.term.head(4), sig.fdr.head(4)))))
    D = pd.DataFrame(drows)
    D.round(4).to_csv(OUT / "differential_summary.csv", index=False)
    DG = pd.DataFrame(dgo)
    DG.to_csv(OUT / "differential_go.csv", index=False)
    md += ["## 4. Differential patterns", "", "z-scored gene scores; top 200 genes with the largest z-difference.", "",
           D[["a", "b", "score", "spearman_between", "delta_vs_gc_rho", "up_pct_vsa", "down_pct_vsa", "up_gc", "down_gc"]].round(3).to_markdown(index=False), "",
           "Paired gene differences (a - b, log2U units) for the VSA families:", "",
           D[["a", "b", "score"] + [f"{f}_{x}" for f in ("VAR", "RIF", "STEVOR") for x in ("mean_diff", "paired_p")]].round(4).to_markdown(index=False), "",
           "GO terms of the top-200 differential genes:", "", DG.to_markdown(index=False), ""]

    head = ["# Gene ranking, GO, feature and differential patterns: raw vs GC-corrected", "",
            "Raw = log2((IP+1)/(INPUT+1)) from RPGC 100 bp tracks; gc_bin = after GC normalisation of IP and INPUT per 100 bp GC class; "
            "gc_gene = raw gene mean minus expectation from the gene's own GC. See run.py docstring for the full method.", "",
            f"_Generated by `run.py`; GO annotations: UniProt UP000001450 ({len(gene_go)} genes with GO terms), ontology {next(l.strip() for l in open(GO_DIR / 'go-basic.obo') if l.startswith('data-version'))}._", ""]
    (HERE / "COMPARISON.md").write_text("\n".join(head + md) + "\n")
    print("\n".join(md)[:20000])


if __name__ == "__main__":
    main()
