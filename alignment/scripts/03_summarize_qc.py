"""Collect per-library QC, compare rebuilt RPGC bigWigs with the repo/GEO tracks, and compute
IP-vs-INPUT fingerprint metrics from the 100 bp bigWigs (no BAMs needed).

Usage: python3 -I 03_summarize_qc.py <WORK/bigwig dir> <repo root>
Writes alignment/qc/qc_summary.tsv, track_comparison.tsv, fingerprint.tsv, fingerprint.png
"""
import json, os, re, sys, csv
import numpy as np, pyBigWig

BW, REPO = sys.argv[1], sys.argv[2]
QC = f"{REPO}/alignment/qc"
PUB = {  # library -> published (repo/GEO) RPGC bigWig
 "3d7_rep1_input": "3d7_rep1.Plasmodium.INPUTcorrect", "3d7_rep1_udnaip": "3d7_rep1.Plasmodium.UDNAIPcorrect",
 "3d7_rep2_input": "3d7_rep2.Plasmodium.INPUT", "3d7_rep2_udnaip": "3d7_rep2.Plasmodium.UDNAIP",
 "DMSO_input": "DMSO.INPUT", "DMSO_udnaip": "DMSO.UDNAIP", "RAPA_input": "RAPA.INPUT", "RAPA_udnaip": "RAPA.UDNAIP"}
libs = [l.split("\t")[0] for l in open(f"{REPO}/alignment/samples.tsv").read().splitlines()[1:]]

def bins(p):
    b = pyBigWig.open(p); out = []
    for c, l in b.chroms().items():
        out.append(np.nan_to_num(np.array(b.stats(c, 0, l, type="mean", nBins=(l + 99) // 100), dtype=float)))
    return np.concatenate(out)

def metrics(path, key):
    t = open(path).read().splitlines(); i = next(k for k, x in enumerate(t) if x.startswith(key))
    return dict(zip(t[i].split("\t"), t[i + 1].split("\t")))

rows, cmp_rows = [], []
for lib in libs:
    q = f"{QC}/{lib}"
    if not os.path.exists(f"{q}/markdup_metrics.txt"): continue
    fp = json.load(open(f"{q}/fastp.json"))
    fs = open(f"{q}/flagstat_aligned.txt").read()
    prim = int(re.search(r"(\d+) \+ 0 primary\n", fs).group(1))
    pm = int(re.search(r"(\d+) \+ 0 primary mapped", fs).group(1))
    md = metrics(f"{q}/markdup_metrics.txt", "LIBRARY")
    ins = metrics(f"{q}/insert_size_metrics.txt", "MEDIAN_INSERT_SIZE")
    gc = metrics(f"{q}/gc_bias_summary.txt", "ACCUMULATION_LEVEL")
    rows.append(dict(library=lib, spots_sra=fp["summary"]["before_filtering"]["total_reads"] // 2,
        pairs_after_trim=fp["summary"]["after_filtering"]["total_reads"] // 2, read_gc_pct=round(100*fp["summary"]["before_filtering"]["gc_content"],1),
        primary_reads=prim, mapped_pct=round(100 * pm / prim, 2),
        dup_pct=round(100 * float(md["PERCENT_DUPLICATION"]), 2), est_library_size=md["ESTIMATED_LIBRARY_SIZE"],
        median_insert=ins["MEDIAN_INSERT_SIZE"], mean_insert=round(float(ins["MEAN_INSERT_SIZE"]), 1),
        gc_at_dropout=gc["AT_DROPOUT"], gc_gc_dropout=gc["GC_DROPOUT"],
        gc_nc_0_19=gc["GC_NC_0_19"], gc_nc_20_39=gc["GC_NC_20_39"], gc_nc_40_59=gc["GC_NC_40_59"], gc_nc_60_79=gc["GC_NC_60_79"],
        reads_mapq1=open(f"{q}/n_reads_mapq1.txt").read().strip()))
    p = f"{REPO}/data/{PUB[lib]}.bin100bp.RPGC.bw"
    R = bins(p)
    for mq in (0, 1, 10, 30):
        f = f"{BW}/{lib}.mapq{mq}.bin100bp.RPGC.bw"
        if not os.path.exists(f): continue
        T = bins(f)
        rk = lambda x: np.argsort(np.argsort(x))
        cmp_rows.append(dict(library=lib, mapq=mq, pearson=round(np.corrcoef(R, T)[0, 1], 4),
            spearman=round(np.corrcoef(rk(R), rk(T))[0, 1], 4), mean_rebuilt=round(T.mean(), 3),
            mean_published=round(R.mean(), 3), sum_ratio=round(T.sum() / R.sum(), 3)))

def write(name, rs):
    if not rs: return
    with open(f"{QC}/{name}", "w", newline="") as fh:
        w = csv.DictWriter(fh, rs[0].keys(), delimiter="\t"); w.writeheader(); w.writerows(rs)
write("qc_summary.tsv", rows); write("track_comparison.tsv", cmp_rows)

# fingerprint (deepTools plotFingerprint-style) from the MAPQ>=1 100 bp tracks
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 4, figsize=(16, 3.8)); fpr = []
for ax, (name, i, p) in zip(axs, [("3d7_rep1", "3d7_rep1_input", "3d7_rep1_udnaip"), ("3d7_rep2", "3d7_rep2_input", "3d7_rep2_udnaip"),
                                  ("DMSO", "DMSO_input", "DMSO_udnaip"), ("RAPA", "RAPA_input", "RAPA_udnaip")]):
    fi, fp_ = f"{BW}/{i}.mapq1.bin100bp.RPGC.bw", f"{BW}/{p}.mapq1.bin100bp.RPGC.bw"
    if not (os.path.exists(fi) and os.path.exists(fp_)): continue
    cur = {}
    for lab, f in (("INPUT", fi), ("IP", fp_)):
        v = np.sort(bins(f)); c = np.cumsum(v) / v.sum(); x = np.arange(1, len(v) + 1) / len(v)
        cur[lab] = (x, c); ax.plot(x, c, label=lab)
        cur[lab + "_top1"] = 1 - np.interp(0.99, x, c); cur[lab + "_top10"] = 1 - np.interp(0.90, x, c)
    d = np.abs(cur["INPUT"][1] - cur["IP"][1]); k = d.argmax()
    fpr.append(dict(pair=name, max_gap=round(float(d[k]), 4), at_rank_fraction=round(float(cur["INPUT"][0][k]), 3),
        input_signal_in_top10pct_bins=round(float(cur["INPUT_top10"]), 4), ip_signal_in_top10pct_bins=round(float(cur["IP_top10"]), 4),
        input_signal_in_top1pct_bins=round(float(cur["INPUT_top1"]), 4), ip_signal_in_top1pct_bins=round(float(cur["IP_top1"]), 4)))
    ax.set_title(name); ax.set_xlabel("bins ranked low to high"); ax.set_ylabel("cumulative fraction of signal"); ax.legend()
plt.tight_layout(); plt.savefig(f"{QC}/fingerprint.png", dpi=110); write("fingerprint.tsv", fpr)
print(open(f"{QC}/qc_summary.tsv").read() if rows else "no libs")
