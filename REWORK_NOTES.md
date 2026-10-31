# Rework notes: reproducing the U-DNA-Seq analyses in Python

Scope: regenerate every computational figure panel from the raw-most data, compare with the legacy files and the
manuscript, and restructure the repository with one folder per figure. Per-check details are in each
`analyses/*/COMPARISON.md`.

## 1. Bottom line

Every quantitative statement in the manuscript that is based on the U-DNA-Seq data is reproduced.

| Claim | Manuscript | Regenerated |
|---|---|---|
| Moderate replicate correlation; DMSO–RAPA highest (S1A/B) | qualitative | Spearman rep1–rep2 0.46, DMSO–RAPA 0.77 (identical to legacy) |
| VSA families > other genes, every dataset (2C) | p < 0.001 | all 12 tests BH p ≤ 1e-11 |
| var genes RAPA > DMSO (2C) | p < 0.001 | Welch BH p = 3e-9; paired t-test p = 1e-18 |
| var enrichment in first ~5–6 kb / exon 1 (3B) | qualitative | mean log2U 0.25 (0–6 kb) vs 0.07 (6–10 kb) |
| rif/stevor enrichment in exon 2 (3B) | qualitative | body peak at 53 % of scaled length; flanks negative |
| var3 and var2csa low (3A) | qualitative | var3 −0.04 to 0.03, var2csa 0.07, other var median 0.21 |
| GC% vs log2U (4B) | r = 0.42 | Pearson 0.436, Spearman 0.415 |
| Hotspot–G4 distance (4C) | median 0, IQR (−206.5, 389.5), G4 size 33.13/29 | identical |

Legacy intermediate files reproduced exactly (float32 precision): all per-sample log2 tracks, the merged 3D7 track,
the per-gene table, the merged gene table/top-200 GO list, all deepTools correlation matrices, the rif/stevor
deepTools matrix, the PCA scores, the GC table and the 100 bp G4 table.

## 2. Points to fix or state in the manuscript

1. **The "merged" 3D7 dataset is mostly replicate 1.** It was built from merged BAMs, so each replicate counts in proportion
   to its reads: IP 89 % rep1, INPUT 57 % rep1 (estimated; fits the legacy track to r = 0.99999999).
   Figures 2B, 3B, 4A–C ("merged") therefore mainly show replicate 1. Consider saying so, or showing rep2 alongside.
2. **Fig. 4B coefficient.** "r = 0.42, Pearson" — Pearson is 0.44, Spearman is 0.42 (0.415).
3. **Fig. 3B legend.** rif/stevor genes were scaled to **1,500 bp**, not 2,000 bp. Both versions are in
   `fig3B_vsa_metaprofiles/output/`. The legacy var heatmap SVG is also titled "Rif/stevor gene U-enrichment".
4. **Fig. 4C filtering order.** The closest G4 was searched among *all* predictions, and pairs were filtered afterwards
   (|score| > 1.2, ≤ 200 bp). This silently drops 2,201 of 5,011 hotspot bins. Filtering G4s first keeps every
   hotspot: median still 0, IQR (−428, 526).
5. **Fig. 4C has no background.** Added a permutation null (random 10 bp bins, 200 permutations): 25.5 % of hotspots
   lie within 100 bp of a strong G4 vs 1.1 % of random bins (24×, p < 0.005); overlap 6.0 % vs 0.2 % (37×).
   Random bins are not GC-matched, so part of this excess follows from the GC correlation.
6. **Statistics wording.** The legacy R code used `pairwise.t.test` with pooled SD; the manuscript says Welch.
   The new code uses Welch + BH (both are reported; every p < 1e-9). RAPA vs DMSO compares the same 61 genes,
   so a paired test is the better design (p even smaller).
7. **Blacklist.** Methods say telomeres were blacklisted "from all downstream analysis". The blacklist does not cover the
   chr05 3′ end, the chr06 5′ end or either end of chr14. The 10 bp analyses (S1C, 3B, 4A, 4C) and the ideogram used no blacklist.
8. **Gene family definition (2C).** Families were assigned from the gene symbol, so var2csa (symbol `VAR2CSA`) and
   three PfEMP1 exon-2 fragments were counted as "Others". Counting them as VAR changes nothing.

## 3. Data / annotation issues found

- `data/BED_files/var_genes/vargene_groups.txt`: `MAL7P1.1` maps to PF3D7_0425800 (duplicate of PFD1235w). As a result
  PF3D7_0700100 is ungrouped. PF3D7_0600600, PF3D7_0833500 and PF3D7_1240900 have no group either.
- var3 IDs used for highlighting (`src/udnaseq/config.py`): PF3D7_0100300, PF3D7_0600400, PF3D7_0937600
  (the three ~4–5 kb Group A var genes). **Please confirm.**
- The legacy GC table stops at chr13:2,249,900 (truncated output). Correlations are unaffected.
- The legacy correlation command lists files in the order RAPA, DMSO but labels them DMSO, RAPA. The stored matrices are
  nevertheless labelled correctly (checked).

## 4. Method details that the legacy outputs depended on

- **deepTools averages are approximate.** `multiBigwigSummary`/`computeMatrix` read bigWig zoom-level summaries. Compared
  with exact per-base means, per-region values differ by up to 0.14 (genes) to 0.31 (10 kb bins) log2 units.
  Correlations change by ≤ 0.12 and no conclusion changes. Every affected analysis outputs both: `deeptools` (reproduces
  the figures) and `exact`.
- **1 bp coordinate offset.** Gene means used the GFF start (1-based) as a 0-based start. The effect is ≤ 0.002 and is
  reproduced on purpose; `*_exact` columns use correct coordinates.
- **Legacy violin EPS p-values** (STEVOR 3.5e-14, 1.3e-15, 6.5e-13) are not reproduced by any dataset or variant.
  They were probably made from an earlier data version. Current values are 3e-14 to 1e-11; the conclusion is unchanged.

## 5. Not reproduced (no data or code in the repository)

Fig. 1 wet-lab panels, the ShinyGO GO analysis itself (input list provided), read mapping, and the exact
coordinates of the Fig. 4A locus (default: chr7 central var cluster, adjustable).
