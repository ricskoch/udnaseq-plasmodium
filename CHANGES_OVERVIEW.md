# What changed against the original analysis (one-page overview)

Results only. Manuscript text is not in this repo. Numbers come from the `COMPARISON.md` next to each analysis in `analyses/`.
Unit: log2(U-enrichment) = log2((IP+1)/(INPUT+1)), 100 bp bins unless stated. "Equal GC" = bins compared inside the same GC class (20 classes), genic = bin centre inside a protein-coding gene.

## 0. Summary table (one screen)

| Original figure / claim | Reanalysis in short | Verdict | Docx section changed | Figure status |
|---|---|---|---|---|
| Fig 1 (KO, IFA, Giemsa, dot blot) | not re-analysed (wet lab) | kept | text cites Fig 1G for "more uracil overall" | original, untouched |
| Fig 2A schematic | n/a | kept | none | original, untouched |
| Fig 2B ideogram: uracil clusters in all chromosomes | holds in rep1, rep2, DMSO, RAPA; tracks reproduced from FASTQ | kept | Results 2 (one sentence on rebuild) | original, untouched |
| Fig 2C: VSA genes higher than other genes | raw: yes in all 4 datasets (+0.09 to +0.38). After GC control: method-dependent (see section 2) | softened | Results 2, Fig 4 section, Discussion | original copy, unedited (RAPA/DMSO mark not removable); regenerated version in repo |
| Fig 2C: var higher in RAPA than DMSO | raw +0.17 (var), +0.10 (rif); after GC control 0.00 and -0.04 | moved to supplement | Results 2, new Fig S4 | new Fig S4 |
| Fig 3A: var groups enriched, var3/var2csa lower | not re-tested at equal GC | kept (raw, descriptive) | none | original, untouched |
| Fig 3B: first var exon / last rif-stevor exon high | holds raw in 4 datasets; follows GC (profile r 0.82-0.92); shrinks after GC normalisation | softened (descriptive) | Results 3 (one sentence + caption) | original, untouched |
| Fig 4A locus (chr7) with G4 bars | window recovered (chr7:507,300-604,000); G4 part dropped | dropped from main | Fig 4 section rewritten | replaced by new Fig 4 |
| Fig 4B: uracil vs GC, r = 0.42 | reproduced r = 0.44 (all bins) | kept (as panel B) | Fig 4 section | new Fig 4, panel B |
| Fig 4C: hotspots near G4 (median distance 0) | raw 25x (rep1), 11x (rep2, RAPA), 0.4x (DMSO); GC-matched 1.25-1.50x (3D7), none in RAPA/DMSO | dropped (1-2 sentences + table) | Fig 4 section, Discussion, Methods | Table S2 (new) |
| S1A/B replicate correlation | moderate agreement confirmed; gene-score rank correlation rep1 vs rep2 0.48 raw, 0.28 after GC control | kept + limitation | Results 2, Discussion | original, untouched |
| S2 / Table S1: GO of top 200 genes | top 200 are 41-73 % VSA raw (4.8 % of all genes), 1-8 % after GC control; list overlap between datasets 0.13-0.36 | kept with caveat | Results 2 | original, untouched |
| New: genic > intergenic | at equal GC +0.02 to +0.08, all 4 datasets, CI excludes 0 | new, robust | Fig 4 section | new Fig 4C |
| New: VSA identity vs subtelomeric position | VSA genes positive in both positions; non-VSA subtelomeric ~0 | new | Fig 4 section | new Fig 4C |
| New: IP vs INPUT GC bias | IP AT dropout 1.6-2.0x INPUT; IP gains coverage with GC | new | Fig 4 section | new Fig 4A |

## 1. Files to look at

Mac: `Plasmodium U-DNA-Seq_working/`
- `Plasmodium_UDNA_cell_working.docx`: tracked changes (author "Claude (reanalysis)").
- `Figures_working/NEW_replaces_Fig4_gcbias.svg`: new Figure 4 (panels A-C); replaces old Fig 4 (`figure_gc.svg`, `figure_g4.svg`, `igv_panel_closeup.svg`). The repo script still writes it as `figure_gcbias.svg`.
- `Figures_working/NEW_adds_FigS4_rapa_dmso.svg`: new Figure S4 (no old counterpart). The repo script writes it as `figure_s4_rapa_dmso.svg`.
- `Figures_working/table_s2_g4.csv`: Table S2.
- `Figures_working/figure_udnaseq_v2.svg`: exact copy of the original Fig 2 SVG (unedited).
Repo `analyses/fig_gcbias`, `analyses/fig_supp` regenerate the new figures from the tracks.

## 2. Important caveat: the VSA excess depends on the method

Three ways to control for GC give different answers for "VSA genes vs other genes", so the manuscript wording "three of four datasets" is conditional on the first method.

| Method | What it compares | VSA vs other (examples) |
|---|---|---|
| A. Equal-GC bin contrast (used in the new text) | 100 bp bins in the same GC class, VSA-gene bins vs other genic bins | +0.075 rep1, +0.126 rep2, +0.144 RAPA, +0.049 merged, +0.019 DMSO (CI includes 0) |
| B. Per-track GC normalisation, then gene means | each IP/INPUT divided by the mean coverage of its GC class | var +0.06 (rep1), +0.10 (rep2), -0.03 (DMSO), +0.01 (RAPA); rif -0.05 (rep1) to +0.06; mixed sign |
| C. Gene-level GC residual | raw gene mean minus a quadratic fit on gene GC of non-VSA genes | about 0 to negative (e.g. var -0.04 rep1, 0.00 rep2, -0.21 DMSO) |

Why they differ: only 29 non-VSA genes reach the GC of VSA genes (28 % or more), so method C extrapolates a curve to a region where it has almost no data, and method A has very few non-VSA bins at high GC. The honest statement is "an excess is visible at equal GC with bin-level matching, but not robust to the way GC is controlled". Suggested wording adjustment in the docx: replace "in three of four datasets" by "with bin-level GC matching in three of four datasets, but not with gene-level GC modelling".

## 3. Per original figure, in detail

**Fig 2C (VSA genes).** Original: significant enrichment of var, rif, stevor vs all other genes (Welch). Reanalysis: raw gene-level differences hold in all 4 datasets (var +0.12 to +0.35, rif +0.09 to +0.38, stevor +0.11 to +0.34). After GC control see section 2. Verdict: softened. Docx: Results 2 (top-200/GO note), Fig 4 section (new paragraphs), Discussion (claims toned down, Limitations paragraph). Figure: original SVG not changed.

**Fig 2C RAPA vs DMSO bracket / S1D.** Original: more uracil in RAPA, significant for var. Reanalysis: raw difference var +0.169, rif +0.100 (p < 1e-4); after per-track GC normalisation -0.0004 and -0.040; gene-GC residual -0.010 and +0.011. Dot blot (wet lab) still shows more uracil overall. Verdict: moved to supplement as "more uracil overall, no detectable change in distribution". Docx: Results 2 rewritten, caption of Fig 2 shortened, Fig S4 added. Figure: new Fig S4; the stray "**" in the original SVG could not be removed safely.

**Fig 3B (intragenic profiles).** Raw shapes hold in all 4 datasets. Last exon of rif/stevor raw +0.21 to +0.54, after GC normalisation +0.02 to +0.21; first exon of var raw +0.25 to +0.50, after +0.08 to +0.13. These exons are GC-rich (33 % and 36 %). Verdict: descriptive only. Docx: one added sentence in Results 3.

**Fig 4B.** Original r = 0.42; reproduced 0.44 (225,867 bins). Kept; now panel B of the new Fig 4.

**Fig 4C (G4).** Original: median distance 0, 25.5 % of hotspots within 100 bp of a G4 vs 1.1 % random. Reproduced exactly in the merged data (25.5 % vs 1.1 %). Per dataset (top 4,034 bins): 26 % rep1, 12 % rep2, 11 % RAPA, 0.4 % DMSO. Against GC-matched nulls fold falls to 1.25 (rep1), 1.50 (rep2), 1.31 (merged), 0.61 (RAPA), 0.07 (DMSO); G4 adds < 0.001 R2 beyond GC. Verdict: dropped from main text. Docx: Fig 4 section (short G4 paragraph), Discussion (mechanism and RPA/G4 hypothesis removed), Methods, Table S2.

**Genic vs intergenic (new).** Genic bins 20.1 % GC, intergenic 14.5 %. Raw difference +0.09 to +0.21; about half survives per-track GC normalisation; at equal GC +0.02 to +0.08.

## 4. What I can and cannot do with the SVGs

| Kind of original SVG | Can do | Cannot do |
|---|---|---|
| Live text labels and shapes (editable in Inkscape/Illustrator) | change text, labels, captions, remove an annotation if it is a separate text/path object | edit pixels inside embedded images |
| Data panels that code generates (2B ideogram, 2C violin, 3A, 3B heatmaps, 4A locus, 4B GC, S1 panels) | regenerate from the tracks with new data or a new layout (the repo scripts write SVG) and swap the panel in | match the original look exactly; the originals were composed by hand in CorelDraw, regenerated panels use a plain style |
| Hand-composed microscopy and schematics (Fig 1 KO figure, Fig 2A, figure_microscope.svg) | leave as is, or edit text labels | redraw or "reanalyse" them; embedded images stay as they are |

Applied so far: Fig 2 SVG copied unchanged (the visible RAPA/DMSO "**" is not an editable text object; removing it needs Inkscape or a raster edit by you). New Fig 4 and S4 are new generated SVGs, not edits of the originals. The docx embeds flattened TIFFs, so the new figures still have to be placed into the document.

## 5. Where each analysis lives (so it is not scattered)

| Question | Folder in `analyses/` |
|---|---|
| Are the published tracks right? | `00_tracks`, `qc_library_bias`, `alignment/` (FASTQ rebuild) |
| Is 1 bp to 10 bp averaging a problem? | `00b_bin_approximation`, `00c_tracks_10bp` (averaging shrinks signal to 0.77-0.82) |
| Do Figs 2C/3B/4B/4C hold per dataset? | `fig3B_per_dataset`, `fig4C_per_dataset`, `replicate_consistency` |
| How strong is the GC bias? | `qc_library_bias`, `fig4B_gc_correlation`, `fig3B_gc_profile` |
| What survives GC control? | `gc_corrected`, `fig4C_gc_matched_null`, `genic_vs_intergenic`, `subtelomere_position` |
| Top genes, GO, differential patterns | `gene_ranking_go` |
| New figures | `fig_gcbias`, `fig_supp` |

`python run_all.py` runs everything in this order and writes `REPRODUCIBILITY_SUMMARY.md`.
