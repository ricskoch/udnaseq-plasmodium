# Genic vs intergenic log2U: raw, GC-normalised and at equal GC

(100 bp bins, nuclear, non-blacklisted; genic = bin centre inside a protein-coding gene span, else intergenic)
    A. Genic minus intergenic mean log2U: raw, and after per-track GC normalisation (gc_bin).
    B. Equal-GC contrast: 20 quantile classes of bin GC; genic-minus-intergenic mean log2U within each class, averaged with weights
       min(n_genic, n_intergenic). Uses raw tracks, so the GC confound is removed without touching real differences at equal GC.
    C. Same equal-GC contrast for VSA-gene bins vs other-gene bins (both genic).
    Bootstrap by 100 kb blocks (200 resamples) for the 95 % CI of B and C.

Bins: genic 181156, intergenic 45702; mean GC genic 20.1 %, intergenic 14.5 %.

| dataset    |   genic_minus_inter_raw |   genic_minus_inter_gcnorm |   retained_pct |   equalGC_genic_minus_inter |   B_lo |   B_hi |   equalGC_VSA_minus_otherGenic |    C_lo |   C_hi |
|:-----------|------------------------:|---------------------------:|---------------:|----------------------------:|-------:|-------:|-------------------------------:|--------:|-------:|
| 3d7_rep1   |                  0.0884 |                     0.0391 |        44.2709 |                      0.0285 | 0.019  | 0.0375 |                         0.0751 |  0.0595 | 0.0937 |
| 3d7_rep2   |                  0.2123 |                     0.1109 |        52.2109 |                      0.0766 | 0.0685 | 0.0832 |                         0.1258 |  0.1123 | 0.1369 |
| DMSO       |                  0.1856 |                     0.0912 |        49.1372 |                      0.0553 | 0.0466 | 0.0617 |                         0.0192 | -0.0042 | 0.045  |
| RAPA       |                  0.1397 |                     0.0662 |        47.4322 |                      0.0191 | 0.0107 | 0.0273 |                         0.1436 |  0.0964 | 0.1827 |
| 3d7_merged |                  0.0937 |                     0.0411 |        43.8472 |                      0.0326 | 0.0221 | 0.0433 |                         0.0492 |  0.032  | 0.0652 |
