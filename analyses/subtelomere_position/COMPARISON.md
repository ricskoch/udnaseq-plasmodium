# Subtelomeric position vs VSA identity at equal GC

Contrast = log2U difference to non-VSA genic bins in internal regions, equal-GC, raw tracks. Subtelomeric = within 100 kb of a chromosome end.

| dataset    | group              |   n_bins |   mean_gc |   equalGC_vs_other_internal |      lo |      hi |
|:-----------|:-------------------|---------:|----------:|----------------------------:|--------:|--------:|
| 3d7_rep1   | VSA/subtelomeric   |     5321 |     0.31  |                      0.0726 |  0.0599 |  0.0867 |
| 3d7_rep1   | VSA/internal       |     2026 |     0.324 |                      0.101  |  0.0628 |  0.1574 |
| 3d7_rep1   | other/subtelomeric |     5836 |     0.198 |                      0.0225 | -0.0009 |  0.0446 |
| 3d7_rep2   | VSA/subtelomeric   |     5321 |     0.31  |                      0.1252 |  0.1158 |  0.1362 |
| 3d7_rep2   | VSA/internal       |     2026 |     0.324 |                      0.1299 |  0.0904 |  0.1863 |
| 3d7_rep2   | other/subtelomeric |     5836 |     0.198 |                     -0.0393 | -0.0558 | -0.0244 |
| DMSO       | VSA/subtelomeric   |     5321 |     0.31  |                      0.0339 |  0.0001 |  0.0656 |
| DMSO       | VSA/internal       |     2026 |     0.324 |                     -0.0139 | -0.0501 |  0.0447 |
| DMSO       | other/subtelomeric |     5836 |     0.198 |                      0.0268 |  0.0036 |  0.0445 |
| RAPA       | VSA/subtelomeric   |     5321 |     0.31  |                      0.1515 |  0.0992 |  0.1975 |
| RAPA       | VSA/internal       |     2026 |     0.324 |                      0.1358 |  0.0706 |  0.238  |
| RAPA       | other/subtelomeric |     5836 |     0.198 |                      0.0362 |  0.0121 |  0.0543 |
| 3d7_merged | VSA/subtelomeric   |     5321 |     0.31  |                      0.0461 |  0.033  |  0.0584 |
| 3d7_merged | VSA/internal       |     2026 |     0.324 |                      0.0684 |  0.0304 |  0.131  |
| 3d7_merged | other/subtelomeric |     5836 |     0.198 |                      0.0067 | -0.0101 |  0.0252 |

Point estimates:

| dataset    |   VSA/internal |   VSA/subtelomeric |   other/subtelomeric |
|:-----------|---------------:|-------------------:|---------------------:|
| 3d7_merged |          0.068 |              0.046 |                0.007 |
| 3d7_rep1   |          0.101 |              0.073 |                0.022 |
| 3d7_rep2   |          0.13  |              0.125 |               -0.039 |
| DMSO       |         -0.014 |              0.034 |                0.027 |
| RAPA       |          0.136 |              0.152 |                0.036 |
