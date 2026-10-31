# Legacy analyses (reference only)

The original R scripts, deepTools commands, notebook and their outputs, moved here unchanged during the
Python rework. The new pipeline in `analyses/` reads these files only to check that it reproduces them.

| Legacy folder | Replaced by |
|---|---|
| `data/log2_bw/` | `analyses/00_tracks` |
| `heatmap_correlation/` | `analyses/figS1AB_sample_correlation` |
| `pca_analysis/` | `analyses/qc_pca` |
| `gene_enrichment/` | `analyses/fig2C_3A_S1D_S2_gene_enrichment` |
| `var_heatmaps/` | `analyses/fig3B_vsa_metaprofiles` |
| `GCcontent/` | `analyses/fig4B_gc_correlation` |
| `g4coordinates/` | `analyses/fig4C_g4_distance`, `analyses/fig2B_S1C_genome_distribution` (S1C) |

Several legacy scripts depend on objects from an interactive R session (`dataset_4`, `combined`,
`merged_blacklisted_GCcontent`, ...) or on files that were never committed, so they cannot be run as they are.
