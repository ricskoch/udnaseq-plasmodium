# U-DNA-Seq analysis: genomic uracil in *Plasmodium falciparum*

Analysis code for **"Genomic uracil within the *Plasmodium falciparum* genome localizes to genes involved in
antigenic variation of the malaria parasite"** (Izrael *et al.*).

Every computational figure panel is regenerated from the raw-most data by one Python script per figure.
Each script compares its output with the previous (legacy) files and the numbers in the manuscript, and
writes the result to a `COMPARISON.md`. See [`REPRODUCIBILITY_SUMMARY.md`](REPRODUCIBILITY_SUMMARY.md)
for the overview and [`REWORK_NOTES.md`](REWORK_NOTES.md) for what was found during the rework.

## Quick start

```bash
conda env create -f environment.yml      # or: pip install -r requirements.txt
conda activate udnaseq-plasmodium
python run_all.py                        # ~7 min plus a one-off ~650 MB download from GEO
```

```bash
python run_all.py fig4C                  # run selected analyses only
```

No R, bedtools or deepTools needed: the analyses use numpy/pandas/scipy/pyBigWig.

## Data flow

```
FASTQ (SRA PRJNA1076816)                                         not in repo
  └─ bwa mem → samtools → picard MarkDuplicates → MAPQ filter
       └─ deepTools bamCoverage --normalizeUsing RPGC
            ├─ data/*.RPGC.bw  (IP + INPUT, 100 bp, per sample)       in repo  ──► 00_tracks
            └─ deepTools bigwigCompare (log2, pseudocount 1)
                 ├─ GEO GSE262467: per-sample 1 bp log2 tracks        fetched  ──► S1A/B, 2C, 3A, S1D
                 └─ GEO GSE262467: merged 3D7 10 bp log2 track        fetched  ──► S1C, 3B, 4A, 4C
```

`scripts/fetch_geo.py` downloads the GEO tracks into `data/geo/` and checks their MD5 sums.
The 1 bp and 10 bp tracks cannot be rebuilt from the 100 bp RPGC tracks, and the legacy figures were made from them.

## Repository layout

```
udnaseq-plasmodium/
├── run_all.py                         runs everything, writes REPRODUCIBILITY_SUMMARY.md
├── src/udnaseq/                       shared code: paths/constants (config.py), bigWig IO,
│                                      genes/blacklist/interval tools, report helpers
├── scripts/fetch_geo.py               GEO download + MD5 check
├── analyses/
│   ├── 00_tracks/                     log2 tracks per sample + merged 3D7 (100 bp) from RPGC coverage
│   ├── figS1AB_sample_correlation/    Fig. S1A/B  correlation heatmaps (100 bp, 10 kb, var genes)
│   ├── qc_pca/                        PCA of the 8 coverage tracks (QC, not in the manuscript)
│   ├── fig2B_S1C_genome_distribution/ Fig. 2B ideogram, Fig. S1C per-chromosome histograms
│   ├── fig2C_3A_S1D_S2_gene_enrichment/ Table S1, top-200 list for GO (Fig. S2), Fig. 2C, S1D, 3A
│   ├── fig3B_vsa_metaprofiles/        Fig. 3B var and rif/stevor scaled heatmaps
│   ├── fig4A_locus_track/             Fig. 4A representative chr7 locus
│   ├── fig4B_gc_correlation/          Fig. 4B GC% vs log2(U-enrichment)
│   └── fig4C_g4_distance/             Fig. 4C distance of uracil hotspots to predicted G4s
│        each: run.py · output/ · COMPARISON.md
├── data/
│   ├── *.RPGC.bw                      IP/INPUT coverage, 4 samples (raw-most data in the repo)
│   ├── Refgen/                        PlasmoDB-57 genome, GFF, chromosome sizes
│   ├── BED_files/                     blacklist, centromeres, G4 predictions, var/rif/stevor genes
│   ├── histone_chipseq/               H3K4me3/H3K9ac/H3K9me3 log2 tracks (not used in the manuscript)
│   └── geo/                           GEO tracks (downloaded, git-ignored)
└── legacy/                            previous R/deepTools/notebook analyses and outputs (reference only)
```

Figures that share one underlying table are grouped in one folder (e.g. Fig. 2C, 3A, S1D and Table S1 all come
from the per-gene means).

## Samples

| Label | Description |
|---|---|
| `3d7_rep1`, `3d7_rep2` | wild-type 3D7 schizonts, two biological replicates |
| `DMSO`, `RAPA` | PfUNG-iKO parasites 88 h after DMSO (mock) or rapamycin (knockout) |
| `3d7_merged` | rep1 + rep2 merged at the BAM level (IP ≈ 89 % rep1 reads, INPUT ≈ 57 %) |

## Not covered by code

- Fig. 1 (growth curves, dot blot, western blot, microscopy): wet-lab data are not in the repository.
- Fig. S2 GO enrichment runs on the ShinyGO web server; `fig2C_3A_S1D_S2_gene_enrichment/output/top200_for_GO.txt` is its input.
- Read mapping (FASTQ → RPGC bigWig) is described above but not scripted here.
- Fig. S3 (protein alignment) is not a data analysis.
