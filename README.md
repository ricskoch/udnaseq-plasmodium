# PCA Analysis: INPUT vs UDNAIP Conditions

This analysis performs Principal Component Analysis (PCA) on bigwig files to compare INPUT and UDNAIP conditions across different experiments (3d7_rep1, 3d7_rep2, DMSO, RAPA).

## Prerequisites

- Conda/Miniconda installed
- Bigwig files in the current directory:
  - `3d7_rep1.Plasmodium.INPUTcorrect.bin100bp.RPGC.bw`
  - `3d7_rep1.Plasmodium.UDNAIPcorrect.bin100bp.RPGC.bw`
  - `3d7_rep2.Plasmodium.INPUT.bin100bp.RPGC.bw`
  - `3d7_rep2.Plasmodium.UDNAIP.bin100bp.RPGC.bw`
  - `DMSO.INPUT.bin100bp.RPGC.bw`
  - `DMSO.UDNAIP.bin100bp.RPGC.bw`
  - `RAPA.INPUT.bin100bp.RPGC.bw`
  - `RAPA.UDNAIP.bin100bp.RPGC.bw`

## Setup

### 1. Create the conda environment

```bash
conda create -n reproducibility-analysis python=3.11 -y
```

### 2. Install required packages

```bash
conda run -n reproducibility-analysis pip install pyBigWig numpy scipy scikit-learn matplotlib seaborn pandas
```

## Running the Analysis

Execute the analysis script:

```bash
conda run -n reproducibility-analysis python pca_analysis.py
```

## Output Files

The script generates the following outputs:

- **`pca_plot.png`** - High-resolution PCA scatter plot (300 DPI)
- **`pca_plot.pdf`** - Vector format PCA plot for publications
- **`pca_results.csv`** - PC scores for all samples

## Results Interpretation

The PCA analysis extracts signal values from bigwig files in 100kb bins and performs dimensionality reduction. The results show:

- **PC1** explains ~69% of variance - represents overall signal distribution differences between experiments
- **PC2** explains ~15% of variance - captures differences between INPUT and UDNAIP conditions

### Plot Legend

**Colors:**
- Blue: INPUT condition
- Orange: UDNAIP condition

**Markers:**
- Circle (○): 3d7_rep1
- Square (□): 3d7_rep2
- Triangle (△): DMSO
- Diamond (◇): RAPA

## Quick Reproduction

```bash
# One-liner to set up and run everything
conda create -n reproducibility-analysis python=3.11 -y && \
conda run -n reproducibility-analysis pip install pyBigWig numpy scipy scikit-learn matplotlib seaborn pandas -q && \
conda run -n reproducibility-analysis python pca_analysis.py
```
