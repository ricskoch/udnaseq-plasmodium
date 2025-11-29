# PCA Analysis: INPUT vs UDNAIP Conditions

This analysis performs Principal Component Analysis (PCA) on bigwig files to compare INPUT and UDNAIP conditions across different experiments (3d7_rep1, 3d7_rep2, DMSO, RAPA).

For detailed methodology, design decisions, and development history, see [`METHODOLOGY.md`](METHODOLOGY.md).

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

### 1. Create the conda environment from the environment file

```bash
conda env create -f environment.yml
```

This will create the `reproducibility-analysis` environment with all required dependencies specified in `environment.yml`.

### 2. Verify the environment

```bash
conda activate reproducibility-analysis
python --version  # Should show Python 3.11.x
```

## Running the Analysis

Execute the analysis in Jupyter notebook:

```bash
conda run -n reproducibility-analysis jupyter notebook pca_analysis.ipynb
```

Or run all cells to generate outputs automatically.

## Output Files

The notebook generates the following outputs:

- **`pca_plot.png`** - 2D PCA scatter plot (PC1 vs PC2, 300 DPI)
- **`pca_plot.pdf`** - Vector format 2D PCA plot for publications
- **3D visualization** - Interactive 3D scatter plot (PC1 vs PC2 vs PC3) displayed in notebook
- **`pca_results.csv`** - PC1, PC2, PC3 scores for all 8 samples

## Results Interpretation

### Data Processing

- **Resolution**: 100bp bin size (233,337 bins per sample)
- **Samples**: 8 total (2 conditions × 4 experiments)
- **Standardization**: StandardScaler normalization applied before PCA

### Optimal Component Selection

Elbow analysis determined **3 principal components** as optimal, explaining **81.13% total variance**:

- **PC1**: 46.63% variance - captures overall U-DNA enrichment levels
- **PC2**: 21.29% variance - represents condition-specific differences between INPUT/UDNAIP
- **PC3**: 13.22% variance - reveals biological heterogeneity across replicates

### Key Findings

1. **Clear INPUT vs UDNAIP separation** visible in 2D and 3D space
2. **Distinct treatment clustering**: DMSO and RAPA conditions show separate clustering patterns
3. **Replicate heterogeneity**: 3d7_rep2_UDNAIP shows extreme PC3 value (+436.68), revealing differential uracil enrichment patterns between biological replicates
4. **Reproducibility**: Validates experimental design with consistent sample separation

### Plot Legend

**Colors:**
- Blue: INPUT condition
- Orange: UDNAIP condition

**Markers:**
- Circle (○): 3d7_rep1
- Square (□): 3d7_rep2
- Triangle (△): DMSO
- Diamond (◇): RAPA

### Biological Context

This analysis examines genomic uracil (U-DNA) distribution in *Plasmodium falciparum* across:
- **3d7_rep1 & 3d7_rep2**: Wild-type 3D7 reference strain (biological replicates)
- **DMSO**: Vehicle control
- **RAPA**: Rapamycin-induced PfUNG knockout

The structured separation supports the hypothesis that uracil distribution is mechanistically controlled, particularly in VSA (Variant Surface Antigen) gene regions involved in antigenic variation.

## Quick Reproduction

```bash
# Set up environment from file and run notebook
conda env create -f environment.yml && \
conda activate reproducibility-analysis && \
jupyter notebook pca_analysis.ipynb
```

Alternatively, without activating the environment:

```bash
conda run -n reproducibility-analysis jupyter notebook pca_analysis.ipynb
```
