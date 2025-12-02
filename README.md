# U-DNA-Seq Analysis: Plasmodium falciparum Genomic Uracil Landscape

This repository contains the bioinformatics analysis pipeline and data for the study:

**"Genomic uracil within the Plasmodium falciparum genome localizes to genes involved in antigenic variation of the malaria parasite"**

The manuscript investigates the distribution and functional impact of uracil in the P. falciparum genome, with a focus on variant surface antigen (VSA) gene families, uracil enrichment, G-quadruplex overlap, and the essentiality of uracil DNA glycosylase (PfUNG). For detailed methodology and experimental context, see the attached manuscript.

## Overview

This repository includes:

- **U-DNA-Seq data processing**: Analysis of uracil-enriched DNA using custom pipelines and deepTools.
- **Gene family and region enrichment**: Quantification of uracil in VSA gene families (var, rif, stevor), G4 regions, and more.
- **Heatmaps and correlation plots**: Genome-wide and gene-level visualization of uracil enrichment and sample relationships.
- **Principal Component Analysis (PCA)**: Multivariate analysis of uracil enrichment across samples and conditions.
- **Reference and annotation files**: Genome, BED, and G4 prediction files for region-based analyses.

## Environment Setup

1. **Install Conda/Miniconda** (if not already installed)
2. **Create the analysis environment:**

  ```bash
  conda env create -f environment.yml
  ```

  This will install all required dependencies for the full analysis pipeline (deepTools, numpy, pandas, matplotlib, jupyter, etc.).

3. **Activate the environment:**

  ```bash
  conda activate reproducibility-analysis
  ```

4. **Verify installation:**

  ```bash
  python --version  # Should show Python 3.11.x
  ```

5. **Run analyses and notebooks as described in each subfolder's README or script header.**

## High-Level Directory Structure

```
udnaseq-plasmodium/
├── README.md                # Project overview and setup (this file)
├── METHODOLOGY.md           # Detailed methods and pipeline documentation
├── environment.yml          # Conda environment specification
├── 0425_Plasmodium_UDNA_cell_final.txt  # Manuscript text
├── data/                    # Raw and processed data files
│   ├── *.bw, *.bedGraph     # BigWig and bedGraph files (uracil enrichment, coverage, etc.)
│   ├── BED_files/           # Genomic regions: blacklist, centromeres, G4, gene families
│   ├── log2_bw/             # Log2 enrichment tracks and processed bedGraphs
│   ├── Refgen/              # Reference genome, annotation, and size files
├── gene_enrichment/         # Scripts and outputs for gene-level uracil enrichment
│   ├── gene_enrichment.R    # R script for gene enrichment analysis
│   ├── list_of_average_uenrichment.csv  # Output: average uracil enrichment per gene
├── heatmap_correlation/     # Heatmap and correlation analysis (deepTools, etc.)
│   ├── deepTools_correlation_plotting.md
│   ├── heatmap_data_all_bin100bp_blacklist.npz
│   ├── heatmap_data_all.npz
│   ├── heatmap_data_var.npz
├── pca_analysis/            # PCA analysis and visualization
│   ├── pca_analysis.ipynb   # Main PCA notebook
│   ├── pca_results.csv      # PCA scores
│   ├── pca_plot.png/.pdf    # PCA plots
├── g4coordinates/           # G4 region analysis and documentation
├── GCcontent/               # GC content analysis
├── var_heatmaps/            # Heatmaps for VSA gene families
└── ...                      # Additional scripts, outputs, and documentation
```

## How to Reproduce the Analysis

1. Set up the environment as above.
2. Explore each analysis module (see subfolders) for scripts and notebooks. Most analyses can be reproduced by running the provided scripts or notebooks in the appropriate environment.
3. For PCA, run:

  ```bash
  jupyter notebook pca_analysis/pca_analysis.ipynb
  ```

4. For heatmaps and correlation, see `heatmap_correlation/` and run deepTools commands as described in the markdown files.
5. For gene enrichment, use the R script in `gene_enrichment/`.

## Notes

- For full details on experimental design, data processing, and interpretation, see the manuscript.
- Some supplementary data (e.g., GO enrichment tables, raw ShinyGO outputs, or scripts for U-DNA-Seq processing) may not be present in this repository. If you require these, please contact the corresponding author.

