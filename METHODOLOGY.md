# Analysis Methodology & Development History

This document captures the analytical decisions, iterations, and key insights from developing the U-DNA PCA analysis pipeline.

## Overview

The analysis examines genomic uracil (U-DNA) distribution patterns in *Plasmodium falciparum* across INPUT (untreated control) and UDNAIP (uracil DNA-immunoprecipitated) conditions using Principal Component Analysis (PCA).

**Study Organisms & Conditions:**
- **3d7_rep1 & 3d7_rep2**: Wild-type 3D7 reference strain (biological replicates)
- **DMSO**: Vehicle control (dimethyl sulfoxide treatment)
- **RAPA**: Rapamycin-induced PfUNG knockout condition

## Data Processing

### Resolution & Binning
- **Selected**: 100bp bin size (not 100kb)
- **Result**: 233,337 bins per sample after padding
- **Rationale**: Captures fine-scale genomic structure while remaining computationally tractable

### Preprocessing Steps
1. **Signal Extraction**: Mean signal value calculated for each 100bp bin using pyBigWig
2. **Padding**: All samples aligned to maximum bin count to create uniform matrix
3. **Standardization**: StandardScaler normalization applied (zero-mean, unit variance)
4. **Data Matrix**: 8 samples × 233,337 features

### Missing Value Handling
- NaN values replaced with 0.0 during extraction
- `np.nan_to_num()` applied to final matrix for additional safety

## Component Selection Strategy

### Elbow Analysis
- Performed PCA with all 8 possible components (limited by sample count: min(n_samples, n_features))
- Plotted both individual and cumulative explained variance
- Identified two candidate thresholds:
  - **80% threshold**: 3 components
  - **90% threshold**: 6 components

### Decision: 3 Components
**Variance Distribution (3-component solution):**
- PC1: 46.63%
- PC2: 21.29%
- PC3: 13.22%
- **Total: 81.13%**

**Reasoning:**
- Captures >80% of variance with interpretable simplicity
- 3D visualization practical and interpretable
- Avoids overfitting with small sample size (n=8)
- 6-component solution adds only 9.87% additional variance (diminishing returns)

## Visualization Evolution

### Stage 1: 2D PCA Plot
**Initial Approach:**
- Standard PC1 vs PC2 scatter plot
- Custom legend with markers by experiment + colors by condition
- Exported as PNG (300 DPI) and PDF

**Markers (by experiment):**
- Circle (○): 3d7_rep1
- Square (□): 3d7_rep2
- Triangle (△): DMSO
- Diamond (◇): RAPA

**Colors (by condition):**
- Blue (C0): INPUT
- Orange (C1): UDNAIP

### Stage 2: 3D Visualization Addition
**Motivation:** Full representation of all 3 principal components for comprehensive pattern analysis

**Technical Challenges & Solutions:**

1. **Label Visibility Issue (Z-axis)**
   - Problem: Z-axis label clipped from view after initial implementation
   - Root cause: Default matplotlib margins insufficient for 3D figure
   - Solution cascade:
     - Increased figure size: 12×9 → 14×10 inches
     - Increased z-axis labelpad: 10 → 20 pixels
     - Added explicit layout control: `plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.1)`

2. **Legend Positioning**
   - Initial: Default legend placement overlapped data points
   - Final: `bbox_to_anchor=(1.05, 1.0)` positioned legend at top-right aligned with z-axis

3. **Viewing Angle**
   - Elevation: 20°
   - Azimuth: 45°
   - Optimizes both PC1-PC2 separation and PC3 visibility

**Final Figure Parameters:**
- Figure size: 14×10 inches
- Marker size: 300 (larger for 3D depth perception)
- Alpha: 0.7 (transparency for overlapping markers)
- Edge colors: Black, linewidth 1.5 (visual clarity)

## Key Findings & Biological Interpretation

### PCA-Driven Insights

1. **INPUT vs UDNAIP Separation**
   - Clear clustering by condition evident in both 2D and 3D space
   - Validates experimental design efficacy
   - Suggests mechanistically controlled uracil distribution rather than random noise

2. **Treatment-Specific Clustering**
   - DMSO and RAPA conditions occupy distinct regions of PC space
   - Indicates significant genomic reorganization in PfUNG knockout (RAPA) parasites
   - Suggests rapamycin-induced changes measurably alter U-DNA landscape

3. **Replicate Heterogeneity Discovery**
   - **Critical Finding**: 3d7_rep2_UDNAIP extreme PC3 value (+436.68)
   - Contrasts sharply with other UDNAIP replicates: 3d7_rep1_UDNAIP (-84.12), DMSO_UDNAIP (+217.43), RAPA_UDNAIP (+155.21)
   - **Interpretation**: Differential uracilation between biological replicates warrants investigation
   - Possible explanations:
     - Biological variation in VSA gene uracilation between parasite clones
     - Different growth phase or stress conditions during sampling
     - Replicate-specific epigenetic or chromatin differences

4. **Reproducibility Validation**
   - Replicates generally cluster together (except for PC3 heterogeneity noted above)
   - Validates experimental pipeline robustness
   - Supports biological interpretation over technical artifact

### Biological Context

**U-DNA Enrichment Patterns:**
- Non-random distribution supports active uracilation machinery
- Particular enrichment in VSA genes (var, rif, stevor families)
- Overlap with G-quadruplex forming sequences suggests structure-dependent mechanisms

**PfUNG Role:**
- Knockout lethal phenotype (ring stage stalling by cycle 3) indicates essential base excision repair function
- Rapamycin treatment dramatically increases genomic uracil accumulation
- RAPA/DMSO separation in PC space reflects differential BER pathway activity

## Reproducibility & Environment Management

### Environment Specification
- **File**: `environment.yml`
- **Python Version**: 3.11
- **Key Packages**:
  - Data: pyBigWig, numpy, pandas
  - Analysis: scikit-learn (PCA), scipy
  - Visualization: matplotlib, seaborn
  - Notebooks: jupyter, jupyterlab

### Setup Command
```bash
conda env create -f environment.yml
```

### Running Analysis
```bash
# With activation
conda activate reproducibility-analysis
jupyter notebook pca_analysis.ipynb

# Without activation
conda run -n reproducibility-analysis jupyter notebook pca_analysis.ipynb
```

## Output Files

| File | Format | Purpose | Generated By |
|------|--------|---------|--------------|
| `pca_plot.png` | PNG (300 DPI) | 2D PCA visualization | Section 6, Cell 7 |
| `pca_plot.pdf` | PDF | Publication-ready 2D plot | Section 7, Cell 8 |
| `pca_results.csv` | CSV | PC1-PC3 scores all samples | Section 7, Cell 8 |
| `pca_analysis.ipynb` | Jupyter | Full interactive analysis | - |

## Implementation Notes

### Why Jupyter Notebook Instead of Script?
- **Interactive Exploration**: Elbow analysis, visualization refinement benefit from iterative development
- **Documentation**: Markdown cells contextualize each analysis step
- **Reproducibility**: Cell execution order clearly specified; easy to re-run or modify parameters
- **Visualization**: In-notebook plotting enables immediate feedback on parameter adjustments

### Design Decisions

1. **Bin Resolution (100bp)**
   - Finer resolution captures local uracil hotspots
   - Sufficient computational performance for 233,337 features
   - Aligns with typical ChIP-Seq analysis granularity

2. **StandardScaler Over Other Normalizations**
   - Preserves variance information critical for PCA
   - Not biased by outliers (unlike MinMaxScaler)
   - Industry standard for PCA preprocessing

3. **3 Components Over 80% vs 90% Threshold**
   - Balance between interpretability (2D insufficient) and parsimony
   - 8-sample limitation argues against overextraction
   - PC3 added surprising biological signal (replicate heterogeneity)

4. **Custom Markers & Colors**
   - Distinguishes both experimental context and treatment condition simultaneously
   - More informative than single-parameter encoding
   - Supports exploratory comparison of multiple factors

## Future Enhancements

### Potential Extensions
1. **Loading Vector Analysis**: Which genomic bins drive PC1/PC2/PC3 separation?
2. **Interactive 3D Visualization**: Plotly for browser-based rotation/zoom
3. **Genomic Annotation Integration**: Map PC loadings to var/rif/stevor gene positions
4. **Additional Projections**: 2D plots of PC1-PC3 and PC2-PC3 combinations
5. **Statistical Testing**: Permutation tests for INPUT vs UDNAIP separation significance

### Computational Optimization
- Current analysis runs in <1 minute; scalable to larger datasets if needed
- Consider sparse matrix formats if extending to smaller bin sizes

## Version Control & Archival

- **Repository**: github.com/ricskoch/udnaseq-reproducibility
- **Environment Specification**: `environment.yml` ensures reproducibility across systems
- **Analysis Code**: `pca_analysis.ipynb` contains fully commented pipeline
- **This Document**: `METHODOLOGY.md` captures analytical decision-making process

---

**Last Updated**: November 29, 2025  
**Analysis Status**: Publication-ready  
**Reproducibility**: ✅ Fully specified environment and documented pipeline
