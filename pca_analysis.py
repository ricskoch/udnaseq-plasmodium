#!/usr/bin/env python3
"""
PCA Analysis of UDNAseq BigWig Files

This script performs PCA analysis on bigwig files comparing INPUT vs UDNAIP
conditions across different experiments (3d7_rep1, 3d7_rep2, DMSO, RAPA).
"""

import pyBigWig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 8)

# Define the bigwig files
bigwig_files = [
    "3d7_rep1.Plasmodium.INPUTcorrect.bin100bp.RPGC.bw",
    "3d7_rep1.Plasmodium.UDNAIPcorrect.bin100bp.RPGC.bw",
    "3d7_rep2.Plasmodium.INPUT.bin100bp.RPGC.bw",
    "3d7_rep2.Plasmodium.UDNAIP.bin100bp.RPGC.bw",
    "DMSO.INPUT.bin100bp.RPGC.bw",
    "DMSO.UDNAIP.bin100bp.RPGC.bw",
    "RAPA.INPUT.bin100bp.RPGC.bw",
    "RAPA.UDNAIP.bin100bp.RPGC.bw",
]

# Parse experiment and condition from filename
def parse_filename(filename):
    """Extract experiment and condition from filename."""
    base = Path(filename).stem
    
    if "3d7_rep1" in base:
        experiment = "3d7_rep1"
    elif "3d7_rep2" in base:
        experiment = "3d7_rep2"
    elif "DMSO" in base:
        experiment = "DMSO"
    elif "RAPA" in base:
        experiment = "RAPA"
    else:
        experiment = "Unknown"
    
    if "UDNAIP" in base:
        condition = "UDNAIP"
    else:
        condition = "INPUT"
    
    return experiment, condition


def extract_signals(bw_file, bin_size=100000):
    """Extract signal values from bigwig file in bins."""
    bw = pyBigWig.open(bw_file)
    signals = []
    
    for chrom in sorted(bw.chroms().keys()):
        chrom_length = bw.chroms()[chrom]
        for start in range(0, chrom_length, bin_size):
            end = min(start + bin_size, chrom_length)
            try:
                val = bw.stats(chrom, start, end, type="mean")
                signals.append(val[0] if val and val[0] is not None else 0.0)
            except:
                signals.append(0.0)
    
    bw.close()
    return np.array(signals)


def main():
    """Main analysis function."""
    print("Loading bigwig files and extracting signals...")
    
    data = {'filename': [], 'experiment': [], 'condition': [], 'signals': []}
    
    for bw_file in bigwig_files:
        try:
            print(f"Processing: {bw_file}")
            experiment, condition = parse_filename(bw_file)
            signals = extract_signals(bw_file)
            
            data['filename'].append(bw_file)
            data['experiment'].append(experiment)
            data['condition'].append(condition)
            data['signals'].append(signals)
            print(f"  - Extracted {len(signals)} bins")
        except Exception as e:
            print(f"  ERROR: {e}")
    
    # Pad signals to max length
    max_len = max(len(s) for s in data['signals'])
    print(f"\nMax signal length: {max_len} bins")
    
    X = np.array([np.pad(s, (0, max_len - len(s)), constant_values=0) if len(s) < max_len else s 
                   for s in data['signals']])
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    print(f"Feature matrix shape: {X.shape}")
    
    # Standardize and perform PCA
    print("\nStandardizing features and performing PCA...")
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    print(f"PC1: {pca.explained_variance_ratio_[0]:.2%}")
    print(f"PC2: {pca.explained_variance_ratio_[1]:.2%}")
    
    # Create results dataframe
    results_df = pd.DataFrame({
        'filename': data['filename'],
        'experiment': data['experiment'],
        'condition': data['condition'],
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1]
    })
    
    print("\nPCA Results:")
    print(results_df)
    
    # Create plot
    from matplotlib.lines import Line2D
    
    fig, ax = plt.subplots(figsize=(12, 10))
    colors = {'INPUT': '#1f77b4', 'UDNAIP': '#ff7f0e'}
    markers = {'3d7_rep1': 'o', '3d7_rep2': 's', 'DMSO': '^', 'RAPA': 'D'}
    
    for _, row in results_df.iterrows():
        ax.scatter(row['PC1'], row['PC2'], c=colors[row['condition']], 
                  marker=markers[row['experiment']], s=300, alpha=0.7, 
                  edgecolors='black', linewidth=2)
        ax.annotate(f"{row['experiment']}\n{row['condition']}", 
                   (row['PC1'], row['PC2']), xytext=(5, 5), 
                   textcoords='offset points', fontsize=8, alpha=0.8)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=12)
    ax.set_title('PCA: INPUT vs UDNAIP Conditions Across Experiments', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=10, label='INPUT', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e', markersize=10, label='UDNAIP', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='3d7_rep1', markeredgecolor='black'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', markersize=10, label='3d7_rep2', markeredgecolor='black'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', markersize=10, label='DMSO', markeredgecolor='black'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='gray', markersize=10, label='RAPA', markeredgecolor='black'),
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=10, framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig('pca_plot.png', dpi=300, bbox_inches='tight')
    plt.savefig('pca_plot.pdf', bbox_inches='tight')
    plt.close()
    
    results_df.to_csv('pca_results.csv', index=False)
    print("\nPlots saved to: pca_plot.png, pca_plot.pdf")
    print("Results saved to: pca_results.csv")


if __name__ == "__main__":
    main()
