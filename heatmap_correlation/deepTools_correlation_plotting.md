# deepTools correlation plotting

## Heatmap using `plotCorrelation`:

1. First we need to create the correct form of data inputs for `plotCorrelation` by `multiBigwigSummary` from .bw (bigWig) files to .npz:
    
    ```bash
    multiBigwigSummary bins \
     -b ./0421.Plasmodium.MAPQfiltered.IPvsINPUT.correct.log2.bw \
    ./0815.Plasmodium.MAPQfiltered.IPvsINPUT.log2.bw \
    ./PlasmodiumH3K4me3.IPvsINPUT.log2.bw \
    ./PlasmodiumH3K9ac.IPvsINPUT.log2.bw \
    ./PlasmodiumH3K9me3.IPvsINPUT.log2.bw \
    --labels 1st_UDNASeq 2nd_UDNASeq H3K4me3 H3K9ac H3K4me3 \
    --binSize 100 \
    --numberOfProcessors "max" -o heatmap_data.npz
    ```
    
    or
    
    ```bash
    multiBigwigSummary BED-file \
     -b bw_1bp/ \
    --labels 3d7_rep1 3d7_rep2 DMSO RAPA \
    --binSize 10 \
    --numberOfProcessors "max" -o var_bin10bp_heatmap.npz \
    --BED ../../BED_files/var_genes/var_genes.bed
    ```
    
    or
    
    ```bash
    multiBigwigSummary bins \
     -b bw_1bp/*.bw \
    --labels 3d7_rep1 3d7_rep2 DMSO RAPA \
    --binSize 10000 \
    --numberOfProcessors "max" -o heatmap_data_all.npz
    ```
    
    or along Bed files
    
    ```bash
    multiBigwigSummary BED-file --BED ../Refgen/PlasmoDB-57_Pfalciparum3D7.gff \
    -b ./0421.Plasmodium.MAPQfiltered.IPvsINPUT.correct.log2.bw \
    ./0815.Plasmodium.MAPQfiltered.IPvsINPUT.log2.bw \
    ./PlasmodiumH3K4me3.IPvsINPUT.log2.bw \
    ./PlasmodiumH3K9ac.IPvsINPUT.log2.bw \
    ./PlasmodiumH3K9me3.IPvsINPUT.log2.bw \
    --labels 1st_UDNASeq 2nd_UDNASeq H3K4me3 H3K9ac H3K4me3 \
    --binSize 100 \
    --numberOfProcessors "max" -o heatmap_data.npz
    ```
    
2. By using the generated .npz format file we can create heatmap:
    
    ```bash
    plotCorrelation -in heatmap_data_all.npz -c spearman \
    -p heatmap --plotTitle "Spearman Correlation of Read Counts" \
    --colorMap RdYlBu_r --plotNumbers \
    -o heatmap_SpearmanCorr_readCounts_all_10000bp.svg
    ```
    
3. To plot scatterplot, use the same .npz file:

```bash
plotCorrelation -in allbw_correlation.npz \
--corMethod spearman -p scatterplot \
--colorMap RdYlBu_r \
--plotNumbers \
--xRange 0 1 \
--plotFile correlation_scatterplot.pdf
```