# computeMatrix 

bash
computeMatrix scale-regions -S IP.vs.INPUT.merged.log2.bw \
-R ../CHIP-seq_data/BED_files/rifinstrevor_genes.bed \
-a 1000 -m 2000 -b 1000 --binSize 10 -o rifinstevor_scaled_matrix.tab


# plotHeatmap

bash
plotHeatmap -m var_scaled_matrix.tab  \
-out UDNA_heatmap_over_rifstevor_scaled.svg \
--averageTypeSummaryPlot mean --plotType std \
--colorMap jet --heatmapWidth 12 \
--yAxisLabel 'log2 U-score' \
--plotTitle "Var gene U-enrichment with 10bp binsize"\
--outFileSortedRegions var_sortedregions.tab \
--missingDataColor white --sortUsing sum \
--linesAtTickMarks --alpha 0.9 
