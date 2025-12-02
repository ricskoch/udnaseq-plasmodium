library('circlize')
bed <- as.data.frame(read.table("0815.Plasmodium.MAPQfiltered.IPvsINPUT.log2.bedGraph",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed_min = min(bed$V4)
bed_max = max(bed$V4)
circos.clear()
circos.par('start.degree' = 90)
circos.par('gap.degree' = 2)
circos.initializeWithIdeogram(cytoband = pfa3d7_cytoband, chromosome.index = paste('Pf3D7_', sprintf("%02d", c(1:14)), '_v3', sep = ''), plotType = NULL)
circos.genomicTrack(`0815.Plasmodium.MAPQfiltered.IPvsINPUT.log2`, ylim=c(0.4, 1.9), panel.fun = function(region, value, ...){
  #COMMENT     Getting rid of unnecessery values, threshold is the minimum value of the defined ylim in genomicTrack
  adj_val = ifelse(value$V4 <= CELL_META$ylim[1], CELL_META$ylim[1], value$V4)
  
  #COMMENT     Min-max normalization for the whole dataframe, so uracil content of chromosomes ARE proportional to one another
  adj_col_to_df = (adj_val - bed_min) / (bed_max - bed_min)
  #COMMENT     The problem is that transition is not emphasized...
  circos.genomicPoints(region, adj_val, pch=16, cex=0.6, col=rgb(adj_col_to_df, 1-adj_col_to_df, 0, 1))
  #circos.genomicPoints(region, adj_val, pch=16, cex=0.7, col=rgb(adj_col_to_df, (1-adj_col_to_df)*0.8, (1-adj_col_to_df)*0.8, 1))
  
  #COMMENT     Min-max normalization for each track, so uracil content of chromosomes are NOT proportional to one another
  #adj_col_to_track = (adj_val - min(adj_val)) / (max(adj_val) - min(adj_val))
  #COMMENT     Although transition in each track is emphasized, there is NOT a strict proportion among chromosomes...
  #circos.genomicPoints(region, adj_val, pch=16, cex=0.6, col=rgb(adj_col_to_track, 1-adj_col_to_track, 0, 1))
  
  #Draws an opacity of 50% black baseline, adjusted to the defined ylim in genomicTrack
  #circos.genomicLines(region, CELL_META$ylim[1], col='#00000050')
  
  #Draws an opacity of 50% black line to the maximum value of the defined ylim in genomicTrack, it slows down graphing significantly
  #circos.genomicLines(region, max(CELL_META$ylim[2]), col='#00000050')
  
}, bg.border = NA)
circos.genomicTrack(vargenes, panel.fun = function(region, value, ...){
  circos.genomicRect(region, value)
},track.height = 0.075, bg.border = NA)


## To present interchromosomal rearrangements use the following:

chrStart = NULL
chrEnd = NULL
posStart = NULL
posEnd = NULL

for (i in 1:1050) {
  traInfo = strsplit(as.character(tra[i,8]), split = ";")
  traQual = as.character(tra[i,7])
  if ((traInfo[[1]][1] == 'PRECISE') & (traQual == 'PASS')) {
    chrInfEnd = strsplit(as.character(traInfo[[1]][5]), split="=")
    chrEnd = c(chrEnd, chrInfEnd[[1]][2])
    posInfEnd = strsplit(as.character(traInfo[[1]][6]), split="=")
    posEnd = c(posEnd, posInfEnd[[1]][2])
    chrStart = c(chrStart, tra[i,1])
    posStart = c(posStart, tra[i,2])
  }
}
traStart = data.frame(chrStart, as.numeric(posStart))
traEnd = data.frame(chrEnd, as.numeric(posEnd))
circos.genomicLink(traStart, traEnd)


#That is also quite good
#circos.genomicPoints(region, adj_val, pch=16, cex=0.7, col=rgb(adj_col_to_df, (1-adj_col_to_df)*0.8, (1-adj_col_to_df)*0.8, 1))
