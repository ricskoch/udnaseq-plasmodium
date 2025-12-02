## Package

- bedtools

## Steps

1. Generate the blacklisted bed file from the merged dataset bedGraph file
    The blacklist was manually chosen to get rid of the telomeres and the Pf11-1 gene, that are overrepresented and probably have a lot of copy number variation

```bash
bedtools subtract -a /data/log2_bw/3d7_merged.IP.vs.INPUT.bin100bp.log2.bedGraph -b /data/BED_files/blacklist_telomer_pf11.txt > /data/log2_bw/3d7_merged.blacklisted.IP.vs.INPUT.bin100bp.log2.bedGraph
```

2. Calculate nucleotide composition of the bed file using bedtools nuc and print the coordinates, log2 U-enrichment and pct_gc values into a bed file that is used in R for correlation plot

```bash
bedtools nuc -fi data/Refgen/PlasmoDB-57_Pfalciparum3D7_Genome.fasta -bed data/log2_bw/3d7_merged.blacklisted.IP.vs.INPUT.bin100bp.log2.bedGraph | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6}' > GCcontent/GCcorrelation/averageGCcontent_100bin_merged_blacklisted.bed
```