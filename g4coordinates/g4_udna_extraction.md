## Packages

- bedtools
- ggplot2
1. Extract the 100 bp bins that have an average log2 enrichment of above 0.5 (that corresponds to around 40% increase in uracil)
    
    ```bash
    awk '$4 > 0.5' data/log2_bw/3d7_merged.IP.vs.INPUT.bin100bp.log2.bedGraph > data/log2_bw/3d7_merged.blacklisted.above0p5.bedGraph
    ```
    
2. Merging “overlapping” intervals:
    
    ```bash
    bedtools merge -i data/log2_bw/3d7_merged.blacklisted.above0p5.bedGraph -c 4 -o mean > data/log2_bw/3d7_merged.blacklisted.above0p5.bedGraph
    ```

3. Check the distance to the closest interval from A to B and look for both downstream and upstream, without overlaps for better representation on a histogram.
bedtools closest -a Merged/merged_above0p5_blacklisted_merged_intervals.bedGraph -b CHIP-seq_data/BED_files/G4_coordinates.bed -D ref > 0p5_to_G4_coordinates.bed
    
    ```bash
    bedtools closest -a data/log2_bw/3d7_merged.blacklisted.above0p5.bedGraph -b data/BED_files/G4_coordinates.bed -D ref > g4coordinates/g4coordinates/0p5_to_G4_coordinates.bed
    ```
    
    Note: this reports all overlapping G4 coordinates as 0, even if there are several G4s in a given interval above 0.5 log2