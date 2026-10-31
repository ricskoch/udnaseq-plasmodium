#!/usr/bin/env bash
# Usage: 02_process_library.sh <library_name>      (names in alignment/samples.tsv)
# SRA -> (stream) adapter trim -> bwa mem -> sort -> picard MarkDuplicates (removed)
#     -> MAPQ-filtered BAMs -> RPGC bigWigs (100 bp) -> per-library QC metrics.
# Large intermediates go to $WORK; only small QC files are written under alignment/qc/.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
LIB="$1"
SRR=$(awk -F'\t' -v l="$LIB" '$1==l{print $2}' "$REPO/alignment/samples.tsv")
[ -n "$SRR" ] || { echo "unknown library $LIB"; exit 1; }
D="$WORK/$LIB"; Q="$REPO/alignment/qc/$LIB"; mkdir -p "$D" "$Q"
GENOME="$WORK/ref/genome.fa"
EGS=$(awk '!/^>/{n+=gsub(/[ACGTacgt]/,"")}END{print n}' "$REF_FA")   # non-N bases

# 1. download (SRA normalized format, includes original base qualities)
[ -d "$WORK/sra/$SRR" ] || prefetch "$SRR" -O "$WORK/sra" --max-size 20G
# 2+3. stream SRA -> adapter trimming -> bwa mem -> coordinate sort, with no FASTQ on disk
# (raw FASTQ would be ~20 GB for the largest library). fastp does adapter trimming only (no
# quality filter, min length 30). bwa mem with default parameters, as in the GEO methods.
# fastq-dump -Z streams without temp files (fasterq-dump needs ~FASTQ-sized scratch, even with --stdout)
fastq-dump -Z --split-spot "$WORK/sra/$SRR/$SRR.sra" 2> "$Q/fastqdump.log" \
 | fastp --stdin --interleaved_in --stdout --detect_adapter_for_pe --disable_quality_filtering \
         --length_required 30 -w "$THREADS" -j "$Q/fastp.json" -h /dev/null 2> "$Q/fastp.log" \
 | bwa mem -p -t "$THREADS" -R "@RG\tID:$LIB\tSM:$LIB\tPL:ILLUMINA" "$GENOME" - 2> "$Q/bwa.log" \
 | samtools sort -@ "$THREADS" -T "$D/sorttmp" -o "$D/aligned.bam" -
samtools index "$D/aligned.bam"
samtools flagstat -@ "$THREADS" "$D/aligned.bam" > "$Q/flagstat_aligned.txt"

# 4. duplicates removed with picard (as in GEO methods), before any MAPQ filter
picard MarkDuplicates -I "$D/aligned.bam" -O "$D/dedup.bam" -M "$Q/markdup_metrics.txt" \
    --REMOVE_DUPLICATES true --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true \
    2> "$Q/markdup.log"
samtools flagstat -@ "$THREADS" "$D/dedup.bam" > "$Q/flagstat_dedup.txt"
rm -f "$D/aligned.bam" "$D/aligned.bam.bai"
rm -rf "$WORK/sra/$SRR"   # free disk; re-downloadable with prefetch

# 5. MAPQ-filtered BAMs + RPGC bigWigs (bin 100 bp, effective genome size = non-N bases)
for MQ in $MAPQS; do
  samtools view -b -q "$MQ" -@ "$THREADS" "$D/dedup.bam" -o "$D/mapq$MQ.bam"; samtools index "$D/mapq$MQ.bam"
  bamCoverage -b "$D/mapq$MQ.bam" -o "$D/$LIB.mapq$MQ.bin${BINSIZE}bp.RPGC.bw" \
     --binSize "$BINSIZE" --normalizeUsing RPGC --effectiveGenomeSize "$EGS" -p "$THREADS" \
     > "$Q/bamCoverage_mapq$MQ.log" 2>&1
  samtools view -c "$D/mapq$MQ.bam" > "$Q/n_reads_mapq$MQ.txt"
  [ "$MQ" = "0" ] || rm -f "$D/mapq$MQ.bam" "$D/mapq$MQ.bam.bai"
done

# 6. library-level QC
picard CollectInsertSizeMetrics -I "$D/dedup.bam" -O "$Q/insert_size_metrics.txt" -H "$Q/insert_size_hist.pdf" \
    --VALIDATION_STRINGENCY SILENT 2> "$Q/insert.log"
picard CollectGcBiasMetrics -I "$D/dedup.bam" -R "$GENOME" -O "$Q/gc_bias_metrics.txt" -CHART "$Q/gc_bias.pdf" \
    -S "$Q/gc_bias_summary.txt" --VALIDATION_STRINGENCY SILENT 2> "$Q/gcbias.log"
samtools idxstats "$D/dedup.bam" > "$Q/idxstats_dedup.txt"
echo "$LIB done; EGS=$EGS"

# 7. keep only bigWigs + QC; BAMs are large and regenerable (set KEEP_BAM=1 to retain)
mkdir -p "$WORK/bigwig"; cp "$D"/*.bw "$WORK/bigwig/"
if [ "${KEEP_BAM:-0}" != "1" ]; then rm -f "$D"/*.bam "$D"/*.bai "$D"/*.bam.bai; fi
