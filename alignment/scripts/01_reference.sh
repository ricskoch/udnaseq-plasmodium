#!/usr/bin/env bash
# Verify the reference FASTA and build the bwa index.
# 1. md5 of the repo FASTA must equal the recorded value.
# 2. Independent check: every sequence must be identical (case-insensitive) to the
#    NCBI GenBank assembly GCA_000002765.3 (P. falciparum 3D7), which PlasmoDB v57 uses.
#    (The PlasmoDB release-57 download URL no longer resolves on the live site.)
# 3. Compute effective genome size (non-N bases) for deepTools RPGC.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$WORK/ref"

echo "$REF_MD5  $REF_FA" | md5sum -c -

NCBI_ZIP="$WORK/ref/GCA_000002765.3.zip"
if [ ! -s "$NCBI_ZIP" ]; then
  curl -sSL -o "$NCBI_ZIP" \
   "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_000002765.3/download?include_annotation_type=GENOME_FASTA"
fi
unzip -o -q "$NCBI_ZIP" -d "$WORK/ref/ncbi"
python3 -I "$(dirname "$0")/verify_reference.py" "$REF_FA" \
  "$(ls "$WORK"/ref/ncbi/ncbi_dataset/data/GCA_000002765.3/*.fna)" \
  | tee "$REPO/alignment/qc/reference_check.txt"

cp "$REF_FA" "$WORK/ref/genome.fa"
bwa index "$WORK/ref/genome.fa" 2> "$WORK/ref/bwa_index.log"
samtools faidx "$WORK/ref/genome.fa"
