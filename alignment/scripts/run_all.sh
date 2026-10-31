#!/usr/bin/env bash
# Run every library in samples.tsv sequentially (one at a time: disk and RAM friendly).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
bash "$(dirname "$0")/01_reference.sh"
tail -n +2 "$REPO/alignment/samples.tsv" | cut -f1 | while read -r LIB; do
  [ -s "$WORK/bigwig/$LIB.mapq1.bin${BINSIZE}bp.RPGC.bw" ] && { echo "skip $LIB"; continue; }
  bash "$(dirname "$0")/02_process_library.sh" "$LIB"
done
