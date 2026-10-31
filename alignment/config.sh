# Shared settings, sourced by all scripts. Override by exporting before running.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${WORK:-/tmp/udna_work}"        # large intermediates (SRA, FASTQ, BAM) - never committed
REF_FA="$REPO/data/Refgen/PlasmoDB-57_Pfalciparum3D7_Genome.fasta"
REF_MD5="e98490036aa0a04b60b7972b3f1f95b9"   # md5 of the repo FASTA (PlasmoDB v57 / GCA_000002765.3)
THREADS="${THREADS:-4}"
BINSIZE=100
# MAPQ thresholds for which bigWigs are produced (0 = no MAPQ filter)
MAPQS="${MAPQS:-0 1 10 30}"
