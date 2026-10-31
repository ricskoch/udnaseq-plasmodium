# Alignment pipeline: GEO FASTQ to RPGC bigWig

Rebuilds the U-DNA-Seq coverage tracks from the raw reads of GEO GSE262467
(SRA project PRJNA1076816, 8 paired-end 2x150 NovaSeq libraries) and checks them against the
bigWigs in `data/`. BAMs and FASTQs are not stored (several GB each); everything is regenerated
from SRA with the scripts here. Only code and small QC files are committed.

## What it does

| Step | Tool | Notes |
|---|---|---|
| Reference | repo `data/Refgen/PlasmoDB-57_Pfalciparum3D7_Genome.fasta` | md5 checked; all 16 sequences verified identical to NCBI GCA_000002765.3 (`qc/reference_check.txt`). The PlasmoDB release-57 download URL no longer resolves. Effective genome size = 23,332,839 (no N bases). |
| Download | `prefetch` | SRA normalized format with original base qualities. ENA holds no FASTQs for this study. |
| Trim | fastp | Adapter trimming only (auto-detected), no quality filter, min length 30. Streamed, no FASTQ on disk. |
| Align | `bwa mem` default parameters | As in the GEO methods. |
| Duplicates | picard MarkDuplicates, `REMOVE_DUPLICATES=true` | As in the GEO methods. Run before any MAPQ filter. |
| MAPQ filter | `samtools view -q` | bigWigs made at MAPQ >= 0, 1, 10, 30 to find what the published tracks used. |
| Coverage | deepTools `bamCoverage --binSize 100 --normalizeUsing RPGC` | No read extension, effective genome size as above. |
| QC | samtools flagstat, picard insert size and GC bias | Per library under `qc/<library>/`. |

Library list and accessions: `samples.tsv`. Run order: `scripts/run_all.sh` (calls
`01_reference.sh`, then `02_process_library.sh` per library), then
`scripts/03_summarize_qc.py <dir with bigWigs> <repo root>` for the tables below.
Environment: `environment.yml`. Runtime on 4 cores: about 20 min per 10M read pairs, about 6.5 h total.

## Result 1: the published tracks are reproduced

The repo/GEO RPGC tracks match the rebuild best at **MAPQ >= 1** (bwa mem MAPQ 0 means
multi-mapped, so those reads were excluded). Correlation of 100 bp bins (`qc/track_comparison.tsv`):

| Library | Pearson r (MAPQ>=1) | Rank r | Mean signal, rebuilt vs published |
|---|---|---|---|
| 3d7_rep1 INPUT | 0.9987 | 0.9965 | 1.607 vs 1.620 |
| 3d7_rep1 IP | 0.9952 | 0.9915 | 1.585 vs 1.627 |
| 3d7_rep2 INPUT | 0.9985 | 0.9887 | 1.580 vs 1.630 |
| 3d7_rep2 IP | 0.9972 | 0.9838 | 1.532 vs 1.598 |
| DMSO INPUT | 0.9985 | 1.0000 | 1.515 vs 1.516 |
| DMSO IP | 0.9999 | 1.0000 | 1.603 vs 1.603 |
| RAPA INPUT | 0.9989 | 1.0000 | 1.587 vs 1.588 |
| RAPA IP | 0.9999 | 1.0000 | 1.588 vs 1.589 |

DMSO and RAPA are reproduced essentially exactly (rank r = 1.0), so the procedure above is what was
used for them. rep1 and rep2 are close (r >= 0.995) but not exact: the published rep1/rep2 tracks have
3-4 % more total signal per genome than my rebuild, and slightly different rank order. The most
likely cause is a different read set or processing for those two older libraries (for example no
trimming, or a different filtering order). I could not pin it down from the GEO methods text. The
differences are small and do not change any track-level conclusion, but the exact rep1/rep2
commands are unknown. Richard may recall them.

No published track is wrong. Nothing needs replacing.

## Result 2: library QC (`qc/qc_summary.tsv`)

| Library | Read pairs (M) | Mapped to Pf | Duplicates | Median insert (bp) | AT dropout | Normalized coverage at 40-59 % GC |
|---|---|---|---|---|---|---|
| rep1 INPUT | 23.9 | 97.7 % | 11.7 % | 234 | 8.6 | 1.07 |
| rep1 IP | 31.6 | 99.4 % | 13.5 % | 186 | 13.5 | 1.68 |
| rep2 INPUT | 21.1 | 98.7 % | 4.7 % | 198 | 8.8 | 1.31 |
| **rep2 IP** | **9.1** | **41.7 %** | 4.7 % | 185 | 17.7 | 2.58 |
| DMSO INPUT | 10.7 | 80.5 % | 8.4 % | 200 | 16.1 | 3.47 |
| DMSO IP | 19.1 | 94.9 % | 7.8 % | 208 | 25.6 | 2.68 |
| RAPA INPUT | 14.0 | 88.6 % | 7.5 % | 258 | 18.0 | 2.66 |
| RAPA IP | 20.3 | 92.0 % | 7.6 % | 192 | 28.8 | 3.10 |

What stands out:

1. **rep2 IP has only 42 % of reads on the P. falciparum genome**, and its raw read GC is 45 %
   (the genome is about 19 % GC; human is about 41 %). This points to host DNA, but I did not align
   the unmapped reads to the human genome, so that is an inference. After MAPQ filtering rep2 IP has
   7.0M reads against 33-53M for the other IPs, so it is about 5-8 times shallower.
2. **DMSO and RAPA libraries have much stronger GC bias than rep1/rep2**, in INPUT and IP alike:
   AT dropout 16-29 against 8.6-13.5, and normalized coverage of 40-59 % GC windows 2.7-3.5
   against 1.1-1.3 in the rep1 INPUT. This agrees with the sibling-thread finding of a steeper GC
   slope in DMSO/RAPA INPUT. The pattern (GC-rich windows over-represented, AT-rich under-represented)
   is typical of extra PCR cycles or a different library prep on low input, but the PCR cycle number
   is not in the data I have.
3. DMSO INPUT also has host contamination (80 % mapped, read GC 31 %), RAPA INPUT less (89 %).
4. rep1 libraries have the highest duplication (12-13 %), rep2 the lowest (4.7 %).
5. RAPA INPUT median insert is 258 bp against 185-234 bp for the rest, so sonication differed.

## Result 3: IP vs INPUT fingerprint (`qc/fingerprint.png`, `qc/fingerprint.tsv`)

Computed from the MAPQ >= 1 100 bp tracks (cumulative signal over bins ranked low to high; the
maximum gap between the IP and INPUT curves is the usual enrichment measure):

| Pair | Max gap | Signal in top 10 % bins, INPUT | IP |
|---|---|---|---|
| rep1 | 0.039 | 0.161 | 0.173 |
| rep2 | 0.091 | 0.164 | 0.218 |
| DMSO | 0.059 | 0.258 | 0.263 |
| RAPA | 0.068 | 0.253 | 0.286 |

Global IP enrichment over INPUT is weak in every pair (a strong ChIP gives a gap above about 0.3).
rep2 is the most enriched pair and rep1 the least. In DMSO and RAPA the top-10 % signal of the INPUT
is already about 0.25 against 0.16 in rep1/rep2 INPUT, which fits the GC bias above: part of what
looks like structure in those tracks is GC-driven coverage unevenness, present in INPUT too.

## Implications for the open G4 question

These are observations from the QC, not conclusions about the manuscript:

- The weak G4 signal in DMSO/RAPA cannot be blamed on a processing difference, because the tracks
  reproduce exactly. It sits in the libraries.
- DMSO/RAPA INPUT and IP carry 2-3 times the GC bias of rep1/rep2. Any GC-correlated feature
  (G4 motifs are GC-rich) is affected more in those libraries. The GC-matched null from the sibling
  thread is the right control.
- rep2 IP is shallow and contaminated with host reads, so its signal-to-noise per bin is lower than
  its enrichment curve suggests.
- Before dropping or keeping any claim, the next checks would be: IP/INPUT ratio after explicit GC
  correction (deepTools `correctGCBias` or a GC-binned ratio), and a down-sampling of rep1 IP to
  rep2 IP depth to see how much of the rep1/rep2 difference is depth.
