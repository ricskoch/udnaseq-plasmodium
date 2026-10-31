"""Central configuration: paths, sample sheet and analysis constants.

Every analysis imports its inputs from here, so a path or threshold only has
to be changed in one place.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
ANALYSES = ROOT / "analyses"
LEGACY = ROOT / "legacy"

# --- Reference genome & annotation (PlasmoDB release 57) -------------------
FASTA = DATA / "Refgen" / "PlasmoDB-57_Pfalciparum3D7_Genome.fasta"
GFF = DATA / "Refgen" / "PlasmoDB-57_Pfalciparum3D7.gff"
CHROM_SIZES = DATA / "Refgen" / "genome.size.chromosomes.tab"

# --- Region annotations ----------------------------------------------------
BED = DATA / "BED_files"
BLACKLIST = BED / "blacklist_telomer_pf11.bed"
CENTROMERES = BED / "CenH3_centromeres.bed"
G4_SITES = BED / "G4_coordinates.bed"  # G4Hunter predictions, Gazanion et al. 2020
VAR_BED = BED / "var_genes" / "var_genes.bed"
VAR_GROUPS_BED = BED / "var_genes" / "vargenes_groups.bed"
RIFSTEVOR_BED = BED / "rifinstrevor_genes.bed"

# --- Raw-most data in the repo: RPGC-normalised coverage, 100 bp bins -------
# (output of bwa mem -> samtools -> picard MarkDuplicates -> bamCoverage RPGC)
SAMPLES = {
    "3d7_rep1": {"ip": DATA / "3d7_rep1.Plasmodium.UDNAIPcorrect.bin100bp.RPGC.bw",
                 "input": DATA / "3d7_rep1.Plasmodium.INPUTcorrect.bin100bp.RPGC.bw",
                 "geo_1bp": "GSM8170137_3d7_rep1.IP.vs.INPUT.bin1bp.log2.bw"},
    "3d7_rep2": {"ip": DATA / "3d7_rep2.Plasmodium.UDNAIP.bin100bp.RPGC.bw",
                 "input": DATA / "3d7_rep2.Plasmodium.INPUT.bin100bp.RPGC.bw",
                 "geo_1bp": "GSM8170139_3d7_rep2.IP.vs.INPUT.bin1bp.log2.bw"},
    "DMSO": {"ip": DATA / "DMSO.UDNAIP.bin100bp.RPGC.bw",
             "input": DATA / "DMSO.INPUT.bin100bp.RPGC.bw",
             "geo_1bp": "GSM8170141_DMSO.IP.vs.INPUT.bin1bp.log2.bw"},
    "RAPA": {"ip": DATA / "RAPA.UDNAIP.bin100bp.RPGC.bw",
             "input": DATA / "RAPA.INPUT.bin100bp.RPGC.bw",
             "geo_1bp": "GSM8170143_RAPA.IP.vs.INPUT.bin1bp.log2.bw"},
}
SAMPLE_ORDER = ["3d7_rep1", "3d7_rep2", "DMSO", "RAPA"]

# --- GEO GSE262467 processed tracks (downloaded by scripts/fetch_geo.py) ----
GEO_DIR = DATA / "geo"
GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262467/suppl/"
GEO_MERGED_10BP = "GSE262467_3d7_merged.MAPQfiltered.IPvsINPUT.bin10bp.log2.bw"
GEO_MD5 = {
    GEO_MERGED_10BP: "f337b7634d060af6814c8f34ea3aa673",
    "GSM8170137_3d7_rep1.IP.vs.INPUT.bin1bp.log2.bw": "3474f743607c8ddc375919c32eb45536",
    "GSM8170139_3d7_rep2.IP.vs.INPUT.bin1bp.log2.bw": "43323be76b637a7195f2f80831ee04e8",
    "GSM8170141_DMSO.IP.vs.INPUT.bin1bp.log2.bw": "8c27a588f44e11e53a4030e1165e740c",
    "GSM8170143_RAPA.IP.vs.INPUT.bin1bp.log2.bw": "2efb27b827ca39971e2f00388401dfc0",
}

# --- Processed tracks produced by analyses/00_tracks ------------------------
TRACKS = ANALYSES / "00_tracks" / "output"

# --- Constants --------------------------------------------------------------
PSEUDOCOUNT = 1.0  # deepTools bigwigCompare default
BIN = 100          # bin size of the RPGC coverage tracks
# Merged 3D7 track = log2 ratio of the RPGC coverage of the *merged BAMs*.
# Merging BAMs before RPGC scaling weights each replicate by its read share.
# BAM read counts are not in the repo, so the shares are estimated in
# analyses/00_tracks (least-squares fit against the legacy merged bedGraph)
# and written to TRACKS/merge_weights.json.
NUCLEAR_CHROMS = [f"Pf3D7_{i:02d}_v3" for i in range(1, 15)]
U_THRESHOLD = 0.5      # log2(U-enrichment) cut-off for "highly uracilated"
G4_MIN_SCORE = 1.2     # |G4Hunter| threshold
G4_MAX_SIZE = 200      # bp

# var gene sub-classes highlighted in Fig. 3A
VAR2CSA = ["PF3D7_1200600"]
# The three short Group A var genes of 3D7 (mini-PfEMP1, ~4-5 kb):
# PFA0015c, MAL6P1.314, PFI1820w.  Please confirm against Zhang et al. 2014.
VAR3 = ["PF3D7_0100300", "PF3D7_0600400", "PF3D7_0937600"]
