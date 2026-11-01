#!/usr/bin/env python3
"""Download public 3D7 WT H3K9me3 ChIP-seq and INPUT bigWigs (GEO GSE157709, trophozoite 30-40 hpi, rep2) to data/h3k9me3/.
GSM4773592 (rep1 IP) returned HTTP 500 from GEO when this was written, so only rep2 is used."""
import sys, urllib.request
from pathlib import Path
D = Path(__file__).resolve().parents[1] / "data" / "h3k9me3"; D.mkdir(parents=True, exist_ok=True)
FILES = ["GSM5353nnn/GSM5353855/suppl/GSM5353855_WT_H3K9me3_T_rep2_HM_ChIP-seq.bw",
         "GSM5353nnn/GSM5353856/suppl/GSM5353856_WT_input_T_rep2_HM_ChIP-seq.bw"]
for f in FILES:
    out = D / Path(f).name
    if out.exists() and out.stat().st_size > 1e6:
        continue
    req = urllib.request.Request("https://ftp.ncbi.nlm.nih.gov/geo/samples/" + f, headers={"User-Agent": "curl/8"})
    out.write_bytes(urllib.request.urlopen(req, timeout=300).read()); print("fetched", out.name, file=sys.stderr)
