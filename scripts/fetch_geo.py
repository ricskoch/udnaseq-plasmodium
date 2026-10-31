#!/usr/bin/env python3
"""Download the processed GEO GSE262467 tracks (~650 MB) into data/geo/ and verify MD5.

  * GSE262467_3d7_merged.MAPQfiltered.IPvsINPUT.bin10bp.log2.bw - merged 3D7, 10 bp bins
  * GSM81701{37,39,41,43}_*.IP.vs.INPUT.bin1bp.log2.bw           - per-sample, 1 bp bins
    (shipped inside GSE262467_RAW.tar)

These are needed for the 1 bp gene-level means and the 10 bp metaprofile/G4
analyses; the 100 bp RPGC coverage tracks in data/ cannot provide that resolution.
"""
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from udnaseq import config  # noqa: E402


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(name: str, dest: Path):
    print(f"downloading {name} ...", flush=True)
    urllib.request.urlretrieve(config.GEO_BASE + name, dest)


def main():
    config.GEO_DIR.mkdir(parents=True, exist_ok=True)
    missing = [f for f in config.GEO_MD5 if not (config.GEO_DIR / f).exists()]
    if config.GEO_MERGED_10BP in missing:
        download(config.GEO_MERGED_10BP, config.GEO_DIR / config.GEO_MERGED_10BP)
    if any(f != config.GEO_MERGED_10BP for f in missing):
        tar = config.GEO_DIR / "GSE262467_RAW.tar"
        download(tar.name, tar)
        with tarfile.open(tar) as t:
            for m in t.getmembers():
                if m.name in config.GEO_MD5 and m.isfile():
                    t.extract(m, config.GEO_DIR, filter="data")
        tar.unlink()
    bad = [f for f, h in config.GEO_MD5.items() if md5(config.GEO_DIR / f) != h]
    if bad:
        sys.exit(f"MD5 mismatch: {bad}")
    print("GEO files present and verified:", *config.GEO_MD5, sep="\n  ")


if __name__ == "__main__":
    main()
