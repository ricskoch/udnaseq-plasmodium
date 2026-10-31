"""Reading/writing genome tracks and interval files (pure Python, no bedtools/deepTools)."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

from . import config


def chrom_sizes() -> dict[str, int]:
    t = pd.read_csv(config.CHROM_SIZES, sep="\t", header=None)
    return dict(zip(t[0], t[1].astype(int)))


def read_bigwig(path, chroms=None) -> dict[str, np.ndarray]:
    """Per-base values of a bigWig as float32 arrays (NaN where no data)."""
    bw = pyBigWig.open(str(path))
    out = {}
    for c, L in bw.chroms().items():
        if chroms is not None and c not in chroms:
            continue
        out[c] = np.asarray(bw.values(c, 0, L, numpy=True), dtype=np.float32)
    bw.close()
    return out


def bin_mean(per_base: np.ndarray, binsize: int) -> np.ndarray:
    """Mean of non-NaN values in consecutive bins (last bin may be partial)."""
    n = len(per_base)
    nb = -(-n // binsize)
    pad = np.full(nb * binsize, np.nan, dtype=np.float64)
    pad[:n] = per_base
    m = pad.reshape(nb, binsize)
    cnt = np.sum(~np.isnan(m), axis=1)
    s = np.nansum(m, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, s / cnt, np.nan)


def read_binned_bigwig(path, binsize=config.BIN, chroms=None) -> dict[str, np.ndarray]:
    return {c: bin_mean(v, binsize) for c, v in read_bigwig(path, chroms).items()}


def bigwig_stats(path, regions: pd.DataFrame, exact=False) -> np.ndarray:
    """Mean signal per region via pyBigWig.stats.

    exact=False reproduces deepTools (multiBigwigSummary/computeMatrix), which
    reads pre-computed zoom-level summaries; these are approximations of the
    true per-base mean.  exact=True returns the true mean over covered bases.
    """
    bw = pyBigWig.open(str(path))
    out = np.full(len(regions), np.nan)
    for i, (c, s, e) in enumerate(zip(regions.chrom, regions.start, regions.end)):
        v = bw.stats(c, int(s), int(e), type="mean", exact=exact)[0]
        out[i] = np.nan if v is None else v
    bw.close()
    return out


def bins_frame(binsize: int, chroms=None, sizes=None) -> pd.DataFrame:
    sizes = sizes or chrom_sizes()
    rows = [(c, s, min(s + binsize, L)) for c, L in sorted(sizes.items())
            if chroms is None or c in chroms for s in range(0, L, binsize)]
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def tracks_to_frame(tracks: dict[str, np.ndarray], binsize: int, sizes=None) -> pd.DataFrame:
    """dict chrom->binned array  ->  bedGraph-like DataFrame (chrom, start, end, value)."""
    sizes = sizes or chrom_sizes()
    parts = []
    for c in sorted(tracks):
        v = tracks[c]
        s = np.arange(len(v)) * binsize
        e = np.minimum(s + binsize, sizes[c])
        parts.append(pd.DataFrame({"chrom": c, "start": s, "end": e, "value": v}))
    return pd.concat(parts, ignore_index=True)


def write_bedgraph(df: pd.DataFrame, path, dropna=True, decimals=6):
    d = df.dropna(subset=["value"]) if dropna else df
    d = d.assign(value=d["value"].round(decimals))
    d[["chrom", "start", "end", "value"]].to_csv(path, sep="\t", header=False, index=False)


def write_bigwig(df: pd.DataFrame, path, sizes=None):
    sizes = sizes or chrom_sizes()
    d = df.dropna(subset=["value"]).sort_values(["chrom", "start"])
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(sorted(sizes.items()))
    for c, g in d.groupby("chrom", sort=True):
        bw.addEntries([c] * len(g), g["start"].astype(int).tolist(),
                      ends=g["end"].astype(int).tolist(), values=g["value"].astype(float).tolist())
    bw.close()


def read_bedgraph(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, comment="#",
                       names=["chrom", "start", "end", "value"], usecols=[0, 1, 2, 3])


def read_bed(path, names=None) -> pd.DataFrame:
    # ragged files (e.g. one extra annotation column on a single line) are tolerated
    d = pd.read_csv(path, sep="\t", header=None, comment="#", names=range(13)).dropna(axis=1, how="all")
    base = ["chrom", "start", "end", "name", "score", "strand"]
    d.columns = (names or base + [f"c{i}" for i in range(6, d.shape[1])])[: d.shape[1]]
    return d


def read_deeptools_matrix(path):
    """Parse a computeMatrix --outFileName (gzipped) matrix -> (params, regions, values)."""
    with gzip.open(path, "rt") as fh:
        params = json.loads(fh.readline()[1:])
        d = pd.read_csv(fh, sep="\t", header=None)
    regions = d.iloc[:, :6]
    regions.columns = ["chrom", "start", "end", "name", "score", "strand"]
    return params, regions, d.iloc[:, 6:].to_numpy(dtype=float)


def load_tracks_npz(name: str):
    """Load a processed track set written by analyses/00_tracks."""
    p = config.TRACKS / f"{name}.npz"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing - run analyses/00_tracks/run.py first")
    z = np.load(p, allow_pickle=False)
    return {k: z[k] for k in z.files}


def geo_path(fname: str) -> Path:
    p = config.GEO_DIR / fname
    if not p.exists():
        raise FileNotFoundError(f"{p} missing - run scripts/fetch_geo.py first")
    return p
