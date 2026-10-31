"""Genomic annotation helpers: genes, blacklist masks, interval statistics, closest-feature distances."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import config


# ---------------------------------------------------------------- genes -----
def protein_coding_genes() -> pd.DataFrame:
    """All protein_coding_gene features of the PlasmoDB-57 GFF (n = 5318).

    Columns: chrom, start (0-based), end, id, symbol, description, strand, category.
    category: VAR / RIF / STEVOR / Others, assigned from the GFF annotation.
    """
    g = pd.read_csv(config.GFF, sep="\t", comment="#", header=None)
    g = g[g[2] == "protein_coding_gene"].copy()
    attr = g[8]
    get = lambda key: attr.str.extract(rf"(?:^|;){key}=([^;]*)")[0]
    out = pd.DataFrame({
        "chrom": g[0].values, "start": (g[3] - 1).values, "end": g[4].values,
        "id": get("ID").values, "symbol": get("Name").values,
        "description": get("description").str.replace("%2C", ",").values,
        "strand": g[6].values,
    })
    out["category"] = out.apply(_category, axis=1)
    # Legacy definition used for the published Fig. 2C: gene symbol only.
    # It leaves var2csa (symbol VAR2CSA) and 3 PfEMP1 exon-2 fragments in "Others".
    sym = out.symbol.fillna("").str.upper()
    out["category_legacy"] = np.where(sym.isin(["VAR", "RIF", "STEVOR"]), sym, "Others")
    return out.reset_index(drop=True)


def _category(r) -> str:
    sym = str(r.symbol).upper() if pd.notna(r.symbol) else ""
    desc = str(r.description).lower() if pd.notna(r.description) else ""
    if sym == "VAR" or "pfemp1" in desc:
        return "VAR"
    if sym == "RIF" or re.search(r"\brifin\b", desc):
        return "RIF"
    if sym == "STEVOR" or "stevor" in desc:
        return "STEVOR"
    return "Others"


# ------------------------------------------------------------ blacklist -----
def blacklist() -> pd.DataFrame:
    return pd.read_csv(config.BLACKLIST, sep="\t", header=None, usecols=[0, 1, 2],
                       names=["chrom", "start", "end"])


def blacklist_mask(chrom: str, starts: np.ndarray, ends: np.ndarray, bl=None) -> np.ndarray:
    """True for intervals overlapping any blacklisted region (bedtools subtract -A semantics)."""
    bl = blacklist() if bl is None else bl
    m = np.zeros(len(starts), bool)
    for _, r in bl[bl.chrom == chrom].iterrows():
        m |= (ends > r.start) & (starts < r.end)
    return m


def add_blacklist_flag(df: pd.DataFrame) -> pd.DataFrame:
    bl = blacklist()
    flag = np.zeros(len(df), bool)
    for c, idx in df.groupby("chrom").groups.items():
        sub = df.loc[idx]
        flag[df.index.get_indexer(idx)] = blacklist_mask(c, sub.start.values, sub.end.values, bl)
    return df.assign(blacklisted=flag)


# --------------------------------------------------- interval statistics -----
def interval_means(per_base: dict[str, np.ndarray], regions: pd.DataFrame,
                   start_col="start", empty_as_zero=False) -> np.ndarray:
    """Mean of per-base values over [start, end) for each region.

    empty_as_zero=False : mean over covered bases only (NaN if none)  - pyBigWig stats 'mean'
    empty_as_zero=True  : sum / region length                        - bigWigAverageOverBed 'mean0'
    """
    out = np.full(len(regions), np.nan)
    for i, (c, s, e) in enumerate(zip(regions.chrom, regions[start_col], regions.end)):
        v = per_base[c][int(s):int(e)]
        ok = ~np.isnan(v)
        if empty_as_zero:
            out[i] = v[ok].sum() / max(len(v), 1)
        elif ok.any():
            out[i] = v[ok].mean()
    return out


def merge_intervals(df: pd.DataFrame, value_col="value") -> pd.DataFrame:
    """bedtools merge -c value -o mean (book-ended intervals are merged)."""
    rows = []
    for c, g in df.sort_values(["chrom", "start"]).groupby("chrom", sort=True):
        s, e, v = g.start.values, g.end.values, g[value_col].values
        brk = np.r_[True, s[1:] > np.maximum.accumulate(e)[:-1]]
        grp = np.cumsum(brk)
        for k in np.unique(grp):
            m = grp == k
            rows.append((c, s[m].min(), e[m].max(), v[m].mean()))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", value_col])


def closest_signed(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """bedtools closest -D ref -t all.

    For every interval in `a`, report the closest interval(s) in `b` and the
    signed distance relative to the reference genome (negative = b upstream
    of a, 0 = overlap).  Ties are all reported, as bedtools does by default.
    Distance follows bedtools convention for half-open coordinates:
    b upstream:  -(a.start - b.end + 1);  b downstream: b.start - a.end + 1.
    """
    rows = []
    for c, ga in a.groupby("chrom", sort=False):
        gb = b[b.chrom == c].reset_index(drop=True)
        if gb.empty:
            continue
        bs, be = gb.start.values, gb.end.values
        for ra in ga.itertuples(index=False):
            up = ra.start - be + 1          # >0 when b ends before a starts
            down = bs - ra.end + 1          # >0 when b starts after a ends
            overlap = (bs < ra.end) & (be > ra.start)
            d = np.where(overlap, 0, np.where(be <= ra.start, -up, down))
            best = np.abs(d).min()
            for j in np.where(np.abs(d) == best)[0]:
                rows.append((*ra, *gb.iloc[j].tolist(), int(d[j])))
    cols = [f"{x}_a" for x in a.columns] + [f"{x}_b" for x in b.columns] + ["distance"]
    return pd.DataFrame(rows, columns=cols)
