#!/usr/bin/env python3
"""Fetch GO annotations for P. falciparum 3D7 (UniProt proteome UP000001450, ORF IDs PF3D7_*) and the GO ontology.

The ShinyGO analysis of the manuscript used an online annotation set; PlasmoDB/GO GAF downloads are not reachable
from every network, so UniProt (direct annotations; ancestors are added in the analysis from go-basic.obo) is used.
Writes data/go/pf3d7_uniprot_go.tsv (committed, small) and data/go/go-basic.obo (git-ignored, 30 MB) and prints
the retrieval date and the go-basic data-version.
"""
import datetime
import sys
import urllib.request
from pathlib import Path

D = Path(__file__).resolve().parents[1] / "data" / "go"
URLS = {"pf3d7_uniprot_go.tsv": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001450&format=tsv&fields=accession,gene_orf,go_id",
        "go-basic.obo": "https://purl.obolibrary.org/obo/go/go-basic.obo"}


def main():
    D.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        p = D / name
        if p.exists() and "--force" not in sys.argv:
            continue
        print("downloading", name, flush=True)
        req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req) as r, open(p, "wb") as fh:
            fh.write(r.read())
    ver = next((l.strip() for l in open(D / "go-basic.obo") if l.startswith("data-version")), "?")
    print(f"retrieved {datetime.date.today()}; go-basic {ver}; annotation rows: {sum(1 for _ in open(D / 'pf3d7_uniprot_go.tsv')) - 1}")


if __name__ == "__main__":
    main()
