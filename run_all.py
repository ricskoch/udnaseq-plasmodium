#!/usr/bin/env python3
"""Run the whole reanalysis: fetch GEO tracks, build tracks, then every figure analysis.

    python run_all.py              # everything (~5 min after the one-off ~650 MB GEO download)
    python run_all.py fig4C        # only analyses whose folder name contains 'fig4C'

Each analysis writes output/ and COMPARISON.md in its own folder; this script
collects the status column of every COMPARISON.md into REPRODUCIBILITY_SUMMARY.md.
"""
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ORDER = [
    "00_tracks",
    "00b_bin_approximation",
    "00c_tracks_10bp",
    "figS1AB_sample_correlation",
    "qc_pca",
    "qc_library_bias",
    "fig2B_S1C_genome_distribution",
    "fig2C_3A_S1D_S2_gene_enrichment",
    "fig3B_vsa_metaprofiles",
    "fig3B_per_dataset",
    "fig3B_gc_profile",
    "gene_ranking_go",
    "fig4A_locus_track",
    "fig4B_gc_correlation",
    "fig4B_per_dataset",
    "fig4C_g4_distance",
    "fig4C_per_dataset",
    "fig4C_gc_matched_null",
    "gc_corrected",
    "replicate_consistency",
]


def run(cmd):
    t = time.time()
    print(f"\n=== {' '.join(map(str, cmd))}", flush=True)
    r = subprocess.run([sys.executable, *map(str, cmd)], cwd=ROOT, capture_output=True, text=True)
    if r.returncode:
        print(r.stdout[-3000:], r.stderr[-3000:])
        sys.exit(f"FAILED: {cmd}")
    print(f"    done in {time.time() - t:.0f}s")


def summary():
    lines = ["# Reproducibility summary", "",
             "Status of every legacy/manuscript check (details in each analysis' COMPARISON.md).", "",
             "| Analysis | MATCH | CLOSE | DIFFERS | UNRESOLVED | INFO |", "|---|---|---|---|---|---|"]
    for name in ORDER:
        f = ROOT / "analyses" / name / "COMPARISON.md"
        if not f.exists():
            continue
        st = re.findall(r"\| (MATCH|CLOSE|DIFFERS|UNRESOLVED|INFO) \|", f.read_text())
        c = {k: st.count(k) for k in ("MATCH", "CLOSE", "DIFFERS", "UNRESOLVED", "INFO")}
        lines.append(f"| [{name}](analyses/{name}/COMPARISON.md) | " + " | ".join(str(v) for v in c.values()) + " |")
    (ROOT / "REPRODUCIBILITY_SUMMARY.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    pick = sys.argv[1:]
    run(["scripts/fetch_geo.py"])
    run(["scripts/fetch_go.py"])
    for name in ORDER:
        if not pick or any(p in name for p in pick) or (name == "00_tracks" and not (ROOT / "analyses/00_tracks/output/tracks_100bp.npz").exists()):
            run([f"analyses/{name}/run.py"])
    summary()


if __name__ == "__main__":
    main()
