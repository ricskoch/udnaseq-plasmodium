#!/usr/bin/env python3
"""Replicate consistency table - does each finding hold in rep1, rep2, DMSO, RAPA (and the merged track)?

Purpose
    One table for the central question of the rework: are the headline findings carried by rep1 only (the merged track
    is ~89 % rep1) or do they reproduce in the independent datasets? Collects effect sizes from the per-dataset analyses.
Method
    Reads the CSVs written by fig2C_3A_S1D_S2_gene_enrichment, fig4B_per_dataset, fig3B_per_dataset, fig4C_per_dataset and
    fig4C_gc_matched_null; nothing is recomputed. Run after those analyses.
Outputs  output/consistency_table.csv, consistency_table.md (also embedded in COMPARISON.md)
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
A = HERE.parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import pandas as pd  # noqa: E402

from udnaseq import config  # noqa: E402

SETS = config.SAMPLE_ORDER + ["3d7_merged"]


def main():
    (HERE / "output").mkdir(exist_ok=True)
    rows = []

    def add(finding, metric, vals, rule):
        rows.append(dict(finding=finding, metric=metric, **{s: vals.get(s, "") for s in SETS}, holds_in=rule))

    # --- Fig 2C: VSA families vs other genes (mean log2U difference, Welch BH p)
    d = pd.read_csv(A / "fig2C_3A_S1D_S2_gene_enrichment/output/stats_fig2C_VSA_vs_others_legacy.csv")
    for fam in ("VAR", "RIF", "STEVOR"):
        x = d[d.family == fam].set_index("dataset")
        add(f"Fig 2C: {fam} genes > other genes", "mean log2U difference (BH p max)",
            {s: f"{x.mean_diff[s]:+.3f}" for s in x.index if s in SETS},
            f"{int(((x.mean_diff > 0) & (x.p_BH < .001)).sum())}/{len(x)} datasets (BH p<0.001)")
    # --- Fig 4B
    d = pd.read_csv(A / "fig4B_per_dataset/output/per_dataset_stats.csv").set_index("dataset")
    add("Fig 4B: GC% vs log2U", "Pearson r (100 bp bins)", {s: f"{d.pearson[s]:.2f}" for s in d.index},
        f"{int((d.pearson > 0.3).sum())}/{len(d)} datasets (r > 0.3)")
    # --- Fig 3B
    d = pd.read_csv(A / "fig3B_per_dataset/output/per_dataset_effects.csv")
    for nm, lab in (("var: 0-6 kb minus 6-10 kb", "Fig 3B: var enriched in first 6 kb vs 6-10 kb"),
                    ("rif/stevor: body 30-90 % minus flanks", "Fig 3B: rif/stevor body > flanks"),
                    ("rif/stevor: body 30-90 % minus first 30 %", "Fig 3B: rif/stevor exon-2 side > first 30 % of body")):
        x = d[d.name == nm].set_index("dataset")
        ok = ((x.diff_ci_lo > 0) & (x.wilcoxon_p < .05)).sum()
        add(lab, "mean log2U difference (95 % CI over genes)", {s: f"{x['diff'][s]:+.2f} ({x.diff_ci_lo[s]:+.2f}, {x.diff_ci_hi[s]:+.2f})" for s in x.index},
            f"{int(ok)}/{len(x)} datasets")
    x = d[d.name == "rif/stevor: body 30-90 % minus flanks"].set_index("dataset")
    add("Fig 3B: rif/stevor profile peak position", "peak of mean profile, % of scaled body", {s: f"{x.peak_pos_frac_of_body[s]:.0%}" for s in x.index},
        "position not stable (51-53 % merged/rep1/rep2, 72-83 % RAPA/DMSO)")
    # --- Fig 4C
    d = pd.read_csv(A / "fig4C_per_dataset/output/per_dataset_summary.csv")
    x = d[d.definition == "matched"].set_index("dataset")
    add("Fig 4C: hotspots within 100 bp of strong G4 (random null)", "% obs / % null (fold), N = 4,034 hotspots",
        {s: f"{x.obs_frac_within_100bp[s]:.1%} / {x.null_frac_within_100bp[s]:.1%} ({x.fold_frac_within_100bp[s]:.1f}x)" for s in x.index},
        f"{int(((x.fold_frac_within_100bp > 1) & (x.p_frac_within_100bp < .05)).sum())}/{len(x)} datasets (DMSO: depleted)")
    d = pd.read_csv(A / "fig4C_gc_matched_null/output/matched_null_summary.csv")
    for m, lab in (("gc10", "10 bp GC"), ("gc10+500", "10 bp GC + 500 bp-window GC")):
        x = d[d.matching == m].set_index("dataset")
        add(f"Fig 4C: same, GC-matched null ({lab})", "% obs / % null (fold)",
            {s: f"{x.obs_frac_within_100bp[s]:.1%} / {x.null_frac_within_100bp[s]:.1%} ({x.fold_frac_within_100bp[s]:.1f}x)" for s in x.index},
            f"{int(((x.fold_frac_within_100bp > 1) & (x.p_frac_within_100bp < .05)).sum())}/{len(x)} datasets")
    mdl = pd.read_csv(A / "fig4C_gc_matched_null/output/model_coefficients.csv")
    g = mdl[mdl.model == "gc+G4"].set_index("dataset")
    add("Fig 4C: G4 effect beyond GC (100 bp OLS, GC deciles)", "log2U coefficient (SE)",
        {s: f"{g.G4_within_100bp_coef[s]:+.3f} ({g.G4_within_100bp_se[s]:.3f})" for s in g.index},
        f"{int((g.G4_within_100bp_coef / g.G4_within_100bp_se > 3).sum())}/{len(g)} datasets with coef > 3 SE")
    T = pd.DataFrame(rows)
    T.to_csv(HERE / "output" / "consistency_table.csv", index=False)
    md = ["| Finding | Metric | " + " | ".join(SETS) + " | Holds in |", "|---|---|" + "---|" * (len(SETS) + 1)]
    md += ["| " + " | ".join([r["finding"], r["metric"]] + [str(r[s]) for s in SETS] + [r["holds_in"]]) + " |"
           for _, r in T.iterrows()]
    (HERE / "output" / "consistency_table.md").write_text("\n".join(md) + "\n")
    (HERE / "COMPARISON.md").write_text(
        "# Replicate consistency table\n\nDoes each finding hold per dataset? Effect sizes collected from the per-dataset analyses "
        "(nothing recomputed here). 3d7_merged is ~89 % rep1 (IP) / 57 % rep1 (INPUT). All datasets are single libraries.\n\n"
        "_Generated by `run.py` on " + __import__("datetime").date.today().isoformat() + ".\n\n" + "\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
