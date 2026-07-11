"""
GSE137398 ONC RGC - RGC Subtype Vulnerability Analysis
Tracks subtype proportion changes across time points post-ONC.
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import chi2_contingency

warnings.filterwarnings("ignore")

CACHE = Path(r"F:\ONC\output\preprocessed.h5ad")
RAW_CACHE = Path(r"F:\ONC\output\qc_filtered.h5ad")
OUTPUT = Path(r"F:\ONC\output\subtype")
OUTPUT.mkdir(parents=True, exist_ok=True)

TIME_ORDER = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]


def main():
    print("Loading preprocessed AnnData...")
    adata = ad.read_h5ad(CACHE)
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"  Subtypes: {adata.obs['rna_subtype'].nunique()}")

    # --- Subtype composition per time point ---
    print("\n=== Subtype Composition ===")
    comp = pd.crosstab(adata.obs["rna_subtype"], adata.obs["time_point"])
    comp_pct = comp.div(comp.sum(axis=0), axis=1) * 100
    comp_pct.to_csv(OUTPUT / "subtype_composition_pct.csv")
    print(f"  Saved: subtype_composition_pct.csv")

    # --- Top subtypes by abundance ---
    top_subtypes = comp.sum(axis=1).nlargest(15).index
    print(f"\n  Top 15 subtypes (total cells):")
    for st in top_subtypes:
        total = comp.loc[st].sum()
        pcts = ", ".join(f"{tp}={comp_pct.loc[st, tp]:.1f}%" for tp in TIME_ORDER)
        print(f"    {st:20s}: {total:5d} cells | {pcts}")

    # --- Vulnerability scores ---
    print("\n=== Subtype Vulnerability ===")
    ctrl_pct = comp_pct["Ctrl"]
    late_pct = comp_pct["2w"]
    # Vulnerability: fold change from Ctrl to 2w (negative = decline = vulnerable)
    vuln = pd.DataFrame({
        "ctrl_pct": ctrl_pct,
        "late_pct": late_pct,
        "fold_change": (late_pct + 0.01) / (ctrl_pct + 0.01),
        "total_cells": comp.sum(axis=1),
    })
    vuln = vuln[vuln["total_cells"] >= 10].copy()
    vuln["vulnerability"] = vuln["fold_change"].apply(
        lambda x: "vulnerable" if x < 0.5 else ("resilient" if x > 1.5 else "stable")
    )

    n_vuln = (vuln["vulnerability"] == "vulnerable").sum()
    n_res = (vuln["vulnerability"] == "resilient").sum()
    n_stable = (vuln["vulnerability"] == "stable").sum()
    print(f"  Vulnerable: {n_vuln}, Resilient: {n_res}, Stable: {n_stable}")

    print("\n  Most vulnerable (declining at 2w):")
    vuln_sorted = vuln.sort_values("fold_change")
    for st in vuln_sorted.head(8).index:
        row = vuln_sorted.loc[st]
        print(f"    {st:20s}: Ctrl={row['ctrl_pct']:.1f}% -> 2w={row['late_pct']:.1f}% (FC={row['fold_change']:.2f})")

    print("\n  Most resilient (increasing at 2w):")
    res_sorted = vuln.sort_values("fold_change", ascending=False)
    for st in res_sorted.head(8).index:
        row = res_sorted.loc[st]
        print(f"    {st:20s}: Ctrl={row['ctrl_pct']:.1f}% -> 2w={row['late_pct']:.1f}% (FC={row['fold_change']:.2f})")

    vuln.to_csv(OUTPUT / "subtype_vulnerability.csv")
    print(f"\n  Saved: subtype_vulnerability.csv")

    # --- Full time-course proportion trajectory ---
    print("\n=== Proportion Trajectories ===")
    trajectory = comp_pct.loc[vuln_sorted.head(6).index.tolist() + res_sorted.head(6).index.tolist()]
    trajectory.to_csv(OUTPUT / "subtype_trajectories.csv")
    print(f"  Saved: subtype_trajectories.csv")

    # --- Chi-squared test for composition shift ---
    print("\n=== Statistical Test ===")
    # Test if Ctrl vs 2w composition is significantly different
    comp_ctrl_2w = comp[["Ctrl", "2w"]]
    # Remove rare subtypes
    comp_filtered = comp_ctrl_2w[comp_ctrl_2w.sum(axis=1) >= 20]
    chi2, p, dof, expected = chi2_contingency(comp_filtered)
    print(f"  Ctrl vs 2w composition: chi2={chi2:.1f}, p={p:.2e}, dof={dof}")
    if p < 0.001:
        print("  *** Highly significant composition shift ***")

    print("\nSubtype analysis complete.")


if __name__ == "__main__":
    main()
