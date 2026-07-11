"""
GSE137398 ONC RGC - Pseudobulk DEG analysis.
Uses pseudobulk matrix (21 samples x 40790 genes) for sample-level stats.
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

PB_PATH = Path(r"F:\ONC\output\deg\pseudobulk.h5ad")
OUTPUT = Path(r"F:\ONC\output\deg")
TIME_ORDER = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]


def normalize_cpm(adata):
    """Normalize to counts per million."""
    lib_sizes = adata.X.sum(axis=1)
    cpm = adata.X / lib_sizes[:, np.newaxis] * 1e6
    adata.X = cpm
    return adata


def run_pb_deg(adata, case_tp, ctrl_tp="Ctrl"):
    """Run t-test DEG on pseudobulk data."""
    key = f"pb_{case_tp}_vs_{ctrl_tp}"
    print(f"\n=== {key} ===")

    ctrl_mask = adata.obs["time_point"] == ctrl_tp
    case_mask = adata.obs["time_point"] == case_tp
    ctrl_samples = adata.obs.loc[ctrl_mask, "sample_id"].tolist()
    case_samples = adata.obs.loc[case_mask, "sample_id"].tolist()

    ctrl_data = adata[ctrl_mask].X.toarray() if hasattr(adata[ctrl_mask].X, "toarray") else adata[ctrl_mask].X
    case_data = adata[case_mask].X.toarray() if hasattr(adata[case_mask].X, "toarray") else adata[case_mask].X

    n_ctrl, n_case = ctrl_data.shape[0], case_data.shape[0]
    print(f"  {ctrl_tp}: {n_ctrl} samples ({', '.join(ctrl_samples)})")
    print(f"  {case_tp}: {n_case} samples ({', '.join(case_samples)})")

    n_genes = adata.n_vars
    log2fc = np.zeros(n_genes)
    pvals = np.ones(n_genes)
    t_stats = np.zeros(n_genes)

    for g in range(n_genes):
        ctrl_vec = ctrl_data[:, g]
        case_vec = case_data[:, g]
        # log2 fold change
        ctrl_mean = ctrl_vec.mean()
        case_mean = case_vec.mean()
        log2fc[g] = np.log2((case_mean + 1) / (ctrl_mean + 1))
        if ctrl_vec.std() == 0 and case_vec.std() == 0:
            pvals[g] = 1.0
        else:
            t_stat, p_val = ttest_ind(case_vec, ctrl_vec, equal_var=False)
            t_stats[g] = t_stat
            pvals[g] = p_val

    # BH correction
    _, pvals_adj, _, _ = multipletests(pvals, method="fdr_bh")

    result = pd.DataFrame({
        "gene": adata.var_names,
        "log2FC": log2fc,
        "t_stat": t_stats,
        "pval": pvals,
        "pval_adj": pvals_adj,
    }).set_index("gene")

    n_sig = (result["pval_adj"] < 0.05).sum()
    n_up = ((result["pval_adj"] < 0.05) & (result["log2FC"] > 0)).sum()
    n_down = ((result["pval_adj"] < 0.05) & (result["log2FC"] < 0)).sum()
    print(f"  DEGs (FDR<0.05): {n_sig} total, {n_up} up, {n_down} down")

    top_up = result[result["pval_adj"] < 0.05].nlargest(10, "log2FC").index.tolist()
    top_down = result[result["pval_adj"] < 0.05].nsmallest(10, "log2FC").index.tolist()
    print(f"  Top up: {', '.join(top_up[:8])}")
    print(f"  Top down: {', '.join(top_down[:8])}")

    out_path = OUTPUT / f"{key}.csv"
    result.sort_values("pval_adj").to_csv(out_path)
    print(f"  Saved: {out_path}")
    return result


def main():
    print("Loading pseudobulk matrix...")
    adata = ad.read_h5ad(PB_PATH)
    print(f"  {adata.n_obs} samples x {adata.n_vars} genes")

    # Per-sample cell counts
    print(f"  Samples per time point: {dict(adata.obs['time_point'].value_counts())}")

    # Normalize to CPM
    adata = normalize_cpm(adata)

    # Run all comparisons
    all_results = {}
    for tp in TIME_ORDER[1:]:
        result = run_pb_deg(adata, tp)
        all_results[tp] = result

    # Summary
    print("\n=== Pseudobulk DEG Summary ===")
    summary_rows = []
    for tp in TIME_ORDER[1:]:
        r = all_results[tp]
        n_sig = (r["pval_adj"] < 0.05).sum()
        n_up = ((r["pval_adj"] < 0.05) & (r["log2FC"] > 0)).sum()
        n_down = ((r["pval_adj"] < 0.05) & (r["log2FC"] < 0)).sum()
        top3_up = r[r["pval_adj"] < 0.05].nlargest(3, "log2FC").index.tolist()
        top3_down = r[r["pval_adj"] < 0.05].nsmallest(3, "log2FC").index.tolist()
        summary_rows.append({
            "time_point": tp, "n_sig": n_sig, "n_up": n_up, "n_down": n_down,
            "top_up": ", ".join(top3_up), "top_down": ", ".join(top3_down),
        })
        print(f"  {tp}: {n_sig} sig ({n_up} up, {n_down} down) | up: {', '.join(top3_up[:3])} | down: {', '.join(top3_down[:3])}")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT / "pb_deg_summary.csv", index=False)
    print(f"\n  Summary saved: pb_deg_summary.csv")
    print("\nPseudobulk DEG complete.")


if __name__ == "__main__":
    main()
