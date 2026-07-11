"""
GSE137398 ONC RGC - Pathway enrichment for DEGs and gene clusters.
Uses gseapy enrichr or fallback to pre-computed GO sets.
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

OUTPUT = Path(r"F:\ONC\output\pathway")
DEG_DIR = Path(r"F:\ONC\output\deg")
TRAJ_DIR = Path(r"F:\ONC\output\trajectory")
OUTPUT.mkdir(parents=True, exist_ok=True)

TIME_ORDER = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]


def run_enrichr(gene_list, description, top_n=15):
    """Run Enrichr enrichment via gseapy."""
    try:
        import gseapy as gp
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Mouse"],
            organism="mouse",
            outdir=None,
            no_plot=True,
        )
        if enr.results is None or len(enr.results) == 0:
            return None
        res = enr.results.copy()
        res["description"] = description
        # Keep top by adjusted p-value
        res = res.sort_values("Adjusted P-value").head(top_n)
        return res
    except Exception as e:
        print(f"  Enrichr failed: {e}")
        return None


def enrich_degs(top_n_genes=200):
    """Run pathway enrichment for top DEGs at each time point."""
    print("=== DEG Pathway Enrichment ===")
    all_results = []

    for tp in TIME_ORDER[1:]:
        deg_path = DEG_DIR / f"{tp}_vs_Ctrl_deg.csv"
        if not deg_path.exists():
            continue
        deg = pd.read_csv(deg_path, index_col=0)
        # Take top genes by absolute logfoldchanges (significant ones)
        sig = deg[deg["pvals_adj"] < 0.05].copy()
        sig = sig.reindex(sig["logfoldchanges"].abs().sort_values(ascending=False).index)
        top_genes = sig.head(top_n_genes).index.tolist()

        print(f"\n  {tp}: {len(sig)} sig DEGs, enriching top {len(top_genes)}...")
        result = run_enrichr(top_genes, f"{tp}_vs_Ctrl")
        if result is not None:
            all_results.append(result)
            # Print top pathways
            for _, row in result.head(8).iterrows():
                overlap = row.get("Overlap", "?")
                print(f"    {row['Term'][:55]:55s} p={row['Adjusted P-value']:.1e} [{overlap}]")
            result.to_csv(OUTPUT / f"pathway_{tp}_vs_Ctrl.csv")

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(OUTPUT / "pathway_all_deg.csv", index=False)
        print(f"\n  Combined results saved.")

    return all_results


def enrich_clusters():
    """Run pathway enrichment for each gene cluster."""
    print("\n=== Cluster Pathway Enrichment ===")
    clusters = pd.read_csv(TRAJ_DIR / "gene_clusters.csv", index_col=0)

    for c in sorted(clusters["cluster"].unique()):
        genes = clusters[clusters["cluster"] == c].index.tolist()
        peak = clusters[clusters["cluster"] == c][TIME_ORDER].mean(axis=0).idxmax()
        print(f"\n  Cluster {c} ({len(genes)} genes, peak={peak}):")
        result = run_enrichr(genes[:500], f"Cluster_{c}")
        if result is not None:
            for _, row in result.head(5).iterrows():
                print(f"    {row['Term'][:55]:55s} p={row['Adjusted P-value']:.1e}")
            result.to_csv(OUTPUT / f"pathway_cluster_{c}.csv")


def simple_gsea_per_time():
    """Simple gene-set scoring using time-course mean expression."""
    print("\n=== Time-course GSEA-like ranking ===")
    clusters = pd.read_csv(TRAJ_DIR / "gene_clusters.csv", index_col=0)

    # For each gene, compute trend: late - early expression
    clusters["trend"] = clusters["2w"] - clusters["Ctrl"]
    clusters["peak_category"] = clusters["cluster"].astype(str)

    # Top increasing and decreasing genes
    increasing = clusters.nlargest(100, "trend")
    decreasing = clusters.nsmallest(100, "trend")

    print(f"  Top increasing (2w > Ctrl): {', '.join(increasing.index[:10])}")
    print(f"  Top decreasing (Ctrl > 2w): {', '.join(decreasing.index[:10])}")

    # Save ranked gene list for GSEA tools
    ranked = clusters[["trend"]].sort_values("trend", ascending=False)
    ranked.to_csv(OUTPUT / "gene_rank_trend.csv")
    print(f"  Saved gene rank: gene_rank_trend.csv")


if __name__ == "__main__":
    print("Starting pathway enrichment analysis...\n")
    enrich_degs(top_n_genes=200)
    enrich_clusters()
    simple_gsea_per_time()
    print("\nPathway analysis complete.")
