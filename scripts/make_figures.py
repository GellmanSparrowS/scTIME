"""
Generate Nature-style figures for GSE137398 ONC RGC analysis.
All figures saved to output/figures/ with underlying data.
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

warnings.filterwarnings("ignore")

OUTPUT = Path(r"F:\ONC\output\figures")
OUTPUT.mkdir(parents=True, exist_ok=True)
DATA_OUT = Path(r"F:\ONC\output\figure_data")
DATA_OUT.mkdir(parents=True, exist_ok=True)

# Nature style settings
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "figure.dpi": 300,
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

TIME_ORDER = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]
TIME_LABELS = ["Control", "12h", "1d", "2d", "4d", "1w", "2w"]
TP_COLORS = {
    "Ctrl": "#2166AC", "12h": "#92C5DE", "1d": "#D1E5F0",
    "2d": "#F4A582", "4d": "#CA0020", "1w": "#B2182B", "2w": "#67001F",
}
STAGE_COLORS = {"Control": "#2166AC", "Early": "#92C5DE", "Middle": "#D6604D", "Late": "#67001F"}

CACHE_PREPROC = Path(r"F:\ONC\output\preprocessed.h5ad")
TRAJ_DIR = Path(r"F:\ONC\output\trajectory")
DEG_DIR = Path(r"F:\ONC\output\deg")
SUBTYPE_DIR = Path(r"F:\ONC\output\subtype")
PATHWAY_DIR = Path(r"F:\ONC\output\pathway")


def fig_umap(adata):
    """Fig a: UMAP colored by time point and stage."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.2))

    # Check if UMAP exists
    if "X_umap" not in adata.obsm:
        print("  No UMAP coordinates, skipping")
        return

    coords = adata.obsm["X_umap"]
    tp_cat = adata.obs["time_point"].cat.codes

    # Subsample for plotting speed
    n = min(len(coords), 20000)
    idx = np.random.choice(len(coords), n, replace=False)

    for ax, color_by in [(ax1, "time_point"), (ax2, "stage")]:
        cats = adata.obs[color_by].cat.categories
        colors = TP_COLORS if color_by == "time_point" else STAGE_COLORS
        for cat in cats:
            mask = adata.obs[color_by].iloc[idx] == cat
            ax.scatter(
                coords[idx[mask.values], 0], coords[idx[mask.values], 1],
                s=0.3, c=colors.get(cat, "#999999"), label=cat, rasterized=True, alpha=0.7,
            )
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.legend(markerscale=3, frameon=False, loc="lower left" if color_by == "time_point" else "lower right")

    ax1.set_title("a  Time point"); ax2.set_title("b  Stage (Control/Early/Middle/Late)")
    fig.tight_layout()
    fig.savefig(OUTPUT / "fig_umap.png")
    fig.savefig(OUTPUT / "fig_umap.pdf")
    plt.close(fig)
    print("  Saved: fig_umap")


def fig_deg_heatmap():
    """Fig b: Heatmap of top pb-DEGs across time points."""
    # Load pb DEG summaries
    pb_summary = pd.read_csv(DEG_DIR / "pb_deg_summary.csv")

    # Collect top DEGs from each time point
    top_genes = set()
    for tp in TIME_ORDER[1:]:
        try:
            pb = pd.read_csv(DEG_DIR / f"pb_{tp}_vs_Ctrl.csv", index_col=0)
            sig = pb[pb["pval_adj"] < 0.05]
            top_genes.update(sig.nlargest(8, "log2FC").index)
            top_genes.update(sig.nsmallest(8, "log2FC").index)
        except FileNotFoundError:
            continue

    # Get time-course mean expression for these genes
    tc = pd.read_csv(TRAJ_DIR / "gene_clusters.csv", index_col=0)
    top_genes = [g for g in top_genes if g in tc.index]
    if not top_genes:
        print("  No overlapping genes for heatmap")
        return

    heatmap_data = tc.loc[top_genes, TIME_ORDER]
    # Z-score
    heatmap_z = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0).divide(heatmap_data.std(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(5, max(3, len(top_genes) * 0.2)))
    sns.heatmap(heatmap_z, cmap="RdBu_r", center=0, xticklabels=TIME_LABELS,
                yticklabels=True, ax=ax, cbar_kws={"label": "Z-score", "shrink": 0.5},
                linewidths=0.3, linecolor="white")
    ax.set_title("Top pseudobulk DEGs across time points")
    ax.set_xlabel(""); ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(OUTPUT / "fig_deg_heatmap.png")
    fig.savefig(OUTPUT / "fig_deg_heatmap.pdf")
    plt.close(fig)
    # Save data
    heatmap_data.to_csv(DATA_OUT / "deg_heatmap_data.csv")
    print("  Saved: fig_deg_heatmap")


def fig_clusters():
    """Fig c: Gene cluster expression patterns."""
    tc = pd.read_csv(TRAJ_DIR / "gene_clusters.csv", index_col=0)
    clusters = sorted(tc["cluster"].unique())

    n_clusters = len(clusters)
    ncols = 3
    nrows = (n_clusters + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7, 2.2 * nrows))
    axes = axes.flatten()

    for i, c in enumerate(clusters):
        ax = axes[i]
        c_genes = tc[tc["cluster"] == c]
        n_genes = len(c_genes)
        mean_expr = c_genes[TIME_ORDER].mean()
        sem_expr = c_genes[TIME_ORDER].std() / np.sqrt(n_genes)

        ax.fill_between(range(len(TIME_ORDER)), mean_expr - sem_expr, mean_expr + sem_expr,
                        alpha=0.2, color="#2166AC")
        ax.plot(range(len(TIME_ORDER)), mean_expr, "o-", color="#2166AC", ms=3, lw=1)
        ax.set_xticks(range(len(TIME_ORDER)))
        ax.set_xticklabels(TIME_LABELS, rotation=45, ha="right", fontsize=5)
        peak = mean_expr.idxmax()
        ax.set_title(f"Cluster {c} (n={n_genes})", fontsize=7)
        ax.set_ylabel("Mean log-expr", fontsize=6)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Gene expression clusters across ONC time course", fontsize=8, y=1.01)
    fig.tight_layout()
    fig.savefig(OUTPUT / "fig_clusters.png")
    fig.savefig(OUTPUT / "fig_clusters.pdf")
    plt.close(fig)
    # Data
    cluster_summary = pd.DataFrame({
        f"Cluster_{c}": tc[tc["cluster"] == c][TIME_ORDER].mean() for c in clusters
    }).T
    cluster_summary.to_csv(DATA_OUT / "cluster_means.csv")
    print("  Saved: fig_clusters")


def fig_subtype():
    """Fig d: Subtype composition shift (top vulnerable + resilient)."""
    vuln = pd.read_csv(SUBTYPE_DIR / "subtype_vulnerability.csv", index_col=0)
    comp = pd.read_csv(SUBTYPE_DIR / "subtype_composition_pct.csv", index_col=0)

    # Select top 6 vulnerable + top 6 resilient (by fold change)
    top_vuln = vuln[vuln["vulnerability"] == "vulnerable"].nsmallest(6, "fold_change").index
    top_res = vuln[vuln["vulnerability"] == "resilient"].nlargest(6, "fold_change").index
    selected = list(top_vuln) + list(top_res)

    comp_sel = comp.loc[[s for s in selected if s in comp.index]]

    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = np.arange(len(TIME_ORDER))
    width = 0.12
    cmap = plt.cm.RdBu_r

    for i, st in enumerate(comp_sel.index):
        offset = (i - len(comp_sel) / 2 + 0.5) * width
        vals = comp_sel.loc[st, TIME_ORDER].values
        color = cmap(i / len(comp_sel))
        ax.bar(x + offset, vals, width, label=st, color=color, edgecolor="white", linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(TIME_LABELS, rotation=45, ha="right")
    ax.set_ylabel("Proportion (%)")
    ax.set_title("Subtype proportion shifts: vulnerable vs resilient RGC subtypes")
    ax.legend(fontsize=4.5, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT / "fig_subtype.png")
    fig.savefig(OUTPUT / "fig_subtype.pdf")
    plt.close(fig)
    comp_sel.to_csv(DATA_OUT / "subtype_figure_data.csv")
    print("  Saved: fig_subtype")


def fig_pathway():
    """Fig e: Pathway enrichment dotplot across time points."""
    pathway_files = sorted(Path(PATHWAY_DIR).glob("pathway_*_vs_Ctrl.csv"))
    if not pathway_files:
        print("  No pathway files")
        return

    all_terms = []
    for f in pathway_files:
        tp = f.stem.replace("pathway_", "").replace("_vs_Ctrl", "")
        df = pd.read_csv(f, index_col=0)
        if "Adjusted P-value" not in df.columns:
            continue
        df = df[df["Adjusted P-value"] < 0.1].head(4)
        df["time_point"] = tp
        all_terms.append(df)

    if not all_terms:
        return

    combined = pd.concat(all_terms, ignore_index=True)
    # Shorten term names
    combined["Term_short"] = combined["Term"].apply(lambda x: x[:60])
    combined["-log10(p)"] = -np.log10(combined["Adjusted P-value"].clip(lower=1e-10))

    fig, ax = plt.subplots(figsize=(7.5, max(3, len(combined) * 0.25)))
    scatter = ax.scatter(
        combined["time_point"], combined["Term_short"],
        s=combined["-log10(p)"] * 15, c=combined["-log10(p)"],
        cmap="Reds", edgecolors="grey", linewidth=0.3,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("GO Biological Process enrichment across ONC time course")
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5)
    cbar.set_label("-log10(adj. p)")
    fig.tight_layout()
    fig.savefig(OUTPUT / "fig_pathway.png")
    fig.savefig(OUTPUT / "fig_pathway.pdf")
    plt.close(fig)
    combined.to_csv(DATA_OUT / "pathway_figure_data.csv", index=False)
    print("  Saved: fig_pathway")


def fig_pb_summary():
    """Fig f: Pseudobulk DEG summary bar chart."""
    pb_summary = pd.read_csv(DEG_DIR / "pb_deg_summary.csv")
    x = np.arange(len(pb_summary))
    width = 0.35

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(x - width/2, pb_summary["n_up"], width, label="Upregulated", color="#CA0020")
    ax.bar(x + width/2, pb_summary["n_down"], width, label="Downregulated", color="#2166AC")
    ax.set_xticks(x)
    ax.set_xticklabels(pb_summary["time_point"])
    ax.set_ylabel("Number of DEGs (FDR < 0.05)")
    ax.set_title("Pseudobulk DEGs per time point vs Control")
    ax.legend(frameon=False)
    # Add count labels
    for i, row in pb_summary.iterrows():
        if row["n_up"] > 0:
            ax.text(i - width/2, row["n_up"] + 5, str(row["n_up"]), ha="center", fontsize=5)
        if row["n_down"] > 0:
            ax.text(i + width/2, row["n_down"] + 5, str(row["n_down"]), ha="center", fontsize=5)
    fig.tight_layout()
    fig.savefig(OUTPUT / "fig_pb_summary.png")
    fig.savefig(OUTPUT / "fig_pb_summary.pdf")
    plt.close(fig)
    pb_summary.to_csv(DATA_OUT / "pb_summary_figure_data.csv", index=False)
    print("  Saved: fig_pb_summary")


def main():
    print("Generating Nature-style figures...\n")

    try:
        adata = ad.read_h5ad(CACHE_PREPROC)
        print(f"Loaded: {adata.n_obs} cells")
    except Exception:
        print("No preprocessed AnnData, skipping UMAP")
        adata = None

    if adata is not None:
        fig_umap(adata)
    fig_deg_heatmap()
    fig_clusters()
    fig_subtype()
    fig_pathway()
    fig_pb_summary()

    print(f"\nFigures saved to: {OUTPUT}")
    print(f"Figure data saved to: {DATA_OUT}")


if __name__ == "__main__":
    main()
