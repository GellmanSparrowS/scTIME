# scTIME

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A modular, checkpoint-resumable framework for single-cell RNA-seq time-course analysis.

`mermaid
flowchart TD
    A[Raw Counts] --> B[QC + Normalize]
    B --> C[HVG + PCA + UMAP]
    C --> D1[sc DEG]
    C --> D2[Pseudobulk DEG]
    C --> D3[Gene Clusters]
    C --> D4[Subtype Composition]
    C --> D5[Pathway Enrichment]
    D1 --> E[Figures + Tables]
    D2 --> E
    D3 --> E
    D4 --> E
    D5 --> E
`

## Why scTIME

Most scRNA-seq pipelines treat each time point independently. scTIME is built for time-course experiments: it computes gene trajectories, identifies transition nodes between time points, quantifies cell-type vulnerability, and enriches pathways at each temporal stage.

- **Dual DEG strategy** — single-cell Wilcoxon for sensitivity; pseudobulk t-test with BH correction for sample-level statistical rigor
- **Time-course aware** — k-means clustering of gene expression trajectories; transition gene identification between consecutive time points
- **Cell-type vulnerability scoring** — fold-change-based classification with chi-squared testing
- **Time-resolved GO enrichment** — Enrichr per time point, not just pooled conditions
- **Checkpoint-resumable** — every step saves intermediate .h5ad files; rerunning skips completed work
- **Memory-efficient** — chunked CSV loading, automatic sparse matrix conversion, float32 storage
- **Publication-ready** — Nature-style figures (PNG + PDF) with source CSV data

## Quick Start

`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm

python scripts/pipeline.py          # Merge + QC + normalize + HVG + PCA + UMAP
python scripts/deg_analysis.py      # Single-cell DEG + pseudobulk construction
python scripts/pb_deg.py            # Pseudobulk t-test DEG
python scripts/trajectory.py        # Gene co-expression clustering + transitions
python scripts/subtype_analysis.py  # Cell-type composition + vulnerability scoring
python scripts/pathway_analysis.py  # GO enrichment per time point
python scripts/make_figures.py      # Publication figures
`

## Configuration

All parameters are class-level attributes. Adjust before running:

`python
MIN_GENES = 200        # Min genes per cell
MAX_GENES = 6000       # Max genes per cell
MIN_COUNTS = 500       # Min UMI count
MAX_MITO_PCT = 20.0    # Max mitochondrial fraction
n_top_genes = 3000     # Number of HVGs
n_comps = 50           # PCA components
n_clusters = 6         # Number of gene trajectory clusters
`

## License

MIT
