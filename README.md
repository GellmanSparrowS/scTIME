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

Most scRNA-seq pipelines treat each time point independently. scTIME is built for time-course experiments: it computes gene expression trajectories across time, identifies transition nodes between consecutive time points, quantifies cell-type vulnerability shifts, and enriches pathways at each temporal stage.

| Capability | Implementation |
|---|---|
| Differential expression | Single-cell Wilcoxon for sensitivity; pseudobulk t-test with BH correction for sample-level statistical rigor |
| Time-course clustering | k-means clustering of genome-wide expression trajectories; transition gene identification |
| Cell-type vulnerability | Fold-change-based classification (vulnerable / stable / resilient) with chi-squared testing |
| Pathway enrichment | Enrichr GO Biological Process per time point, not just pooled conditions |
| Checkpoint-resumable | Every step saves intermediate .h5ad files; rerunning skips completed work |
| Memory-efficient | Chunked CSV loading, automatic sparse matrix conversion, float32 storage |
| Figures | Nature-style PNG + PDF with source CSV data |

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

All scripts include `if __name__ == "__main__"` entry points. Data paths and QC thresholds are configured as class-level attributes.

## Configuration

Edit class-level attributes in each script before running:

`python
# pipeline.py -- QC thresholds
MIN_GENES = 200          # Min genes per cell
MAX_GENES = 6000         # Max genes per cell
MIN_COUNTS = 500         # Min UMI count
MAX_MITO_PCT = 20.0      # Max mitochondrial fraction

# trajectory.py -- clustering
n_clusters = 6           # Number of gene trajectory clusters

# pb_deg.py -- statistical test
# Uses Welch t-test with Benjamini-Hochberg correction
`

## License

MIT
