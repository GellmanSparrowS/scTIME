# scTIME: Single-Cell Time-course Analysis Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A modular, checkpoint-resumable pipeline for analyzing single-cell RNA-seq time-course experiments. Designed for multi-condition, multi-timepoint studies with built-in support for quality control, differential expression, gene trajectory clustering, cell-type composition analysis, and pathway enrichment.

## Features

- **Checkpoint-resumable**: every major step saves intermediate .h5ad files; re-running skips completed work
- **Memory-efficient**: chunked CSV loading with automatic conversion to sparse matrices
- **Dual DEG strategy**: single-cell-level Wilcoxon for sensitivity + pseudobulk-level t-test for statistical rigor
- **Time-course aware**: gene clustering by expression trajectory, transition gene identification
- **Cell-type vulnerability**: composition shift analysis with chi-squared testing
- **Pathway enrichment**: integrated Enrichr support via gseapy
- **Publication figures**: Nature-style PNG + PDF with underlying CSV data

## Pipeline Steps

`
Raw counts → QC filter → Normalize → HVG selection → PCA → UMAP
    │
    ├── Single-cell DEG (Wilcoxon)
    ├── Pseudobulk DEG (t-test)
    ├── Gene time-course clustering (KMeans)
    ├── Cell-type proportion analysis
    └── GO pathway enrichment (Enrichr)
`

## Quick Start

### Requirements
`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm
`

### Usage
`ash
python scripts/pipeline.py          # QC + preprocessing (resumable)
python scripts/deg_analysis.py      # Single-cell DEG + pseudobulk
python scripts/pb_deg.py            # Pseudobulk t-test DEG
python scripts/trajectory.py        # Gene clustering + transitions
python scripts/subtype_analysis.py  # Cell-type composition
python scripts/pathway_analysis.py  # GO enrichment
python scripts/make_figures.py      # Publication figures
`

All scripts include a `if __name__ == "__main__"` entry point and can be run independently.

## Repository Structure

`
├── scripts/                  # Analysis scripts (independently runnable)
├── manuscript/               # Manuscript draft (local only)
├── SKILL.md                  # Methods reference & best practices
├── PROGRESS.md               # Analysis progress log
└── README.md / README_zh.md  # Documentation
`

## Output Structure

`
output/
├── merged.h5ad               # Raw merged counts
├── qc_filtered.h5ad          # After QC
├── normalized_hvg.h5ad       # Normalized + HVG
├── preprocessed.h5ad         # PCA + UMAP, HVG subset
├── deg/                      # Differential expression tables
├── trajectory/               # Gene clusters & transitions
├── subtype/                  # Cell-type composition analysis
├── pathway/                  # GO enrichment results
├── figures/                  # PNG + PDF figures
└── figure_data/              # Underlying data for figures
`

## Configuration

Key parameters are defined as class-level attributes in `ONCPipeline` and can be adjusted before running:

- `MIN_GENES` / `MAX_GENES` — gene count thresholds
- `MAX_MITO_PCT` — mitochondrial percentage cutoff
- `n_top_genes` — number of HVGs
- `n_clusters` — number of gene trajectory clusters

## License

MIT
