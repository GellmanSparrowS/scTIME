# scTIME &middot; [中文](#中文文档)

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![scanpy](https://img.shields.io/badge/scanpy-1.12-red)](https://scanpy.readthedocs.io/)

> A modular, checkpoint-resumable framework for single-cell RNA-seq time-course analysis. Goes from raw counts to publication figures in one pipeline.

`mermaid
flowchart LR
    A[Raw Counts] --> B[QC Filter]
    B --> C[Normalize + HVG]
    C --> D[PCA + UMAP]
    D --> E1[sc DEG<br/>Wilcoxon]
    D --> E2[pb DEG<br/>t-test]
    D --> E3[Gene Clusters<br/>KMeans]
    D --> E4[Subtype<br/>Composition]
    D --> E5[Pathway<br/>Enrichr]
    E1 --> F[Figures + Tables]
    E2 --> F
    E3 --> F
    E4 --> F
    E5 --> F
`

---

## Why scTIME?

Most scRNA-seq pipelines treat each time point independently. **scTIME is built for time-course experiments**: it computes gene trajectories, identifies transition nodes between time points, quantifies cell-type vulnerability across conditions, and enriches pathways at each temporal stage.

- **Dual DEG strategy**: single-cell Wilcoxon for sensitivity; pseudobulk t-test with BH correction for sample-level statistical rigor — no pseudo-replication
- **Time-course aware**: k-means clustering of 40K+ gene expression trajectories; transition gene identification between consecutive time points
- **Cell-type vulnerability scoring**: fold-change-based classification (vulnerable / stable / resilient) with chi-squared testing
- **Time-resolved GO enrichment**: Enrichr per time point, not just pooled conditions
- **Checkpoint-resumable**: every step saves intermediate .h5ad files; rerunning skips completed work — saves hours on large datasets
- **Memory-efficient**: chunked CSV loading → automatic sparse matrix conversion → loat32 storage
- **Publication-ready**: Nature-style figures (PNG + PDF) with corresponding CSV source data

## Quick Start

`ash
pip install pandas numpy scipy anndata scanpy \
            statsmodels gseapy scikit-learn \
            seaborn matplotlib tqdm

python scripts/pipeline.py          # QC + preprocessing (resumable)
python scripts/deg_analysis.py      # Single-cell DEG + pseudobulk
python scripts/pb_deg.py            # Pseudobulk t-test DEG
python scripts/trajectory.py        # Gene co-expression clustering
python scripts/subtype_analysis.py  # Cell-type vulnerability scoring
python scripts/pathway_analysis.py  # GO enrichment per time point
python scripts/make_figures.py      # Publication figures
`

All scripts include `if __name__ == "__main__"` entry points and can run independently or chained.

## Pipeline Architecture

| Step | Script | Output |
|------|---------|--------|
| Merge + QC | pipeline.py | merged.h5ad, qc_filtered.h5ad |
| Normalize + HVG + PCA + UMAP | pipeline.py | 
ormalized_hvg.h5ad, preprocessed.h5ad |
| Single-cell DEG | deg_analysis.py | deg/*_vs_Ctrl_deg.csv |
| Pseudobulk DEG | pb_deg.py | deg/pb_*_vs_Ctrl.csv |
| Gene clustering | 	rajectory.py | 	rajectory/gene_clusters.csv |
| Subtype analysis | subtype_analysis.py | subtype/subtype_vulnerability.csv |
| Pathway enrichment | pathway_analysis.py | pathway/pathway_*_vs_Ctrl.csv |
| Figures | make_figures.py | igures/*.png, igure_data/*.csv |

## Configuration

All parameters are class-level attributes. Adjust before running:

`python
# pipeline.py
class ONCPipeline:
    MIN_GENES = 200       # Min genes per cell
    MAX_GENES = 6000      # Max genes per cell
    MIN_COUNTS = 500      # Min UMI count
    MAX_MITO_PCT = 20.0   # Max mitochondrial %
    MITO_PREFIX = mt-   # Mouse MT gene prefix

# trajectory.py
    n_clusters = 6        # Number of gene trajectory clusters

# deg_analysis.py / pb_deg.py
    method = wilcoxon   # sc DEG method
    # pb DEG uses Welch t-test with BH correction
`

## Output Structure

`
output/
├── merged.h5ad              # Raw merged counts
├── qc_filtered.h5ad         # Post-QC
├── normalized_hvg.h5ad      # Normalized, all genes
├── preprocessed.h5ad        # PCA + UMAP, HVG subset only
├── deg/                     # Differential expression tables
├── trajectory/              # Gene clusters + transitions
├── subtype/                 # Composition + vulnerability scores
├── pathway/                 # GO enrichment results
├── figures/                 # Publication figures (PNG + PDF)
└── figure_data/             # Underlying CSV per figure
`

## Documentation

See [SKILL.md](SKILL.md) for detailed methods, memory management patterns, and best practices developed during this project.

---

## 中文文档

scTIME 是一套面向单细胞 RNA-seq 时间序列实验的模块化分析框架。与大多数将每个时间点独立处理的流程不同，scTIME 专为时间序列设计：计算基因轨迹、识别时间点间转换节点、量化跨条件的细胞类型脆弱性，并在每个时间阶段进行通路富集。

**核心能力**：单细胞与伪Bulk双重差异表达策略 &middot; 全基因组表达轨迹 k 均值聚类 &middot; 基于倍变化的细胞类型脆弱性评分 &middot; Enrichr 时间分辨 GO 富集 &middot; 断点续跑 &middot; 内存高效的稀疏矩阵处理 &middot; Nature 风格出版图表 + CSV 源数据。

**快速开始**：
`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm
python scripts/pipeline.py     # 预处理（支持断点续跑）
python scripts/pb_deg.py       # 伪Bulk差异分析
python scripts/trajectory.py   # 基因聚类
`

所有脚本均包含独立入口点。详细方法和经验总结见 [SKILL.md](SKILL.md)。

## License

MIT &mdash; see [LICENSE](LICENSE)
