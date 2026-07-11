# scTIME &middot; Single-Cell Time-course Analysis Pipeline &middot; [中文](#中文)

A modular, checkpoint-resumable computational framework for analyzing single-cell RNA-seq time-course experiments. Designed for multi-condition, multi-timepoint studies with built-in support for quality control, differential expression at both single-cell and pseudobulk resolution, gene trajectory clustering, cell-type composition analysis, and pathway enrichment.

## Features

- **Dual DEG strategy**: single-cell Wilcoxon for sensitivity; pseudobulk t-test with BH correction for sample-level statistical rigor
- **Time-course aware**: k-means clustering of genome-wide expression trajectories; transition gene identification between consecutive time points
- **Cell-type vulnerability scoring**: quantitative composition shift analysis with chi-squared testing and fold-change-based classification
- **Time-resolved pathway enrichment**: GO Biological Process enrichment per time point via Enrichr
- **Checkpoint-resumable**: every major step saves intermediate AnnData objects; rerunning skips completed work
- **Memory-efficient**: chunked CSV loading with automatic sparse matrix conversion
- **Publication-ready**: Nature-style figures (PNG + PDF) with underlying CSV data

## Installation

`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm
`

## Usage

`ash
python scripts/pipeline.py          # Merge + QC + normalize + HVG + PCA + UMAP
python scripts/deg_analysis.py      # Single-cell DEG + pseudobulk construction
python scripts/pb_deg.py            # Pseudobulk t-test DEG (sample-level)
python scripts/trajectory.py        # Gene co-expression clustering + transitions
python scripts/subtype_analysis.py  # Cell-type composition + vulnerability scoring
python scripts/pathway_analysis.py  # GO enrichment per time point
python scripts/make_figures.py      # Publication figures
`

All scripts include `if __name__ == "__main__"` entry points and can be run independently or as a chained pipeline. Data paths are configured through class-level attributes.

## Repository Structure

`
scripts/          Analysis modules (7 scripts, independently runnable)
output/           Generated results (excluded from version control)
  deg/            Differential expression tables
  trajectory/     Gene clusters and transition genes
  subtype/        Cell-type composition and vulnerability scores
  pathway/        GO enrichment results
  figures/        Publication figures (PNG + PDF)
  figure_data/    Underlying CSV data for each figure
`

## Configuration

Key parameters are class-level attributes and can be adjusted before execution:

| Parameter | Default | Description |
|-----------|---------|-------------|
| MIN_GENES / MAX_GENES | 200 / 6000 | Gene count thresholds for QC |
| MAX_MITO_PCT | 20.0 | Mitochondrial percentage cutoff |
| n_top_genes | 3000 | Number of highly variable genes |
| n_comps | 50 | PCA components |
| n_clusters | 6 | Gene trajectory clusters |

## License

MIT

---

## 中文 {#中文}

scTIME 是一套模块化、支持断点续跑的单细胞 RNA-seq 时间序列分析框架。适用于多条件、多时间点的实验设计，内置质量控制、单细胞与伪Bulk双重差异表达、基因轨迹聚类、细胞类型组成分析和通路富集功能。

### 核心特性

- **双重差异策略**：单细胞级 Wilcoxon（高灵敏度）+ 伪Bulk级 t 检验及 BH 校正（样本级统计严谨性）
- **时间序列感知**：全基因组表达轨迹的 k 均值聚类；相邻时间点间转换基因鉴定
- **细胞类型脆弱性评分**：定量组成偏移分析，结合卡方检验和倍变化分类
- **时间分辨通路富集**：通过 Enrichr 逐时间点进行 GO Biological Process 富集
- **断点续跑**：每个主要步骤保存中间 AnnData 对象，重复运行自动跳过已完成步骤
- **内存高效**：分块 CSV 读取并自动转换为稀疏矩阵
- **出版级图表**：Nature 风格图表（PNG + PDF）及对应的 CSV 源数据

### 安装与运行

`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm
python scripts/pipeline.py     # 全流程预处理
python scripts/pb_deg.py       # 伪Bulk差异分析
python scripts/trajectory.py   # 基因聚类
`

所有脚本均包含入口点，可独立运行或链式执行。数据路径通过类属性配置。
